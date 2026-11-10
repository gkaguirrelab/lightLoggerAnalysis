"""Browse FLIC eye recordings and preview the existing Python cleaner locally.

Run this file with the same Python environment used for cleanEyeTracking.py.
The browser sends parameter values; the existing cleaner does all numerical
processing and plotting. Recording files are only ever opened for reading.
"""

import argparse
import inspect
import json
import math
import pickle
import re
import sys
import subprocess
from pathlib import Path

import matplotlib

# WebAgg serves the original Matplotlib canvas, including its navigation tools.
# Choose the backend before importing the cleaner, which imports pyplot.
matplotlib.use("WebAgg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_webagg import WebAggApplication
import tornado.ioloop
import tornado.web

WEB_DIRECTORY: Path = Path(__file__).resolve().parent
sys.path.insert(0, str(WEB_DIRECTORY.parent / "pythonCode"))
import cleanEyeTracking


# These labels describe the metadata fields, not a second implementation of
# their behavior. Parameter defaults and validation still come from Python.
PARAMETER_LABELS: dict[str, str] = {
    "pupil_diameter_thresholds": "Pupil diameter bounds (mm)",
    "pupil_change_threshold": "Maximum pupil change per sample (mm)",
    "aperture_change_threshold": "Maximum aperture change per sample (mm)",
    "blink_pre_ns": "Padding before a blink (ns)",
    "blink_post_ns": "Padding after a blink (ns)",
    "gaze_x_upper_bound": "Maximum gaze X (pixels)",
    "gaze_y_upper_bound": "Maximum gaze Y (pixels)",
    "gaze_spike_threshold": "Maximum gaze change per sample (pixels)",
    "dilation_samples": "Mask dilation window (samples)",
    "general_max_gap": "General fill gap (sample-index distance)",
    "blink_max_gap": "Blink fill gap (sample-index distance)",
    "smoothing_samples": "Smoothing window (positive odd samples)",
}
PARAMETER_DEFAULTS: dict = {
    name: inspect.signature(cleanEyeTracking.clean_eye_data).parameters[name].default
    for name in PARAMETER_LABELS
}


class MetadataUnpickler(pickle.Unpickler):
    """Read the cleaner's plain dictionary without executing pickle globals."""

    def find_class(self, module: str, name: str):
        raise ValueError("Cleaning metadata must contain only plain Python values.")


def validate_parameters(parameters: dict) -> dict:
    """Check browser/metadata values before passing them to the cleaner.

    This validates input types and ranges only. The cleaner remains the sole
    source of sample selection, filtering, interpolation, and smoothing.
    """
    if not isinstance(parameters, dict) or set(parameters) != set(PARAMETER_DEFAULTS):
        raise ValueError("Cleaning metadata must contain all 13 supported parameter fields.")

    validated: dict = {}
    for name, default in PARAMETER_DEFAULTS.items():
        value = parameters[name]
        if isinstance(default, tuple):
            if not isinstance(value, (list, tuple)) or len(value) != 2:
                raise ValueError(f"{name} needs a lower and upper bound.")
            values = list(value)
        else:
            values = [value]

        for number in values:
            if isinstance(number, bool) or not isinstance(number, (int, float)):
                raise ValueError(f"{name} must be numeric.")
            if not math.isfinite(number) or number < 0:
                raise ValueError(f"{name} must be finite and nonnegative.")

        if isinstance(default, int):
            if int(value) != value:
                raise ValueError(f"{name} must be a whole number.")
            value = int(value)
        if isinstance(default, tuple):
            if value[0] > value[1]:
                raise ValueError("Pupil diameter lower bound must not exceed the upper bound.")
            value = tuple(value)
        validated[name] = value
    return validated


def read_parameters(eye_data_directory: Path) -> dict:
    """Load the settings saved for this recording without changing the file."""
    with (eye_data_directory / "eye_cleaning_metadata.pkl").open("rb") as metadata_file:
        parameters = MetadataUnpickler(metadata_file).load()
    return validate_parameters(parameters)


def list_recordings(data_directory: Path) -> dict:
    """Return subject/activity/attempt choices beneath an experiment directory."""
    if not data_directory.is_dir():
        raise ValueError("Choose the experiment directory containing the FLIC subject folders.")
    subjects: dict = {}
    for subject_directory in sorted(data_directory.iterdir()):
        if not subject_directory.is_dir() or not re.fullmatch(r"FLIC_\d+", subject_directory.name):
            continue
        activities: dict = {}
        for activity_directory in sorted(subject_directory.iterdir()):
            if not activity_directory.is_dir():
                continue
            recordings = sorted(
                (directory.name for directory in activity_directory.iterdir()
                 if directory.is_dir() and directory.name.isdigit()
                 and (directory / "eyeData").is_dir()),
                key=int,
            )
            if recordings:
                activities[activity_directory.name] = recordings
        if activities:
            subjects[subject_directory.name] = activities
    if not subjects:
        raise ValueError("No FLIC_<subject>/<activity>/<recording>/eyeData directories found.")
    return subjects


def resolve_recording(data_directory: Path, selection: dict) -> Path:
    """Resolve one selected eyeData folder, keeping reads inside the chosen root."""
    parts = [selection.get(name, "") for name in ("subject", "activity", "recording")]
    if any(not isinstance(part, str) or not part or Path(part).name != part
           or part in (".", "..") for part in parts):
        raise ValueError("Select a subject, activity, and recording.")
    eye_directory = data_directory.joinpath(*parts, "eyeData").resolve()
    if not eye_directory.is_relative_to(data_directory) or not eye_directory.is_dir():
        raise ValueError("The selected recording is not inside the loaded directory.")
    # A symlinked file must not redirect a recording read outside this root.
    for filename in ("3d_eye_states.csv", "blinks.csv", "alternative_camera_gaze.csv",
                     "eye_cleaning_metadata.pkl"):
        file_path = eye_directory / filename
        if not file_path.resolve().is_relative_to(data_directory) or not file_path.is_file():
            raise ValueError(f"Recording is missing a local {filename} file.")
    return eye_directory


def choose_data_directory(initial_directory: str) -> str:
    """Run the native picker in a short-lived process so Tk cannot linger."""
    # Keep Tk out of the long-running web server. When selection or cancellation
    # finishes, this process exits and the OS removes all of its native windows.
    result = subprocess.run(
        [sys.executable, str(WEB_DIRECTORY / "folder_picker.py"), initial_directory],
        capture_output=True, text=True, check=True,
    )
    return json.loads(result.stdout)


def zoom_time_axis(event) -> None:
    """Zoom all six shared time axes around the cursor when the wheel scrolls."""
    if event.inaxes is None or event.xdata is None:
        return
    axis = event.inaxes
    left, right = axis.get_xlim()
    scale = 0.8 if event.button == "up" else 1.25
    axis.set_xlim(event.xdata + (left - event.xdata) * scale,
                  event.xdata + (right - event.xdata) * scale)
    axis.figure.canvas.draw_idle()


class ApiHandler(tornado.web.RequestHandler):
    """Handle local, same-origin requests for directory selection and previews."""

    def prepare(self) -> None:
        # Bind to loopback and require same-origin JSON requests so an unrelated
        # website cannot ask the local server to read directories on this machine.
        if self.request.host not in self.application.allowed_hosts:
            raise tornado.web.HTTPError(403)
        if self.request.method == "POST":
            expected_origin = f"http://{self.request.host}"
            if self.request.headers.get("Origin") != expected_origin:
                raise tornado.web.HTTPError(403)
            if not self.request.headers.get("Content-Type", "").startswith("application/json"):
                raise tornado.web.HTTPError(415)

    def get(self) -> None:
        self.write({"directory": self.application.initial_directory,
                    "labels": PARAMETER_LABELS, "defaults": PARAMETER_DEFAULTS})

    def post(self) -> None:
        # Requests run sequentially on the local event loop. This deliberately
        # avoids concurrent calls into Matplotlib's process-wide figure state.
        previous_figures = set(plt.get_fignums())
        try:
            request = json.loads(self.request.body)
            if not isinstance(request, dict):
                raise ValueError("Expected a JSON object.")
            if self.request.path in ("/api/directory", "/api/choose-directory"):
                if self.request.path == "/api/choose-directory":
                    starting_directory = str(self.application.data_directory or self.application.initial_directory)
                    chosen_directory = choose_data_directory(starting_directory)
                    if not chosen_directory:
                        self.write({"cancelled": True})
                        return
                else:
                    chosen_directory = request["directory"]
                directory = Path(chosen_directory).expanduser().resolve()
                subjects = list_recordings(directory)
                self.application.data_directory = directory
                self.write({"directory": str(directory), "subjects": subjects})
                return

            if self.application.data_directory is None:
                raise ValueError("Load a data directory first.")
            eye_directory = resolve_recording(self.application.data_directory, request)
            saved_parameters = read_parameters(eye_directory)
            parameters = saved_parameters if request.get("parameters") is None else validate_parameters(request["parameters"])

            # Reuse both the existing cleaner AND its plotting helper. The two
            # preview options suppress PDF writes and desktop windows; no CSV,
            # metadata, or other recording output is written by this call.
            cleaned_data, effective_parameters = cleanEyeTracking.clean_eye_data(
                eye_directory, eye_directory, verbose=False,
                plot_output_path=None, **parameters,
            )
            figure = plt.gcf()
            # The page already identifies the recording; omit the repeated
            # figure heading in this browser preview only.
            figure.suptitle("")
            figure.set_label("")
            # Let long panel titles wrap when the sidebar leaves a narrow canvas.
            for axis in figure.axes:
                axis.title.set_wrap(True)
            # Match the page's 70% scale without changing the plotted samples
            # or the original cleaner's PDF/desktop figure styling.
            for text in figure.findobj(match=matplotlib.text.Text):
                text.set_fontsize(text.get_fontsize() * 0.7)
            figure.canvas.mpl_connect("scroll_event", zoom_time_axis)
            figure.canvas.manager.toolbar.push_current()
            self.write({"figure_url": f"/plot/{figure.number}",
                        "parameters": effective_parameters,
                        "saved_parameters": saved_parameters,
                        "sample_count": len(cleaned_data)})

            # Retain just the displayed preview. A failed recalculation leaves
            # the previous figure available rather than replacing it with blank.
            for figure_number in previous_figures:
                plt.close(figure_number)
        except (ValueError, KeyError, TypeError, OSError, pickle.UnpicklingError, EOFError, subprocess.SubprocessError) as error:
            for figure_number in set(plt.get_fignums()) - previous_figures:
                plt.close(figure_number)
            self.set_status(400)
            self.write({"error": str(error)})


class ViewerApplication(WebAggApplication):
    """Serve local interface files alongside Matplotlib's own browser viewer."""

    def __init__(self, port: int, initial_directory: str = ""):
        super().__init__(url_prefix="/plot")
        self.data_directory: Path | None = None
        self.initial_directory = initial_directory
        self.allowed_hosts = {f"127.0.0.1:{port}", f"localhost:{port}"}
        self.add_handlers(r".*", [
            (r"/api/(?:config|directory|choose-directory|plot)", ApiHandler),
            (r"/plot/([0-9]+)", FigurePage, {"url_prefix": "/plot"}),
            (r"/(|index.html|app.js|style.css)", tornado.web.StaticFileHandler,
             {"path": str(WEB_DIRECTORY), "default_filename": "index.html"}),
        ])


class FigurePage(WebAggApplication.SingleFigurePage):
    """Use a responsive shell around Matplotlib's unmodified figure canvas."""

    def get_template_path(self) -> str:
        return str(WEB_DIRECTORY)


def main() -> None:
    """Start the local viewer; no data is copied or saved on startup."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--port", type=int, default=8765)
    parser.add_argument("--data-directory", default="", help="Experiment folder containing FLIC subject folders.")
    arguments = parser.parse_args()
    application = ViewerApplication(arguments.port, arguments.data_directory)
    application.listen(arguments.port, address="127.0.0.1", max_buffer_size=1024 * 1024)
    print(f"FLIC eye viewer: http://127.0.0.1:{arguments.port}", flush=True)
    print("Press Ctrl+C to stop. Recording files are read-only.", flush=True)
    try:
        tornado.ioloop.IOLoop.current().start()
    except KeyboardInterrupt:
        plt.close("all")


if __name__ == "__main__":
    main()
