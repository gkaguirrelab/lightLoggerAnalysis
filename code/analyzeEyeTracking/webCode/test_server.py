"""Check preview parity, directory navigation, and preservation of source files."""

import hashlib
import json
import pickle
import tempfile
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd
from tornado.testing import AsyncHTTPTestCase

import server


def create_example_recording(eye_directory: Path) -> None:
    """Create a small complete eyeData fixture without using participant data."""
    eye_directory.mkdir(parents=True)
    sample_numbers = np.arange(240)
    timestamps = 1_790_000_000_000_000_000 + sample_numbers * 5_000_000
    pupil = 4 + 0.1 * np.sin(sample_numbers / 10)
    # An abrupt but physiological excursion makes threshold changes measurable.
    pupil[100:110] = 6
    pd.DataFrame({
        "timestamp [ns]": timestamps,
        "pupil diameter left [mm]": pupil,
        "pupil diameter right [mm]": pupil + 0.1,
        "eyelid aperture left [mm]": 8 + 0.1 * np.sin(sample_numbers / 12),
        "eyelid aperture right [mm]": 7.8 + 0.1 * np.sin(sample_numbers / 12),
    }).to_csv(eye_directory / "3d_eye_states.csv", index=False)
    gaze = pd.DataFrame({"timestamp [ns]": timestamps,
                         "gaze x [px]": 300 + 10 * np.sin(sample_numbers / 15),
                         "gaze y [px]": 200 + 10 * np.cos(sample_numbers / 15)})
    gaze.to_csv(eye_directory / "alternative_camera_gaze.csv", index=False)
    gaze.to_csv(eye_directory / "gaze.csv", index=False)
    pd.DataFrame({"start timestamp [ns]": [timestamps[40]],
                  "end timestamp [ns]": [timestamps[45]]}).to_csv(eye_directory / "blinks.csv", index=False)
    with (eye_directory / "eye_cleaning_metadata.pkl").open("wb") as output:
        pickle.dump(server.PARAMETER_DEFAULTS, output)
    # Sentinels must remain byte-for-byte unchanged during every preview.
    (eye_directory / "cleaned_eye_data.csv").write_text("existing cleaned output\n")
    (eye_directory / "eye_cleaned_vs_raw.pdf").write_bytes(b"existing PDF")
    (eye_directory / "egocentric_mapper_metadata.pkl").write_bytes(b"existing mapper metadata")


def snapshot_files(directory: Path) -> dict:
    """Record contents and modification times, including the list of files."""
    return {str(path.relative_to(directory)): (hashlib.sha256(path.read_bytes()).hexdigest(), path.stat().st_mtime_ns)
            for path in directory.rglob("*") if path.is_file()}


class ViewerTests(AsyncHTTPTestCase):
    """Exercise the same HTTP endpoints used by the browser."""

    def get_app(self):
        return server.ViewerApplication(8765)

    def setUp(self):
        super().setUp()
        self._app.allowed_hosts.add(f"127.0.0.1:{self.get_http_port()}")
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary_directory.name)
        self.eye_directory = self.root / "FLIC_0021" / "chat" / "1" / "eyeData"
        create_example_recording(self.eye_directory)
        self.selection = {"subject": "FLIC_0021", "activity": "chat", "recording": "1"}

    def tearDown(self):
        server.plt.close("all")
        self.temporary_directory.cleanup()
        super().tearDown()

    def post(self, path: str, payload: dict):
        """Send a same-origin request as the local page does."""
        return self.fetch(path, method="POST", body=json.dumps(payload), headers={
            "Content-Type": "application/json", "Origin": self.get_url("").rstrip("/")})

    def test_preview_matches_existing_figure_without_writes(self):
        """All plotted samples/styles match Python, including after adjustment."""
        original_files = snapshot_files(self.root)
        response = self.post("/api/directory", {"directory": str(self.root)})
        self.assertEqual(response.code, 200)
        self.assertEqual(json.loads(response.body)["subjects"], {"FLIC_0021": {"chat": ["1"]}})

        server.cleanEyeTracking.clean_eye_data(self.eye_directory, self.eye_directory,
                                               verbose=False, plot_output_path=None)
        original_figure = server.plt.gcf()
        expected_lines = [[(line.get_xdata().copy(), line.get_ydata().copy(), line.get_color(), line.get_marker())
                           for line in axis.lines] for axis in original_figure.axes]
        expected_titles = [axis.get_title() for axis in original_figure.axes]
        server.plt.close(original_figure)

        response = self.post("/api/plot", self.selection)
        self.assertEqual(response.code, 200, response.body)
        result = json.loads(response.body)
        figure = server.plt.gcf()
        self.assertEqual(len(figure.axes), 6)
        self.assertEqual([axis.get_title() for axis in figure.axes], expected_titles)
        for axis, lines in zip(figure.axes, expected_lines):
            self.assertEqual(len(axis.lines), len(lines))
            for line, (x_values, y_values, color, marker) in zip(axis.lines, lines):
                np.testing.assert_array_equal(line.get_xdata(), x_values)
                np.testing.assert_array_equal(line.get_ydata(), y_values)
                self.assertEqual((line.get_color(), line.get_marker()), (color, marker))
        self.assertEqual(self.fetch(result["figure_url"]).code, 200)
        self.assertEqual(self.fetch("/plot/js/mpl.js").code, 200)
        self.assertEqual(self.fetch("/").code, 200)
        self.assertEqual(self.fetch("/server.py").code, 404)

        original_cleaned_pupil = figure.axes[0].lines[1].get_ydata().copy()
        parameters = dict(server.PARAMETER_DEFAULTS, pupil_diameter_thresholds=[1.5, 4.05])
        response = self.post("/api/plot", {**self.selection, "parameters": parameters})
        self.assertEqual(response.code, 200, response.body)
        self.assertFalse(np.allclose(original_cleaned_pupil, server.plt.gcf().axes[0].lines[1].get_ydata(), equal_nan=True))
        self.assertEqual(len(server.plt.get_fignums()), 1)
        self.assertEqual(snapshot_files(self.root), original_files)

    def test_invalid_parameters_and_paths_preserve_files(self):
        """Reject invalid values and traversal without changing the recording."""
        before = snapshot_files(self.root)
        self.post("/api/directory", {"directory": str(self.root)})
        for parameters in [dict(server.PARAMETER_DEFAULTS, smoothing_samples=2),
                           dict(server.PARAMETER_DEFAULTS, blink_pre_ns=-1),
                           dict(server.PARAMETER_DEFAULTS, pupil_diameter_thresholds=[8, 2])]:
            self.assertEqual(self.post("/api/plot", {**self.selection, "parameters": parameters}).code, 400)
        self.assertEqual(self.post("/api/plot", {**self.selection, "activity": ".."}).code, 400)
        (self.eye_directory / "eye_cleaning_metadata.pkl").unlink()
        self.assertEqual(self.post("/api/plot", self.selection).code, 400)
        before.pop(str((self.eye_directory / "eye_cleaning_metadata.pkl").relative_to(self.root)))
        self.assertEqual(snapshot_files(self.root), before)

    def test_display_and_save_are_independent(self):
        """None never saves; explicit paths save even when display is disabled."""
        original_files = snapshot_files(self.root)
        with patch.object(server.plt, "show") as show_figure:
            server.cleanEyeTracking.clean_eye_data(self.eye_directory, self.eye_directory)
        show_figure.assert_called_once_with(block=False)
        self.assertEqual(snapshot_files(self.root), original_files)

        output_path = self.root / "preview" / "comparison.pdf"
        with patch.object(server.plt, "show") as show_figure:
            server.cleanEyeTracking.clean_eye_data(
                self.eye_directory, self.eye_directory,
                verbose=False, plot_output_path=output_path,
            )
        show_figure.assert_not_called()
        self.assertTrue(output_path.read_bytes().startswith(b"%PDF"))

    def test_folder_picker_selection_and_cancel(self):
        """A selected folder loads choices; cancelling preserves that folder."""
        before = snapshot_files(self.root)
        with patch.object(server, "choose_data_directory", return_value=str(self.root)):
            response = self.post("/api/choose-directory", {})
        self.assertEqual(response.code, 200)
        self.assertIn("FLIC_0021", json.loads(response.body)["subjects"])
        selected_directory = self._app.data_directory
        with patch.object(server, "choose_data_directory", return_value=""):
            response = self.post("/api/choose-directory", {})
        self.assertEqual(json.loads(response.body), {"cancelled": True})
        self.assertEqual(self._app.data_directory, selected_directory)
        self.assertEqual(snapshot_files(self.root), before)


if __name__ == "__main__":
    import unittest
    unittest.main()
