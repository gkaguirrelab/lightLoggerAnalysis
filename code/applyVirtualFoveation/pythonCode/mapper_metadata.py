"""Capture the mapper's Python console output and extract alignment metadata."""

from contextlib import contextmanager, redirect_stderr, redirect_stdout
import re
import sys
import tempfile
from typing import Any, Iterator, TextIO, TypedDict


class TemporalAlignment(TypedDict):
    offset_seconds: float
    pearson_correlation: float


class ParsedMapperMetadata(TypedDict):
    schema_version: int
    temporal_alignment: TemporalAlignment
    temporal_alignment_estimates: list[TemporalAlignment]
    metric_description: str


class MapperInputs(TypedDict):
    neon_timeseries_dir: str
    alternative_vid_path: str


class MapperSettings(TypedDict):
    mapping_choice: str
    refresh_time_threshold_sec: float
    render_video: bool
    render_video_comparison: bool
    optic_flow_algorithm: str
    image_matcher: str


class MapperMetadata(ParsedMapperMetadata):
    started_at_utc: str
    elapsed_seconds: float
    status: str
    inputs: MapperInputs
    settings: MapperSettings


class _Tee:
    """Forward each text write to both the live console and the transcript."""

    def __init__(self, console: TextIO, transcript: TextIO) -> None:
        self.console: TextIO = console
        self.transcript: TextIO = transcript

    def write(self, text: str) -> int:
        # Persist each fragment immediately, then preserve the console write result.
        self.transcript.write(text)
        self.transcript.flush()
        return self.console.write(text)

    def flush(self) -> None:
        self.transcript.flush()
        self.console.flush()

    def __getattr__(self, name: str) -> Any:
        # Stream-specific attributes (e.g. encoding/isatty) belong to the console.
        # Their types vary, so this compatibility fallback deliberately uses Any.
        return getattr(self.console, name)


@contextmanager
def capture_mapper_output() -> Iterator[TextIO]:
    """Tee stdout/stderr to a temporary file, deleted even if mapping fails.

    This captures Python writes (including progress bars), not native code or
    child-process file descriptors, or logging handlers bound to older streams.
    """
    
    # Create a temporary file that will be deleted per video that is processed
    # that contains its output, but also let that output continue going to the cell 
    with tempfile.TemporaryFile(mode="w+", encoding="utf-8") as transcript:

        # Redirect the standard out and the stadard error to these Tee objects 
        # That is so that any time that a write happens, it happens both to 
        # the console and the text file at the same time
        with redirect_stdout(_Tee(sys.stdout, transcript)), redirect_stderr(
            _Tee(sys.stderr, transcript)
        ):
            yield transcript


# Match the mapper's message, capturing the two values as text.
# \s* means optional whitespace; (\S+) captures a non-whitespace value.
# Parentheses escaped as \( and \) match literal parentheses in the message.
_ALIGNMENT: re.Pattern[str] = re.compile(
    r"Estimated offset of alternative egocentric video with respect to Neon "
    r"scene video:\s*(\S+)\s*seconds\s*"
    r"\(Pearson correlation:\s*([^\s)]+)\s*\)",
    re.IGNORECASE,  # Accept uppercase or lowercase letters in the message.
)

# Match ANSI terminal controls, such as color changes and cursor movements.
# \x1b\[ starts the control; the next groups match its parameters and final command.
_ANSI: re.Pattern[str] = re.compile(r"\x1b\[[0-?]*[ -/]*[@-~]")


def parse_mapper_output(transcript: TextIO) -> ParsedMapperMetadata:
    """Read captured output and extract the reported time offsets/correlations."""
    # Finish writing any text still waiting in the file's buffer.
    transcript.flush()
    # Rewind the file so reading starts at its first character.
    transcript.seek(0)
    # Collect each alignment estimate in the order it appeared.
    estimates: list[TemporalAlignment] = []

    # Read one line at a time instead of loading the whole file into memory.
    for line in transcript:
        # Remove invisible terminal formatting before looking for the message.
        clean_line: str = _ANSI.sub("", line)
        # Find every matching message on this line (usually there is just one).
        for match in _ALIGNMENT.finditer(clean_line):
            # float accepts negatives, scientific notation, nan, and inf.
            # Invalid numeric text raises ValueError rather than being silently skipped.
            # Group 1 is the text between "scene video:" and "seconds".
            offset: float = float(match.group(1))
            # Group 2 is the text after "Pearson correlation:".
            correlation: float = float(match.group(2))
            # Store the two numbers together as one estimate.
            estimates.append({
                "offset_seconds": offset,
                "pearson_correlation": correlation,
            })

    # Every successful mapper run must report an alignment estimate.
    if not estimates:
        raise ValueError("No estimated offset/Pearson correlation found in mapper output.")

    return {
        "schema_version": 1,
        "temporal_alignment": estimates[-1],  # Last reported estimate.
        "temporal_alignment_estimates": estimates,
        "metric_description": (
            "Pearson correlation of video motion signals at the estimated time "
            "offset; measures temporal alignment, not spatial gaze accuracy."
        ),
    }
