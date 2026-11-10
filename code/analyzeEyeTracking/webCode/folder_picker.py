"""Select a local folder, then exit so no Tk window stays with the web server."""

import json
import sys
from pathlib import Path
import tkinter as tk
from tkinter import filedialog


def main() -> None:
    """Return the chosen path as JSON; an empty string means cancellation."""
    root = tk.Tk()
    root.withdraw()
    try:
        starting_directory = sys.argv[1] if len(sys.argv) > 1 else ""
        selected_directory = filedialog.askdirectory(
            parent=root,
            title="Select your FLIC_data folder",
            initialdir=starting_directory or str(Path.home()),
            mustexist=True,
        )
        print(json.dumps(selected_directory))
    finally:
        root.destroy()


if __name__ == "__main__":
    main()
