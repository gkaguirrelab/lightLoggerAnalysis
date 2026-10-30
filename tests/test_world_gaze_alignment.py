"""Exercise timestamp alignment without importing camera reconstruction dependencies."""
import ast
import os
from pathlib import Path
import tempfile
import types
import unittest

import numpy as np
import pandas as pd

SOURCE = Path(__file__).resolve().parents[1] / 'code/preprocessRecordingData/preprocessing_pipeline.py'
tree = ast.parse(SOURCE.read_text())
tree.body = [node for node in tree.body if isinstance(node, ast.FunctionDef)
             and node.name == '_align_gaze_to_raw_world_frames']


class WorldGazeAlignmentTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.files = []
        # A gap across chunk boundaries; the video contains blank frames at 20/30.
        for index, times in enumerate(([100, 110], [140, 150], [])):
            metadata = self.root / f'metadata{index}.npy'
            frames = self.root / f'frames{index}.npy'
            matrix = np.zeros((len(times), 4), dtype=np.int64)
            matrix[:, 0] = times
            np.save(metadata, matrix)
            np.save(frames, np.zeros((len(times), 1, 1)))
            self.files.append((str(metadata), str(frames)))
        namespace = dict(np=np, pd=pd, os=os, sensors_utility=types.SimpleNamespace(
            group_sensors_files=lambda path: {'W': self.files}))
        exec(compile(tree, str(SOURCE), 'exec'), namespace)
        self.align = namespace['_align_gaze_to_raw_world_frames']
        self.origin = 1_700_000_000_000_000_000
        pd.DataFrame({'timestamp [ns]': self.origin + np.arange(0, 51, 10)}).to_csv(
            self.root / 'alternative_camera_timestamps.csv', index=False)

    def write_gaze(self, offsets, x):
        pd.DataFrame({'timestamp [ns]': self.origin + np.array(offsets),
                      'gaze x [px]': x, 'gaze y [px]': x}).to_csv(
                          self.root / 'cleaned_eye_data.csv', index=False)

    def test_gap_nearest_tie_and_nan(self):
        self.write_gaze([0, 9, 11, 20, 30, 39, 41, 50],
                        [1, 2, 3, 4, 5, np.nan, 7, 8])
        result = self.align(str(self.root), str(self.root))
        np.testing.assert_allclose(result[:, 0], [1, 2, np.nan, 8], equal_nan=True)
        self.assertEqual(result.shape, (4, 2))

    def test_outside_coverage(self):
        self.write_gaze([9, 11, 39, 41], [1, 2, 3, 4])
        result = self.align(str(self.root), str(self.root))
        np.testing.assert_allclose(result[:, 0], [np.nan, 1, 3, np.nan], equal_nan=True)

    def test_unsorted_gaze_rejected(self):
        self.write_gaze([0, 20, 10], [1, 2, 3])
        with self.assertRaisesRegex(ValueError, 'Cleaned gaze timestamps'):
            self.align(str(self.root), str(self.root))

    def test_frame_metadata_mismatch_rejected(self):
        self.write_gaze([0, 50], [1, 2])
        np.save(self.files[0][1], np.zeros((1, 1, 1)))
        with self.assertRaisesRegex(ValueError, 'Frame/metadata length mismatch'):
            self.align(str(self.root), str(self.root))

    def test_relative_cleaned_times_recover_trimmed_origin(self):
        # Use realistic intervals, larger than float precision at epoch scale.
        self.origin = 1_700_000_000_000_000_000
        for (metadata, _), times in zip(self.files, ([100, 110], [140, 150], [])):
            matrix = np.zeros((len(times), 4), dtype=np.int64)
            matrix[:, 0] = np.array(times, dtype=np.int64) * 1_000_000
            np.save(metadata, matrix)
        pd.DataFrame({'timestamp [ns]': self.origin + np.arange(0, 51, 10) * 1_000_000}).to_csv(
            self.root / 'alternative_camera_timestamps.csv', index=False)
        offsets = np.array([0, 5, 10, 20, 30, 40, 45, 50]) * 1_000_000
        pd.DataFrame({'timestamp [ns]': self.origin + offsets}).to_csv(
            self.root / 'alternative_camera_gaze.csv', index=False)
        pd.DataFrame({'timestamp [ns]': self.origin + np.array([4, 46]) * 1_000_000}).to_csv(
            self.root / '3d_eye_states.csv', index=False)
        retained = (self.origin + offsets[1:-1]).astype(float)
        cleaned = pd.DataFrame({'time [s]': (retained - retained[0]) / 1e9,
                                'gaze x [px]': [1, 2, 3, 4, 5, 6],
                                'gaze y [px]': [1, 2, 3, 4, 5, 6]})
        path = self.root / 'cleaned_eye_data.csv'
        cleaned.to_csv(path, index=False)
        result = self.align(str(self.root), str(self.root), str(self.root))
        np.testing.assert_allclose(result[:, 0], [np.nan, 2, 5, np.nan], equal_nan=True)
        # A same-length file from another cleaning timeline must not silently align.
        cleaned.loc[2, 'time [s]'] += .001
        cleaned.to_csv(path, index=False)
        with self.assertRaisesRegex(ValueError, 'does not match'):
            self.align(str(self.root), str(self.root), str(self.root))


if __name__ == '__main__':
    unittest.main()
