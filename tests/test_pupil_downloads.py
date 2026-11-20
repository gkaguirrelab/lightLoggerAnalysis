"""Exercise download and migration without scientific imports or Cloud access."""
import ast
from collections import defaultdict
import io
import json
import os
import pathlib
import re
import shutil
import tempfile
import time
import types
import unittest
import warnings
import zipfile
from typing import Iterable
from unittest.mock import Mock, patch

SOURCE = pathlib.Path(__file__).resolve().parents[1] / 'code/preprocessRecordingData/preprocessing_pipeline.py'


def load_functions():
    names = {'download_pupil_cloud_recordings', '_get_desired_subjects_and_activities'}
    tree = ast.parse(SOURCE.read_text())
    tree.body = [n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name in names]
    request_mock = Mock()
    request_mock.ConnectionError = ConnectionError
    request_mock.Timeout = TimeoutError
    request_mock.exceptions.ChunkedEncodingError = EOFError
    namespace = dict(globals(), requests=request_mock, time=Mock(),
                     get_subject_ids=lambda: {2001: None}, get_activity_names=lambda: {"walk"},
                     natsorted=sorted)
    exec(compile(tree, str(SOURCE), 'exec'), namespace)
    notebook = json.loads(SOURCE.with_name('migrate_neon_raw_exports.ipynb').read_text())
    exec(''.join(notebook['cells'][2]['source']), namespace)
    return types.SimpleNamespace(**namespace)


def response(files):
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, 'w') as archive:
        for name, content in files.items():
            archive.writestr(name, content)
    result = Mock()
    result.__enter__ = Mock(return_value=result)
    result.__exit__ = Mock(return_value=False)
    result.iter_content.return_value = [buffer.getvalue()]
    return result


class DownloadTests(unittest.TestCase):
    def setUp(self):
        self.api = load_functions()
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = pathlib.Path(self.temp.name)

    def test_two_formats_and_independent_resume(self):
        listing = Mock()
        listing.json.return_value = {'result': [{'name': 'FLIC_2001_walk_1', 'id': 'abc'}]}
        export = response({'recording/world_timestamps.csv': 'timestamp', 'recording/scene.mp4': 'video'})
        native = response({'abc/info.json': '{}'})
        self.api.requests.get.side_effect = [listing, export, native]
        self.api.download_pupil_cloud_recordings('test', dst_dir=str(self.root))
        attempt = self.root / 'FLIC_2001/walk/Neon/1'
        self.assertTrue((attempt / 'timeseries_and_video_export/recording/scene.mp4').is_file())
        self.assertTrue((attempt / 'native_recording_data/abc/info.json').is_file())
        calls = self.api.requests.get.call_args_list
        self.assertTrue(calls[1].args[0].endswith('recordings:raw-data-export'))
        self.assertTrue(calls[2].args[0].endswith('recordings.zip'))
        self.assertEqual(calls[2].kwargs['params'], {'ids': ['abc']})
        shutil.rmtree(attempt / 'native_recording_data')
        self.api.requests.get.reset_mock()
        self.api.requests.get.side_effect = [listing, response({'abc/info.json': '{}'})]
        self.api.download_pupil_cloud_recordings('test', dst_dir=str(self.root))
        self.assertEqual(self.api.requests.get.call_count, 2)
        self.assertTrue(self.api.requests.get.call_args.args[0].endswith('recordings.zip'))

    def test_only_selected_subject_activity_pairs_are_downloaded(self):
        listing = Mock()
        listing.json.return_value = {'result': [
            {'name': f'FLIC_{subject}_{activity}_1', 'id': f'{subject}-{activity}'}
            for subject in (2001, 2002, 2003) for activity in ('walk', 'work')
        ]}
        self.api.requests.get.side_effect = [listing] + [response({'data': 'ok'}) for _ in range(4)]
        with self.assertWarns(UserWarning):
            self.api.download_pupil_cloud_recordings(
                'test', dst_dir=str(self.root),
                subjects_and_activities_to_process={2001: {'walk', 'work'}, 2002: {'walk'}},
                subjects_and_activities_to_skip={2001: {'walk'}},
            )
        self.assertFalse((self.root / 'FLIC_2001/walk').exists())
        self.assertTrue((self.root / 'FLIC_2001/work/Neon/1').is_dir())
        self.assertTrue((self.root / 'FLIC_2002/walk/Neon/1').is_dir())
        self.assertFalse((self.root / 'FLIC_2002/work').exists())
        self.assertFalse((self.root / 'FLIC_2003').exists())
        self.assertEqual(self.api.requests.get.call_count, 5)

    def test_migration_and_native_exclusion(self):
        old = self.root / '1/export'
        old.mkdir(parents=True)
        (old / 'world_timestamps.csv').write_text('original')
        native = self.root / '1/native_recording_data'
        native.mkdir()
        (native / 'world_timestamps.csv').write_text('different')
        self.api._migrate_neon_export(str(self.root / '1'))
        expected = self.root / '1/timeseries_and_video_export/export'
        self.assertEqual((expected / 'world_timestamps.csv').read_text(), 'original')
        self.api._migrate_neon_export(str(self.root / '1'))
        self.assertTrue((native / 'world_timestamps.csv').is_file())

    def test_failed_extraction_is_removed(self):
        listing = Mock()
        listing.json.return_value = {'result': [{'name': 'FLIC_2001_walk_1', 'id': 'abc'}]}
        invalid = response({})
        invalid.iter_content.return_value = [b'not a zip']
        self.api.requests.get.side_effect = [listing, invalid]
        with self.assertRaises(zipfile.BadZipFile):
            self.api.download_pupil_cloud_recordings('test', dst_dir=str(self.root))
        self.assertFalse((self.root / 'FLIC_2001/walk/Neon/1/timeseries_and_video_export').exists())

    def test_leftover_zip_does_not_skip_download(self):
        dst = self.root / 'FLIC_2001/walk/Neon/1/timeseries_and_video_export'
        dst.mkdir(parents=True)
        (dst / 'download.zip').write_bytes(b'partial')
        listing = Mock()
        listing.json.return_value = {'result': [{'name': 'FLIC_2001_walk_1', 'id': 'abc'}]}
        self.api.requests.get.side_effect = [listing, response({'gaze.csv': 'data'}), response({'info.json': '{}'})]
        self.api.download_pupil_cloud_recordings('test', dst_dir=str(self.root))
        self.assertEqual((dst / 'gaze.csv').read_text(), 'data')
        self.assertFalse((dst / 'download.zip').exists())

    def test_network_error_is_not_retried(self):
        listing = Mock()
        listing.json.return_value = {'result': [{'name': 'FLIC_2001_walk_1', 'id': 'abc'}]}
        self.api.requests.get.side_effect = [listing, ConnectionError('offline')]
        with self.assertRaises(ConnectionError):
            self.api.download_pupil_cloud_recordings('test', dst_dir=str(self.root))
        self.assertEqual(self.api.requests.get.call_count, 2)
        self.api.time.sleep.assert_not_called()


if __name__ == '__main__':
    unittest.main()
