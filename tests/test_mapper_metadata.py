"""Lightweight tests; no GPU, video data, or Pupil Labs install required."""
import ast
from contextlib import redirect_stdout, redirect_stderr
from datetime import datetime, timezone
import io
import os
from pathlib import Path
import pickle
import sys
import tempfile
import time
import types
import unittest
from unittest.mock import patch
from typing import Any, Literal, TextIO, cast

SOURCE: Path = Path(__file__).resolve().parents[1] / 'code/applyVirtualFoveation/pythonCode'
sys.path.insert(0, str(SOURCE))
from mapper_metadata import (
    MapperMetadata, ParsedMapperMetadata, capture_mapper_output, parse_mapper_output,
)

MESSAGE: str = ('Estimated offset of alternative egocentric video with respect to Neon '
           'scene video: -4.952 seconds (Pearson correlation: 0.9717797291080065)')


class MapperMetadataTests(unittest.TestCase):
    def test_tee_and_parse(self) -> None:
        # Stand in for the notebook streams and verify both destinations.
        out: io.StringIO = io.StringIO()
        err: io.StringIO = io.StringIO()
        with redirect_stdout(out), redirect_stderr(err):
            with capture_mapper_output() as transcript:
                print('\x1b[32m' + MESSAGE + '\x1b[0m')
                print('progress\r100%', file=sys.stderr)
                result: ParsedMapperMetadata = parse_mapper_output(transcript)
        self.assertTrue(transcript.closed)
        self.assertIn(MESSAGE, out.getvalue())
        self.assertEqual(err.getvalue(), 'progress\r100%\n')
        self.assertEqual(result['temporal_alignment']['offset_seconds'], -4.952)
        self.assertEqual(result['temporal_alignment']['pearson_correlation'], 0.9717797291080065)

    def test_missing_repeated_and_scientific_notation(self) -> None:
        with self.assertRaisesRegex(ValueError, 'No estimated offset/Pearson correlation'):
            parse_mapper_output(io.StringIO('100% complete'))
        result: ParsedMapperMetadata = parse_mapper_output(io.StringIO(MESSAGE + '\n' + MESSAGE.replace('-4.952', '+1.2e-3')))
        self.assertEqual(len(result['temporal_alignment_estimates']), 2)
        self.assertEqual(result['temporal_alignment']['offset_seconds'], .0012)

    def test_failure_restores_streams_and_closes_file(self) -> None:
        # Exceptions must unwind both redirections and the temporary file.
        original: tuple[TextIO, TextIO] = (sys.stdout, sys.stderr)
        with self.assertRaisesRegex(RuntimeError, 'mapper failed'):
            with capture_mapper_output() as transcript:
                raise RuntimeError('mapper failed')
        self.assertTrue(transcript.closed)
        self.assertEqual((sys.stdout, sys.stderr), original)

    def test_wrapper_writes_pickle(self) -> None:
        # Load only the wrapper and Args to avoid unrelated scientific imports.
        tree: ast.Module = ast.parse((SOURCE / 'virtual_foveation.py').read_text())
        tree.body = [node for node in tree.body if isinstance(node, (ast.FunctionDef, ast.ClassDef))
                     and node.name in ('Args', 'run_egocentric_video_mapper')]
        # Typed stubs exercise the wrapper without running video conversion.
        def convert_video(src: str, dst: str) -> None:
            Path(dst).touch()

        def run_mapper(args: object) -> None:
            print(MESSAGE)

        def preview(*args: object) -> None:
            pass

        # exec requires heterogeneous globals, so Any is intentional here.
        namespace: dict[str, Any] = dict(os=os, pickle=pickle, time=time, datetime=datetime,
                         timezone=timezone, Literal=Literal,
                         MapperMetadata=MapperMetadata, ParsedMapperMetadata=ParsedMapperMetadata,
                         capture_mapper_output=capture_mapper_output,
                         parse_mapper_output=parse_mapper_output,
                         convert_avi_to_mp4=convert_video)
        exec(compile(tree, str(SOURCE / 'virtual_foveation.py'), 'exec'), namespace)
        modules: dict[str, types.ModuleType] = {name: types.ModuleType(name) for name in (
            'pupil_labs', 'pupil_labs.egocentric_video_mapper',
            'pupil_labs.egocentric_video_mapper.__main__',
            'pupil_labs.egocentric_video_mapper.utils')}
        setattr(modules['pupil_labs.egocentric_video_mapper.__main__'], 'main', run_mapper)
        setattr(modules['pupil_labs.egocentric_video_mapper.utils'], 'show_videos_preview', preview)
        # Inspect the persisted artifact, not just the parser's return value.
        with tempfile.TemporaryDirectory() as output, patch.dict(sys.modules, modules), redirect_stdout(io.StringIO()):
            namespace['run_egocentric_video_mapper']('neon', 'world.avi', output)
            with open(Path(output) / 'egocentric_mapper_metadata.pkl', 'rb') as stream:
                metadata: MapperMetadata = cast(MapperMetadata, pickle.load(stream))
            self.assertEqual(metadata['status'], 'completed')
            self.assertEqual(metadata['inputs']['alternative_vid_path'], 'world.avi')
            self.assertEqual(metadata['temporal_alignment']['offset_seconds'], -4.952)
            self.assertEqual(os.listdir(output), ['egocentric_mapper_metadata.pkl'])


if __name__ == '__main__':
    unittest.main()
