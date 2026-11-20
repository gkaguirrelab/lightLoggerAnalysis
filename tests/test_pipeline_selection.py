"""Check subject/activity selection without importing camera or MATLAB libraries."""
import ast
from collections import defaultdict
import copy
import os
import pathlib
import re
import tempfile
import types
import unittest
import warnings
from typing import Iterable, Literal
from unittest.mock import Mock, patch

import numpy as np

SOURCE = pathlib.Path(__file__).resolve().parents[1] / 'code/preprocessRecordingData/preprocessing_pipeline.py'


def load_functions(*names):
    tree = ast.parse(SOURCE.read_text())
    names = set(names) | {'_get_desired_subjects_and_activities',
                          '_extract_num_from_id'}
    tree.body = [node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name in names]
    namespace = dict(globals(), natsorted=sorted,
                     get_subject_ids=Mock(return_value={21: 30, 22: None}),
                     get_activity_names=Mock(return_value={'walk', 'work'}))
    exec(compile(tree, str(SOURCE), 'exec'), namespace)
    return namespace


class SelectionTests(unittest.TestCase):
    def setUp(self):
        self.api = load_functions('verify_data_integrity')
        self.select = self.api['_get_desired_subjects_and_activities']

    def test_defaults_and_independent_sets(self):
        result = self.select(None, None)
        self.assertEqual(result, {21: {'walk', 'work'}, 22: {'walk', 'work'}})
        result[21].clear()
        self.assertEqual(result[22], {'walk', 'work'})
        self.assertEqual(self.select(None, None)[21], {'walk', 'work'})

    def test_only_requested_pairs_and_skip_wins_without_mutation(self):
        requested = {21: {'walk', 'work'}, 22: {'walk'}, 23: set()}
        skipped = {21: {'walk'}, 22: {'work'}, 99: {'walk'}}
        original = copy.deepcopy((requested, skipped))
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always')
            result = self.select(requested, skipped)
        self.assertEqual(result, {21: {'work'}, 22: {'walk'}})
        self.assertEqual(len(caught), 1)
        self.assertIn("(21, 'walk')", str(caught[0].message))
        self.assertIn('both', str(caught[0].message))
        self.assertEqual((requested, skipped), original)
        self.api['get_subject_ids'].assert_not_called()
        self.api['get_activity_names'].assert_not_called()

    def test_empty_means_nothing_and_all_skipped_is_empty(self):
        self.assertEqual(self.select({}, None), {})
        with self.assertWarns(UserWarning):
            self.assertEqual(self.select({21: {'walk'}}, {21: {'walk'}}), {})
        self.api['get_subject_ids'].assert_not_called()

    def test_integrity_checks_only_selected_activities(self):
        with tempfile.TemporaryDirectory() as directory:
            pathlib.Path(directory, 'FLIC_21').mkdir()
            pathlib.Path(directory, 'FLIC_22').mkdir()
            with self.assertWarns(UserWarning):
                result = self.api['verify_data_integrity'](
                    src_dir=directory, output_path=None,
                    subjects_and_activities_to_process={21: {'walk', 'work'}, 22: {'walk'}},
                    subjects_and_activities_to_skip={21: {'walk'}},
                )
            self.assertEqual(set(result[21]['activities']), {'work'})
            self.assertEqual(set(result[22]['activities']), {'walk'})
            self.assertEqual(self.api['verify_data_integrity'](
                src_dir=directory, output_path=None, subjects_and_activities_to_process={}), {})

    def test_tag_task_generation_uses_selected_pairs(self):
        api = load_functions('generate_tag_task_start_ends')
        api['video_io'] = Mock()
        api['video_io'].inspect_video_FPS.return_value = 10
        api['video_io'].inspect_video_frame_count.return_value = 100
        api['scipy'] = types.SimpleNamespace(io=Mock())

        with tempfile.TemporaryDirectory() as directory:
            for subject in ('FLIC_0021', 'FLIC_22', 'FLIC_99'):
                for activity in ('walk', 'work'):
                    folder = pathlib.Path(directory, subject, activity, 'GKA')
                    folder.mkdir(parents=True)
                    (folder / 'W.avi').touch()

            with self.assertWarns(UserWarning):
                api['generate_tag_task_start_ends'](
                    src_dir=directory, dst_dir=directory, task_length_seconds=4,
                    subjects_and_activities_to_process={21: {'walk', 'work'}, 22: {'walk'}},
                    subjects_and_activities_to_skip={21: {'walk'}},
                )

            calls = api['scipy'].io.savemat.call_args_list
            self.assertEqual({str(pathlib.Path(call.args[0]).relative_to(directory)) for call in calls}, {
                'FLIC_0021/work/tag_task_start_end.mat',
                'FLIC_22/walk/tag_task_start_end.mat',
            })
            for call in calls:
                result = call.args[1]['tag_task_start_end']
                np.testing.assert_array_equal(result['tag'], [1, 60])
                np.testing.assert_array_equal(result['task'], [61, 100])

    def test_combining_preserves_different_activity_sets(self):
        api = load_functions('combine_spds')
        engine = Mock()
        matlab = types.ModuleType('matlab')
        matlab.engine = types.ModuleType('matlab.engine')
        matlab.engine.start_matlab = Mock(return_value=engine)
        api['scipy'] = types.SimpleNamespace(io=Mock())
        with tempfile.TemporaryDirectory() as directory:
            for subject, activity in ((21, 'work'), (22, 'walk')):
                folder = pathlib.Path(directory, 'GRAY', f'FLIC_{subject:04d}', activity)
                folder.mkdir(parents=True)
                (folder / f'FLIC_{subject:04d}_{activity}_virtuallyFoveated_SPDResults.mat').touch()
            pathlib.Path(directory, 'GRAY', 'FLIC_99').mkdir()
            with patch.dict('sys.modules', {'matlab': matlab, 'matlab.engine': matlab.engine}), patch.object(os, 'remove'):
                api['combine_spds'](
                    src_dir=directory, dst_dir=directory,
                    subjects_and_activities_to_process={21: {'work'}, 22: {'walk'}},
                    projection_types={'virtuallyFoveated'},
                )
            result = api['scipy'].io.savemat.call_args.args[1]['spds']
            self.assertEqual(set(result), {'FLIC_0021', 'FLIC_0022'})
            self.assertEqual(set(result['FLIC_0021']), {'work'})
            self.assertEqual(set(result['FLIC_0022']), {'walk'})

    def test_mean_uses_only_subjects_selected_for_each_activity(self):
        api = load_functions('generate_mean_spds')
        engine = Mock()
        matlab = types.ModuleType('matlab')
        matlab.engine = types.ModuleType('matlab.engine')
        matlab.engine.start_matlab = Mock(return_value=engine)
        api['spd_util'] = Mock()
        api['spd_util'].load_spds.side_effect = [
            {'GRAY': {'FLIC_21': {'work': {'projection': 'work21'}}}},
            {'GRAY': {'FLIC_22': {'walk': {'projection': 'walk22'}}}},
        ]
        with tempfile.TemporaryDirectory() as directory, patch.dict(
            'sys.modules', {'matlab': matlab, 'matlab.engine': matlab.engine}
        ):
            api['generate_mean_spds'](
                src_dir=directory, dst_dir=directory, dimension='subject',
                subjects_and_activities_to_process={21: {'work'}, 22: {'walk'}},
            )
        self.assertEqual(engine.meanSPDs.call_count, 2)
        self.assertEqual(engine.meanSPDs.call_args_list[0].args[0], [{'projection': 'walk22'}])
        self.assertEqual(engine.meanSPDs.call_args_list[1].args[0], [{'projection': 'work21'}])


if __name__ == '__main__':
    unittest.main()
