"""End-to-end integrity tests on explicitly synthetic saved predictions only."""
from pathlib import Path
import ast
import copy
import csv
import json
import os
import shutil
import sys
import tempfile
import unittest
from concurrent.futures import ThreadPoolExecutor
from unittest.mock import patch

CODE=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(CODE));sys.path.insert(0,str(CODE/'tests'))
from validation_integrity import (audit_validation, checkpoint_issues,
    ensure_validation_plan, record_task_event, ValidationIntegrityError)
from validation_execution import execute_validation_tasks, ValidationTasksFailed
from analysis_input_guard import write_completion_marker, run_guarded_training
from validation_fixture import create_complete_validation_fixture, load_checkpoint, save_checkpoint


class ValidationIntegrityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.base_temp=tempfile.TemporaryDirectory(prefix='complete_synthetic_base_',dir=CODE/'tests')
        cls.base=Path(cls.base_temp.name)/'revision_R1_20261014'
        cls.fixture=create_complete_validation_fixture(cls.base)
        cls.plan=cls.fixture['plan']
        # Validate the positive control once before exercising damaged copies.
        report=audit_validation(cls.base)
        if report['status']!='COMPLETE':raise RuntimeError('Invalid synthetic positive control')

    @classmethod
    def tearDownClass(cls):
        cls.base_temp.cleanup()

    def setUp(self):
        self.temp=tempfile.TemporaryDirectory(prefix='validation_case_',dir=CODE/'tests')
        self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name)
        self.rev=self.root/'revision_R1_20261014'
        shutil.copytree(self.base,self.rev)
        self.audit=copy.deepcopy(self.fixture['input_audit'])
        for key in ['analysis_matrix','feature_table']:
            self.audit[key]=str(self.rev/Path(self.audit[key]).relative_to(self.base))
        self.cp=self.rev/'04_validation/checkpoints/bootstrap_0499.json.gz'

    def mutate(self,callback):
        value=load_checkpoint(self.cp);callback(value);save_checkpoint(self.cp,value)

    def assert_blocked(self,code=None):
        with self.assertRaises(ValidationIntegrityError) as caught:
            write_completion_marker(self.root,self.audit)
        report=json.loads((self.rev/'04_validation/VALIDATION_INTEGRITY.json').read_text(encoding='utf-8'))
        self.assertIn(report['status'],['PARTIAL','FAILED'])
        self.assertFalse((self.root/'ANALYSIS_COMPLETE_AWAITING_DOCUMENT_SYNC.json').exists())
        if code:self.assertIn(code,[x['code'] for x in report['issues']])
        return report

    def test_complete_510_checkpoints_and_actual_ledger_allow_completion(self):
        marker=write_completion_marker(self.root,self.audit)
        self.assertEqual(marker['status'],'COMPLETE')
        self.assertEqual(marker['final_patient_count'],25)
        self.assertEqual(marker['validation_counts']['point']['checkpoint_success_counts'],10)
        self.assertEqual(marker['validation_counts']['bootstrap']['checkpoint_success_counts'],500)
        self.assertEqual(marker['validation_counts']['bootstrap']['attempted'],500)
        self.assertEqual(marker['validation_counts']['bootstrap']['primary_paired_folds'],2500)

    def test_missing_checkpoint_blocks_completion(self):
        self.cp.unlink();self.assert_blocked('checkpoints_missing_or_invalid')

    def test_duplicate_internal_checkpoint_identity_blocks_completion(self):
        shutil.copy2(self.cp,self.cp.with_name('duplicate_bootstrap.json.gz'))
        self.assert_blocked('duplicate_checkpoint_identity')

    def test_wrong_checkpoint_filename_blocks_completion(self):
        self.cp.rename(self.cp.with_name('bootstrap_0500.json.gz'))
        self.assert_blocked('checkpoint_filename')

    def test_wrong_internal_run_number_blocks_completion(self):
        self.mutate(lambda cp:cp.update(run_id=500));self.assert_blocked('checkpoint_identity')

    def test_wrong_original_seed_blocks_completion(self):
        self.mutate(lambda cp:cp.update(seed=12345));self.assert_blocked('checkpoint_seed')

    def test_missing_applicable_model_blocks_completion(self):
        self.mutate(lambda cp:cp['folds'][0]['models'].pop('clinical'))
        self.assert_blocked('model_missing_or_failed')

    def test_missing_outer_fold_blocks_completion(self):
        self.mutate(lambda cp:cp['folds'].pop());self.assert_blocked('fold_identity')

    def test_duplicate_outer_fold_number_blocks_completion(self):
        self.mutate(lambda cp:cp['folds'][4].update(fold=0));self.assert_blocked('fold_identity')

    def test_missing_primary_paired_model_blocks_completion(self):
        self.mutate(lambda cp:cp['folds'][2]['models'].pop('clinical_diameter_mean_hu_radiomics'))
        self.assert_blocked('primary_pair_missing')

    def test_missing_formal_paired_results_blocks_completion_even_with_510_valid_runs(self):
        (self.rev/'04_validation/paired_incremental_performance.csv').unlink()
        self.assert_blocked('paired_summary_missing')

    def test_wrong_formal_paired_delta_blocks_completion(self):
        path=self.rev/'04_validation/paired_incremental_performance.csv'
        with path.open(encoding='utf-8-sig',newline='') as f:rows=list(csv.DictReader(f))
        rows[0]['paired_delta_C_index']='0.123456'
        with path.open('w',encoding='utf-8-sig',newline='') as f:
            writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
        self.assert_blocked('paired_summary_value')

    def test_unpaired_primary_validation_outcomes_block_completion(self):
        def change(cp):cp['folds'][1]['models']['clinical_diameter_mean_hu_radiomics']['validation_times'][0]+=1
        self.mutate(change);self.assert_blocked('primary_pair_missing')

    def test_wrong_matrix_fingerprint_blocks_completion(self):
        self.mutate(lambda cp:cp.update(analysis_matrix_sha256='0'*64))
        self.assert_blocked('checkpoint_input_fingerprint')

    def test_wrong_canonical_configuration_fingerprint_blocks_completion(self):
        self.mutate(lambda cp:cp.update(configuration_sha256='0'*64))
        self.assert_blocked('checkpoint_input_fingerprint')

    def test_wrong_configuration_file_fingerprint_blocks_completion(self):
        self.mutate(lambda cp:cp.update(configuration_file_sha256='0'*64))
        self.assert_blocked('checkpoint_input_fingerprint')

    def test_actual_configuration_file_drift_blocks_completion(self):
        p=self.rev/'00_protocol_and_provenance/revision_config.json'
        p.write_text(p.read_text(encoding='utf-8')+'\n',encoding='utf-8')
        self.assert_blocked('validation_plan_error')

    def test_missing_attempt_ledger_cannot_be_reconstructed_from_500_files(self):
        (self.rev/'04_validation/VALIDATION_TASK_EVENTS.jsonl').unlink()
        report=self.assert_blocked('task_ledger_missing')
        self.assertIsNone(report['counts']['bootstrap']['attempted'])
        self.assertEqual(report['counts']['bootstrap']['checkpoint_success_counts'],500)

    def test_missing_inner_tuning_candidate_blocks_completion(self):
        self.mutate(lambda cp:cp['folds'][0]['radiomics_development']['tuning'].pop())
        self.assert_blocked('development_tuning_grid')

    def test_model_specification_transform_drift_blocks_completion(self):
        self.mutate(lambda cp:cp['folds'][0]['models']['clinical']['transform']['columns'].pop())
        self.assert_blocked('transform_columns')

    def test_train_validation_patient_leakage_blocks_completion(self):
        def change(cp):cp['folds'][0]['train_ids'][0]=cp['folds'][0]['validation_ids'][0]
        self.mutate(change);self.assert_blocked('fold_patient_leakage')

    def test_task_failure_even_with_existing_valid_checkpoint_blocks_completion(self):
        record_task_event(self.rev,'bootstrap',499,'attempted',attempt_id='retry_failure')
        record_task_event(self.rev,'bootstrap',499,'failed',attempt_id='retry_failure',exception={'type':'RuntimeError','message':'Synthetic failure'})
        report=self.assert_blocked('task_success_missing')
        self.assertEqual(report['counts']['bootstrap']['attempted'],501)
        self.assertEqual(report['counts']['bootstrap']['failed'],1)

    def test_successful_retry_retains_old_failure_and_original_seed(self):
        for event,attempt,exception in [('attempted','prior_failure',None),('failed','prior_failure',{'message':'Synthetic prior failure'}),('attempted','same_seed_retry',None),('reused','same_seed_retry',None)]:
            record_task_event(self.rev,'bootstrap',499,event,attempt_id=attempt,exception=exception)
        report=audit_validation(self.rev)
        self.assertEqual(report['status'],'COMPLETE')
        self.assertEqual(report['counts']['bootstrap']['attempted'],502)
        self.assertEqual(report['counts']['bootstrap']['failed'],1)
        self.assertEqual(report['counts']['bootstrap']['reused'],1)
        self.assertEqual(report['prior_failures'][0]['seed'],20261015+499)

    def test_stale_completion_is_preserved_but_invalidated_after_checkpoint_loss(self):
        write_completion_marker(self.root,self.audit)
        self.cp.unlink()
        self.assert_blocked('checkpoints_missing_or_invalid')
        archived=list((self.rev/'04_validation/invalidated_completion_markers').glob('*.json'))
        self.assertEqual(len(archived),1)
        self.assertEqual(json.loads(archived[0].read_text(encoding='utf-8'))['status'],'COMPLETE')

    def test_summary_and_publication_functions_gate_before_any_output(self):
        self.cp.unlink()
        for filename in ['validate_pipeline.py','publication_data.py','publication_figures.py']:
            with self.subTest(filename=filename):
                tree=ast.parse((CODE/filename).read_text(encoding='utf-8'))
                function=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name==('summarize' if filename=='validate_pipeline.py' else 'main'))
                namespace={'audit_validation':audit_validation,'ROOT':self.rev}
                module=ast.Module(body=[function],type_ignores=[])
                exec(compile(ast.fix_missing_locations(module),filename,'exec'),namespace)
                with self.assertRaises(ValidationIntegrityError):namespace[function.name]()
        self.assertFalse((self.rev/'04_validation/full_pipeline_performance.csv').exists())
        self.assertFalse((self.rev/'05_figures_and_tables').exists())

    def test_guarded_pipeline_stops_before_summary_on_partial_validation(self):
        self.cp.unlink();calls=[]
        with patch('analysis_input_guard.verify_analysis_inputs',return_value=self.audit):
            with self.assertRaises(ValidationIntegrityError):
                run_guarded_training(self.root,lambda *args:calls.append(args))
        self.assertEqual(calls,[('diagnostics_and_tests.py',),('validate_pipeline.py','--mode','point'),('validate_pipeline.py','--mode','bootstrap')])

    def test_rsf_applicability_uses_actual_base_matrix_thresholds(self):
        for n,events,expected in [(100,70,21),(100,69,19),(99,70,19)]:
            with self.subTest(n=n,events=events),tempfile.TemporaryDirectory(dir=CODE/'tests') as folder:
                fixture=create_complete_validation_fixture(Path(folder)/'revision_R1_20261014',n_patients=n,n_events=events,modes=('point',),run_ids={'point':[0]})
                plan=fixture['plan'];self.assertEqual(len(plan['applicable_models']),expected)
                cp=load_checkpoint(fixture['checkpoint_paths'][0])
                self.assertEqual(checkpoint_issues(cp,plan,fixture['checkpoint_paths'][0]),[])
                if expected==21:
                    cp['folds'][0]['models'].pop('clinical_size_intensity_RSF')
                    self.assertIn('model_missing_or_failed',[i['code'] for i in checkpoint_issues(cp,plan)])

    def test_499_of_500_real_executor_exceptions_stop_and_count_all_attempts(self):
        # Preserve point10 but remove all prior bootstrap files/events.
        for path in (self.rev/'04_validation/checkpoints').glob('bootstrap_*.json.gz'):
            if path.name!='bootstrap_0000.json.gz':path.unlink()
        ledger=self.rev/'04_validation/VALIDATION_TASK_EVENTS.jsonl'
        events=[json.loads(line) for line in ledger.read_text(encoding='utf-8').splitlines()]
        ledger.write_text(''.join(json.dumps(e)+'\n' for e in events if e['mode']=='point'),encoding='utf-8')
        def task(run_id,mode):
            if run_id==0:return {'path':str(self.rev/'04_validation/checkpoints/bootstrap_0000.json.gz'),'reused':False}
            raise RuntimeError('Synthetic task exception '+str(run_id))
        with self.assertRaises(ValidationTasksFailed):
            execute_validation_tasks(self.rev,'bootstrap',range(500),task,workers=4,executor_factory=ThreadPoolExecutor)
        report=self.assert_blocked('checkpoints_missing_or_invalid')
        counts=report['counts']['bootstrap']
        self.assertEqual((counts['planned'],counts['attempted'],counts['succeeded'],counts['failed']),(500,500,1,499))
        failures=[json.loads(line) for line in (self.rev/'04_validation/VALIDATION_TASK_FAILURES.jsonl').read_text(encoding='utf-8').splitlines()]
        self.assertEqual(len(failures),499)
        self.assertEqual({x['original_seed'] for x in failures},{20261015+i for i in range(1,500)})
        self.assertTrue(all(x['exception'] and len(x['input_fingerprints']['analysis_matrix_sha256'])==64 for x in failures))


if __name__=='__main__':unittest.main()
