"""Release-input unit and synthetic-entry tests, never scientific model evidence."""
from pathlib import Path
import csv,json,shutil,subprocess,sys,tempfile,unittest
from unittest.mock import patch
CODE=Path(__file__).resolve().parents[1];sys.path.insert(0,str(CODE))
from legacy_audit_input import ensure_legacy_pipeline_audit
from release_identity import verify_release
from frozen_input_pipeline import verify_rosters,rows,write,compare_csv,verify_bindings,prepare,verify_outputs
PACKAGE=CODE.parent
class PublicInputTests(unittest.TestCase):
 def setUp(self):
  self.temp=tempfile.TemporaryDirectory(prefix='public input test ');self.addCleanup(self.temp.cleanup);self.root=Path(self.temp.name);self.copy=self.root/'package copy'
  shutil.copytree(PACKAGE,self.copy,ignore=shutil.ignore_patterns('__pycache__'))
 def test_reviewed_legacy_csv_valid(self):
  p=ensure_legacy_pipeline_audit(self.root/'run',self.copy/'code');self.assertEqual(p.read_bytes(),(self.copy/'code/reporting_inputs/01_data_audit/legacy_pipeline_audit.csv').read_bytes())
 def test_legacy_tamper_rejected(self):
  p=self.copy/'code/reporting_inputs/01_data_audit/legacy_pipeline_audit.csv';p.write_bytes(p.read_bytes()+b'changed')
  with self.assertRaises(RuntimeError):ensure_legacy_pipeline_audit(self.root/'run',self.copy/'code')
 def test_legacy_missing_rejected(self):
  (self.copy/'code/reporting_inputs/01_data_audit/legacy_pipeline_audit.csv').unlink()
  with self.assertRaises(FileNotFoundError):ensure_legacy_pipeline_audit(self.root/'run',self.copy/'code')
 def test_manifest_version_or_hash_rejected(self):
  p=self.copy/'RELEASE_INPUTS.json';d=json.loads(p.read_text());d['version']='wrong';p.write_text(json.dumps(d))
  with self.assertRaises(RuntimeError):verify_release(self.copy)
 def test_matrix_entry_from_other_cwd_and_spaces(self):
  cwd=self.root/'different working directory';cwd.mkdir();dest=self.root/'actual matrix input run'
  r=subprocess.run([sys.executable,'-B',str(self.copy/'code/run_revision.py'),'--from-matrix','--preflight-only','--project-root',str(dest)],cwd=cwd,capture_output=True,text=True)
  self.assertEqual(r.returncode,0,r.stderr);d=json.loads((dest/'INPUT_PREFLIGHT.json').read_text());self.assertEqual((d['patients'],d['active_ROIs']),(120,210));self.assertFalse((dest/'ANALYSIS_COMPLETE_AWAITING_DOCUMENT_SYNC.json').exists())
 def test_rosters_and_17_final_masks(self):
  im,rp,cp,mp=verify_rosters(self.copy);self.assertEqual((len(im),len(rp),len(cp)),(230,120,126));self.assertEqual(sum(x['final_mask_kind']=='packaged_final' for x in mp),17)
 def test_missing_final_mask_rejected(self):
  next((self.copy/'final_masks').glob('*.nrrd')).unlink()
  with self.assertRaises(ValueError):verify_rosters(self.copy)
 def test_corrupt_final_mask_rejected(self):
  next((self.copy/'final_masks').glob('*.nrrd')).write_bytes(b'changed')
  with self.assertRaises(ValueError):verify_rosters(self.copy)
 def test_wrong_CT_mapping_rejected(self):
  p=self.copy/'frozen_inputs/SOURCE_FILE_MAP.csv';a=rows(p);a[0]['ct_sha256']='0'*64;write(p,a)
  with self.assertRaises(ValueError):verify_rosters(self.copy)
 def test_changed_representative_rejected(self):
  p=self.copy/'frozen_inputs/FINAL_PATIENT_AND_REPRESENTATIVE_ROSTER.csv';a=rows(p);a[0]['sample_id']=a[1]['sample_id'];write(p,a)
  with self.assertRaises(ValueError):verify_rosters(self.copy)
 def test_changed_clinical_count_rejected(self):
  p=self.copy/'frozen_inputs/FINAL_PATIENT_AND_REPRESENTATIVE_ROSTER.csv';a=rows(p);a[0]['lesions_number']='999';write(p,a)
  with self.assertRaises(ValueError):verify_rosters(self.copy)
 def test_duplicate_or_extra_ROI_rejected(self):
  p=self.copy/'frozen_inputs/FINAL_IMAGE_AND_ALIAS_ROSTER.csv';a=rows(p);a.append(dict(a[0]));write(p,a)
  with self.assertRaises(ValueError):verify_rosters(self.copy)
 def test_missing_patient_rejected(self):
  p=self.copy/'frozen_inputs/FINAL_PATIENT_AND_REPRESENTATIVE_ROSTER.csv';a=rows(p);write(p,a[1:])
  with self.assertRaises(ValueError):verify_rosters(self.copy)
 def test_nonfinite_feature_rejected(self):
  ref=self.copy/'code/data/analysis_matrix.csv';a=rows(ref);key=next(k for k in a[0] if k.startswith('original_'));a[0][key]='nan';p=self.root/'bad.csv';write(p,a)
  with self.assertRaises(ValueError):compare_csv(ref,p,1e-8,1e-10)
 def test_changed_feature_value_rejected(self):
  ref=self.copy/'code/data/analysis_matrix.csv';a=rows(ref);key=next(k for k in a[0] if k.startswith('original_'));a[0][key]='999999999999999999';p=self.root/'bad.csv';write(p,a)
  with self.assertRaises(ValueError):compare_csv(ref,p,1e-8,1e-10)
 def test_input_change_after_preflight_rejected(self):
  p=self.root/'input.csv';p.write_text('original');d=verify_release();import hashlib
  binding={'release_manifest_sha256':d['manifest_sha256'],'files':{str(p):hashlib.sha256(p.read_bytes()).hexdigest()},'source_files':{}};p.write_text('changed')
  with self.assertRaises(ValueError):verify_bindings(binding)
 def test_recorded_extraction_failure_blocks_training_gate(self):
  root,rev,binding=prepare(self.root/'failure gate')
  write(rev/'02_extraction/extraction_failures.csv',[{'sample_id':'101_1_0','reason':'synthetic extraction failure'}])
  with self.assertRaises(Exception):verify_outputs(root,False)
  self.assertFalse((root/'ANALYSIS_COMPLETE_AWAITING_DOCUMENT_SYNC.json').exists())
 def test_decreased_feature_rows_block_training_gate(self):
  root,rev,binding=prepare(self.root/'missing feature gate');p=rev/'02_extraction/arterial_lesion_features.csv';a=rows(p);write(p,a[1:])
  with self.assertRaises(ValueError):verify_outputs(root,False)
 def test_synthetic_child_failure_returns_nonzero_without_complete(self):
  import frozen_run_revision as driver
  root=self.root/'synthetic driver failure'
  argv=['run_revision.py','--from-matrix','--train','--rscript',sys.executable,'--project-root',str(root)]
  with patch.object(sys,'argv',argv),patch.object(driver,'check_environment',return_value={'scope':'SYNTHETIC driver test'}),patch.object(driver.subprocess,'run',return_value=subprocess.CompletedProcess([],9)):
   status=driver.main()
  self.assertEqual(status,1);d=json.loads((root/'EXECUTION_FAILURE.json').read_text());self.assertEqual(d['status'],'FAILED');self.assertFalse((root/'ANALYSIS_COMPLETE_AWAITING_DOCUMENT_SYNC.json').exists())
if __name__=='__main__':unittest.main(verbosity=2)
