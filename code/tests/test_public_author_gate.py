"""Public author-input gate tests; zero scientific model fitting."""
from pathlib import Path
import json,sys,tempfile,unittest
CODE=Path(__file__).resolve().parents[1];sys.path.insert(0,str(CODE))
from frozen_input_pipeline import prepare,rows,verify_rosters
from public_execution_gate import prepare_execution_gate
from run_new_cohort import require_frozen_model_inputs
PACKAGE=CODE.parent

class AuthorInputGateTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory(prefix='TACE author gate ');self.addCleanup(self.tmp.cleanup)
        self.root=Path(self.tmp.name)/'new run'
    def test_author_adoption_does_not_certify_physician_evidence(self):
        prepare(self.root)
        audit=prepare_execution_gate(self.root,False)
        marker=json.loads((self.root/'COHORT_READY_FOR_MODELING.json').read_text())
        self.assertEqual(audit['actual_patient_count'],120)
        self.assertEqual(audit['actual_active_unique_mask_count'],210)
        self.assertTrue(marker['input_decisions_complete'])
        self.assertFalse(marker['medical_decisions_complete'])
        self.assertFalse(marker['medical_evidence_requirements_satisfied'])
        self.assertFalse(marker['submission_ready'])
        self.assertEqual(marker['cohort_analysis_scope'],'author_specified_inputs')
        self.assertFalse((self.root/'MODEL_ANALYSIS_COMPLETE.json').exists())
    def test_changed_adoption_record_blocks_low_level_model_entry(self):
        prepare(self.root);prepare_execution_gate(self.root,False)
        path=self.root/'PUBLIC_INPUT_ADOPTION.csv'
        path.write_bytes(path.read_bytes()+b'\nchanged')
        with self.assertRaises(Exception):require_frozen_model_inputs(self.root)
        self.assertFalse((self.root/'MODEL_ANALYSIS_COMPLETE.json').exists())
    def test_eight_original_masks_no_alias_and_clinical_counts_preserved(self):
        images,reps,candidates,mapping=verify_rosters()
        selected=[r for r in images if r['patient_id'] in ['83','94']]
        self.assertEqual(len(selected),8)
        self.assertTrue(all(r['action']=='retain' and r['active_analysis_record']=='True' and not r['alias_of'] for r in selected))
        self.assertTrue(all(m['final_mask_kind']=='official_original' for m in mapping if m['patient_id'] in ['83','94']))
        cs={r['patient_id']:r for r in candidates}
        self.assertEqual((cs['83']['source_clinical_lesion_count'],cs['94']['source_clinical_lesion_count']),('4','4'))
        rs={r['patient_id']:r for r in reps}
        self.assertEqual((rs['83']['sample_id'],rs['94']['sample_id']),('83_1_2','94_1_1'))
        self.assertEqual(sum(m['final_mask_kind']=='packaged_final' for m in mapping),17)
    def test_author_gate_marker_cannot_claim_medical_completion(self):
        prepare(self.root);prepare_execution_gate(self.root,False)
        p=self.root/'COHORT_READY_FOR_MODELING.json';d=json.loads(p.read_text())
        d['medical_evidence_requirements_satisfied']=True;p.write_text(json.dumps(d))
        with self.assertRaises(RuntimeError):require_frozen_model_inputs(self.root)

if __name__=='__main__':unittest.main(verbosity=2)
