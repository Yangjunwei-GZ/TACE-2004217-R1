"""Attribution-only compatibility and fail-closed evidence-scope checks."""
from pathlib import Path
import copy,sys,unittest
CODE=Path(__file__).resolve().parents[1];sys.path.insert(0,str(CODE))
from decision_provenance import REQUIRED_SAMPLES,LEGACY_KIND,PHYSICIAN_KIND,provenance_from_rows

def corrected_rows():
    return [dict(sample_id=sid,action='retain',alias_of='',decision_kind=PHYSICIAN_KIND,
                 decision_maker_role='physician',communication_source_role='author',
                 physician_retention_decision_relayed='True',physician_decision_date='not_provided',
                 physician_signature_provided='False',independent_review_evidence_provided='False',
                 decision_blinding_provided='False',independent_lesion_identity_confirmed='False',
                 decision_date='2026-10-26',
                 decision_date_scope='author_communication_record_not_physician_decision_date')
            for sid in REQUIRED_SAMPLES]

class DecisionProvenanceTests(unittest.TestCase):
    def test_physician_decision_and_author_relay_are_distinct_from_review_evidence(self):
        before=corrected_rows();original=copy.deepcopy(before);result=provenance_from_rows(before)
        self.assertTrue(result['physician_retention_decision_relayed'])
        self.assertEqual(result['physician_decision_date'],'not_provided')
        self.assertEqual(result['cohort_analysis_scope'],'physician_decision_relayed_by_author_inputs')
        for key in ['physician_signature_provided','independent_review_evidence_provided','decision_blinding_provided','independent_lesion_identity_confirmed']:
            self.assertFalse(result[key])
        self.assertEqual(before,original,'Reading provenance must not mutate the source rows')
    def test_historical_adoption_is_not_silently_promoted(self):
        old=[dict(sample_id=s,decision_kind=LEGACY_KIND) for s in REQUIRED_SAMPLES]
        result=provenance_from_rows(old)
        self.assertEqual(result['cohort_analysis_scope'],'author_specified_inputs')
        self.assertFalse(result['physician_retention_decision_relayed'])
    def test_mixed_or_missing_kind_fails(self):
        for value in ['',LEGACY_KIND,'PHYSICIAN_SIGNED']:
            with self.subTest(value=value):
                rows=corrected_rows();rows[0]['decision_kind']=value
                with self.assertRaises(ValueError):provenance_from_rows(rows)
    def test_missing_and_duplicate_decisions_fail(self):
        rows=corrected_rows()
        for value in [rows[:-1],rows+[rows[0]],rows[:-1]+[rows[0]]]:
            with self.assertRaises(ValueError):provenance_from_rows(value)
    def test_communication_date_cannot_become_physician_date(self):
        rows=corrected_rows();rows[0]['physician_decision_date']=rows[0]['decision_date']
        with self.assertRaises(ValueError):provenance_from_rows(rows)
    def test_no_new_signature_independence_blinding_or_identity_can_be_inferred(self):
        for key in ['physician_signature_provided','independent_review_evidence_provided','decision_blinding_provided','independent_lesion_identity_confirmed']:
            for value in ['True','']:
                with self.subTest(key=key,value=value):
                    rows=corrected_rows();rows[0][key]=value
                    with self.assertRaises(ValueError):provenance_from_rows(rows)
    def test_changed_action_alias_or_actor_cannot_pass_this_correction_scope(self):
        for key,value in [('action','exclude'),('alias_of','83_1_2'),('decision_maker_role','author'),('communication_source_role','physician'),('decision_date_scope','physician_review_date'),('physician_retention_decision_relayed','False')]:
            with self.subTest(key=key):
                rows=corrected_rows();rows[0][key]=value
                with self.assertRaises(ValueError):provenance_from_rows(rows)
    def test_legacy_kind_cannot_hide_corrected_physician_actor(self):
        rows=corrected_rows()
        for row in rows:row['decision_kind']=LEGACY_KIND
        with self.assertRaises(ValueError):provenance_from_rows(rows)

if __name__=='__main__':unittest.main(verbosity=2)
