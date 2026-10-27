"""Decision attribution is distinct from signed, independent or complete review evidence.

The v1.1.0 author-adoption label remains readable only with its historical
meaning. A corrected physician decision relayed by the author is represented
explicitly; no missing medical-review date or evidence is inferred by code.
"""
import csv
from pathlib import Path

REQUIRED_SAMPLES = ['83_1_' + str(i) for i in range(5)] + ['94_1_' + str(i) for i in range(3)]
LEGACY_KIND = 'AUTHOR_RESEARCH_INPUT_ADOPTION'
PHYSICIAN_KIND = 'PHYSICIAN_DECISION_RELAYED_BY_AUTHOR'
UNKNOWN_VALUES = {'', 'unknown', 'not_provided', '未提供'}


def provenance_from_rows(rows):
    selected = [row for row in rows if row.get('sample_id') in REQUIRED_SAMPLES]
    if len(selected) != len(REQUIRED_SAMPLES) or {r.get('sample_id') for r in selected} != set(REQUIRED_SAMPLES):
        raise ValueError('All eight decision records must occur exactly once')
    kinds = {r.get('decision_kind') for r in selected}
    if kinds == {LEGACY_KIND}:
        if any(row.get('decision_maker_role') == 'physician' or str(row.get('physician_retention_decision_relayed', '')).lower() == 'true' for row in selected):
            raise ValueError('Historical author-adoption kind conflicts with corrected physician attribution')
        return dict(
            cohort_analysis_scope='author_specified_inputs',
            analysis_execution_basis='author_specified_inputs',
            decision_provenance_kind=LEGACY_KIND,
            physician_retention_decision_relayed=False,
            provenance_interpretation='Historical author-adoption metadata; not reinterpreted as a physician decision. The v1.1.1 correction is recorded separately.',
            physician_signature_provided=False, independent_review_evidence_provided=False,
            decision_blinding_provided=False, independent_lesion_identity_confirmed=False,
            physician_decision_date='not_provided')
    if kinds != {PHYSICIAN_KIND}:
        raise ValueError('Decision provenance must be explicit and consistent; unknown or mixed kinds cannot pass')
    for row in selected:
        sid = row['sample_id']
        required = {
            'decision_maker_role':'physician', 'communication_source_role':'author',
            'decision_date_scope':'author_communication_record_not_physician_decision_date'}
        for key, expected in required.items():
            if row.get(key) != expected:
                raise ValueError(sid + ': invalid or missing decision provenance field ' + key)
        if str(row.get('physician_retention_decision_relayed', '')).lower() != 'true':
            raise ValueError(sid + ': physician retention decision relay must be explicit')
        for key in ['physician_signature_provided','independent_review_evidence_provided',
                    'decision_blinding_provided','independent_lesion_identity_confirmed']:
            if str(row.get(key, '')).lower() != 'false':
                raise ValueError(sid + ': unsupported or missing evidence scope field ' + key)
        if 'physician_decision_date' not in row or str(row['physician_decision_date']).lower() not in UNKNOWN_VALUES:
            raise ValueError(sid + ': physician decision date was not supplied for this correction')
        if row.get('action') != 'retain' or row.get('alias_of', ''):
            raise ValueError(sid + ': this correction must retain the original mask without aliasing')
    return dict(
        cohort_analysis_scope='physician_decision_relayed_by_author_inputs',
        analysis_execution_basis='physician_decision_relayed_by_author_inputs',
        decision_provenance_kind=PHYSICIAN_KIND,
        physician_retention_decision_relayed=True,
        provenance_interpretation='The eight original masks from P120/P125 were retained without contour changes following a physician decision communicated by the author.',
        physician_signature_provided=False, independent_review_evidence_provided=False,
        decision_blinding_provided=False, independent_lesion_identity_confirmed=False,
        physician_decision_date='not_provided')


def decision_provenance(path):
    with Path(path).open(encoding='utf-8-sig', newline='') as handle:
        return provenance_from_rows(list(csv.DictReader(handle)))
