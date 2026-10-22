"""Fail-closed identity audit between extraction and every model stage.

Only reads frozen rosters, features, source images and the actual model matrix.
Never repairs an extraction failure by changing the declared cohort.
"""
from pathlib import Path
from collections import Counter
import csv
import hashlib
import json
import math
import datetime


class InputConsistencyError(RuntimeError):
    pass


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def _read(path):
    with Path(path).open(encoding='utf-8-sig', newline='') as handle:
        reader = csv.DictReader(handle)
        values = list(reader)
        if not reader.fieldnames or len(set(reader.fieldnames)) != len(reader.fieldnames):
            raise ValueError('Missing or duplicated CSV column names: ' + str(path))
        return values, reader.fieldnames


def _true(value):
    return str(value).lower() == 'true' or value is True


def _write(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2), encoding='utf-8')


def verify_analysis_inputs(run_root, *, revision=None, feature_names=None,
                           expected_patients=None, expected_images=None,
                           expected_representatives=None, report_path=None,
                           verify_image_hashes=True):
    """Require exact patient, canonical lesion and representative identities.

    Expected inputs are read from FINAL_* rosters by default. Active images are
    identified by active_analysis_record, not the number of mask files. All
    extraction failures stop the run, even if a representative succeeded.
    """
    root = Path(run_root)
    rev = Path(revision) if revision is not None else root / 'revision_R1_20261014'
    report_path = Path(report_path) if report_path else root / 'POST_EXTRACTION_INPUT_AUDIT.json'
    issues = []
    report = {'status': 'FAIL', 'checked_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
              'analysis_matrix': str(rev / '03_model_development/analysis_matrix.csv'),
              'feature_table': str(rev / '02_extraction/arterial_lesion_features.csv'),
              'issues': issues, 'model_training_authorized': False}

    def issue(kind, **details):
        issues.append({'type': kind, **details})

    def load(path, label):
        try:
            return _read(path)
        except (OSError, ValueError, csv.Error) as exc:
            issue('unreadable_required_input', input=label, path=str(path), reason=str(exc))
            return [], []

    representative_path = root / 'FINAL_PATIENT_AND_REPRESENTATIVE_ROSTER.csv'
    image_path = root / 'FINAL_IMAGE_AND_ALIAS_ROSTER.csv'
    if expected_representatives is None:
        representatives, _ = load(representative_path, 'final_representatives')
    else:
        representatives = list(expected_representatives)
    if expected_images is None:
        image_roster, _ = load(image_path, 'final_images_and_aliases')
    else:
        image_roster = list(expected_images)
    active = [r for r in image_roster
              if 'active_analysis_record' not in r or _true(r['active_analysis_record'])]
    if expected_patients is None:
        disposition_path = root / 'FINAL_CANDIDATE_DISPOSITION.csv'
        if disposition_path.is_file():
            disposition, _ = load(disposition_path, 'final_candidate_disposition')
            patients = [str(r['patient_id']) for r in disposition if r.get('eligibility') == 'include']
            report['candidate_disposition_sha256'] = sha256(disposition_path)
        else:
            patients = [str(r.get('patient_id', '')) for r in representatives]
    else:
        patients = [str(r.get('patient_id', '')) if isinstance(r, dict) else str(r)
                    for r in expected_patients]
    matrix_path = rev / '03_model_development/analysis_matrix.csv'
    feature_path = rev / '02_extraction/arterial_lesion_features.csv'
    failure_path = rev / '02_extraction/extraction_failures.csv'
    matrix, matrix_columns = load(matrix_path, 'actual_analysis_matrix')
    features, feature_columns = load(feature_path, 'actual_lesion_features')
    failures, _ = load(failure_path, 'extraction_failures')
    for records in [representatives, active, matrix, features]:
        for r in records:
            for key in ['patient_id', 'sample_id']:
                r[key] = str(r.get(key) or '')
    if failures:
        issue('extraction_failures', count=len(failures), records=failures)
    for label, records, key in [('expected_patients', [{'patient_id': p} for p in patients], 'patient_id'),
                                 ('final_representatives', representatives, 'patient_id'),
                                 ('active_images', active, 'sample_id'),
                                 ('actual_analysis_matrix', matrix, 'patient_id'),
                                 ('actual_lesion_features', features, 'sample_id')]:
        if not records:
            issue('empty_required_input', input=label)
        identifiers = [str(r.get(key, '')) for r in records]
        duplicated = sorted(x for x, n in Counter(identifiers).items() if n > 1)
        if duplicated:
            issue('duplicate_identifiers', input=label, key=key, identifiers=duplicated)
        if any(not x or x in ['None', 'nan'] for x in identifiers):
            issue('missing_identifier', input=label, key=key)
        if label != 'expected_patients' and any(not r.get('patient_id') or not r.get('sample_id') for r in records):
            issue('missing_identity_columns', input=label)

    def compare_sets(label, expected, actual):
        missing, unexpected = sorted(expected - actual), sorted(actual - expected)
        if missing or unexpected:
            issue('identity_set_mismatch', input=label, missing=missing, unexpected=unexpected)

    expected_pids = set(patients)
    image_ids = [str(r.get('sample_id') or '') for r in image_roster]
    duplicate_source_ids = sorted(sid for sid, n in Counter(image_ids).items() if n > 1)
    if duplicate_source_ids:
        issue('duplicate_final_image_roster_identifiers', sample_ids=duplicate_source_ids)
    for r in image_roster:
        if str(r.get('patient_id')) in expected_pids and r.get('action') in ['retain', 'replace', 'voxel_edit'] and 'active_analysis_record' in r and not _true(r['active_analysis_record']):
            issue('included_canonical_marked_inactive', sample_id=r.get('sample_id'), patient_id=r.get('patient_id'))
    compare_sets('final_patient_and_representative_roster', expected_pids,
                 {str(r.get('patient_id', '')) for r in representatives})
    compare_sets('actual_analysis_matrix_patients', expected_pids,
                 {str(r.get('patient_id', '')) for r in matrix})
    compare_sets('active_lesion_patients', expected_pids,
                 {str(r.get('patient_id', '')) for r in active})
    compare_sets('actual_lesion_features', {r.get('sample_id') for r in active},
                 {r.get('sample_id') for r in features})
    expected_by_sid = {r.get('sample_id'): r for r in active}
    actual_by_sid = {r.get('sample_id'): r for r in features}
    matrix_by_pid = {str(r.get('patient_id')): r for r in matrix}
    rep_by_pid = {str(r.get('patient_id')): r for r in representatives}
    expected_counts = Counter(str(r.get('patient_id')) for r in active)
    for r in active:
        if ('action' in r and r['action'] not in ['retain', 'replace', 'voxel_edit']) or r.get('alias_of') not in [None, '', '未提供']:
            issue('noncanonical_record_marked_active', sample_id=r.get('sample_id'), action=r.get('action'), alias_of=r.get('alias_of'))
    for pid in expected_pids:
        ct_hashes = {r.get('ct_sha256') for r in active if str(r.get('patient_id')) == pid}
        if len(ct_hashes) > 1:
            issue('patient_uses_multiple_CT_series', patient_id=pid)
    image_hash_cache = {}
    for sid in sorted(set(expected_by_sid) & set(actual_by_sid)):
        e, a = expected_by_sid[sid], actual_by_sid[sid]
        if str(e.get('patient_id')) != str(a.get('patient_id')):
            issue('lesion_patient_mismatch', sample_id=sid, expected=e.get('patient_id'), actual=a.get('patient_id'))
        for key, target in [('ct_path', 'ct_path'), ('final_mask_path', 'mask_path')]:
            expected_path = e.get(key) or (e.get('mask_path') if key == 'final_mask_path' else None)
            if expected_path and Path(expected_path).resolve() != Path(a.get(target, '')).resolve():
                issue('image_path_mismatch', sample_id=sid, field=target, expected=expected_path, actual=a.get(target))
        if verify_image_hashes:
            for field, path_field in [('ct_sha256', 'ct_path'), ('final_mask_sha256', 'mask_path')]:
                expected_hash = e.get(field) or (e.get('mask_sha256') if field == 'final_mask_sha256' else None)
                path = a.get(path_field)
                if not expected_hash or not path:
                    issue('missing_image_hash_binding', sample_id=sid, field=field)
                    continue
                try:
                    resolved = str(Path(path).resolve())
                    if resolved not in image_hash_cache:
                        image_hash_cache[resolved] = sha256(path)
                    actual_hash = image_hash_cache[resolved]
                    if actual_hash != expected_hash:
                        issue('image_hash_mismatch', sample_id=sid, field=field, expected=expected_hash, actual=actual_hash)
                except OSError as exc:
                    issue('unreadable_source_image', sample_id=sid, field=path_field, reason=str(exc))
    for r in representatives:
        pid, sid = str(r.get('patient_id')), r.get('sample_id')
        if sid not in expected_by_sid or str(expected_by_sid[sid].get('patient_id')) != pid:
            issue('representative_not_active_canonical', patient_id=pid, sample_id=sid)
    # The declared representative is also checked against actual extracted
    # original-grid physical volumes, including a deterministic tie breaker.
    ranked = {}
    for r in features:
        pid, sid = str(r.get('patient_id')), r.get('sample_id')
        try:
            volume = float(r.get('physical_volume_mm3', ''))
            if not math.isfinite(volume) or volume <= 0:
                raise ValueError('Non-positive or non-finite volume')
            if pid not in ranked or (-volume, sid) < ranked[pid][0]:
                ranked[pid] = ((-volume, sid), sid)
            e = expected_by_sid.get(sid, {})
            if e.get('physical_volume_mm3') not in [None, ''] and not math.isclose(volume, float(e['physical_volume_mm3']), rel_tol=1e-12, abs_tol=1e-8):
                issue('lesion_volume_mismatch', sample_id=sid, expected=e['physical_volume_mm3'], actual=volume)
        except (TypeError, ValueError) as exc:
            issue('invalid_lesion_volume', sample_id=sid, reason=str(exc))
    for pid in sorted(expected_pids):
        e, a = rep_by_pid.get(pid, {}), matrix_by_pid.get(pid, {})
        if e.get('sample_id') != a.get('sample_id'):
            issue('representative_mismatch', patient_id=pid, expected=e.get('sample_id'), actual=a.get('sample_id'))
        ranked_sid = ranked.get(pid, (None, None))[1]
        if e.get('sample_id') != ranked_sid:
            issue('representative_volume_rule_mismatch', patient_id=pid, expected=e.get('sample_id'), ranked=ranked_sid)
        try:
            if float(a.get('n_segmented_arterial_lesions', '')) != expected_counts[pid]:
                issue('active_lesion_count_mismatch', patient_id=pid, expected=expected_counts[pid], actual=a.get('n_segmented_arterial_lesions'))
        except (TypeError, ValueError):
            issue('missing_active_lesion_count', patient_id=pid)
    if feature_names is None:
        _, reference_columns = load(Path(__file__).resolve().parent / 'data/analysis_matrix.csv', 'protocol_feature_schema')
        feature_names = [x for x in reference_columns if x.startswith('original_')]
    feature_names = list(feature_names)
    if len(feature_names) != 107 or len(set(feature_names)) != 107:
        issue('invalid_protocol_feature_schema', feature_count=len(feature_names))
    for label, records, columns in [('actual_analysis_matrix', matrix, matrix_columns),
                                     ('actual_lesion_features', features, feature_columns)]:
        found = {x for x in columns if x.startswith('original_')}
        compare_sets(label + '_107_feature_columns', set(feature_names), found)
        bad = []
        for r in records:
            for key in feature_names:
                try:
                    if not math.isfinite(float(r.get(key, ''))):
                        raise ValueError('non-finite')
                except (TypeError, ValueError):
                    bad.append({'sample_id': r.get('sample_id'), 'feature': key})
        if bad:
            issue('missing_or_nonfinite_features', input=label, count=len(bad), records=bad)
        if 'ct_phase' in columns:
            wrong_phase = [r.get('sample_id') for r in records if not _same_number(r.get('ct_phase'), 1)]
            if wrong_phase:
                issue('non_target_phase_record', input=label, sample_ids=wrong_phase)
    for pid, a in matrix_by_pid.items():
        f = actual_by_sid.get(a.get('sample_id'))
        if f is None:
            continue
        changed = [key for key in feature_names if a.get(key) != f.get(key)
                   and not _same_number(a.get(key), f.get(key))]
        if changed:
            issue('representative_feature_values_mismatch', patient_id=pid, sample_id=a.get('sample_id'), features=changed)
        for key in ['ct_path', 'mask_path', 'physical_volume_mm3']:
            if key in f and key in a and a[key] != f[key] and not (key == 'physical_volume_mm3' and _same_number(a[key], f[key])):
                issue('representative_source_mismatch', patient_id=pid, field=key, expected=f[key], actual=a[key])
    clinical_path = rev / '01_data_audit/source_clinical_233.csv'
    if clinical_path.is_file():
        clinical, clinical_columns = load(clinical_path, 'final_source_clinical_values')
        clinical_by_pid = {str(r.get('patient_id')): r for r in clinical}
        for pid in expected_pids:
            baseline, actual = clinical_by_pid.get(pid), matrix_by_pid.get(pid, {})
            if baseline is None:
                issue('missing_source_clinical_patient', patient_id=pid)
                continue
            changed = [k for k in clinical_columns if k not in actual or
                       (baseline[k] != actual[k] and not _same_number(baseline[k], actual[k]))]
            if changed:
                issue('source_clinical_value_mismatch', patient_id=pid, fields=changed)
    for key, path in [('analysis_matrix_sha256', matrix_path), ('feature_table_sha256', feature_path),
                      ('extraction_failures_sha256', failure_path)]:
        if path.is_file():
            report[key] = sha256(path)
    report.update(actual_patient_count=len(matrix), actual_unique_patient_count=len(matrix_by_pid),
                  actual_active_unique_mask_count=len(features), actual_unique_sample_count=len(actual_by_sid),
                  expected_patient_count=len(patients), expected_active_unique_mask_count=len(active),
                  expected_feature_count=len(feature_names), actual_patient_ids=sorted(matrix_by_pid),
                  actual_representatives={pid: r.get('sample_id') for pid, r in sorted(matrix_by_pid.items())},
                  source_image_hashes_verified=verify_image_hashes)
    report['bound_files'] = {str(path): sha256(path) for path in
                             [matrix_path, feature_path, failure_path, representative_path, image_path,
                              root / 'FINAL_CANDIDATE_DISPOSITION.csv'] if path.is_file()}
    for path in [clinical_path, rev / '01_data_audit/phase_pairing_audit.csv',
                 rev / '01_data_audit/phase_analysis_eligibility.csv']:
        if path.is_file():
            report['bound_files'][str(path)] = sha256(path)
    report['status'] = 'FAIL' if issues else 'PASS'
    report['model_training_authorized'] = not issues
    _write(report_path, report)
    if issues:
        raise InputConsistencyError(f'Post-extraction input audit failed ({len(issues)} issues). No training may start. Report: {report_path}')
    return report


def _same_number(a, b):
    try:
        return math.isfinite(float(a)) and math.isfinite(float(b)) and math.isclose(float(a), float(b), rel_tol=1e-12, abs_tol=1e-12)
    except (TypeError, ValueError):
        return False


def assert_inputs_unchanged(run_root, audit, stage):
    """Bind each model stage to exactly the audited matrix and feature bytes."""
    changed = []
    for path_key, hash_key in [('analysis_matrix', 'analysis_matrix_sha256'), ('feature_table', 'feature_table_sha256')]:
        path = Path(audit[path_key])
        actual = sha256(path) if path.is_file() else None
        if actual != audit[hash_key]:
            changed.append({'input': path_key, 'expected_sha256': audit[hash_key], 'actual_sha256': actual})
    for file, expected in audit.get('bound_files', {}).items():
        path = Path(file)
        actual = sha256(path) if path.is_file() else None
        if actual != expected:
            changed.append({'input': file, 'expected_sha256': expected, 'actual_sha256': actual})
    if changed:
        report = {'status': 'FAIL', 'stage': stage, 'issues': changed, 'completion_authorized': False}
        _write(Path(run_root) / 'MODEL_INPUT_CHANGED.json', report)
        raise InputConsistencyError('Model input changed after successful audit: ' + stage)


def run_guarded_training(run_root, run_stage, *, code=None, verify_image_hashes=True, feature_names=None):
    """Single route from extracted matrix to diagnostics and full validation."""
    audit = verify_analysis_inputs(run_root, verify_image_hashes=verify_image_hashes, feature_names=feature_names)
    from validation_integrity import audit_validation
    for stage, args in [('diagnostics_and_tests.py', ()),
                        ('validate_pipeline.py', ('--mode', 'point')),
                        ('validate_pipeline.py', ('--mode', 'bootstrap')),
                        ('validate_pipeline.py', ('--mode', 'summary'))]:
        assert_inputs_unchanged(run_root, audit, stage)
        if args == ('--mode', 'summary'):
            audit_validation(Path(run_root) / 'revision_R1_20261014')
        run_stage(stage, *args)
    assert_inputs_unchanged(run_root, audit, 'after_full_validation')
    audit_validation(Path(run_root) / 'revision_R1_20261014')
    return audit


def verify_model_output(run_root, audit):
    """Check final development scores contain exactly the audited model inputs."""
    root = Path(run_root)
    assert_inputs_unchanged(root, audit, 'completion')
    path = Path(audit['analysis_matrix']).parent / 'FINAL_REVISION_PATIENT_SCORES.csv'
    matrix, _ = _read(audit['analysis_matrix'])
    try:
        scores, _ = _read(path)
    except (OSError, ValueError) as exc:
        _write(root / 'MODEL_OUTPUT_INPUT_AUDIT.json', {'status': 'FAIL', 'reason': str(exc)})
        raise InputConsistencyError('Actual model output is unavailable') from exc
    fields = set(matrix[0]) if matrix else set()
    # Scores are created by appending a score to the complete actual matrix.
    score_by_pid = {str(r.get('patient_id')): r for r in scores}
    problems = []
    if len(scores) != len(matrix) or len(score_by_pid) != len(scores) or set(score_by_pid) != set(audit['actual_patient_ids']):
        problems.append('Model score patient identities differ from the audited matrix')
    for row in matrix:
        scored = score_by_pid.get(str(row['patient_id']), {})
        mismatch = [k for k in fields if k not in scored or (row[k] != scored[k] and not _same_number(row[k], scored[k]))]
        if mismatch:
            problems.append({'patient_id': row['patient_id'], 'changed_or_missing_input_columns': sorted(mismatch)})
    result = {'status': 'FAIL' if problems else 'PASS', 'issues': problems,
              'actual_model_patient_count': len(scores), 'actual_model_unique_patient_count': len(score_by_pid),
              'model_score_sha256': sha256(path), 'analysis_matrix_sha256': audit['analysis_matrix_sha256']}
    _write(root / 'MODEL_OUTPUT_INPUT_AUDIT.json', result)
    if problems:
        raise InputConsistencyError('Actual model output disagrees with the audited input matrix')
    return result


def write_completion_marker(run_root, audit):
    """Only actual audited matrix/feature rows and verified scores supply counts."""
    root = Path(run_root)
    actual_model = verify_model_output(root, audit)
    from validation_integrity import audit_validation, verify_paired_summary
    validation=audit_validation(root / 'revision_R1_20261014')
    paired=verify_paired_summary(root / 'revision_R1_20261014',validation)
    matrix, _ = _read(audit['analysis_matrix'])
    features, _ = _read(audit['feature_table'])
    log_path=root/'ACTUAL_EXECUTION_LOG.json'
    provenance_path=root/'ACTUAL_RUN_PROVENANCE.json'
    executed=json.loads(log_path.read_text(encoding='utf-8')) if log_path.is_file() else []
    provenance=json.loads(provenance_path.read_text(encoding='utf-8')) if provenance_path.is_file() else {}
    reextracted=any(s.get('stage')=='extract_matrix.py' for s in executed) and bool(features) and len([k for k in features[0] if k.startswith('original_')])==107 and all(
        r.get('feature_provenance')=='newly extracted with locked corrected settings' for r in features)
    marker = {
        'status': 'COMPLETE',
        'data_scope': audit.get('data_scope','REAL_ANALYSIS'),
        'actual_new_results': str(root / 'revision_R1_20261014/04_validation'),
        'final_patient_count': len(matrix),
        'actual_model_patient_count': actual_model['actual_model_patient_count'],
        'final_active_unique_mask_count': len(features),
        'analysis_matrix_sha256': audit['analysis_matrix_sha256'],
        'feature_table_sha256': audit['feature_table_sha256'],
        'model_score_sha256': actual_model['model_score_sha256'],
        'validation_counts': validation['counts'],
        'paired_summary_sha256': paired['paired_csv_sha256'],
        'validation_integrity_report': str(root / 'revision_R1_20261014/04_validation/VALIDATION_INTEGRITY.json'),
        'validation_integrity_sha256': sha256(root / 'revision_R1_20261014/04_validation/VALIDATION_INTEGRITY.json'),
        'all_107_features_reextracted': reextracted,
        'old_checkpoints_copied': provenance.get('copied_old_checkpoints'),
        'remaining': 'Sync manuscripts, all tables/figures, 18 responses and release archive to these actual results; inspect rendered documents before submission.'}
    _write(root / 'ANALYSIS_COMPLETE_AWAITING_DOCUMENT_SYNC.json', marker)
    return marker


def verify_source_audit_inputs(run_root, *, archived_matrix=False):
    """Guard the original reproducer too; it cannot bypass the extraction gate.

    Post-review rosters take priority whenever present. Without them this entry
    only reproduces the source audit's declared technical cohort; it does not
    label that historical cohort as a newly medically finalized cohort.
    """
    root = Path(run_root)
    if any((root / name).exists() for name in ['FINAL_PATIENT_AND_REPRESENTATIVE_ROSTER.csv', 'FINAL_IMAGE_AND_ALIAS_ROSTER.csv']):
        return verify_analysis_inputs(root, verify_image_hashes=not archived_matrix)
    rev = root / 'revision_R1_20261014'
    try:
        pair, _ = _read(rev / '01_data_audit/phase_pairing_audit.csv')
        eligibility, _ = _read(rev / '01_data_audit/phase_analysis_eligibility.csv')
        patients = [r['patient_id'] for r in eligibility if _true(r.get('phase_analysis_eligible'))]
        active = [dict(r, active_analysis_record=True, final_mask_path=r.get('mask_path'),
                       final_mask_sha256=r.get('mask_sha256')) for r in pair
                  if str(r.get('ct_phase')) == '1' and _true(r.get('eligible')) and r['patient_id'] in patients]
        ranked = sorted(active, key=lambda r: (r['patient_id'], -float(r['physical_volume_mm3']), r['sample_id']))
        representatives = {}
        for r in ranked:
            representatives.setdefault(r['patient_id'], r)
    except (OSError, ValueError, KeyError) as exc:
        _write(root / 'POST_EXTRACTION_INPUT_AUDIT.json', {'status': 'FAIL', 'scope': 'source_audit_reproduction', 'reason': str(exc), 'model_training_authorized': False})
        raise InputConsistencyError('Required source audit is unavailable; no training') from exc
    audit = verify_analysis_inputs(root, expected_patients=patients, expected_images=active,
                                   expected_representatives=list(representatives.values()),
                                   verify_image_hashes=not archived_matrix)
    audit['scope'] = 'historical_source_audit_reproduction_not_new_medical_finalization'
    _write(root / 'POST_EXTRACTION_INPUT_AUDIT.json', audit)
    return audit
