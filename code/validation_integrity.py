"""Fail-closed audit of the locked TACE validation protocol.

Original identity and task audits are extended by stored-survival evidence checks.
Checkpoint presence never implies that a task was attempted or completed successfully.
"""
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
import csv
import gzip
import hashlib
import json
import math
import os
import re
import uuid
import shutil
import statistics
from survival_evidence import source_fingerprint, evidence_issues, SCHEMA_VERSION


class ValidationIntegrityError(RuntimeError):
    def __init__(self, message, report=None):
        super().__init__(message)
        self.report = report or {}


CLINICAL = ['age', 'sex_female', 'child_b', 'log_afp', 'log_tumour_count', 'bclc_b']
BASE_SPECS = {
    'clinical': CLINICAL,
    'clinical_diameter': CLINICAL + ['selected_diameter_mm'],
    'clinical_selected_volume': CLINICAL + ['log_selected_volume'],
    'clinical_total_volume': CLINICAL + ['log_total_volume'],
    'clinical_mean_hu': CLINICAL + ['mean_hu'],
    'clinical_median_hu': CLINICAL + ['median_hu'],
    'clinical_p10_hu': CLINICAL + ['p10_hu'],
    'clinical_diameter_mean_hu': CLINICAL + ['selected_diameter_mm', 'mean_hu'],
    'clinical_diameter_p10_hu': CLINICAL + ['selected_diameter_mm', 'p10_hu'],
    'clinical_axial_diameter': CLINICAL + ['clinical_max_axial_diameter_mm'],
    'clinical_legacy_seven': ['age', 'sex_female', 'child_b', 'log_afp', 'afp_high_400', 'tumour_multiple', 'bclc_b'],
}
PAIRED_BASELINES = ('clinical', 'clinical_diameter', 'clinical_selected_volume', 'clinical_total_volume', 'clinical_diameter_mean_hu', 'clinical_diameter_p10_hu')
for _name in PAIRED_BASELINES:
    BASE_SPECS[_name + '_radiomics'] = BASE_SPECS[_name] + ['radiomics_score']
BASE_SPECS['radiomics_only'] = ['radiomics_score']
BASE_SPECS['clinical_age_spline'] = CLINICAL
RSF_SPECS = {
    'clinical_size_intensity_RSF': CLINICAL + ['selected_diameter_mm', 'mean_hu'],
    'clinical_size_intensity_radiomics_RSF': CLINICAL + ['selected_diameter_mm', 'mean_hu', 'radiomics_score'],
}
LOCKED_COUNTS = {'point': 10, 'bootstrap': 500}
LOCKED_SCIENTIFIC_CONFIG = {
    'primary_phase': 1, 'geometry_numeric_tolerance': 0.0001,
    'features_original': 107, 'correlation_threshold': 0.85,
    'max_filtered_features': 12, 'seed_point': 20261014,
    'seed_bootstrap': 20261015, 'baseline_penalizer': 0.1,
    'candidates': ['ridge', 'top_two_cox', 'top_two_elastic_net'],
    'ridge_penalties': [0.1, 1.0, 10.0],
    'two_feature_penalties': [0.01, 0.1, 1.0],
    'prediction_months': [12, 24, 36], 'ibs_monthly_grid': [1, 36],
}
INPUT_KEYS = ('analysis_matrix_sha256', 'configuration_sha256', 'configuration_file_sha256')
CHECKPOINT_RE = re.compile(r'^(point|bootstrap)_([0-9]{4})\.json\.gz$')


def _utc():
    return datetime.now(timezone.utc).isoformat()


def _sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def _canonical_sha(value):
    # Matches the historical validator's json.dumps(CFG, sort_keys=True).
    return hashlib.sha256(json.dumps(value, sort_keys=True).encode()).hexdigest()


def _paths(revision):
    root = Path(revision)
    return root, root / '04_validation', root / '00_protocol_and_provenance/revision_config.json', root / '03_model_development/analysis_matrix.csv'


def _read_json(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))


def _write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + '.writing')
    tmp.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False), encoding='utf-8')
    os.replace(str(tmp), str(path))


def _issue(issues, code, message, **location):
    issues.append(dict(code=code, message=message, **location))


def _publish_integrity_report(root, out, report):
    """Persist current status and preserve, then revoke, any stale active marker."""
    if report['status'] != 'COMPLETE':
        marker = root.parent / 'ANALYSIS_COMPLETE_AWAITING_DOCUMENT_SYNC.json'
        if marker.is_file():
            archive = out / 'invalidated_completion_markers'
            archive.mkdir(parents=True,exist_ok=True)
            target = archive / (datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S')+'_'+uuid.uuid4().hex+'.json')
            os.replace(str(marker),str(target))
            report['previous_completion_marker_preserved_and_invalidated'] = str(target)
        _write_json(root.parent/'ANALYSIS_STATUS.json',{
            'status':report['status'],'validation_status':report['status'],
            'formal_outputs_authorized':False,'audit_utc':report['audited_utc'],
            'integrity_report':str(out/'VALIDATION_INTEGRITY.json')})
    _write_json(out/'VALIDATION_INTEGRITY.json',report)


def _number(value):
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value)


def _integer(value):
    return isinstance(value, int) and not isinstance(value, bool)


def _fingerprints(revision):
    _, _, cfg_path, matrix_path = _paths(revision)
    observed = {}
    for key, path in [('analysis_matrix_sha256', matrix_path), ('configuration_file_sha256', cfg_path)]:
        try:
            observed[key] = _sha(path)
        except Exception as exc:
            observed[key] = None
            observed[key + '_error'] = str(exc)
    try:
        observed['configuration_sha256'] = _canonical_sha(_read_json(cfg_path))
    except Exception as exc:
        observed['configuration_sha256'] = None
        observed['configuration_sha256_error'] = str(exc)
    return observed


def _build_plan(revision, specs=None):
    _, _, cfg_path, matrix_path = _paths(revision)
    cfg = _read_json(cfg_path)
    for key, expected in LOCKED_SCIENTIFIC_CONFIG.items():
        if cfg.get(key) != expected:
            raise ValidationIntegrityError('Locked scientific protocol mismatch: ' + key)
    for key, expected in [('point_repeats', 10), ('bootstrap_replicates', 500), ('outer_folds', 5), ('inner_folds', 3)]:
        if not _integer(cfg.get(key)) or cfg[key] != expected:
            raise ValidationIntegrityError('Locked protocol mismatch: %s must be %s' % (key, expected))
    primary = ['clinical_diameter_mean_hu', 'clinical_diameter_mean_hu_radiomics']
    if cfg.get('primary_comparison') != primary:
        raise ValidationIntegrityError('Locked primary comparison mismatch')
    for seed_key in ('seed_point', 'seed_bootstrap'):
        if not _integer(cfg.get(seed_key)):
            raise ValidationIntegrityError('Missing integer protocol seed: ' + seed_key)
    if specs is not None and specs != BASE_SPECS:
        raise ValidationIntegrityError('Model specifications differ from the locked 19-model registry')
    with matrix_path.open(encoding='utf-8-sig', newline='') as handle:
        matrix = list(csv.DictReader(handle))
    if not matrix:
        raise ValidationIntegrityError('Actual model input matrix is empty')
    ids, outcomes = [], {}
    for row in matrix:
        pid = str(row.get('patient_id', '')).strip()
        if not pid or pid in outcomes:
            raise ValidationIntegrityError('Actual model input has missing/duplicate patient identity')
        try:
            event = float(row['death_event'])
            time = float(row['os_months'])
        except (KeyError, TypeError, ValueError) as exc:
            raise ValidationIntegrityError('Invalid matrix survival outcome: ' + str(exc)) from exc
        if event not in (0, 1) or not math.isfinite(time) or time <= 0:
            raise ValidationIntegrityError('Invalid matrix survival outcome for patient ' + pid)
        ids.append(pid)
        outcomes[pid] = {'event': int(event), 'time': time}
    n, deaths = len(ids), sum(x['event'] for x in outcomes.values())
    feature_names = sorted(k for k in matrix[0] if k.startswith('original_'))
    if len(feature_names) != 107:
        raise ValidationIntegrityError('Actual matrix must contain all 107 protocol features')
    applicable = dict(BASE_SPECS)
    if n >= 100 and deaths >= 70:
        applicable.update(RSF_SPECS)
    return dict(
        schema_version=2, protocol='TACE_locked_10x5_plus_500_patient_bootstrap',
        prediction_evidence_schema_version=SCHEMA_VERSION,
        execution_source_sha256=source_fingerprint(Path(__file__).resolve().parent),
        inputs={k: v for k, v in _fingerprints(revision).items() if k in INPUT_KEYS},
        configuration=cfg, planned=dict(LOCKED_COUNTS), outer_folds=5, inner_folds=3,
        applicable_models=list(applicable), model_specs=applicable,
        primary_comparison=primary,
        matrix={'n': n, 'events': deaths, 'patient_ids': ids, 'outcomes': outcomes,
                'feature_names':feature_names},
        tasks=[dict(mode=mode, run_id=i, seed=cfg['seed_' + mode] + i)
               for mode, count in LOCKED_COUNTS.items() for i in range(count)],
    )


def ensure_validation_plan(revision, specs=None):
    """Create an immutable plan, or reject any drift from its locked input files."""
    _, out, _, _ = _paths(revision)
    plan_path = out / 'VALIDATION_PLAN.json'
    candidate = _build_plan(revision, specs)
    if plan_path.exists():
        saved = _read_json(plan_path)
        comparison = {k: v for k, v in saved.items() if k not in ('created_utc', 'plan_sha256')}
        if comparison != candidate or saved.get('plan_sha256') != _canonical_sha(comparison):
            raise ValidationIntegrityError('Validation plan/input/configuration drift; use a new run directory')
        return saved
    plan = dict(candidate, created_utc=_utc(), plan_sha256=_canonical_sha(candidate))
    _write_json(plan_path, plan)
    return plan


def record_task_event(revision, mode, run_id, event, *, attempt_id, exception=None, details=None):
    """Append one task event even when inputs have changed or become unreadable.

    Called by the single parent executor, so parallel workers never race to append
    the ledger.  All earlier failures are retained when a task is retried.
    """
    if event not in ('attempted', 'succeeded', 'failed', 'reused'):
        raise ValueError('Unsupported task event: ' + str(event))
    _, out, _, _ = _paths(revision)
    out.mkdir(parents=True, exist_ok=True)
    try:
        plan = _read_json(out / 'VALIDATION_PLAN.json')
    except Exception:
        plan = {}
    cfg = plan.get('configuration', {})
    seed = cfg.get('seed_' + str(mode))
    if _integer(seed) and _integer(run_id):
        seed += run_id
    else:
        seed = None
    if exception is not None and not isinstance(exception, (str, dict)):
        exception = {'type': type(exception).__name__, 'message': str(exception)}
    value = dict(schema_version=1, utc=_utc(), mode=mode, run_id=run_id,
                 attempt_id=str(attempt_id), event=event, seed=seed,
                 original_seed=seed, exception=exception,
                 input_fingerprints=plan.get('inputs', {}),
                 observed_input_fingerprints=_fingerprints(revision),
                 plan_sha256=plan.get('plan_sha256'), details=details or {})
    data = json.dumps(value, ensure_ascii=False, allow_nan=False) + '\n'
    with (out / 'VALIDATION_TASK_EVENTS.jsonl').open('a', encoding='utf-8', newline='') as handle:
        handle.write(data)
        handle.flush()
        os.fsync(handle.fileno())
    if event == 'failed':
        with (out / 'VALIDATION_TASK_FAILURES.jsonl').open('a', encoding='utf-8', newline='') as handle:
            handle.write(data)
            handle.flush()
            os.fsync(handle.fileno())
    return value


def _transform_issues(value, expected_columns, issues, location):
    if not isinstance(value, dict):
        _issue(issues, 'transform_missing', 'Training transform is missing', **location)
        return
    columns = value.get('columns')
    if columns != expected_columns:
        _issue(issues, 'transform_columns', 'Transform columns do not match locked model/development features', **location)
    length = len(expected_columns)
    for key in ('median', 'mean', 'sd', 'keep'):
        values = value.get(key)
        if not isinstance(values, list) or len(values) != length:
            _issue(issues, 'transform_shape', 'Transform %s has an incorrect length' % key, **location)
            continue
        if key == 'keep':
            if any(not isinstance(x, bool) for x in values):
                _issue(issues, 'transform_keep', 'Transform keep flags must be boolean', **location)
        elif any(not _number(x) or (key == 'sd' and x <= 0) for x in values):
            _issue(issues, 'transform_nonfinite', 'Transform %s contains invalid values' % key, **location)


def checkpoint_issues(checkpoint_dict, plan, path=None):
    """Audit a single loaded checkpoint without scanning the other 509 runs."""
    issues = []
    run = checkpoint_dict
    location = {'checkpoint': str(path) if path is not None else None}
    if not isinstance(run, dict):
        _issue(issues, 'checkpoint_schema', 'Checkpoint root must be an object', **location)
        return issues
    mode, run_id = run.get('mode'), run.get('run_id')
    location.update(mode=mode, run_id=run_id)
    valid_id = mode in LOCKED_COUNTS and _integer(run_id) and 0 <= run_id < LOCKED_COUNTS.get(mode, 0)
    if not valid_id:
        _issue(issues, 'checkpoint_identity', 'Invalid checkpoint mode/run identifier', **location)
    if path is not None:
        match = CHECKPOINT_RE.fullmatch(Path(path).name)
        if not match or match.group(1) != mode or int(match.group(2)) != run_id:
            _issue(issues, 'checkpoint_filename', 'Filename and checkpoint identity differ', **location)
    for key in INPUT_KEYS:
        if run.get(key) != plan['inputs'].get(key):
            _issue(issues, 'checkpoint_input_fingerprint', 'Checkpoint %s differs from locked actual input' % key, **location)
    if valid_id and run.get('seed') != plan['configuration']['seed_' + mode] + run_id:
        _issue(issues, 'checkpoint_seed', 'Checkpoint seed differs from its original protocol seed', **location)
    if run.get('n') != plan['matrix']['n']:
        _issue(issues, 'checkpoint_population', 'Run row count differs from actual model input', **location)
    if run.get('failures') != []:
        _issue(issues, 'checkpoint_task_failures', 'Checkpoint contains task/fold failures or lacks a failure list', failures=run.get('failures'), **location)
    folds = run.get('folds')
    if not isinstance(folds, list):
        _issue(issues, 'folds_missing', 'Outer-fold records are missing', **location)
        return issues
    fold_ids = [f.get('fold') if isinstance(f, dict) else None for f in folds]
    if len(folds) != 5 or any(not _integer(i) for i in fold_ids) or sorted(i for i in fold_ids if _integer(i)) != list(range(5)):
        _issue(issues, 'fold_identity', 'Expected exactly outer folds 0..4 once each', **location)
    applicable = plan['applicable_models']
    expected_ids = set(plan['matrix']['patient_ids'])
    outcomes = plan['matrix']['outcomes']
    validation_counter, sample_counters = Counter(), []
    statistics = {name: [] for name in applicable}
    for record in folds:
        if not isinstance(record, dict):
            _issue(issues, 'fold_schema', 'Fold must be an object', **location)
            continue
        loc = dict(location, fold=record.get('fold'))
        train, val = record.get('train_ids'), record.get('validation_ids')
        ids_ok = isinstance(train, list) and isinstance(val, list) and bool(train) and bool(val)
        if not ids_ok or any(not isinstance(x, str) or x not in expected_ids for x in (train or []) + (val or [])):
            _issue(issues, 'fold_population', 'Fold patient identifiers are missing or outside the actual matrix', **loc)
            train, val = [], []
        if set(train) & set(val):
            _issue(issues, 'fold_patient_leakage', 'Training and validation patient identities overlap', **loc)
        counter = Counter(train + val)
        sample_counters.append(counter)
        validation_counter.update(val)
        if sum(counter.values()) != plan['matrix']['n']:
            _issue(issues, 'fold_sample_count', 'Train plus validation rows differ from actual sample size', **loc)
        if mode == 'point' and counter != Counter(plan['matrix']['patient_ids']):
            _issue(issues, 'point_population', 'Point run must contain each input patient exactly once', **loc)
        development = record.get('radiomics_development')
        if not isinstance(development, dict):
            _issue(issues, 'development_missing', 'Training-only radiomics development is missing', **loc)
        else:
            chosen, coefficients = development.get('features'), development.get('coefficients')
            if not isinstance(chosen, list) or not chosen or len(chosen) != len(set(chosen)) or any(not isinstance(x, str) or not x.startswith('original_') for x in chosen):
                _issue(issues, 'development_features', 'Development features are missing or invalid', **loc)
                chosen = []
            if any(x not in plan['matrix']['feature_names'] for x in chosen):
                _issue(issues, 'development_feature_not_in_matrix', 'Selected radiomics feature is absent from actual matrix', **loc)
            if not isinstance(coefficients, list) or len(coefficients) != len(chosen) or any(not _number(x) for x in coefficients):
                _issue(issues, 'development_coefficients', 'Selected feature coefficients are invalid', **loc)
            _transform_issues(development.get('transform'), chosen, issues, loc)
            tuning = development.get('tuning')
            matching = []
            if isinstance(tuning, list):
                matching = [x for x in tuning if isinstance(x, dict) and x.get('kind') == development.get('kind') and x.get('penalty') == development.get('penalty') and x.get('success_folds') == 3 and x.get('failed_folds') == 0 and x.get('failure') == [] and _number(x.get('C_index'))]
            if not matching or development.get('kind') not in plan['configuration'].get('candidates', []) or not _number(development.get('penalty')):
                _issue(issues, 'development_tuning', 'No valid three-fold tuning record supports the fitted radiomics specification', **loc)
            cfg = plan['configuration']
            expected_grid = {(kind, float(p)) for kind in cfg['candidates'] for p in
                             (cfg['ridge_penalties'] if kind == 'ridge' else cfg['two_feature_penalties'])}
            tuning_records = tuning if isinstance(tuning,list) else []
            actual_grid = [(t.get('kind'),t.get('penalty')) for t in tuning_records if isinstance(t,dict)]
            if len(actual_grid) != len(expected_grid) or set(actual_grid) != expected_grid:
                _issue(issues, 'development_tuning_grid', 'All locked inner tuning candidates must have explicit records', **loc)
            for t in tuning_records:
                if not isinstance(t,dict):continue
                success, failed = t.get('success_folds'),t.get('failed_folds')
                if not _integer(success) or not _integer(failed) or success < 0 or failed < 0 or success+failed != 3 or not isinstance(t.get('failure'),list) or len(t['failure']) != failed:
                    _issue(issues, 'development_inner_fold_counts', 'Inner tuning attempted/succeeded/failed folds are incomplete', **loc)
            valid_tuning = [t for t in tuning_records if isinstance(t,dict) and t.get('success_folds') == 3 and t.get('failed_folds') == 0 and _number(t.get('C_index'))]
            if valid_tuning:
                best=max(t['C_index'] for t in valid_tuning)
                shortlist=[t for t in valid_tuning if t['C_index'] >= best-.01]
                expected_choice=sorted(shortlist,key=lambda t:(t['kind']=='ridge',-t['penalty'],-t['C_index'],t['kind']))[0]
                if (development.get('kind'),development.get('penalty')) != (expected_choice['kind'],expected_choice['penalty']):
                    _issue(issues, 'development_tuning_selection', 'Fitted radiomics setting differs from locked training-only selection rule', **loc)
        models = record.get('models')
        if not isinstance(models, dict):
            _issue(issues, 'models_missing', 'Fold model records are missing', **loc)
            models = {}
        unknown = sorted(set(models) - set(applicable))
        if unknown:
            _issue(issues, 'unexpected_models', 'Models outside the applicable locked specification', models=unknown, **loc)
        for name in applicable:
            stat = models.get(name)
            ml = dict(loc, model=name)
            if not isinstance(stat, dict) or 'error' in stat:
                _issue(issues, 'model_missing_or_failed', 'Applicable model is missing or failed', error=stat.get('error') if isinstance(stat, dict) else None, **ml)
                continue
            counts = [stat.get(k) for k in ('concordant', 'discordant', 'tied_risk', 'comparable_pairs')]
            valid_stats = all(_integer(x) and x >= 0 for x in counts) and counts[-1] > 0 and sum(counts[:3]) == counts[-1]
            ci = stat.get('C_index')
            if not valid_stats or not _number(ci) or not 0 <= ci <= 1:
                _issue(issues, 'model_statistics', 'C-index or comparable-pair counts are incomplete/invalid', **ml)
            elif not math.isclose(ci, (counts[0] + .5 * counts[2]) / counts[3], rel_tol=1e-10, abs_tol=1e-10):
                _issue(issues, 'model_cindex_counts', 'C-index disagrees with concordance counts', **ml)
                valid_stats = False
            if valid_stats and _number(ci):
                statistics[name].append(stat)
            for key in ('risk', 'validation_events', 'validation_times'):
                values = stat.get(key)
                if not isinstance(values, list) or len(values) != len(val) or any(not _number(x) and not (key == 'validation_events' and isinstance(x, bool)) for x in values):
                    _issue(issues, 'model_predictions', 'Prediction/outcome %s is invalid or differs from fold validation length' % key, **ml)
            expected_events = [outcomes[x]['event'] for x in val]
            expected_times = [outcomes[x]['time'] for x in val]
            if stat.get('validation_events') != expected_events:
                _issue(issues, 'model_event_alignment', 'Model event values differ from matrix patients in fold order', **ml)
            times = stat.get('validation_times')
            if not isinstance(times, list) or len(times) != len(expected_times) or any(not _number(a) or not math.isclose(a, b, rel_tol=1e-12, abs_tol=1e-12) for a, b in zip(times or [], expected_times)):
                _issue(issues, 'model_time_alignment', 'Model survival times differ from matrix patients in fold order', **ml)
            _transform_issues(stat.get('transform'), plan['model_specs'][name], issues, ml)
            if name not in RSF_SPECS:
                coefficients = stat.get('coefficients')
                tr = stat.get('transform', {})
                expected_coefficients = sum(x is True for x in tr.get('keep', [])) if isinstance(tr, dict) else -1
                # The age spline replaces one standardized age column by the
                # spline basis; its actual coefficient count is intentionally
                # not inferred from the scalar transform.
                if not isinstance(coefficients, list) or any(not _number(x) for x in coefficients) or (name != 'clinical_age_spline' and len(coefficients) != expected_coefficients):
                    _issue(issues, 'model_coefficients', 'Fitted model coefficients are invalid', **ml)
        first, second = plan['primary_comparison']
        a, b = models.get(first, {}), models.get(second, {})
        if not isinstance(a, dict) or not isinstance(b, dict) or 'C_index' not in a or 'C_index' not in b or any(a.get(k) != b.get(k) for k in ('validation_events', 'validation_times')):
            _issue(issues, 'primary_pair_missing', 'Primary comparison lacks matched fold predictions and outcomes', **loc)
    if sample_counters and any(x != sample_counters[0] for x in sample_counters):
        _issue(issues, 'fold_resample_mismatch', 'Outer folds do not partition one common bootstrap/point sample', **location)
    if sample_counters and validation_counter != sample_counters[0]:
        _issue(issues, 'validation_fold_partition', 'Validation folds do not cover the run sample exactly once', **location)
    if sample_counters and run.get('unique_patients') != len(sample_counters[0]):
        _issue(issues, 'checkpoint_unique_population', 'unique_patients differs from actual run patient identities', **location)
    summaries = run.get('summary')
    if not isinstance(summaries, dict):
        _issue(issues, 'summary_missing', 'Run model summaries are missing', **location)
        summaries = {}
    unknown_summaries = set(summaries) - set(applicable) - set(RSF_SPECS)
    if unknown_summaries:
        _issue(issues, 'unexpected_summary_models', 'Summary has unplanned model names', models=sorted(unknown_summaries), **location)
    for name in applicable:
        summary = summaries.get(name)
        valid = statistics[name]
        if not isinstance(summary, dict) or len(valid) != 5 or summary.get('successful_folds') != 5 or summary.get('attempted_folds') != 5:
            _issue(issues, 'model_fold_summary', 'Applicable model does not have exactly five valid outer folds', model=name, **location)
            continue
        denominator = sum(x['comparable_pairs'] for x in valid)
        expected = sum(x['C_index'] * x['comparable_pairs'] for x in valid) / denominator
        if not _number(summary.get('C_index')) or not math.isclose(summary['C_index'], expected, rel_tol=1e-12, abs_tol=1e-12):
            _issue(issues, 'summary_recomputation', 'Run C-index differs from recomputed pair-weighted fold values', model=name, **location)
    for name in set(RSF_SPECS) - set(applicable):
        if name in summaries and (summaries[name].get('C_index') is not None or summaries[name].get('successful_folds') != 0):
            _issue(issues, 'inapplicable_rsf_summary', 'Inapplicable RSF summary must not claim successful validation', model=name, **location)
    if not issues:
        issues.extend(evidence_issues(run,plan))
    return issues


def audit_validation(revision, *, modes=('point', 'bootstrap'), require_complete=True, allow_legacy=False):
    """Audit actual checkpoints, original task events and locked fingerprints.

    Strict default requires all 510 planned tasks.  Historical evidence may be
    inspected with allow_legacy=True; unknown historical attempts are reported as
    null and can never yield a new COMPLETE claim.
    """
    root, out, _, _ = _paths(revision)
    report = dict(schema_version=1, audited_utc=_utc(), revision=str(root.resolve()),
                  status='FAILED', scope=list(modes), legacy_inspection=bool(allow_legacy),
                  input_fingerprints=_fingerprints(revision), issues=[], modes={},
                  checkpoint_count=0, valid_checkpoint_count=0, prior_failures=[])
    # Preserve observed event counts even when the input/config plan itself is
    # invalid; never infer attempted executions from checkpoint presence.
    raw_events=[]
    early_ledger=out/'VALIDATION_TASK_EVENTS.jsonl'
    if early_ledger.is_file():
        for line in early_ledger.read_text(encoding='utf-8').splitlines():
            try:
                event=json.loads(line)
                if isinstance(event,dict):raw_events.append(event)
            except (ValueError,TypeError):pass
    report['counts']={mode:{'planned':count,**{key:(sum(e.get('mode')==mode and e.get('event')==key for e in raw_events) if early_ledger.is_file() else None) for key in ['attempted','succeeded','failed','reused']}} for mode,count in LOCKED_COUNTS.items()}
    report['counts_interpretation']='Observed task events; checkpoint counts are separately audited after plan verification.'
    issues = report['issues']
    try:
        plan = ensure_validation_plan(revision)
    except Exception as exc:
        _issue(issues, 'validation_plan_error', str(exc), exception_type=type(exc).__name__)
        _publish_integrity_report(root,out,report)
        if require_complete:
            raise ValidationIntegrityError('Validation plan/input verification failed', report) from exc
        return report
    report.update(plan_sha256=plan['plan_sha256'], applicable_models=plan['applicable_models'],
                  plan={'model_names':plan['applicable_models'],'planned':plan['planned'],'outer_folds':5,'inner_folds':3},
                  matrix=plan['matrix'], original_seeds={m: plan['configuration']['seed_' + m] for m in LOCKED_COUNTS})
    if set(modes) != set(LOCKED_COUNTS) or len(modes) != 2:
        _issue(issues, 'incomplete_audit_scope', 'Global COMPLETE requires both point and bootstrap modes')
    requested = [m for m in modes if m in LOCKED_COUNTS]
    if len(requested) != len(modes):
        _issue(issues, 'invalid_audit_scope', 'Unknown validation mode in audit scope')
    ledger_path = out / 'VALIDATION_TASK_EVENTS.jsonl'
    ledger = []
    if ledger_path.exists():
        for line_no, line in enumerate(ledger_path.read_text(encoding='utf-8').splitlines(), 1):
            try:
                value = json.loads(line)
                if not isinstance(value, dict):
                    raise ValueError('Task event must be an object')
                ledger.append(value)
            except Exception as exc:
                _issue(issues, 'task_ledger_invalid', str(exc), line=line_no)
    elif not allow_legacy:
        _issue(issues, 'task_ledger_missing', 'No actual task-attempt ledger exists; checkpoints do not establish attempts')
    attempts, terminal, attempt_owners = {}, {}, {}
    for index, entry in enumerate(ledger, 1):
        mode, run_id, attempt = entry.get('mode'), entry.get('run_id'), entry.get('attempt_id')
        loc = dict(mode=mode, run_id=run_id, attempt_id=attempt, ledger_line=index)
        if mode not in LOCKED_COUNTS or not _integer(run_id) or not 0 <= run_id < LOCKED_COUNTS.get(mode, 0) or not isinstance(attempt, str) or not attempt:
            _issue(issues, 'task_ledger_identity', 'Invalid task/attempt identity', **loc)
            continue
        key = (mode, run_id, attempt)
        if attempt in attempt_owners and attempt_owners[attempt] != (mode,run_id):
            _issue(issues,'task_attempt_identity_reused','Attempt ID is shared by different planned task identities',**loc)
        attempt_owners[attempt]=(mode,run_id)
        if entry.get('event') == 'failed':
            report['prior_failures'].append(entry)
        if entry.get('event') not in ('attempted', 'succeeded', 'failed', 'reused'):
            _issue(issues, 'task_ledger_event', 'Invalid task event', **loc)
            continue
        if entry.get('input_fingerprints') != plan['inputs'] or any(entry.get('observed_input_fingerprints', {}).get(k) != plan['inputs'][k] for k in INPUT_KEYS) or entry.get('plan_sha256') != plan['plan_sha256']:
            _issue(issues, 'task_ledger_fingerprint', 'Task attempt used or observed a different input/configuration', **loc)
        expected_seed = plan['configuration']['seed_' + mode] + run_id
        if entry.get('seed') != expected_seed or entry.get('original_seed') != expected_seed:
            _issue(issues, 'task_ledger_seed', 'Task event does not preserve the original protocol seed', **loc)
        if entry['event'] == 'attempted':
            if key in attempts:
                _issue(issues, 'task_ledger_duplicate_attempt', 'Attempt identifier was submitted more than once', **loc)
            attempts[key] = entry
        else:
            if key not in attempts:
                _issue(issues, 'task_terminal_without_attempt', 'Terminal task event has no preceding attempted event', **loc)
            if key in terminal:
                _issue(issues, 'task_duplicate_terminal', 'Attempt has multiple terminal events', **loc)
            if entry['event'] == 'failed' and not entry.get('exception'):
                _issue(issues, 'task_failure_exception_missing', 'Failed task lacks an exception record', **loc)
            terminal[key] = entry
    paths = sorted((out / 'checkpoints').glob('*')) if (out / 'checkpoints').exists() else []
    actual, identities = {}, {}
    report['checkpoint_manifest']=[]
    for path in paths:
        if not path.is_file():
            continue
        report['checkpoint_count'] += 1
        try:
            with gzip.open(path, 'rt', encoding='utf-8') as handle:
                checkpoint = json.load(handle)
            cp_issues = checkpoint_issues(checkpoint, plan, path)
        except Exception as exc:
            _issue(issues, 'checkpoint_unreadable', str(exc), checkpoint=str(path))
            continue
        issues.extend(cp_issues)
        report['checkpoint_manifest'].append({'file':path.name,'sha256':_sha(path),
            'mode':checkpoint.get('mode') if isinstance(checkpoint,dict) else None,
            'run_id':checkpoint.get('run_id') if isinstance(checkpoint,dict) else None,
            'valid':not cp_issues})
        if isinstance(checkpoint, dict):
            identity = (checkpoint.get('mode'), checkpoint.get('run_id'))
            try:
                if identity in identities:
                    _issue(issues, 'duplicate_checkpoint_identity', 'Multiple files claim the same validation task', checkpoint=str(path), other_checkpoint=identities[identity])
                identities[identity] = str(path)
                if not cp_issues:
                    actual[identity] = checkpoint
            except TypeError:
                _issue(issues, 'checkpoint_identity_unhashable', 'Checkpoint identity is not a scalar', checkpoint=str(path))
    report['valid_checkpoint_count'] = len(actual)
    for mode in requested:
        planned = LOCKED_COUNTS[mode]
        mode_attempts = [key for key in attempts if key[0] == mode]
        events = [value for key, value in terminal.items() if key[0] == mode]
        successes = [x for x in events if x['event'] == 'succeeded']
        failures = [x for x in events if x['event'] == 'failed']
        reused = [x for x in events if x['event'] == 'reused']
        valid_ids = sorted(run_id for m, run_id in actual if m == mode)
        missing = sorted(set(range(planned)) - set(valid_ids))
        latest = {}
        for entry in events:
            latest[entry['run_id']] = entry
        incomplete_terminal = sorted(i for i in range(planned) if latest.get(i, {}).get('event') not in ('succeeded', 'reused'))
        open_attempts = [list(key) for key in mode_attempts if key not in terminal]
        if missing:
            _issue(issues, 'checkpoints_missing_or_invalid', 'Planned task checkpoints are missing or fail integrity', mode=mode, run_ids=missing)
        if incomplete_terminal and not allow_legacy:
            _issue(issues, 'task_success_missing', 'No successful/reused terminal attempt exists for all planned tasks', mode=mode, run_ids=incomplete_terminal)
        if open_attempts:
            _issue(issues, 'task_attempt_unfinished', 'Submitted attempts have no terminal record', mode=mode, attempts=open_attempts)
        primary_pairs = sum(1 for key in actual if key[0] == mode)
        report['modes'][mode] = dict(
            planned=planned, attempted=len(mode_attempts) if ledger_path.exists() else None,
            succeeded=len(successes) if ledger_path.exists() else None,
            failed=len(failures) if ledger_path.exists() else None,
            reused=len(reused) if ledger_path.exists() else None,
            unique_tasks_attempted=len({key[1] for key in mode_attempts}) if ledger_path.exists() else None,
            actual_checkpoint_files=sum(1 for key in identities if key[0] == mode),
            checkpoint_success_counts=len(valid_ids), completed_run_ids=valid_ids,
            missing_or_invalid_run_ids=missing, primary_paired_runs=primary_pairs,
            primary_paired_folds=primary_pairs * 5,
        )
    # A failed last attempt cannot be hidden by a previously valid checkpoint.
    # Prior failed attempts may remain in a successfully retried task's audit.
    if allow_legacy:
        _issue(issues, 'legacy_attempts_unverifiable', 'Historical inspection does not certify newly attempted execution; no COMPLETE claim is issued')
    report['prior_failure_count'] = len(report['prior_failures'])
    report['counts'] = report['modes']
    first,second=plan['primary_comparison']
    report['primary_paired_results']=[{'mode':mode,'run_id':run_id,'seed':cp['seed'],
        'baseline':first,'combined':second,
        'baseline_C_index':cp['summary'][first]['C_index'],
        'combined_C_index':cp['summary'][second]['C_index'],
        'paired_delta_C_index':cp['summary'][second]['C_index']-cp['summary'][first]['C_index'],
        'paired_outer_folds':5}
        for (mode,run_id),cp in sorted(actual.items())]
    report['task_ledger_sha256']=_sha(ledger_path) if ledger_path.is_file() else None
    report['status'] = 'COMPLETE' if not issues else ('PARTIAL' if actual else 'FAILED')
    report['complete'] = report['status'] == 'COMPLETE'
    report['claims'] = {
        'point_repeats_completed': report['modes'].get('point', {}).get('checkpoint_success_counts', 0),
        'bootstrap_replicates_completed': report['modes'].get('bootstrap', {}).get('checkpoint_success_counts', 0),
        'all_510_planned_tasks_complete': report['status'] == 'COMPLETE',
    }
    _publish_integrity_report(root,out,report)
    if require_complete and report['status'] != 'COMPLETE':
        raise ValidationIntegrityError('Validation incomplete: %s (%d integrity issue groups); see %s' % (report['status'], len(issues), out / 'VALIDATION_INTEGRITY.json'), report)
    return report


def audit_validation_reference(source_revision, evidence_directory):
    """Inspect an immutable historical reference in a separate evidence snapshot.

    Never writes to the source revision or copies its checkpoints into a new
    scientific run's 04_validation directory. Missing historical attempt records
    remain missing; successful files do not manufacture attempted-task counts.
    """
    source=Path(source_revision)
    destination=Path(evidence_directory)
    if destination.exists() and any(destination.iterdir()):
        raise ValidationIntegrityError('Historical evidence destination must be new and empty')
    destination.mkdir(parents=True,exist_ok=True)
    snapshot=destination/'reference_revision'
    manifest=[]
    names=['00_protocol_and_provenance/revision_config.json','03_model_development/analysis_matrix.csv',
           '04_validation/VALIDATION_PLAN.json','04_validation/VALIDATION_TASK_EVENTS.jsonl',
           '04_validation/VALIDATION_TASK_FAILURES.jsonl']
    names += [str(p.relative_to(source)) for p in sorted((source/'04_validation/checkpoints').glob('*')) if p.is_file()]
    for name in names:
        original=source/name
        if not original.is_file():continue
        target=snapshot/name;target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(original,target)
        digest=_sha(original)
        if _sha(target)!=digest:raise ValidationIntegrityError('Historical evidence copy/hash mismatch')
        manifest.append({'source_path':str(original),'snapshot_path':str(target),'sha256':digest})
    _write_json(destination/'REFERENCE_PROVENANCE.json',{'scope':'HISTORICAL_REFERENCE_ONLY_NOT_NEW_MODEL_EXECUTION',
        'new_model_runs':0,'files':manifest})
    report=audit_validation(snapshot,require_complete=False)
    if report['status']!='COMPLETE':
        raise ValidationIntegrityError('Historical validation cannot be certified for reuse under the full protocol',report)
    return report


def verify_paired_summary(revision, integrity=None):
    """Check formal paired results against the actual complete run checkpoints."""
    root,out,_,_=_paths(revision)
    integrity=integrity if integrity is not None else audit_validation(revision)
    if integrity.get('status')!='COMPLETE':
        raise ValidationIntegrityError('Cannot verify formal pairs from incomplete validation',integrity)
    path=out/'paired_incremental_performance.csv'
    issues=[];values=[]
    try:
        with path.open(encoding='utf-8-sig',newline='') as f:rows=list(csv.DictReader(f))
    except (OSError,ValueError) as exc:
        rows=[];_issue(issues,'paired_summary_missing',str(exc))
    expected={(base,base+'_radiomics') for base in PAIRED_BASELINES}
    found=[(r.get('baseline'),r.get('combined')) for r in rows]
    if len(found)!=len(expected) or set(found)!=expected:
        _issue(issues,'paired_summary_specifications','Formal paired comparisons are missing, duplicated or unexpected')
    deltas={pair:{'point':{},'bootstrap':{}} for pair in expected}
    for entry in integrity['checkpoint_manifest']:
        if not entry['valid']:continue
        cp_path=out/'checkpoints'/entry['file']
        if _sha(cp_path)!=entry['sha256']:
            _issue(issues,'paired_checkpoint_changed','Checkpoint changed after validation audit',checkpoint=str(cp_path))
            continue
        with gzip.open(cp_path,'rt',encoding='utf-8') as f:cp=json.load(f)
        for base,combined in expected:
            deltas[(base,combined)][cp['mode']][cp['run_id']]=cp['summary'][combined]['C_index']-cp['summary'][base]['C_index']
    def quantile(v,q):
        ordered=sorted(v);position=(len(ordered)-1)*q;low=int(math.floor(position));high=int(math.ceil(position))
        return ordered[low]+(ordered[high]-ordered[low])*(position-low)
    for row in rows:
        pair=(row.get('baseline'),row.get('combined'))
        if pair not in deltas:continue
        point=deltas[pair]['point'];boot=deltas[pair]['bootstrap']
        if len(point)!=10 or len(boot)!=500:
            _issue(issues,'paired_checkpoint_coverage','Paired checkpoint results are incomplete',baseline=pair[0]);continue
        expected_values={'paired_delta_C_index':point[0],
            'repeat_mean_delta':statistics.mean(point.values()),'repeat_SD_delta':statistics.stdev(point.values()),
            'CI_lower':quantile(list(boot.values()),.025),'CI_upper':quantile(list(boot.values()),.975),
            'bootstrap_joint_complete':len(boot),'point_joint_complete':len(point)}
        for key,value in expected_values.items():
            try:actual=float(row.get(key,''))
            except (TypeError,ValueError):actual=float('nan')
            if not math.isfinite(actual) or not math.isclose(actual,value,rel_tol=1e-10,abs_tol=1e-12):
                _issue(issues,'paired_summary_value','Formal paired value differs from actual checkpoints',baseline=pair[0],field=key,expected=value,actual=row.get(key))
        crossing=expected_values['CI_lower']<=0<=expected_values['CI_upper']
        if str(row.get('ci_crosses_zero')).lower()!=str(crossing).lower():
            _issue(issues,'paired_summary_interval','Paired CI crossing flag differs from actual bootstrap values',baseline=pair[0])
        values.append({'baseline':pair[0],'combined':pair[1],**expected_values})
    report={'status':'FAILED' if issues else 'COMPLETE','scope':'FORMAL_PAIRED_SUMMARY',
            'issues':issues,'recomputed_from_actual_checkpoints':values,
            'paired_csv_sha256':_sha(path) if path.is_file() else None,
            'validation_plan_sha256':integrity.get('plan_sha256')}
    _write_json(out/'PAIRED_SUMMARY_INTEGRITY.json',report)
    if issues:
        _publish_integrity_report(root,out,{**integrity,'status':'PARTIAL','complete':False,
            'issues':integrity.get('issues',[])+issues,'formal_paired_summary_complete':False})
        raise ValidationIntegrityError('Formal paired comparison results are incomplete or inconsistent',report)
    return report
