"""Record fitted survival objects without refitting or altering model predictions.

The right-continuous Cox baseline is zero before its first observed time.
No extrapolation is reported beyond the last training follow-up time. The
associated probability evidence is not itself a censoring-adjusted metric.
"""
import hashlib
import json
from pathlib import Path
import numpy as np

SCHEMA_VERSION = 2
PRIMARY_MODELS = ('clinical_diameter_mean_hu', 'clinical_diameter_mean_hu_radiomics')
MONTHS = list(range(1, 37))


def source_fingerprint(directory):
    files = {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
             for p in sorted(Path(directory).glob('*.py'))}
    return hashlib.sha256(json.dumps(files, sort_keys=True).encode()).hexdigest()


def fitted_baseline(model, training_times):
    """Save the baseline already fitted by sksurv, including Coxnet centering."""
    if hasattr(model, 'cum_baseline_hazard_'):
        base = model
    elif hasattr(model, '_get_baseline_model'):
        # sksurv 0.23.1 stores the actual fitted alpha-specific Breslow object.
        base = model._get_baseline_model(None)
    else:
        # The original pipeline's no-varying-predictor fallback is training KM.
        curve = model.predict_survival_function(np.empty((1, 0)))[0]
        return dict(kind='training_KM_null_model', times=curve.x.tolist(),
                    baseline_survival=curve.y.tolist(),
                    before_first_time_survival=1.0,
                    max_training_followup=float(np.max(training_times)),
                    prediction_offset=0.0)
    hazard = base.cum_baseline_hazard_
    survival = base.baseline_survival_
    offset = float(np.asarray(getattr(model, 'offset_', [0.0])).ravel()[-1])
    return dict(kind='fitted_Breslow', times=hazard.x.tolist(),
                cumulative_hazard=hazard.y.tolist(),
                baseline_survival=survival.y.tolist(),
                before_first_time_cumulative_hazard=0.0,
                max_training_followup=float(np.max(training_times)),
                prediction_offset=offset,
                risk_definition='model.predict(X); Coxnet subtracts stored prediction_offset')


def monthly_survival(baseline, linear_predictor):
    grid = np.asarray(MONTHS, dtype=float)
    times = np.asarray(baseline['times'], dtype=float)
    idx = np.searchsorted(times, grid, side='right') - 1
    if baseline['kind'] == 'fitted_Breslow':
        values = np.asarray(baseline['cumulative_hazard'], dtype=float)
        h0 = np.zeros(len(grid)); keep = idx >= 0; h0[keep] = values[idx[keep]]
        probabilities = np.exp(-np.outer(np.exp(np.asarray(linear_predictor, dtype=float)), h0))
    else:
        values = np.asarray(baseline['baseline_survival'], dtype=float)
        surv = np.ones(len(grid)); keep = idx >= 0; surv[keep] = values[idx[keep]]
        probabilities = np.tile(surv, (len(linear_predictor), 1))
    supported = grid <= baseline['max_training_followup']
    if (not np.isfinite(probabilities).all() or (probabilities < 0).any()
            or (probabilities > 1).any() or (np.diff(probabilities, axis=1) > 1e-12).any()):
        raise ValueError('Invalid validation survival probability evidence')
    return dict(months=MONTHS, supported_by_training_followup=supported.tolist(),
                unsupported_reason=[None if ok else 'beyond_max_training_followup' for ok in supported],
                survival_probability=[[float(v) if supported[j] else None
                                       for j, v in enumerate(row)] for row in probabilities],
                row_order='fold.validation_ids; includes bootstrap multiplicities',
                estimator='exp(-fitted_Breslow_H0(t) * exp(model.predict(X))); right-continuous; H0=0 before first time',
                extrapolation='None beyond maximum training follow-up')


def evidence_issues(run, plan):
    """Additional fail-closed evidence gate; original scientific audit remains."""
    problems = []
    def bad(code, **where):
        problems.append(dict(code='survival_evidence_' + code, **where))
    if run.get('prediction_evidence_schema_version') != SCHEMA_VERSION:
        bad('missing_or_wrong_schema'); return problems
    if run.get('execution_source_sha256') != plan.get('execution_source_sha256'):
        bad('execution_source_mismatch')
    n = plan['matrix']['n']; ids = plan['matrix']['patient_ids']
    draw = run.get('draw_source_row_indices')
    if not isinstance(draw, list) or len(draw) != n or any(type(x) is not int or not 0 <= x < n for x in draw):
        bad('draw_indices'); return problems
    expected = (np.random.default_rng(run['seed']).integers(0, n, n).tolist()
                if run['mode'] == 'bootstrap' else list(range(n)))
    if draw != expected: bad('seed_draw_mismatch')
    draw_ids = [ids[i] for i in draw]
    if run.get('draw_patient_ids') != draw_ids: bad('draw_patient_id_alignment')
    from sklearn.model_selection import KFold
    unique = np.asarray(sorted(set(draw_ids)))
    expected_splits = []
    for tr, va in KFold(plan['outer_folds'], shuffle=True, random_state=run['seed']).split(unique):
        tr_ids, va_ids = set(unique[tr]), set(unique[va])
        expected_splits.append(([i for i, pid in enumerate(draw_ids) if pid in tr_ids],
                                [i for i, pid in enumerate(draw_ids) if pid in va_ids]))
    for fold in run.get('folds', []):
        fid = fold.get('fold')
        if type(fid) is not int or not 0 <= fid < len(expected_splits): continue
        tr, va = expected_splits[fid]
        if fold.get('train_draw_row_indices') != tr or fold.get('validation_draw_row_indices') != va:
            bad('fold_draw_alignment', fold=fid)
        if fold.get('train_ids') != [draw_ids[i] for i in tr] or fold.get('validation_ids') != [draw_ids[i] for i in va]:
            bad('fold_patient_alignment', fold=fid)
        if fold.get('training_events') != [plan['matrix']['outcomes'][draw_ids[i]]['event'] for i in tr]:
            bad('training_events', fold=fid)
        # CSV float() and pandas' parser can differ by one binary64 ULP.
        # Preserve exact row/ID/event checks; permit only an absolute 1e-12-month
        # parser tolerance, not changes to the stored clinical outcome data.
        expected_times = np.asarray([plan['matrix']['outcomes'][draw_ids[i]]['time'] for i in tr], dtype=float)
        actual_times = np.asarray(fold.get('training_times', []), dtype=float)
        if (actual_times.shape != expected_times.shape
                or not np.isfinite(actual_times).all()
                or not np.allclose(actual_times, expected_times, rtol=0.0, atol=1e-12)):
            bad('training_times', fold=fid)
        development = fold.get('radiomics_development', {})
        if 'fitted_baseline' not in development: bad('radiomics_baseline_missing', fold=fid)
        inner = fold.get('inner_splits')
        if not isinstance(inner, list) or len(inner) != plan['inner_folds']:
            bad('inner_splits_missing', fold=fid)
        else:
            outer_train_ids = [draw_ids[i] for i in tr]
            inner_unique = np.asarray(sorted(set(outer_train_ids)))
            for k, (ii, jj) in enumerate(KFold(plan['inner_folds'], shuffle=True,
                    random_state=run['seed']+10000+fid).split(inner_unique)):
                train_set, valid_set = set(inner_unique[ii]), set(inner_unique[jj])
                expected_inner = dict(fold=k, seed=run['seed']+10000+fid,
                    train_outer_training_row_indices=[i for i,p in enumerate(outer_train_ids) if p in train_set],
                    validation_outer_training_row_indices=[i for i,p in enumerate(outer_train_ids) if p in valid_set],
                    train_ids=[p for p in outer_train_ids if p in train_set],
                    validation_ids=[p for p in outer_train_ids if p in valid_set])
                if inner[k] != expected_inner: bad('inner_split_alignment', fold=fid, inner_fold=k)
        for name, result in fold.get('models', {}).items():
            if name.endswith('_RSF'): continue
            base = result.get('fitted_baseline')
            if not isinstance(base, dict): bad('baseline_missing', fold=fid, model=name); continue
            if base.get('kind') not in ['fitted_Breslow', 'training_KM_null_model']:
                bad('baseline_kind', fold=fid, model=name); continue
            times = np.asarray(base.get('times', []), float)
            surv = np.asarray(base.get('baseline_survival', []), float)
            if (len(times) == 0 or len(times) != len(surv) or not np.isfinite(times).all()
                    or not np.isfinite(surv).all() or (np.diff(times) <= 0).any()
                    or (surv < 0).any() or (surv > 1).any() or (np.diff(surv) > 1e-12).any()):
                bad('baseline_values', fold=fid, model=name); continue
            if base.get('max_training_followup') != max(fold['training_times']):
                bad('baseline_training_support', fold=fid, model=name)
            if base.get('kind') == 'fitted_Breslow':
                h = np.asarray(base.get('cumulative_hazard', []), float)
                if (len(h) != len(times) or not np.isfinite(h).all() or (h < 0).any()
                        or (np.diff(h) < -1e-12).any() or not np.allclose(surv, np.exp(-h), rtol=1e-12, atol=1e-12)):
                    bad('baseline_hazard', fold=fid, model=name); continue
            if name in PRIMARY_MODELS:
                actual = result.get('monthly_survival')
                expected_prediction = monthly_survival(base, result['risk'])
                if actual != expected_prediction: bad('monthly_prediction_mismatch', fold=fid, model=name)
    return problems
