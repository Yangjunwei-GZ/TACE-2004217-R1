"""Synthetic validation evidence only: no image extraction or model fitting.

The saved predictions are deliberately constructed fixtures. They must never be
used as scientific results. Membership, seeds, counters and survival outcomes
follow the real execution schema so integrity gates can be exercised end to end.
"""
from pathlib import Path
import collections
import csv
import gzip
import hashlib
import json

import numpy as np
from sklearn.model_selection import KFold


DEFAULT_CONFIG = {
    'revision_date': '2026-10-14', 'primary_phase': 1,
    'geometry_numeric_tolerance': 0.0001, 'features_original': 107,
    'correlation_threshold': 0.85, 'max_filtered_features': 12,
    'inner_folds': 3, 'outer_folds': 5, 'point_repeats': 10,
    'bootstrap_replicates': 500, 'seed_point': 20261014,
    'seed_bootstrap': 20261015, 'baseline_penalizer': 0.1,
    'candidates': ['ridge', 'top_two_cox', 'top_two_elastic_net'],
    'ridge_penalties': [0.1, 1.0, 10.0],
    'two_feature_penalties': [0.01, 0.1, 1.0],
    'primary_comparison': ['clinical_diameter_mean_hu',
                           'clinical_diameter_mean_hu_radiomics'],
    'prediction_months': [12, 24, 36], 'ibs_monthly_grid': [1, 36],
    'independent_radiologist_review_required': True,
    'public_code_required': True, 'no_hr_eligibility_rule': True,
    'no_dca': True,
}


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_csv(path, rows):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with path.open('w', encoding='utf-8-sig', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def load_checkpoint(path):
    with gzip.open(path, 'rt', encoding='utf-8') as handle:
        return json.load(handle)


def save_checkpoint(path, checkpoint):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, 'wt', encoding='utf-8', compresslevel=1) as handle:
        json.dump(checkpoint, handle, allow_nan=False, separators=(',', ':'))


def _transform(columns):
    return {'columns': list(columns), 'median': [0.0] * len(columns),
            'mean': [0.0] * len(columns), 'sd': [1.0] * len(columns),
            'keep': [True] * len(columns)}


def _stat(rows):
    """Exact survival comparable-pair counts for a perfect negative-time risk."""
    concordant = discordant = tied_risk = 0
    for index, first in enumerate(rows):
        for second in rows[index + 1:]:
            ta, tb = first['os_months'], second['os_months']
            ea, eb = bool(first['death_event']), bool(second['death_event'])
            if ta == tb:
                # An event and censoring at the same time are comparable.
                if ea != eb:
                    tied_risk += 1
            elif (ta < tb and ea) or (tb < ta and eb):
                concordant += 1
    total = concordant + discordant + tied_risk
    if not total:
        raise ValueError('Synthetic fold unexpectedly has no comparable pairs')
    return {'C_index': (concordant + tied_risk * 0.5) / total,
            'concordant': concordant, 'discordant': discordant,
            'tied_risk': tied_risk, 'comparable_pairs': total}


def create_complete_validation_fixture(revision, *, n_patients=25,
                                       n_events=None,
                                       modes=('point', 'bootstrap'),
                                       run_ids=None, write_ledger=True):
    """Write real-format synthetic checkpoints; return plan, paths and input audit.

    ``run_ids`` may map modes to explicit run IDs for a deliberately partial
    fixture. Protocol counts remain locked to 10/500, including partial cases.
    """
    from validation_integrity import ensure_validation_plan, record_task_event
    revision = Path(revision)
    if n_events is None:
        n_events = n_patients
    config_path = revision / '00_protocol_and_provenance/revision_config.json'
    config_path.parent.mkdir(parents=True, exist_ok=True)
    config_path.write_text(json.dumps(DEFAULT_CONFIG, indent=2), encoding='utf-8')
    rows = []
    clinical = ['age', 'sex_female', 'child_b', 'log_afp', 'log_tumour_count',
                'bclc_b', 'selected_diameter_mm', 'log_selected_volume',
                'log_total_volume', 'mean_hu', 'median_hu', 'p10_hu',
                'clinical_max_axial_diameter_mm', 'afp_high_400',
                'tumour_multiple', 'radiomics_score']
    for index in range(n_patients):
        row = {'patient_id': 'SYN%03d' % index,
               'sample_id': 'SYN%03d_1_0' % index,
               'death_event': int(index < n_events),
               'os_months': float(index + 10)}
        row.update({column: float(index + 1) for column in clinical})
        row.update({'original_fixture_%03d' % feature: float(index + feature)
                    for feature in range(107)})
        rows.append(row)
    matrix_path = revision / '03_model_development/analysis_matrix.csv'
    write_csv(matrix_path, rows)
    feature_path = revision / '02_extraction/arterial_lesion_features.csv'
    write_csv(feature_path, rows)
    write_csv(matrix_path.parent / 'FINAL_REVISION_PATIENT_SCORES.csv',
              [dict(row, synthetic_test_score=0.5) for row in rows])
    plan = ensure_validation_plan(revision)
    paths = []
    for mode in modes:
        count = DEFAULT_CONFIG['point_repeats' if mode == 'point'
                               else 'bootstrap_replicates']
        identifiers = list(range(count)) if run_ids is None else list(run_ids.get(mode, []))
        for run_id in identifiers:
            seed = DEFAULT_CONFIG['seed_point' if mode == 'point' else 'seed_bootstrap'] + run_id
            if mode == 'bootstrap':
                selected = np.random.default_rng(seed).integers(0, n_patients, n_patients)
                run_rows = [rows[int(index)] for index in selected]
            else:
                run_rows = list(rows)
            unique = np.array(sorted({row['patient_id'] for row in run_rows}))
            checkpoint = {'run_id': run_id, 'mode': mode, 'seed': seed,
                          'synthetic_fixture':True,'model_fits_performed':0,
                          'n': len(run_rows), 'unique_patients': len(unique),
                          'folds': [], 'failures': [],
                          **plan['inputs']}
            for fold, (train_indices, val_indices) in enumerate(KFold(
                    5, shuffle=True, random_state=seed).split(unique)):
                train_ids = set(unique[train_indices])
                val_ids = set(unique[val_indices])
                train = [row for row in run_rows if row['patient_id'] in train_ids]
                validation = [row for row in run_rows if row['patient_id'] in val_ids]
                statistic = _stat(validation)
                models = {}
                for name in plan['applicable_models']:
                    columns = plan['model_specs'][name]
                    model = dict(statistic,
                                 risk=[-row['os_months'] for row in validation],
                                 validation_events=[row['death_event'] for row in validation],
                                 validation_times=[row['os_months'] for row in validation],
                                 transform=_transform(columns))
                    if 'RSF' not in name:
                        model['coefficients'] = [0.0] * len(columns)
                    models[name] = model
                development = {
                    'kind': 'ridge', 'penalty': 0.1,
                    'features': ['original_fixture_000'],
                    'coefficients': [1.0],
                    'transform': _transform(['original_fixture_000']),
                    'tuning': [{'kind': kind, 'penalty': penalty,
                                'C_index': 0.7 if kind == 'ridge' and penalty == 0.1 else 0.6,
                                'success_folds': 3, 'failed_folds': 0, 'failure': []}
                               for kind in DEFAULT_CONFIG['candidates'] for penalty in
                               (DEFAULT_CONFIG['ridge_penalties'] if kind == 'ridge' else DEFAULT_CONFIG['two_feature_penalties'])],
                }
                checkpoint['folds'].append({
                    'fold': fold, 'train_ids': [row['patient_id'] for row in train],
                    'validation_ids': [row['patient_id'] for row in validation],
                    'models': models, 'radiomics_development': development})
            summary = {}
            for name in plan['applicable_models']:
                statistics = [fold['models'][name] for fold in checkpoint['folds']]
                total = sum(statistic['comparable_pairs'] for statistic in statistics)
                summary[name] = {
                    'C_index': sum(statistic['C_index'] * statistic['comparable_pairs']
                                   for statistic in statistics) / total,
                    'successful_folds': 5, 'attempted_folds': 5}
            checkpoint['summary'] = summary
            checkpoint['elapsed_seconds'] = 0.0
            path = revision / '04_validation/checkpoints' / ('%s_%04d.json.gz' % (mode, run_id))
            save_checkpoint(path, checkpoint)
            paths.append(path)
            if write_ledger:
                attempt_id = 'synthetic:%s:%d' % (mode, run_id)
                record_task_event(revision, mode, run_id, 'attempted', attempt_id=attempt_id,
                                  details={'synthetic_fixture': True})
                record_task_event(revision, mode, run_id, 'succeeded', attempt_id=attempt_id,
                                  details={'synthetic_fixture': True, 'checkpoint_path': str(path)})
    audit = {'status': 'PASS', 'analysis_matrix': str(matrix_path),
             'data_scope':'SYNTHETIC_TEST_FIXTURE_ONLY',
             'analysis_matrix_sha256': digest(matrix_path),
             'feature_table': str(feature_path),
             'feature_table_sha256': digest(feature_path),
             'actual_patient_ids': [row['patient_id'] for row in rows],
             'actual_patient_count': len(rows),
             'actual_active_unique_mask_count': len(rows)}
    if tuple(modes)==('point','bootstrap') and run_ids is None:
        paired=[{'baseline':base,'combined':base+'_radiomics','paired_delta_C_index':0.0,
                 'repeat_mean_delta':0.0,'repeat_SD_delta':0.0,'CI_lower':0.0,'CI_upper':0.0,
                 'bootstrap_joint_complete':500,'point_joint_complete':10,'ci_crosses_zero':True}
                for base in ['clinical','clinical_diameter','clinical_selected_volume','clinical_total_volume','clinical_diameter_mean_hu','clinical_diameter_p10_hu']]
        write_csv(revision/'04_validation/paired_incremental_performance.csv',paired)
    return {'revision': revision, 'config': DEFAULT_CONFIG.copy(), 'plan': plan,
            'matrix_path': matrix_path, 'feature_path': feature_path,
            'checkpoint_paths': paths, 'input_audit': audit,
            'expected_model_names': plan['applicable_models']}
