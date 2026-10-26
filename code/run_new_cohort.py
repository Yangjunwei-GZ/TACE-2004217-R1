"""Explicit, fail-closed new-cohort entry; defaults to preflight, never trains implicitly."""
from pathlib import Path
import argparse
import csv
import datetime
import hashlib
import json
import os
import shutil
import subprocess
import sys
import traceback
sys.path.insert(0, str(Path(__file__).resolve().parent))
from analysis_input_guard import verify_analysis_inputs, assert_inputs_unchanged, verify_model_output
from survival_evidence import source_fingerprint

CODE = Path(__file__).resolve().parent
REQUIRED_ROI_DECISIONS = ['83_1_' + str(i) for i in range(5)] + ['94_1_' + str(i) for i in range(3)]


def utc(): return datetime.datetime.now(datetime.timezone.utc).isoformat()
def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def read_json(path): return json.loads(Path(path).read_text(encoding='utf-8-sig'))
def write_json(path, value):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False), encoding='utf-8')


def decision_path(run_root):
    return Path(run_root) / 'PUBLIC_INPUT_ADOPTION.csv'


def decision_issues(run_root):
    """Only explicit sourced final actions count. Acceptability/phase alone cannot substitute."""
    path = decision_path(run_root)
    if not path.is_file():
        return [dict(type='missing_medical_decision_file', path=str(path), sample_ids=REQUIRED_ROI_DECISIONS)]
    with path.open(encoding='utf-8-sig', newline='') as handle:
        rows = list(csv.DictReader(handle))
    problems = []
    for sid in REQUIRED_ROI_DECISIONS:
        found = [r for r in rows if r.get('sample_id') == sid]
        if len(found) != 1:
            problems.append(dict(type='missing_or_duplicate_final_medical_decision', sample_id=sid, rows=len(found))); continue
        row = found[0]; action = row.get('action', '')
        if action not in ['retain', 'replace', 'voxel_edit', 'alias', 'exclude']:
            problems.append(dict(type='missing_explicit_final_action', sample_id=sid, recorded_action=action))
        if not row.get('decision_source') or row.get('decision_source') in ['未提供', 'unknown', 'pending']:
            problems.append(dict(type='missing_final_decision_source', sample_id=sid))
        if action in ['replace', 'voxel_edit']:
            mask = row.get('final_mask_path', '')
            if not mask or not Path(mask).is_file():
                problems.append(dict(type='missing_final_modified_mask', sample_id=sid, path=mask))
        if action == 'alias' and (not row.get('alias_of') or row['alias_of'] == sid):
            problems.append(dict(type='missing_valid_alias_target', sample_id=sid))
    return problems


def preflight(run_root, *, verify_image_hashes=True):
    run_root = Path(run_root); rev = run_root / 'revision_R1_20261014'
    # Revoke a stale PASS before touching any input. An unexpected read/hash
    # exception must not leave an older successful preflight usable downstream.
    write_json(run_root / 'MODEL_PREFLIGHT.json', dict(status='CHECKING', checked_utc=utc(),
               model_training_authorized=False, issues=[], scientific_model_runs_started=0))
    issues = decision_issues(run_root)
    required = ['03_model_development/analysis_matrix.csv', '02_extraction/arterial_lesion_features.csv',
                '02_extraction/extraction_failures.csv', '01_data_audit/source_clinical_233.csv',
                '01_data_audit/phase_pairing_audit.csv', '01_data_audit/phase_analysis_eligibility.csv',
                '00_protocol_and_provenance/revision_config.json']
    for rel in required:
        if not (rev / rel).is_file(): issues.append(dict(type='missing_required_input', path=str(rev / rel)))
    for name in ['FINAL_CANDIDATE_DISPOSITION.csv', 'FINAL_IMAGE_AND_ALIAS_ROSTER.csv',
                 'FINAL_PATIENT_AND_REPRESENTATIVE_ROSTER.csv']:
        if not (run_root / name).is_file(): issues.append(dict(type='missing_required_roster', path=str(run_root / name)))
    cfg = rev / '00_protocol_and_provenance/revision_config.json'
    if cfg.is_file() and sha(cfg) != sha(CODE / 'protocol/revision_config.json'):
        issues.append(dict(type='locked_configuration_bytes_changed', path=str(cfg)))
    if issues:
        report = dict(status='BLOCKED', checked_utc=utc(), model_training_authorized=False, issues=issues,
                      scientific_model_runs_started=0, boundary='No phase or mask decision inferred by code')
        write_json(run_root / 'MODEL_PREFLIGHT.json', report)
        raise RuntimeError('Medical/input preflight blocked: ' + str(run_root / 'MODEL_PREFLIGHT.json'))
    schema = read_json(CODE / 'protocol/feature_schema.json')
    audit = verify_analysis_inputs(run_root, feature_names=schema['feature_names'],
                                   verify_image_hashes=verify_image_hashes)
    # Source clinical counts remain clinical; active ROI count is a separate field.
    with (rev / '03_model_development/analysis_matrix.csv').open(encoding='utf-8-sig', newline='') as handle:
        matrix = list(csv.DictReader(handle))
    if any(float(r['lesions_number']) != 4 for r in matrix if r['patient_id'] in ['83', '94']):
        raise RuntimeError('Source clinical lesion count for patients 83/94 must remain 4')
    audit['bound_files'][str(decision_path(run_root))] = sha(decision_path(run_root))
    audit['bound_files'][str(cfg)] = sha(cfg)
    audit.update(execution_source_sha256=source_fingerprint(CODE), checked_utc=utc(),
                 actual_event_count=sum(int(float(r['death_event'])) for r in matrix),
                 scientific_model_runs_started=0,
                 analysis_execution_basis='author_specified_inputs',
                 medical_evidence_requirements_satisfied=False,
                 submission_ready=False,
                 medical_evidence_limitation='Author adoption of the eight original masks does not supply new physician boundary adjudication or independent lesion identity for P120L03/L04; unknown actual review metadata remain unknown.')
    write_json(run_root / 'MODEL_PREFLIGHT.json', audit)
    return audit


def require_execution_release(run_root, audit):
    marker = Path(run_root) / 'COHORT_READY_FOR_MODELING.json'
    value = read_json(marker) if marker.is_file() else {}
    # A sourced author adoption decision can freeze a research input. It is not
    # a substitute for the separately required radiologist medical evidence.
    if (value.get('status') != 'READY' or value.get('input_decisions_complete') is not True
            or value.get('analysis_matrix_sha256') != audit['analysis_matrix_sha256']
            or value.get('medical_decisions_sha256') != sha(decision_path(run_root))):
        raise RuntimeError('Explicit input-adoption freeze is missing or does not match matrix/decision records')
    if (value.get('cohort_analysis_scope') != 'author_specified_inputs'
            or value.get('medical_evidence_requirements_satisfied') is not False
            or value.get('submission_ready') is not False):
        raise RuntimeError('Author input adoption must not be labelled as completed physician evidence or submission approval')
    bound = value.get('bound_files', {})
    for path, digest in audit['bound_files'].items():
        if bound.get(path) != digest:
            raise RuntimeError('Cohort freeze must bind every audited input file: ' + path)


def require_frozen_model_inputs(run_root):
    """Low-level model entries and workers cannot bypass the same cohort gate."""
    run_root = Path(run_root)
    path = run_root / 'MODEL_PREFLIGHT.json'
    audit = read_json(path) if path.is_file() else {}
    if audit.get('status') != 'PASS' or audit.get('model_training_authorized') is not True:
        raise RuntimeError('No passing frozen cohort preflight; model execution forbidden')
    if decision_issues(run_root):
        raise RuntimeError('Final sourced medical decisions remain incomplete; model execution forbidden')
    require_execution_release(run_root, audit)
    assert_inputs_unchanged(run_root, audit, 'low_level_model_entry')
    if source_fingerprint(CODE) != audit.get('execution_source_sha256'):
        raise RuntimeError('Execution source differs from the frozen cohort preflight')
    return audit


def execute_stage(run_root, audit, script, arguments):
    assert_inputs_unchanged(run_root, audit, script)
    if source_fingerprint(CODE) != audit['execution_source_sha256']:
        raise RuntimeError('Execution code changed after preflight')
    env = dict(os.environ, TACE_PROJECT_ROOT=str(run_root), PYTHONDONTWRITEBYTECODE='1')
    cmd = [sys.executable, '-B', '-X', 'utf8', str(CODE / script), *arguments]
    logdir = Path(run_root) / 'execution_logs'; logdir.mkdir(exist_ok=True)
    label = script.replace('.py', '') + '_' + '_'.join(arguments).replace('--', '')
    logfile = logdir / (datetime.datetime.now().strftime('%Y%m%d_%H%M%S_%f') + '_' + label + '.log')
    start = utc()
    with logfile.open('w', encoding='utf-8') as handle:
        completed = subprocess.run(cmd, env=env, stdout=handle, stderr=subprocess.STDOUT, check=False)
    event = dict(stage=script, arguments=arguments, command=cmd, started_utc=start, finished_utc=utc(),
                 exit_code=completed.returncode, log=str(logfile), matrix_sha256=audit['analysis_matrix_sha256'],
                 execution_source_sha256=audit['execution_source_sha256'])
    with (logdir / 'STAGE_EVENTS.jsonl').open('a', encoding='utf-8') as handle:
        handle.write(json.dumps(event, ensure_ascii=False) + '\n')
    if completed.returncode:
        raise RuntimeError('Stage failed (exit %s): %s; %s' % (completed.returncode, script, logfile))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--run-root', type=Path, required=True)
    parser.add_argument('--stage', choices=['preflight', 'development', 'point', 'bootstrap', 'summary', 'all'], default='preflight')
    parser.add_argument('--execute', action='store_true', help='Required in addition to completed frozen cohort for any model fit')
    parser.add_argument('--workers', type=int, default=4)
    args = parser.parse_args()
    if args.stage != 'preflight' and not args.execute:
        raise RuntimeError('Non-preflight stage requires explicit --execute')
    audit = preflight(args.run_root)
    if args.stage == 'preflight': print(json.dumps(audit, ensure_ascii=False)); return
    require_execution_release(args.run_root, audit)
    rev = args.run_root / 'revision_R1_20261014'
    for d in ['03_model_development', '04_validation']:
        (rev / d).mkdir(parents=True, exist_ok=True)
    stages = ['development', 'point', 'bootstrap', 'summary'] if args.stage == 'all' else [args.stage]
    for stage in stages:
        if stage == 'development': execute_stage(args.run_root, audit, 'diagnostics_and_tests.py', [])
        else:
            execute_stage(args.run_root, audit, 'validate_pipeline.py', ['--mode', stage, '--workers', str(args.workers)])
    if args.stage in ['all', 'summary']:
        from validation_integrity import audit_validation, verify_paired_summary
        integrity = audit_validation(rev); verify_paired_summary(rev, integrity)
        actual = verify_model_output(args.run_root, audit)
        write_json(args.run_root / 'MODEL_ANALYSIS_COMPLETE.json', dict(
            status='MODEL_VALIDATION_COMPLETE_AWAITING_SECONDARY_AND_DOCUMENT_SYNC', completed_utc=utc(),
            actual_input_audit=audit, actual_model_output_audit=actual,
            actual_validation_counts=integrity['counts'],
            validation_integrity_sha256=sha(rev / '04_validation/VALIDATION_INTEGRITY.json'),
            model_training_and_validation_complete=True, secondary_uncertainty_complete=False,
            publication_materials_complete=False,
            analysis_execution_basis='author_specified_inputs',
            medical_evidence_requirements_satisfied=False,
            submission_ready=False))


if __name__ == '__main__':
    try: main()
    except Exception as exc:
        destination = None
        if '--run-root' in sys.argv and sys.argv.index('--run-root') + 1 < len(sys.argv):
            destination = Path(sys.argv[sys.argv.index('--run-root') + 1])
        failure = dict(status='FAILED_OR_BLOCKED', utc=utc(), arguments=sys.argv[1:],
                       exception_type=type(exc).__name__, exception=str(exc), traceback=traceback.format_exc())
        if destination is not None:
            destination.mkdir(parents=True, exist_ok=True)
            marker = destination / 'MODEL_ANALYSIS_COMPLETE.json'
            if marker.is_file():
                archive = destination / 'invalidated_completion_markers'
                archive.mkdir(exist_ok=True)
                target = archive / (datetime.datetime.now().strftime('%Y%m%d_%H%M%S_%f') + '_MODEL_ANALYSIS_COMPLETE.json')
                marker.replace(target)
                failure['prior_completion_marker_preserved_and_invalidated'] = str(target)
            write_json(destination / 'MODEL_EXECUTION_STATUS.json', dict(
                status='BLOCKED_OR_FAILED', updated_utc=utc(), completion_authorized=False,
                details=failure))
            with (destination / 'ENTRY_FAILURES.jsonl').open('a', encoding='utf-8') as handle:
                handle.write(json.dumps(failure, ensure_ascii=False) + '\n')
        print(json.dumps(failure, ensure_ascii=False), file=sys.stderr)
        sys.exit(1)
