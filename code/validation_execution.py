"""Audited validation task execution; a failed future is never a successful run."""
from pathlib import Path
import concurrent.futures
import datetime
import json
import traceback
import uuid

from validation_integrity import (ensure_validation_plan, record_task_event,
                                  audit_validation, checkpoint_issues,
                                  ValidationIntegrityError)


class ValidationTasksFailed(RuntimeError):
    pass


def execute_validation_tasks(revision, mode, run_ids, task, *, workers=4,
                             executor_factory=concurrent.futures.ProcessPoolExecutor):
    """Submit explicit planned seeds and persist attempts before dispatch.

    `attempted` means a concrete invocation submitted/attempted for that run ID;
    it is not reconstructed from checkpoint files. Cached checkpoints are
    separately recorded as reused. Failures from prior attempts remain in JSONL.
    """
    revision = Path(revision)
    plan = ensure_validation_plan(revision)
    ids = list(run_ids)
    if mode not in ['point', 'bootstrap'] or len(set(ids)) != len(ids):
        raise ValueError('Validation mode and unique run IDs are required')
    cfg = json.loads((revision / '00_protocol_and_provenance/revision_config.json').read_text(encoding='utf-8-sig'))
    limit = cfg['point_repeats'] if mode == 'point' else cfg['bootstrap_replicates']
    if not ids or any(type(i) is not int or i < 0 or i >= limit for i in ids):
        raise ValueError('Validation IDs must be within the locked protocol range')
    invocation_id = str(uuid.uuid4())
    ledger = []
    submitted = 0
    succeeded = 0
    reused = 0
    failures = []
    futures = {}

    def failed(run_id, attempt_id, exc):
        entry = record_task_event(revision, mode, run_id, 'failed', attempt_id=attempt_id,
                                  exception={'type': type(exc).__name__, 'message': str(exc),
                                             'traceback': traceback.format_exc()},
                                  details={'invocation_id': invocation_id})
        ledger.append(entry); failures.append(entry)

    try:
        with executor_factory(max_workers=workers) as pool:
            for run_id in ids:
                attempt_id = invocation_id + ':' + str(run_id)
                ledger.append(record_task_event(revision, mode, run_id, 'attempted',
                                                attempt_id=attempt_id,
                                                details={'invocation_id': invocation_id, 'attempt_semantics': 'task_invocation_submitted'}))
                submitted += 1
                try:
                    futures[pool.submit(task, run_id, mode)] = (run_id, attempt_id)
                except Exception as exc:
                    failed(run_id, attempt_id, exc)
            for future in concurrent.futures.as_completed(futures):
                run_id, attempt_id = futures[future]
                try:
                    value = future.result()
                    path = Path(value['path'] if isinstance(value, dict) else value)
                    expected = revision / '04_validation/checkpoints' / f'{mode}_{run_id:04d}.json.gz'
                    if path.resolve() != expected.resolve():
                        raise ValueError('Task returned an unexpected checkpoint path: ' + str(path))
                    import gzip
                    with gzip.open(path, 'rt', encoding='utf-8') as handle:
                        checkpoint = json.load(handle)
                    problems = checkpoint_issues(checkpoint, plan, path)
                    if problems:
                        raise ValueError('Checkpoint incomplete/invalid: ' + json.dumps(problems, ensure_ascii=False))
                    cached = bool(value.get('reused')) if isinstance(value, dict) else False
                    event = 'reused' if cached else 'succeeded'
                    ledger.append(record_task_event(revision, mode, run_id, event, attempt_id=attempt_id,
                                                    details={'invocation_id': invocation_id, 'checkpoint_path': str(path)}))
                    if cached: reused += 1
                    else: succeeded += 1
                except Exception as exc:
                    failed(run_id, attempt_id, exc)
    except Exception as exc:
        # Infrastructure exceptions are separate from completed task outcomes;
        # retain all accepted jobs and explicitly mark unfinished ones failed.
        terminal_ids = {entry.get('attempt_id') for entry in ledger
                        if entry.get('event') in ['succeeded', 'failed', 'reused']}
        for run_id in ids[:submitted]:
            attempt_id = invocation_id + ':' + str(run_id)
            if attempt_id not in terminal_ids:
                failed(run_id, attempt_id, exc)
        infrastructure = {'type': type(exc).__name__, 'message': str(exc), 'traceback': traceback.format_exc()}
    else:
        infrastructure = None
    # Always inspect actual checkpoint content. A report of submitted/terminal
    # events is additional evidence, never a substitute for checkpoint audits.
    integrity = audit_validation(revision, require_complete=False)
    result = {'status': 'FAILED' if failures or infrastructure else ('MODE_TASKS_FINISHED' if len(ids) == limit else 'PARTIAL'),
              'invocation_id': invocation_id, 'mode': mode,
              'protocol_planned': limit, 'this_invocation_planned': len(ids),
              'this_invocation_attempted': submitted, 'this_invocation_succeeded': succeeded,
              'this_invocation_reused': reused, 'this_invocation_failed': len(failures),
              'run_ids': ids, 'failures': failures, 'infrastructure_exception': infrastructure,
              'global_validation_status': integrity['status'],
              'finished_utc': datetime.datetime.now(datetime.timezone.utc).isoformat()}
    output = revision / '04_validation/invocations'
    output.mkdir(parents=True, exist_ok=True)
    (output / (invocation_id + '.json')).write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding='utf-8')
    if failures or infrastructure:
        raise ValidationTasksFailed(f'{mode}: {len(failures)} failed of {submitted} attempted; structured records retained: {output}')
    return result
