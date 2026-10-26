"""Bind the public reproduction to adopted research inputs, not medical approval."""
from pathlib import Path
import json, shutil
from frozen_input_pipeline import PACKAGE, sha, dump
from release_identity import verify_release
from run_new_cohort import preflight, require_execution_release, require_frozen_model_inputs

def prepare_execution_gate(run_root, raw):
    """Set up the same fail-closed current-source gate before a new model run.

    The public decision record is an explicit author input-adoption summary.
    It contains neither private physician forms nor fabricated new review data.
    """
    root=Path(run_root)
    release=verify_release()
    source=PACKAGE/'frozen_inputs/AUTHOR_INPUT_ADOPTION.csv'
    target=root/'PUBLIC_INPUT_ADOPTION.csv'
    if target.exists() and target.read_bytes()!=source.read_bytes():
        raise RuntimeError('Existing author-input record differs; cannot overwrite')
    if not target.exists():shutil.copy2(source,target)
    audit=preflight(root,verify_image_hashes=bool(raw))
    marker=dict(status='READY',input_decisions_complete=True,medical_decisions_complete=False,
                medical_evidence_requirements_satisfied=False,all_decisions_physician_confirmed=False,
                medical_evidence_status='INCOMPLETE',submission_ready=False,
                cohort_analysis_scope='author_specified_inputs',
                analysis_matrix_sha256=audit['analysis_matrix_sha256'],
                medical_decisions_sha256=sha(target),bound_files=audit['bound_files'],
                release_manifest_sha256=release['manifest_sha256'],
                explanation='The legacy medical_decisions_sha256 key binds author research-input adoption only; it does not certify independent physician evidence.')
    dump(root/'COHORT_READY_FOR_MODELING.json',marker)
    require_execution_release(root,audit)
    require_frozen_model_inputs(root)
    return audit
