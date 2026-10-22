"""Restore the packaged historical audit, without requiring submitted DOCX/PDF files."""
from pathlib import Path
import hashlib,json
EXPECTED_SHA256 = "991ff138088dc65007d55808f4269621c942c6eff0762d447d913e021194b38f"

def ensure_legacy_pipeline_audit(project_root, code_root=None):
    code=Path(code_root) if code_root is not None else Path(__file__).resolve().parent
    source=code/'reporting_inputs/01_data_audit/legacy_pipeline_audit.csv'
    if not source.is_file():
        raise FileNotFoundError('Missing packaged historical audit: '+str(source)+'. Restore the exact release input; do not invent historical findings.')
    from release_identity import verify_release
    release=verify_release(code.parent)
    if release['files'].get('code/reporting_inputs/01_data_audit/legacy_pipeline_audit.csv')!=EXPECTED_SHA256:
        raise RuntimeError('Historical audit release version mismatch')
    data=source.read_bytes()
    if hashlib.sha256(data).hexdigest()!=EXPECTED_SHA256:
        raise RuntimeError('Packaged historical audit hash mismatch')
    target=Path(project_root)/'revision_R1_20261014/01_data_audit/legacy_pipeline_audit.csv'
    target.parent.mkdir(parents=True,exist_ok=True)
    if target.exists() and target.read_bytes()!=data:
        raise RuntimeError('Existing historical audit differs; preserve and resolve explicitly: '+str(target))
    created=not target.exists()
    if created:target.write_bytes(data)
    provenance=target.with_suffix('.provenance.json')
    provenance.write_text(json.dumps({'role':'Frozen historical pipeline audit, not a newly executed historical document audit','packaged_source':'reporting_inputs/01_data_audit/legacy_pipeline_audit.csv','sha256':EXPECTED_SHA256,'created':created},indent=2),encoding='utf-8')
    return target
