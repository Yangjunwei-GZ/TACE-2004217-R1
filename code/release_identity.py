"""Pinned identity of audited release inputs, never a runtime self-approved hash."""
from pathlib import Path
import hashlib,json
PACKAGE=Path(__file__).resolve().parent.parent
EXPECTED_RELEASE_MANIFEST_SHA256='c1af5ab2699fef1d22bd2b08f5475b4f20a9526a68eac028f72e1bf506c408db'
RELEASE_VERSION='local-input-route-fix-20261022-r2'
def verify_release(package=None):
 root=Path(package) if package is not None else PACKAGE;p=root/'RELEASE_INPUTS.json'
 if not p.is_file():raise FileNotFoundError('Missing pinned release input manifest')
 raw=p.read_bytes();actual=hashlib.sha256(raw).hexdigest()
 if actual!=EXPECTED_RELEASE_MANIFEST_SHA256:raise RuntimeError('Release input manifest hash/version mismatch')
 d=json.loads(raw)
 if d['version']!=RELEASE_VERSION:raise RuntimeError('Wrong release version')
 for name,expected in d['files'].items():
  path=(root/name).resolve()
  if root.resolve() not in path.parents or not path.is_file():raise FileNotFoundError('Missing release input: '+name)
  if hashlib.sha256(path.read_bytes()).hexdigest()!=expected:raise RuntimeError('Release input hash mismatch: '+name)
 return {'manifest_sha256':actual,'version':d['version'],'files':d['files']}
