"""Pinned identity of audited release inputs, never a runtime self-approved hash."""
from pathlib import Path
import hashlib,json
PACKAGE=Path(__file__).resolve().parent.parent
EXPECTED_RELEASE_MANIFEST_SHA256='2af42a767a921ecb92f1390315c91bd9a5e0b7112635c10269f1422b53664ebd'
RELEASE_VERSION='v1.1.0'
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
