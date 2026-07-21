"""Layout-agnostic repo-root resolution for data/ and cache/.

The SAME sede/ package runs in two tree layouts:
  - dev hub :  <root>/sede/          (research/sede-dev)          -> root is sede/'s parent
  - release :  <root>/src/sede/      (papers/N-*)                 -> root is one level higher

Both keep data/ and cache/ at <root>. Resolving them by walking to the repo root (rather than a
hardcoded number of '..'s) lets the identical file be copied verbatim between trees — so the code
never forks over path depth again. See research/hub/tools/sync_public_repos.py.
"""
import os

def _repo_root():
    here = os.path.dirname(os.path.abspath(__file__))      # <...>/sede
    parent = os.path.dirname(here)                          # <root>  OR  <root>/src
    # Release repos nest the package under src/; the dev hub does not.
    return os.path.dirname(parent) if os.path.basename(parent) == 'src' else parent

REPO_ROOT = _repo_root()
DATA_DIR  = os.path.join(REPO_ROOT, 'data')
CACHE_DIR = os.path.join(REPO_ROOT, 'cache')
