import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
REPO = ROOT.parent
for p in (ROOT, REPO):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))
