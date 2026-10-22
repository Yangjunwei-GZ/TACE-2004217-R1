"""Portable frozen input entry; see README for both supported routes."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parent))
from frozen_run_revision import main
if __name__=="__main__":sys.exit(main() or 0)
