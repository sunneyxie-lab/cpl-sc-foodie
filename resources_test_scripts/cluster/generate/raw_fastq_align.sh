#!/bin/bash
#SBATCH --job-name=raw_fastq_align_fixtures
#SBATCH --partition=cpu1,cpu2
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=24:00:00
# Run inside a cluster allocation with the original BAM/FASTQ sources mounted.
# Upload the completed selected set from a cluster login node.
set -euo pipefail

resources_dir=""
if [[ $# -gt 0 ]]; then
  if [[ $# -ne 2 || "$1" != --resources-dir || -z "$2" ]]; then
    echo "usage: $0 [--resources-dir <root>]" >&2
    exit 1
  fi
  resources_dir="$2"
fi
installed_repo=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../../.." && pwd)
if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
  repo="$SLURM_SUBMIT_DIR"
elif [[ -f "$installed_repo/resources_test_scripts/share/raw_fastq_align_extract.py" ]]; then
  repo="$installed_repo"
else
  repo="${SLURM_SUBMIT_DIR:?error: submit this job from the checkout}"
fi
if [[ ! -f "$repo/resources_test_scripts/share/raw_fastq_align_extract.py" ]]; then
  echo "error: checkout not found at $repo; submit this job from the checkout" >&2
  exit 1
fi
command -v python3 >/dev/null || { echo 'error: Python 3 is required' >&2; exit 1; }
command -v samtools >/dev/null || { echo 'error: samtools >=1.12 is required' >&2; exit 1; }
command -v bqtools >/dev/null || { echo 'error: bqtools 0.6.1 is required' >&2; exit 1; }

# Resolve relative roots in the caller's directory before any directory change.
# Python owns preflight, per-sample extraction, failure cleanup and publication.
python3 - "$repo" "$resources_dir" <<'PY'
import shlex
import sys
from pathlib import Path

repo = Path(sys.argv[1]).resolve()
sys.path.insert(0, str(repo / "resources_test_scripts/share"))
from fixture_integrity import resources_root
from raw_fastq_align_extract import generate_fixture_set

try:
    root = resources_root(repo, sys.argv[2] or None)
    result = generate_fixture_set(repo, root)
except Exception as exc:
    raise SystemExit(f"error: {exc}")
print(f"Completed fixture set: {result}")
upload = repo / "resources_test_scripts/cluster/upload.sh"
print(f"On cluster login, upload with: bash {shlex.quote(str(upload))} raw_fastq_align --resources-dir {shlex.quote(str(root))}")
PY
