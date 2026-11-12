#!/bin/bash
# Run on the cluster login node, never through sbatch.
set -euo pipefail

command -v python3 >/dev/null || { echo 'error: Python 3 is required' >&2; exit 1; }
script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
exec python3 "$script_dir/common.py" "$@"
