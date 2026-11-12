#!/bin/bash
#SBATCH --job-name=consensus_dedup_fixtures
#SBATCH --partition=cpu1,cpu2
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=8:00:00
#SBATCH --output=logs/consensus_dedup_fixtures.%j.out
#SBATCH --error=logs/consensus_dedup_fixtures.%j.err
# Submit from the checkout with sbatch, or run with bash inside an allocation.
# common must already exist; upload the completed set separately from login.
set -euo pipefail

command -v python3 >/dev/null || { echo "error: Python 3 is required on PATH" >&2; exit 1; }
script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
  # sbatch executes a spooled copy; its directory is not the checkout.
  driver="$SLURM_SUBMIT_DIR/resources_test_scripts/cluster/generate/consensus_dedup.py"
else
  driver="$script_dir/consensus_dedup.py"
fi
[[ -f "$driver" ]] || { echo "error: cannot find $driver; submit from the checkout" >&2; exit 1; }
# Do not change directory: Python resolves an explicit relative root against the caller.
exec python3 "$driver" "$@"
