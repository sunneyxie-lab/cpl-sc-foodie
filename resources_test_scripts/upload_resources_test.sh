#!/bin/bash
# Uploads the fixtures consensus_dedup.sh built to the lab's OSS bucket. Run it on
# the lab cluster's login node, with a key that can write under the fixtures' prefix.
set -euo pipefail

cd "$(git rev-parse --show-toplevel)"

die() {
  echo "error: $*" >&2
  exit 1
}

out=resources_test/consensus_dedup
[[ -s "$out/samples.csv" ]] || die "$out/samples.csv not found; run consensus_dedup.sh first"
while IFS=, read -r id input; do
  [[ -s "$out/$input" ]] || die "$id: $out/$input is listed in samples.csv but missing"
done < <(tail -n +2 "$out/samples.csv")
[[ -n "$(ls -A "$out/genome" 2> /dev/null)" ]] || die "$out/genome is missing or empty"

source resources_test_scripts/oss_remote.sh

echo "Uploading $out to $remote"
rclone sync --progress "$out" "$remote"
