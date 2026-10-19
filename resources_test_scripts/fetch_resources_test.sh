#!/bin/bash
set -euo pipefail

cd "$(git rev-parse --show-toplevel)"

source resources_test_scripts/oss_remote.sh

rclone sync "$remote" resources_test/consensus_dedup
