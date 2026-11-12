#!/bin/bash
# Run on the cluster login node; never submit transfers to compute nodes.
set +x
set -euo pipefail

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
source "$script_dir/../share/transfer.sh"
transfer_main upload "$HOME/.config/oss/xielab-readwrite.env" "$@"
