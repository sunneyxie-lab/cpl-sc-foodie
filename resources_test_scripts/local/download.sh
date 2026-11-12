#!/bin/bash
set +x
set -euo pipefail

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
source "$script_dir/../share/transfer.sh"
transfer_main download "$HOME/.config/oss/xielab-readonly.env" "$@"
