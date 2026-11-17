#!/bin/bash
# Sourced by the upload/download entrypoints; transfer exactly one fixture set.

transfer_die() {
  printf 'error: %s\n' "$*" >&2
  exit 1
}

transfer_usage() {
  printf 'Usage: bash %s <common|consensus_dedup|raw_fastq_align> [--resources-dir <root>]\n' "$0" >&2
  printf '%s\n' \
    'The set and --resources-dir may appear in either order.' \
    "The default root is this checkout's resources_test/, not the current directory." \
    'Sync replaces the selected destination in place, deleting destination-only files.' \
    'Interruption can leave that destination invalid; rerun the same named transfer.' \
    'Other sets are never transferred or repaired automatically.' >&2
}

transfer_main() {
  local operation=$1
  shift
  local fixture_set= root_argument= root_seen=false
  while (($#)); do
    case "$1" in
      --resources-dir)
        if [[ "$root_seen" == true || $# -lt 2 || -z "$2" || "$2" == --* ]]; then
          transfer_usage
          transfer_die '--resources-dir requires one nonempty root and may appear only once'
        fi
        root_argument=$2
        root_seen=true
        shift 2
        ;;
      common|consensus_dedup|raw_fastq_align)
        if [[ -n "$fixture_set" ]]; then
          transfer_usage
          transfer_die 'select exactly one fixture set'
        fi
        fixture_set=$1
        shift
        ;;
      *)
        transfer_usage
        transfer_die "unsupported fixture set or option: $1"
        ;;
    esac
  done
  if [[ -z "$fixture_set" ]]; then
    transfer_usage
    transfer_die 'select exactly one fixture set'
  fi

  command -v python3 > /dev/null || transfer_die 'python3 not found on PATH'
  local share_dir repo resources_dir directory remote
  share_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
  repo=$(cd -- "$share_dir/../.." && pwd)
  resources_dir=$(python3 - "$repo" "$root_argument" <<'PY'
import sys
from pathlib import Path

repo = Path(sys.argv[1])
sys.path.insert(0, str(repo / "resources_test_scripts" / "share"))
from fixture_integrity import resources_root

print(resources_root(repo, sys.argv[2] or None))
PY
  ) || transfer_die 'cannot resolve the resources root'
  directory="$resources_dir/$fixture_set"
  [[ ! -L "$directory" ]] || transfer_die "selected destination must not be a symlink: $directory"

  if [[ "$operation" == upload ]]; then
    local credential_file="$HOME/.config/oss/xielab-readwrite.env"
    [[ -r "$credential_file" ]] || transfer_die "cannot read required credential file: $credential_file"
    # Only the upload file, not ambient OSS variables, supplies credentials.
    unset OSS_ACCESS_KEY_ID OSS_ACCESS_KEY_SECRET
    set +x
    source "$credential_file"
    set +x
  fi
  source "$share_dir/oss_remote.sh" "$operation"
  remote="oss:xielab-foodie/cpl-sc-foodie/resources_test/$fixture_set/"

  if [[ "$operation" == upload ]]; then
    python3 "$share_dir/fixture_integrity.py" validate "$fixture_set" --resources-dir "$resources_dir" || \
      transfer_die 'selected fixtures failed validation; nothing was uploaded'
    printf 'Uploading %s to %s\n' "$directory" "$remote"
    rclone sync --progress "$directory/" "$remote" || \
      transfer_die "upload failed; the selected remote prefix may be incomplete. Rerun this upload for $fixture_set"
  else
    if [[ "$fixture_set" != common ]]; then
      python3 "$share_dir/fixture_integrity.py" validate common --resources-dir "$resources_dir" || \
        transfer_die 'common must be valid before downloading a workflow set; nothing was downloaded'
    fi
    printf 'Downloading %s to %s\n' "$remote" "$directory"
    rclone sync "$remote" "$directory/" || \
      transfer_die "download failed; the selected local directory may be incomplete. Rerun this download for $fixture_set"
    python3 "$share_dir/fixture_integrity.py" validate "$fixture_set" --resources-dir "$resources_dir" || \
      transfer_die "downloaded fixtures failed validation; do not use them. Rerun this download for $fixture_set, or restore compatible common explicitly"
  fi
  printf 'Completed %s of validated %s fixtures.\n' "$operation" "$fixture_set"
}
