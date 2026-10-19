#!/bin/bash
# The lab cluster's compute nodes can't reach any container registry, so run this
# on the login node: it pulls the Singularity images consensus_dedup.sh needs into
# the cache Nextflow reads. With --check it only reports whether all are cached.
set -euo pipefail

cd "$(git rev-parse --show-toplevel)"

die() {
  echo "error: $*" >&2
  exit 1
}

check_only=false
case "${1:-}" in
  "") ;;
  --check) check_only=true ;;
  *) die "unknown argument $1; the only option is --check" ;;
esac

export NXF_SINGULARITY_CACHEDIR="${NXF_SINGULARITY_CACHEDIR:-$HOME/.singularity/nxf_cache}"

module=target/nextflow/bismark/bismark_prepare/main.nf
[[ -f "$module" ]] || die "$module not found; copy the workstation's target/nextflow/ build here"

tool_image=$(awk -F'"' '
  /"container" : \{/ { c = 1; next }
  c && $2 == "registry" { r = $4 }
  c && $2 == "image" { i = $4 }
  c && $2 == "tag" { print r "/" i ":" $4; exit }
' "$module")
[[ -n "$tool_image" ]] || die "no container image found in $module"

# cluster.config's singularity.registry makes Nextflow look the publish image up as dockerproxy.net/...
images=("$tool_image" dockerproxy.net/nextflow/bash:latest)

missing=()
for image in "${images[@]}"; do
  # Nextflow's cache file name for an image
  sif="$NXF_SINGULARITY_CACHEDIR/$(tr '/:' '--' <<< "$image").img"
  if [[ -s "$sif" ]]; then
    echo "cached: $image"
  else
    missing+=("$image")
  fi
done

if [[ "$check_only" == true ]]; then
  (( ${#missing[@]} == 0 )) || die "not in $NXF_SINGULARITY_CACHEDIR: ${missing[*]}"
  exit 0
fi

command -v singularity > /dev/null || die "singularity not found on PATH; ask the cluster admins"
mkdir -p "$NXF_SINGULARITY_CACHEDIR"
for image in "${missing[@]}"; do
  sif="$NXF_SINGULARITY_CACHEDIR/$(tr '/:' '--' <<< "$image").img"
  echo "pulling: $image"
  # Pull under a temporary name so an interrupted pull never looks cached.
  rm -f "$sif.pulling"
  singularity pull "$sif.pulling" "docker://$image" || die "pulling $image failed"
  mv "$sif.pulling" "$sif"
done
