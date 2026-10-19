#!/bin/bash
#SBATCH --job-name=consensus_dedup_fixtures
#SBATCH --partition=cpu1,cpu2
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=8:00:00
#SBATCH --output=logs/consensus_dedup_fixtures.%j.out
#SBATCH --error=logs/consensus_dedup_fixtures.%j.err
# Builds the consensus_dedup test fixtures on the lab cluster into
# resources_test/consensus_dedup. Submit from the checkout with sbatch, or run it
# with bash inside an interactive allocation. Compute nodes can't reach OSS, so
# upload afterwards from the login node with upload_resources_test.sh.
# --reuse-subsets keeps the BAM subsets of an earlier run and only rechecks them.
set -euo pipefail

cd "$(git rev-parse --show-toplevel)"

SOURCE_BAMS=(
  /lustre/grp/xxllab/fangj/03.project/02.snFOODIE/H0O260120CPL0071F1S/reRun/Align_bam/H0O260120CPL0071F1S_001.bam
  /lustre/grp/xxllab/fangj/03.project/02.snFOODIE/H0O260120CPL0071F1S/reRun/Align_bam/H0O260120CPL0071F1S_002.bam
)
REFERENCE_FASTA=/lustre/grp/xxllab/zhangzhen/FOODIE_Training_CP_20260305/resources/GRCh38/genome.fa
REGION_CHROM=chr14
REGION_START=18200001   # 1-based, inclusive
REGION_END=32900000

die() {
  echo "error: $*" >&2
  exit 1
}

reuse_subsets=false
case "${1:-}" in
  "") ;;
  --reuse-subsets) reuse_subsets=true ;;
  *) die "unknown argument $1; the only option is --reuse-subsets" ;;
esac

for path in "${SOURCE_BAMS[@]}" "$REFERENCE_FASTA"; do
  [[ -r "$path" ]] || die "cannot read $path; run this on the lab cluster"
done

command -v samtools > /dev/null || die "samtools not found on PATH"
command -v singularity > /dev/null || die "singularity not found on PATH; ask the cluster admins"
command -v nextflow > /dev/null || die "nextflow not found on PATH. Install it with: curl -s https://get.nextflow.io | NXF_VER=25.10.6 bash (needs Java 17-24)"

samtools_version=$(samtools --version | head -1 | awk '{ print $2 }')
IFS=. read -r st_major st_minor _ <<< "$samtools_version"
if (( st_major < 1 || (st_major == 1 && st_minor < 12) )); then
  die "samtools $samtools_version is too old; view -e needs 1.12 or newer"
fi

nextflow -version > /dev/null 2>&1 \
  || die "nextflow -version failed. Nextflow 25.10.6 needs Java 17-24; if the cluster's Java is older, unpack a JDK 17 (e.g. Temurin) and set JAVA_HOME"

[[ -f target/nextflow/bismark/bismark_prepare/main.nf ]] \
  || die "target/nextflow/bismark/bismark_prepare/main.nf not found; copy the workstation's target/nextflow/ build here"

export NXF_SINGULARITY_CACHEDIR="${NXF_SINGULARITY_CACHEDIR:-$HOME/.singularity/nxf_cache}"
resources_test_scripts/cluster_images.sh --check \
  || die "compute nodes can't reach the registries; run resources_test_scripts/cluster_images.sh on the login node first"

out=resources_test/consensus_dedup
if [[ "$reuse_subsets" == true ]]; then
  rm -rf "$out/genome" "$out/samples.csv"
else
  rm -rf "$out"
fi
mkdir -p "$out"
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

cpus=${SLURM_CPUS_PER_TASK:-$(nproc)}
if [[ -n "${SLURM_MEM_PER_NODE:-}" ]]; then
  mem_gb=$(( SLURM_MEM_PER_NODE / 1024 ))
else
  mem_gb=$(( $(awk '/^MemAvailable:/ { print $2 }' /proc/meminfo) / 1024 / 1024 ))
fi
# Both subsets run at once; -@ counts threads on top of the main one.
samtools_threads=$(( cpus > 3 ? cpus / 2 - 1 : 1 ))

# Every command gets its own "|| die" so a failure names the sample and the step;
# set -e alone would end the background job without a message.
subset_bam() {
  local bam="$1" id sub
  id=$(basename "$bam" .bam)
  sub="$out/$id.bam"
  if [[ "$reuse_subsets" == true ]]; then
    [[ -f "$sub" ]] || die "$id: --reuse-subsets given, but $sub does not exist"
    echo "Reusing $sub"
  else
    echo "Subsetting $id to $REGION_CHROM:$REGION_START-$REGION_END"
    samtools view -h -b --no-PG -@ "$samtools_threads" \
      -e "rname==\"$REGION_CHROM\" && pos>=$REGION_START && pos<=$REGION_END && mpos>=$REGION_START && mpos<=$REGION_END" \
      -o "$sub" "$bam" \
      || die "$id: samtools view failed while subsetting $bam"
  fi

  samtools quickcheck "$sub" || die "$id: $sub is truncated or not a valid BAM"

  local n_records unpaired n_unbarcoded n_barcodes
  n_records=$(samtools view -c "$sub") || die "$id: counting records failed"
  (( n_records > 0 )) || die "$id: subset is empty"

  unpaired=$(samtools view "$sub" | cut -f1 | sort -T "$tmp" | uniq -c | awk '$1 != 2 && n++ < 5') \
    || die "$id: pairing check failed"
  [[ -z "$unpaired" ]] || die "$id: read names not occurring exactly twice: $unpaired"

  n_unbarcoded=$(samtools view "$sub" | awk -F'\t' '$1 !~ /_CB:Z:[^_]+$/ { n++ } END { print n + 0 }') \
    || die "$id: barcode suffix check failed"
  (( n_unbarcoded == 0 )) || die "$id: $n_unbarcoded records have no _CB:Z:<barcode> read name suffix"

  n_barcodes=$(samtools view "$sub" | cut -f1 | sed 's/.*_CB:Z://' | sort -T "$tmp" -u | wc -l) \
    || die "$id: counting barcodes failed"
  (( n_barcodes >= 2 )) || die "$id: $n_barcodes distinct barcodes, expected at least 2"

  echo "$id: $n_records records, $n_barcodes barcodes"
}

pids=()
for bam in "${SOURCE_BAMS[@]}"; do
  subset_bam "$bam" &
  pids+=($!)
done
failed=false
for i in "${!pids[@]}"; do
  wait "${pids[$i]}" || { echo "error: subset of ${SOURCE_BAMS[$i]} failed; see above" >&2; failed=true; }
done
[[ "$failed" == false ]] || exit 1

echo "Extracting $REGION_CHROM from $REFERENCE_FASTA"
# No early exit from awk: it would SIGPIPE zcat and fail the pipeline under pipefail.
zcat -f "$REFERENCE_FASTA" \
  | awk -v chrom="$REGION_CHROM" '
      /^>/ { split(substr($0, 2), w, /[ \t]/); keep = (w[1] == chrom) }
      keep
    ' > "$tmp/$REGION_CHROM.fa"
[[ -s "$tmp/$REGION_CHROM.fa" ]] || die "$REGION_CHROM not found in $REFERENCE_FASTA"

# Without the cap, bismark_prepare's labels ask for more than the allocation and
# the local executor refuses to run it. The head JVM keeps a few GB.
echo "process.resourceLimits = [cpus: $cpus, memory: $(( mem_gb > 8 ? mem_gb - 4 : mem_gb / 2 )).GB]" > "$tmp/resources.config"

# Everything Nextflow needs is local or cached; offline keeps it from trying the network.
export NXF_OFFLINE=true
nextflow run target/nextflow/bismark/bismark_prepare/main.nf \
  -profile singularity \
  -c resources_test_scripts/cluster.config \
  -c "$tmp/resources.config" \
  -work-dir "$tmp/work" \
  --input "$tmp/$REGION_CHROM.fa" \
  --output genome \
  --publish_dir "$tmp/prepare"
cp -r "$tmp/prepare/genome" "$out/genome"

{
  echo "id,input"
  for bam in "${SOURCE_BAMS[@]}"; do
    id=$(basename "$bam" .bam)
    echo "$id,$id.bam"
  done
} > "$out/samples.csv"

echo "Fixtures ready in $out; upload them from the login node with resources_test_scripts/upload_resources_test.sh"
