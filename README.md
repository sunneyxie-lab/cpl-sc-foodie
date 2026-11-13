# cpl-sc-foodie

Nextflow workflows for single-cell FOODIE data.

## Build the workflows

Use Viash 0.9.7, Nextflow 25.10.6, Java 17 through 24 and Docker.
Set `repo` to the absolute path of your checkout, then build once:

```bash
repo=/absolute/path/to/cpl-sc-foodie
genome=/absolute/path/to/prepared/bismark/genome
(cd "$repo" && viash ns build --parallel --setup cachedbuild)
```

Run the generated workflows under `target/nextflow/workflows/`, rather than the
source files under `src/`. Keep the complete generated tree and its dependencies,
including the configuration files beside each `main.nf`.

Set `genome` to a directory containing reference FASTAs and prepared Bismark
indexes. Use the same genome for both workflows. The settings above apply to
both single-sample and multi-sample runs below.

## CI builds and releases

GitHub Actions builds the workflows and tests every local component with synthetic
data on pull requests, pushes to `main` and manual runs. CI uses Viash 0.9.7 and
Java 17. Real-data integration remains separate; CI downloads no OSS fixtures.

After all component tests pass, pushes to `main` publish production Docker images
tagged `build_main`, then update the private `build/main` branch with the source
and generated `target/` tree. Test-only components and test fixtures are excluded.
Pull requests never publish. Dependencies from biobox retain their upstream images.

For a release, run **Actions → Build → Run workflow** and enter an `X.Y.Z` version,
for example `0.1.0`. This publishes matching image tags and updates
`release/0.1.x`. Leaving the version empty creates a development build for the
selected branch. Development tags and release branches move on later builds;
do not reuse a release version if its images must remain unchanged.

### GHCR setup

Publishing uses the workflow's `GITHUB_TOKEN`, not a developer's local token.
For each existing package under `ghcr.io/sunneyxie-lab/cpl-sc-foodie/`, grant
`sunneyxie-lab/cpl-sc-foodie` **Write** access under **Package settings → Manage
Actions access**. Allow public packages in the organization's package settings.

Set each production package's visibility to **Public**, including `fditools`,
`seqkit`, `trim_galore` and `restore_consensus_tags`. New GHCR packages start
private. On the first publishing run, new images may be pushed but fail the
anonymous-access check. Make those packages public and rerun the failed jobs.
The generated branch is published only after every image passes that check.
Public package visibility cannot be changed back to private.

Publishing is not atomic: an image failure can leave some tags updated while the
generated branch stays unchanged. The workflow never changes repository or
package visibility automatically.

## Choose an analysis folder

We recommend a separate folder outside the project checkout for each analysis
or batch. Run Nextflow from that folder so results, intermediate files and logs
stay separate from the project files.

## Process one sample

Pass one sample's ID and input paths directly on the command line.

### Create an analysis folder

```bash
mkdir -p /absolute/path/to/analysis/sample
cd /absolute/path/to/analysis/sample

sample=sample
```

Replace the example paths and sample ID. IDs must have no whitespace or leading
`@`.

### Align raw FASTQs

`raw_fastq_align` extracts barcodes, trims paired reads, aligns them and writes QC
reports. Supply the original paired gzip FASTQs and a matching barcode whitelist.

```bash
nextflow run "$repo/target/nextflow/workflows/raw_fastq_align/main.nf" \
  -profile docker \
  -work-dir work/raw_fastq_align \
  --id "$sample" \
  --read1 /absolute/path/to/sample.R1.fastq.gz \
  --read2 /absolute/path/to/sample.R2.fastq.gz \
  --whitelist "$repo/assets/barcodes.txt" \
  --genome "$genome" \
  --publish_dir results/raw_fastq_align
```

The first-pass BAM is `results/raw_fastq_align/$sample.raw_fastq_align.bam`.
It retains barcode suffixes in read names and can be passed directly to
`consensus_dedup`.

### Deduplicate and call consensus reads

`consensus_dedup` calls cells, collapses duplicate read pairs per cell and
re-aligns consensus reads. It publishes a coordinate-sorted BAM, its index and
reports.

```bash
nextflow run "$repo/target/nextflow/workflows/consensus_dedup/main.nf" \
  -profile docker \
  -work-dir work/consensus_dedup \
  --id "$sample" \
  --input "results/raw_fastq_align/$sample.raw_fastq_align.bam" \
  --genome "$genome" \
  --publish_dir results/consensus_dedup
```

Run this after alignment finishes. To start from an existing BAM, replace
`--input` with its absolute path. It must be a first-pass, non-directional paired
Bismark BAM whose read names end in `_CB:Z:<barcode>`.

To use an existing cell list instead of calling cells, add
`--cell_list /absolute/path/to/sample.cell_list.tsv`. The list must contain
`barcode<TAB>count` without a header.

## Process multiple samples

Use `--param_list` with a CSV containing one row per sample instead of passing
one sample's ID and input paths on the command line.

### Create an analysis folder

Use a separate analysis folder for the batch. Keep the shared `repo` and
`genome` settings from the build section:

```bash
mkdir -p /absolute/path/to/analysis/batch
cd /absolute/path/to/analysis/batch
```

### Align raw FASTQs

Create `samples.raw.csv` in this folder:

```csv
id,read1,read2
sample1,data/sample1.R1.fastq.gz,data/sample1.R2.fastq.gz
sample2,data/sample2.R1.fastq.gz,data/sample2.R2.fastq.gz
```

Replace the example paths with your FASTQs. Input paths inside a CSV resolve
relative to that CSV's folder. Each sample needs a unique ID so its outputs
don't overwrite another sample's files.

```bash
nextflow run "$repo/target/nextflow/workflows/raw_fastq_align/main.nf" \
  -profile docker \
  -work-dir work/raw_fastq_align \
  --param_list samples.raw.csv \
  --whitelist "$repo/assets/barcodes.txt" \
  --genome "$genome" \
  --publish_dir results/raw_fastq_align
```

### Deduplicate and call consensus reads

After alignment finishes, create `samples.consensus.csv` in the same folder:

```csv
id,input
sample1,results/raw_fastq_align/sample1.raw_fastq_align.bam
sample2,results/raw_fastq_align/sample2.raw_fastq_align.bam
```

```bash
nextflow run "$repo/target/nextflow/workflows/consensus_dedup/main.nf" \
  -profile docker \
  -work-dir work/consensus_dedup \
  --param_list samples.consensus.csv \
  --genome "$genome" \
  --publish_dir results/consensus_dedup
```

Nextflow schedules samples concurrently as resources allow. Command-line inputs
such as `--genome` apply to every row. CSV columns use argument names without
`--`, and values in a row override the corresponding command-line inputs.
For existing cell lists, add a `cell_list` column to `samples.consensus.csv` with
the path for each sample. Omit it to call cells independently for each sample.

## Paths, configuration and resuming

Relative command-line paths resolve from the analysis folder. Results go under
`results/`, intermediate task files under `work/`, and Nextflow keeps `.nextflow/`
and `.nextflow.log` in the analysis folder.

Re-run the same command from the same analysis folder with `-resume` to reuse
completed tasks. Keep both `.nextflow/` and `work/` for resuming.

Add `-c /absolute/path/to/analysis.config` for site-specific execution settings.
For a Singularity runtime, replace `-profile docker` with `-profile singularity`.
The runtime must have access to the required container images. CI requires public
GHCR images, so successfully published builds need no image-pull credentials.
Unpublished local images must be built locally or made available to the runtime.

Show all workflow arguments without running tasks:

```bash
nextflow run "$repo/target/nextflow/workflows/raw_fastq_align/main.nf" --help
nextflow run "$repo/target/nextflow/workflows/consensus_dedup/main.nf" --help
```

See [test fixture instructions](resources_test_scripts/README.md) for test data
and integration runs.
