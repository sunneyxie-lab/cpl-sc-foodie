# cpl-sc-foodie

Nextflow workflows for single-cell FOODIE data:

- `raw_fastq_align` extracts barcodes, trims paired reads, aligns them and writes
  QC reports. Its first-pass BAM is the input to `consensus_dedup`.
- `consensus_dedup` calls cells, collapses duplicate read pairs per cell and
  re-aligns consensus reads. It publishes a coordinate-sorted BAM, its index,
  the cell list and reports.

## Requirements and setup

Use Nextflow 25.10.6 and Java 17 through 24, with a container runtime:

- **Local Docker:** use `-profile docker`, as in the examples below.
- **Cluster Singularity:** replace `-profile docker` with `-profile singularity`
  in both workflow commands. The runtime must be available on the execution
  nodes. Configure your site's scheduler separately; see
  [cluster configuration](docs/usage.md#runtime-and-site-configuration).

Prepare the original paired gzip FASTQs, a matching barcode whitelist, and a
Bismark genome directory containing reference FASTAs and prepared indexes.
Use the same genome for both workflows.

### Get a published build

Clone the `main_build` branch into a new folder. It includes the generated
workflows and uses public container images, so no Viash build or image-pull
credentials are needed:

```bash
git clone --branch main_build --single-branch \
  https://github.com/sunneyxie-lab/cpl-sc-foodie.git /absolute/path/to/cpl-sc-foodie
repo=/absolute/path/to/cpl-sc-foodie
genome=/absolute/path/to/prepared/bismark/genome
```

Replace these paths with your own. Keep `repo` and `genome` set in the shell
used for the commands below.

**Reproducibility warning:** `main_build` and its container tags move with later
builds. Publication alone does not certify that component or integration checks
passed. Pinning a generated commit does not preserve overwritten container image
bytes. See [release tagging and consumption](docs/development.md#consume-and-tag-a-release)
for versioned builds.

Run the generated workflows under `target/nextflow/workflows/`, rather than the
source files under `src/`. Keep the complete generated tree and its dependencies,
including the configuration files beside each `main.nf`. Do not rebuild a
published tree. To build a source checkout instead, see
[building from source](docs/development.md#build-from-source).

## Process multiple samples

Use `--param_list` with a CSV containing one row per sample. For a single sample
without CSV files, see the [single-sample walkthrough](docs/usage.md#process-one-sample).

### 1. Create an analysis folder

Use a separate folder outside the project checkout for each analysis or batch
so results, intermediate files and logs stay separate from the project files:

```bash
mkdir -p /absolute/path/to/analysis/batch
cd /absolute/path/to/analysis/batch
```

### 2. Align raw FASTQs

Create `samples.raw.csv` in this folder:

```csv
id,read1,read2
sample1,data/sample1.R1.fastq.gz,data/sample1.R2.fastq.gz
sample2,data/sample2.R1.fastq.gz,data/sample2.R2.fastq.gz
```

Replace the example paths with your FASTQs. Paths inside a CSV resolve relative
to that CSV's folder. The command below uses `assets/barcodes.txt`; supply a
different whitelist if your data require one.

```bash
nextflow run "$repo/target/nextflow/workflows/raw_fastq_align/main.nf" \
  -profile docker \
  -work-dir work/raw_fastq_align \
  --param_list samples.raw.csv \
  --whitelist "$repo/assets/barcodes.txt" \
  --genome "$genome" \
  --publish_dir results/raw_fastq_align
```

Each first-pass BAM is written to
`results/raw_fastq_align/<id>.raw_fastq_align.bam`, alongside QC reports and
statistics. The BAM retains barcode suffixes in read names and can be passed
directly to the next workflow.

### 3. Deduplicate and call consensus reads

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

Final outputs are under `results/consensus_dedup/`. For each sample, look for
`<id>.consensus_dedup.bam`, its `.bam.bai` index, `<id>.cell_list.tsv` and reports.
When cells are called, outputs also include per-barcode counts, a cell-calling
QC summary and a knee plot.

## Input rules, paths and resuming

- Each sample needs a unique ID so outputs do not overwrite another sample's
  files. IDs must have no whitespace or leading `@`.
- CSV columns use argument names without `--`. Command-line inputs such as
  `--genome` apply to every row; values in a row override the corresponding
  command-line inputs. Nextflow schedules samples concurrently as resources allow.
- To start from existing BAMs, put their paths in `samples.consensus.csv`. Each
  must be a first-pass, non-directional paired Bismark BAM whose read names end
  in `_CB:Z:<barcode>`.
- To use existing cell lists, add a `cell_list` column to
  `samples.consensus.csv` with the path for each sample. Each list must contain
  `barcode<TAB>count` without a header. Omit the column to call cells
  independently for each sample.

Relative command-line paths resolve from the analysis folder. With the commands
above, results go under `results/`, intermediate task files under `work/`, and
Nextflow keeps `.nextflow/` and `.nextflow.log` in the analysis folder.

Re-run the same command from the same analysis folder with `-resume` to reuse
completed tasks. Keep both `.nextflow/` and `work/` for resuming.

## More help

Show all workflow arguments without running tasks:

```bash
nextflow run "$repo/target/nextflow/workflows/raw_fastq_align/main.nf" --help
nextflow run "$repo/target/nextflow/workflows/consensus_dedup/main.nf" --help
```

- [Usage guide](docs/usage.md): single-sample commands, runtime and cluster
  configuration, paths and resuming.
- [Development and publishing](docs/development.md): source builds, CI,
  publishing, releases and publication failures.
- [Test fixtures and integration runs](resources_test_scripts/README.md).
