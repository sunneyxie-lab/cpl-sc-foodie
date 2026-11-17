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

Alternatively, consume a published generated tree without rebuilding it; see
[CI builds and releases](#ci-builds-and-releases). Source branches, even ones
named like versions, do not contain generated `target/`.

Set `genome` to a directory containing reference FASTAs and prepared Bismark
indexes. Use the same genome for both workflows. The settings above apply to
both single-sample and multi-sample runs below.

## CI builds and releases

### Workflows and outputs

CI uses Viash 0.9.7, Nextflow 25.10.6, Java 17 and Docker on standard
`ubuntu-latest` runners.

| Workflow | Trigger and responsibility |
| --- | --- |
| `viash-test.yml` | Pull requests, pushes to `main`, or reusable calls: all local production components' synthetic tests. PRs are read-only and never publish. |
| `main-build.yml` | Pushes to `main`, or reusable integration/release calls: the single production publisher. Main publishing runs independently of component checks, not after them. |
| `integration-test.yml` | Manual dispatch: publish a dedicated build, then run both complete real-data integrations. A reusable release call tests the supplied generated commit without republishing. No scheduled, push or PR integration runs. |
| `release-build.yml` | Manual dispatch with required `vX.Y.Z` version: publish first, then run component checks and both integrations in parallel against the exact generated commit. |

For manual runs, open the corresponding workflow under **Actions → Run workflow**
and select a **source branch**. Its selected immutable source commit is built.
Do not edit or dispatch from generated output branches:

| Build | Generated branch | Local production image tag / generated `_viash.yaml` version |
| --- | --- | --- |
| Main push | `main_build` | `main_build` |
| Standalone integration | `integration_build` | `integration_build` |
| Release | `release` | Supplied `vX.Y.Z`, e.g. `v0.1.0` |

Images use `ghcr.io/sunneyxie-lab/cpl-sc-foodie/<component>:<version>`, preserving
namespaced component paths. Biobox dependencies retain their upstream images.
Generated branches contain source and the complete production `target/`, not
fixtures or test-only targets/images. Integration checkers are built locally,
never pushed; production containers are pulled from published references.
Each integration downloads and validates common fixtures before its workflow
fixtures anonymously. See the existing [fixture and integration instructions](resources_test_scripts/README.md)
for local usage and scientific coverage.

### Consume and tag a release

`main_build`, `integration_build` and `release` move on later builds. Version tags
are manual: no workflow creates Git tags or GitHub Release objects, or removes
historical release branches/tags. After a successful release run, inspect its
recorded generated `published_sha` and tag that exact commit:

```bash
git fetch origin release
git tag v0.1.0 <published_sha>
git push origin refs/tags/v0.1.0
```

Never tag source `main` or a later moving `release` tip instead. To consume it:

```bash
git fetch origin refs/tags/v0.1.0:refs/tags/v0.1.0
git checkout --detach v0.1.0
repo="$(pwd)"
```

Run `"$repo/target/nextflow/workflows/<name>/main.nf"` with the analysis arguments
shown below; do not rebuild the published tree. Keeping the tag preserves that
generated tree, not overwritten Docker image bytes. For stable consumption,
never reuse its release image version.

### Repository, GHCR and runner prerequisites

The public-runner design requires the user to make the repository public.
That disclosure includes existing code, branches and history; it is a user-owned
decision, not an automated visibility change or a sensitive-history audit.

Publishing uses the workflow's `GITHUB_TOKEN`, not a PAT or local token.
For each existing package under `ghcr.io/sunneyxie-lab/cpl-sc-foodie/`, grant
`sunneyxie-lab/cpl-sc-foodie` **Write** access under **Package settings → Manage
Actions access**. Allow public packages in the organization's package settings.
Set every production package's visibility to **Public**, including `fditools`,
`seqkit`, `trim_galore` and `restore_consensus_tags`.

New GHCR packages start private: a first run may push images but fail anonymous
manifest checks. Make those packages public, then **Re-run all jobs**.
The generated branch is published only after every image passes its anonymous
check. Public package visibility cannot be changed back to private. Workflows
never change repository or package visibility automatically.

Standard public Linux runners provide 4 CPUs, 16 GB RAM and 14 GB SSD; private
runners provide only 2 CPUs/8 GB and do not establish intended integration
acceptance. CI retains 4-CPU/13-GB process requests and `process.maxForks = 1`.
Disk, total memory and runtime sufficiency remain unmeasured. Each scientific
integration command has a 120-minute budget, not a measured runtime guarantee.

### Failures, retries and concurrency

Publication is not atomic and has no automatic rollback or retry:

- **Before publication:** invalid inputs, generation or component-discovery
  failures leave remote state unchanged.
- **During image pushes or anonymous checks:** some tags may already be replaced,
  while the generated branch stays unchanged. Cancellation, timeout or a lost
  runner can leave the same partial state. Inspect Actions logs and registry
  state, resolve access/capacity issues, then **Re-run all jobs**.
- **During branch publication:** images already exist. A failed push leaves the
  previous branch, but connection loss after server acceptance may leave the new
  branch despite a failed job. Inspect the remote commit before a full rerun.
  Failure to resolve the generated SHA or write the summary after a push does
  not undo publication.
- **After publication:** failed/cancelled component or integration checks,
  fixture validation failures and runner exhaustion leave branches and images
  published. Publication alone does not certify a release.

Publishing retries, including main, require **Re-run all jobs**; stale partial
retries are rejected before remote mutation. Standalone integration retries
also require a fresh full publish-and-test attempt, not failed-jobs-only retries.
Release test-only retries can reuse the recorded generated SHA and versioned
images. Full reruns can republish tags, including release versions: do not reuse
a version if its image bytes must stay stable.

Publishers serialize by mode; all release versions share one publisher lock.
Standalone integration holds one shared lock across publication and both tests,
regardless of source branch. Releases serialize publication and both check
suites. In-progress runs are not automatically cancelled, but GitHub's default
concurrency queue can replace older pending runs; not every queued run will run.
Serialization prevents workflow publishers overlapping, not external consumers
pulling during publication. Shared main/integration tags remain mutable.
Status and failure notifications use ordinary GitHub Actions statuses and logs.

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
