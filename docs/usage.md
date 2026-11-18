# Usage guide

For setup and the batch walkthrough, start with the [README](../README.md).
The examples below use its `repo` and `genome` variables and a published build.

## Process one sample

Pass one sample's ID and input paths directly on the command line instead of
using `--param_list`.

### Create an analysis folder

Use a separate folder outside the project checkout so results, intermediate
files and logs stay separate from the project files:

```bash
mkdir -p /absolute/path/to/analysis/sample
cd /absolute/path/to/analysis/sample

sample=sample
```

Replace the example paths and sample ID. IDs must have no whitespace or leading
`@`.

### Align raw FASTQs

Supply the original paired gzip FASTQs and a matching barcode whitelist:

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

Run this after alignment finishes:

```bash
nextflow run "$repo/target/nextflow/workflows/consensus_dedup/main.nf" \
  -profile docker \
  -work-dir work/consensus_dedup \
  --id "$sample" \
  --input "results/raw_fastq_align/$sample.raw_fastq_align.bam" \
  --genome "$genome" \
  --publish_dir results/consensus_dedup
```

Outputs include `results/consensus_dedup/$sample.consensus_dedup.bam`, its
`.bai` index, the cell list and reports.

To start from an existing BAM, replace `--input` with its absolute path. It must
be a first-pass, non-directional paired Bismark BAM whose read names end in
`_CB:Z:<barcode>`. Use the same prepared genome as for the first-pass alignment.

To use an existing cell list instead of calling cells, add
`--cell_list /absolute/path/to/sample.cell_list.tsv`. The list must contain
`barcode<TAB>count` without a header.

## Runtime and site configuration

The examples use Docker. For Singularity, replace `-profile docker` with
`-profile singularity` in each command. Both generated workflows provide these
profiles.

Add `-c /absolute/path/to/analysis.config` to either workflow command for
site-specific execution settings. A container profile selects the runtime, not
a cluster scheduler. Configure your site's Nextflow executor, queue and resource
limits in that file according to your cluster's policies.

On a cluster:

- Make the generated checkout, inputs, genome, analysis folder and work directory
  accessible from the execution nodes. Containers must also be able to access
  these paths, including the targets of any input symlinks.
- Put the Singularity image cache on storage shared by the launch host and
  execution nodes. Set it before running Nextflow:

  ```bash
  export NXF_SINGULARITY_CACHEDIR=/absolute/path/to/shared/singularity/cache
  mkdir -p "$NXF_SINGULARITY_CACHEDIR"
  ```

- The launch environment needs access to the image registries when downloading
  uncached images. Follow your site's rules for image pulls and job submission.

See Nextflow's [Singularity documentation](https://docs.seqera.io/nextflow/container/singularity)
for automatic mounts and image caching, and its
[executor reference](https://docs.seqera.io/nextflow/reference/executor) for
scheduler-specific settings.

Successfully published builds use public GHCR images and need no image-pull
credentials. Unpublished local images must be built locally or made available
to the runtime. For source builds, see the [development guide](development.md#build-from-source).

## Paths and resuming

Relative command-line paths resolve from the analysis folder. Input paths in a
CSV resolve relative to that CSV's folder. With the example commands, results
go under `results/`, intermediate task files under `work/`, and Nextflow keeps
`.nextflow/` and `.nextflow.log` in the analysis folder.

Re-run the same command from the same analysis folder with `-resume` to reuse
completed tasks. Keep both `.nextflow/` and `work/` for resuming.

For test data and integration runs, see the
[test fixture instructions](../resources_test_scripts/README.md).
