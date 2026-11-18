nextflow.enable.dsl=2

include { raw_fastq_align } from params.rootDir + '/target/nextflow/workflows/raw_fastq_align/main.nf'
include { fditools } from params.rootDir + '/target/nextflow/fditools/main.nf'
include { trim_galore } from params.rootDir + '/target/nextflow/trim_galore/main.nf'
include { consensus_dedup } from params.rootDir + '/target/nextflow/workflows/consensus_dedup/main.nf'
include { raw_fastq_align_test } from params.rootDir + '/target/_test/nextflow/test_workflows/raw_fastq_align_test/main.nf'
include { trim_galore_parity_test } from params.rootDir + '/target/_test/nextflow/test_workflows/trim_galore_parity_test/main.nf'
include { consensus_dedup_test } from params.rootDir + '/target/_test/nextflow/test_workflows/consensus_dedup_test/main.nf'

params.resources_test = 'resources_test/'

def sha256(path) {
  def digest = java.security.MessageDigest.getInstance('SHA-256')
  path.withInputStream { stream ->
    byte[] buffer = new byte[1024 * 1024]
    int n
    while ((n = stream.read(buffer)) > 0) digest.update(buffer, 0, n)
  }
  digest.digest().encodeHex().toString()
}

process check_consensus_mates {
  label 'singlecpu'
  label 'lowmem'
  container 'quay.io/biocontainers/samtools:1.19.2--h50ea8bc_1'
  input:
  tuple val(id), path(bam), path(bai)
  output:
  val(id)
  script:
  """
  set -eo pipefail
  samtools view '${bam}' | awk -F '\\t' '
    {
      first = int(\$2 / 64) % 2; second = int(\$2 / 128) % 2
      if (first + second != 1) exit 1
      if (first) r1[\$1]++; else r2[\$1]++
      names[\$1] = 1
    }
    END {
      for (name in names) if (r1[name] != 1 || r2[name] != 1) exit 1
    }'
  samtools view -c '${bam}' chr14:23661200-24661200 > indexed_count.txt
  test "\$(cat indexed_count.txt)" -gt 0
  """
}

workflow test_wf {
  def resourcesRoot = file(params.resources_test).toAbsolutePath().normalize()
  def validation = new ProcessBuilder([
    'python3', file(params.rootDir).resolve('resources_test_scripts/share/fixture_integrity.py').toString(),
    'validate', 'raw_fastq_align', '--resources-dir', resourcesRoot.toString()
  ]).inheritIO().start()
  assert validation.waitFor() == 0 : "Raw fixture validation failed in $resourcesRoot"

  res = resourcesRoot.resolve('raw_fastq_align')
  whitelist = file(params.rootDir).resolve('assets/barcodes.txt')
  def expected = ['H0O260120CPL0071F1S_001', 'H0O260120CPL0071F1S_002']
  assert whitelist.exists() && sha256(whitelist) == 'ecc4d583077e5d3b12d0bef70059820bf2579598c4177e8b1f926ca7e9ffee37'
  samples = res.resolve('samples.csv').splitCsv(header: true)
  assert samples.collect { it.id } == expected : 'Expected exactly the two ordered regional real samples'
  def runs = samples.collect { row ->
    assert row.keySet() == ['id', 'read1', 'read2', 'cell_list'].toSet()
    assert row.read1 == "${row.id}.R1.fastq.gz" && row.read2 == "${row.id}.R2.fastq.gz"
    assert row.cell_list == "../common/cell_lists/${row.id}.cell_list.tsv"
    def r = [id: row.id, read1: res.resolve(row.read1), read2: res.resolve(row.read2),
             cell_list: res.resolve(row.cell_list).normalize(),
             genome: resourcesRoot.resolve('common/genome'), whitelist: whitelist]
    ['read1', 'read2', 'cell_list'].each { assert r[it].exists() : "Missing ${r[it]}" }
    r
  }

  checked = Channel.fromList(runs)
    | map { r -> [r.id, r] }
    | raw_fastq_align.run(
        fromState: { id, state -> [id: id, read1: state.read1, read2: state.read2, whitelist: state.whitelist, genome: state.genome] },
        toState: { id, output, state -> state + output + [first_pass_bam: output.output_bam] }
      )
    | raw_fastq_align_test.run(
        fromState: { id, state -> [id: id, read1: state.read1, read2: state.read2] + state.findAll { it.key.startsWith('output_') } },
        toState: { id, output, state -> state }
      )
    | fditools.run(
        key: 'fditools_parity',
        fromState: [read1: 'read1', read2: 'read2', bc_whitelist: 'whitelist'],
        args: [output_read1: 'parity_extracted_R1.fastq.gz', output_read2: 'parity_extracted_R2.fastq.gz'],
        toState: [parity_extracted_read1: 'output_read1', parity_extracted_read2: 'output_read2']
      )
    | trim_galore.run(
        key: 'trim_galore_parity',
        fromState: { id, state -> [input: [state.parity_extracted_read1, state.parity_extracted_read2]] },
        args: [output_read1: 'parity_trimmed_R1.fastq.gz', output_read2: 'parity_trimmed_R2.fastq.gz',
               fastqc_html: 'parity_fastqc.*.html', fastqc_zip: 'parity_fastqc.*.zip'],
        toState: [parity_rust_read1: 'output_read1', parity_rust_read2: 'output_read2']
      )
    | trim_galore_parity_test.run(
        fromState: { id, state -> [id: id, extracted_read1: state.parity_extracted_read1,
          extracted_read2: state.parity_extracted_read2, rust_read1: state.parity_rust_read1, rust_read2: state.parity_rust_read2] },
        args: [output: 'trimming_parity.json'],
        toState: { id, output, state -> state + [parity_evidence: output.output] }
      )
    | consensus_dedup.run(
        fromState: { id, state -> [id: id, input: state.first_pass_bam, genome: state.genome, cell_list: state.cell_list] },
        toState: { id, output, state -> state + [consensus: output, input_bam: state.first_pass_bam, input_cell_list: state.cell_list] }
      )
    | consensus_dedup_test.run(
        fromState: { id, state -> [id: id, output_bam: state.consensus.output_bam, output_bai: state.consensus.output_bai,
          bismark_report: state.consensus.output_bismark_report, family_size_histogram: state.consensus.output_family_size_histogram,
          grouping_metrics: state.consensus.output_grouping_metrics, simplex_stats: state.consensus.output_simplex_stats,
          cell_list: state.consensus.output_cell_list, input_bam: state.input_bam, input_cell_list: state.input_cell_list] },
        toState: { id, output, state -> state }
      )

  checked
    | map { id, state -> [id, state.consensus.output_bam, state.consensus.output_bai] }
    | check_consensus_mates
    | toSortedList()
    | map { ids -> assert ids == expected : "Expected exactly two completed sample events, got $ids" }
}
