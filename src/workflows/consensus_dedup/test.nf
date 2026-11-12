nextflow.enable.dsl=2

include { consensus_dedup } from params.rootDir + "/target/nextflow/workflows/consensus_dedup/main.nf"
include { consensus_dedup_test } from params.rootDir + "/target/_test/nextflow/test_workflows/consensus_dedup_test/main.nf"

params.resources_test = "resources_test/"

workflow test_wf {
  def resourcesRoot = file(params.resources_test).toAbsolutePath().normalize()
  def validation = new ProcessBuilder([
    'python3', file(params.rootDir).resolve('resources_test_scripts/share/fixture_integrity.py').toString(),
    'validate', 'consensus_dedup', '--resources-dir', resourcesRoot.toString()
  ]).inheritIO().start()
  assert validation.waitFor() == 0 : "Consensus fixture validation failed in $resourcesRoot"

  res = resourcesRoot.resolve("consensus_dedup")
  samples = res.resolve("samples.csv").splitCsv(header: true)
  def expected = ['H0O260120CPL0071F1S_001', 'H0O260120CPL0071F1S_002']
  assert samples.collect { it.id } == expected : "Expected exactly the two ordered regional samples in ${res.resolve('samples.csv')}"

  def runs = samples.collect { row ->
    assert row.keySet() == ['id', 'input', 'cell_list'].toSet()
    assert row.input == "${row.id}.bam"
    assert row.cell_list == "../common/cell_lists/${row.id}.cell_list.tsv"
    [id: row.id, input: res.resolve(row.input), genome: resourcesRoot.resolve("common/genome"),
     cell_list: res.resolve(row.cell_list).normalize()]
  }
  // The first sample once more without a cell list, so cell calling runs too.
  runs += [runs[0] + [id: runs[0].id + "_called", cell_list: null]]

  output_ch = Channel.fromList(runs)
    | map { r -> [r.id, r] }
    | consensus_dedup.run(
        fromState: { id, state -> [id: id, input: state.input, genome: state.genome, cell_list: state.cell_list].findAll { it.value != null } },
        toState: { id, output, state -> output + [input_bam: state.input, input_cell_list: state.cell_list] }
      )
    | view { output ->
        assert output.size() == 2 : "outputs should contain two elements: [id, state]"
        "Output: $output"
      }
    | consensus_dedup_test.run(
        fromState: { id, state -> [
          id: id,
          output_bam: state.output_bam,
          output_bai: state.output_bai,
          bismark_report: state.output_bismark_report,
          family_size_histogram: state.output_family_size_histogram,
          grouping_metrics: state.output_grouping_metrics,
          simplex_stats: state.output_simplex_stats,
          cell_list: state.output_cell_list,
          cell_counts: state.output_cell_counts,
          cell_qc: state.output_cell_qc,
          knee_plot: state.output_knee_plot,
          input_bam: state.input_bam,
          input_cell_list: state.input_cell_list,
        ].findAll { it.value != null } }
      )
    | toSortedList()
    | map { output_list ->
        assert output_list.size() == 3 : "output channel should contain three events"
      }
}
