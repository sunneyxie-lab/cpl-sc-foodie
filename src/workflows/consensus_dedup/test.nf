nextflow.enable.dsl=2

include { consensus_dedup } from params.rootDir + "/target/nextflow/workflows/consensus_dedup/main.nf"
include { consensus_dedup_test } from params.rootDir + "/target/_test/nextflow/test_workflows/consensus_dedup_test/main.nf"

params.resources_test = "resources_test/"

workflow test_wf {
  res = file(params.resources_test).resolve("consensus_dedup")
  assert res.resolve("genome").exists() : "Test resources not found in $res; run resources_test_scripts/fetch_resources_test.sh first"
  samples = res.resolve("samples.csv").splitCsv(header: true)
  assert samples.size() == 2 : "expected 2 samples in ${res.resolve('samples.csv')}, got ${samples.size()}"

  output_ch = Channel.fromList(samples)
    | map { row -> [row.id, [id: row.id, input: res.resolve(row.input), genome: res.resolve("genome")]] }
    | consensus_dedup.run(
        fromState: ["id", "input", "genome"],
        toState: { id, output, state -> output + [input_bam: state.input] }
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
          input_bam: state.input_bam,
        ]}
      )
    | toSortedList()
    | map { output_list ->
        assert output_list.size() == 2 : "output channel should contain two events"
      }
}
