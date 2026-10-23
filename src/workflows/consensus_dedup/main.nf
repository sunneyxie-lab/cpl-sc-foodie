workflow run_wf {
  take: input_ch
  main:

  output_ch = input_ch
    | map { id, state ->
        assert id ==~ /[^@\s]\S*/ : "Sample ID '$id' must have no whitespace and no leading '@': it becomes the consensus read name prefix"
        [id, state]
      }
    // fgumi group keeps cells apart by the CB tag, but the input has the barcode
    // only as a read name suffix.
    | tag_cell_barcode.run(
        fromState: ["input": "input"],
        toState: ["cb_bam": "output"]
      )
    // The input is grouped by name, but its header says SO:unsorted and Bismark's
    // name order may differ from samtools' natural order, which fixmate expects.
    | samtools_sort.run(
        key: "samtools_sort_queryname",
        fromState: ["input": "cb_bam"],
        args: [natural_sort: true],
        toState: ["qname_bam": "output"]
      )
    | samtools_fixmate.run(
        fromState: ["input": "qname_bam"],
        toState: ["fixmate_bam": "output"]
      )
    | fgumi_sort.run(
        fromState: ["input": "fixmate_bam"],
        args: [order: "template-coordinate"],
        toState: ["tc_bam": "output"]
      )
    | fgumi_group.run(
        fromState: ["input": "tc_bam"],
        args: [
          strategy: "identity",
          no_umi: true,
          family_size_histogram: "family_sizes.txt",
          grouping_metrics: "grouping_metrics.txt",
        ],
        toState: [
          "grouped_bam": "output",
          "family_size_histogram": "family_size_histogram",
          "grouping_metrics": "grouping_metrics",
        ]
      )
    // A passed-in cell_list takes the state key call_cells would write.
    | call_cells.run(
        runIf: { id, state -> !state.cell_list },
        fromState: { id, state -> [
          input: state.grouped_bam,
          id: id,
          drop_lowest: state.drop_lowest,
          max_cells: state.max_cells,
          fallback_cells: state.fallback_cells,
        ]},
        toState: [
          "cell_list": "output_cell_list",
          "cell_counts": "output_cell_counts",
          "cell_qc": "output_cell_qc",
          "knee_plot": "output_knee_plot",
        ]
      )
    | filter_cells.run(
        fromState: ["input": "grouped_bam", "cell_list": "cell_list"],
        toState: ["cells_bam": "output"]
      )
    | fgumi_simplex.run(
        fromState: { id, state -> [
          input: state.cells_bam,
          min_reads: state.min_reads,
          read_name_prefix: id,
        ]},
        args: [stats: "simplex_stats.txt"],
        toState: ["consensus_bam": "output", "simplex_stats": "stats"]
      )
    // output_reads must stay null: viash would otherwise give it a default path,
    // and its -o overrides -1/-2 so no R1/R2 get written.
    | samtools_fastq.run(
        fromState: ["input": "consensus_bam"],
        args: [
          no_suffix: true,
          read1: "consensus_R1.fq.gz",
          read2: "consensus_R2.fq.gz",
          output_reads_both: "consensus_other.fq.gz",
          output: "consensus_stdout.fq",
          output_reads: null,
        ],
        toState: ["fastq_r1": "read1", "fastq_r2": "read2"]
      )
    | bismark_align.run(
        fromState: [
          "genome": "genome",
          "mates1": "fastq_r1",
          "mates2": "fastq_r2",
          "maxins": "maxins",
        ],
        toState: ["bismark_bam": "output", "bismark_report": "report"]
      )
    | restore_consensus_tags.run(
        fromState: ["input": "bismark_bam", "consensus": "consensus_bam"],
        toState: ["tagged_bam": "output"]
      )
    | samtools_sort.run(
        key: "samtools_sort_coordinate",
        fromState: ["input": "tagged_bam"],
        toState: ["sorted_bam": "output"]
      )
    | samtools_index.run(
        fromState: { id, state -> [
          input: state.sorted_bam,
          output: state.sorted_bam.name + ".bai",
        ]},
        toState: ["bai": "output"]
      )
    | setState([
        "output_bam": "sorted_bam",
        "output_bai": "bai",
        "output_bismark_report": "bismark_report",
        "output_family_size_histogram": "family_size_histogram",
        "output_grouping_metrics": "grouping_metrics",
        "output_simplex_stats": "simplex_stats",
        "output_cell_list": "cell_list",
        "output_cell_counts": "cell_counts",
        "output_cell_qc": "cell_qc",
        "output_knee_plot": "knee_plot",
      ])

  emit:
  output_ch
}
