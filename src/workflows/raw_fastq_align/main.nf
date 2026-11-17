def collectFastqcReports(directory) {
  def files = directory.listFiles().findAll { it.isFile() }
  def html = files.findAll { it.name.endsWith('.html') }
  def zip = files.findAll { it.name.endsWith('.zip') }
  assert html.size() == 1 && zip.size() == 1 : "Expected one top-level HTML and ZIP in $directory"
  [html: html[0], zip: zip[0]]
}

workflow run_wf {
  take: input_ch
  main:
  output_ch = input_ch
    | map { id, state ->
        assert id ==~ /[^@\s]\S*/ : "Sample ID '$id' must have no whitespace and no leading '@'"
        [id, state]
      }
    | fastqc.run(
        key: 'fastqc_raw_r1',
        fromState: [input: 'read1'], args: [outdir: 'raw_r1_qc'],
        toState: [raw_r1_qc_dir: 'outdir']
      )
    | fastqc.run(
        key: 'fastqc_raw_r2',
        fromState: [input: 'read2'], args: [outdir: 'raw_r2_qc'],
        toState: [raw_r2_qc_dir: 'outdir']
      )
    | map { id, state ->
        def r1 = collectFastqcReports(state.raw_r1_qc_dir)
        def r2 = collectFastqcReports(state.raw_r2_qc_dir)
        [id, state + [output_raw_fastqc_html: [r1.html, r2.html], output_raw_fastqc_zip: [r1.zip, r2.zip]]]
      }
    | seqkit.run(
        key: 'seqkit_raw',
        fromState: { id, state -> [input: [state.read1, state.read2]] },
        args: [out_file: 'raw.stats.tsv'], toState: [output_raw_stats: 'out_file']
      )
    | fditools.run(
        fromState: [read1: 'read1', read2: 'read2', bc_whitelist: 'whitelist'],
        args: [output_read1: 'extracted_R1.fastq.gz', output_read2: 'extracted_R2.fastq.gz'],
        toState: [extracted_read1: 'output_read1', extracted_read2: 'output_read2']
      )
    | fastqc.run(
        key: 'fastqc_extracted_r1',
        fromState: [input: 'extracted_read1'], args: [outdir: 'extracted_r1_qc'],
        toState: [extracted_r1_qc_dir: 'outdir']
      )
    | fastqc.run(
        key: 'fastqc_extracted_r2',
        fromState: [input: 'extracted_read2'], args: [outdir: 'extracted_r2_qc'],
        toState: [extracted_r2_qc_dir: 'outdir']
      )
    | map { id, state ->
        def r1 = collectFastqcReports(state.extracted_r1_qc_dir)
        def r2 = collectFastqcReports(state.extracted_r2_qc_dir)
        [id, state + [output_extracted_fastqc_html: [r1.html, r2.html], output_extracted_fastqc_zip: [r1.zip, r2.zip]]]
      }
    | trim_galore.run(
        fromState: { id, state -> [input: [state.extracted_read1, state.extracted_read2]] },
        args: [output_read1: 'trimmed_R1.fastq.gz', output_read2: 'trimmed_R2.fastq.gz',
               fastqc_html: 'trimmed_fastqc.*.html',
               fastqc_zip: 'trimmed_fastqc.*.zip'],
        toState: [trimmed_read1: 'output_read1', trimmed_read2: 'output_read2',
                  output_trimmed_fastqc_html: 'fastqc_html', output_trimmed_fastqc_zip: 'fastqc_zip']
      )
    | seqkit.run(
        key: 'seqkit_extracted_trimmed',
        fromState: { id, state -> [input: [state.extracted_read1, state.extracted_read2, state.trimmed_read1, state.trimmed_read2]] },
        args: [out_file: 'extracted_trimmed.stats.tsv'], toState: [output_extracted_trimmed_stats: 'out_file']
      )
    | bismark_align.run(
        key: 'bismark_first_pass',
        fromState: [genome: 'genome', mates1: 'trimmed_read1', mates2: 'trimmed_read2'],
        args: [non_directional: true, maxins: 2000],
        toState: [output_bam: 'output', output_bismark_report: 'report']
      )
    | setState([
        output_bam: 'output_bam', output_bismark_report: 'output_bismark_report',
        output_raw_stats: 'output_raw_stats', output_extracted_trimmed_stats: 'output_extracted_trimmed_stats',
        output_raw_fastqc_html: 'output_raw_fastqc_html', output_raw_fastqc_zip: 'output_raw_fastqc_zip',
        output_extracted_fastqc_html: 'output_extracted_fastqc_html', output_extracted_fastqc_zip: 'output_extracted_fastqc_zip',
        output_trimmed_fastqc_html: 'output_trimmed_fastqc_html', output_trimmed_fastqc_zip: 'output_trimmed_fastqc_zip'
      ])
  emit:
  output_ch
}
