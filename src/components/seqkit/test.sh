#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail
source "$meta_resources_dir/test_helpers.sh"
setup_test_env

printf '@a\nACGT\n+\nIIII\n@b\nGGCCAA\n+\nIIIIII\n' | gzip > "$meta_temp_dir/input.1.fastq.gz"
printf '@c\nATATATAT\n+\n!!!!!!!!\n' | gzip > "$meta_temp_dir/input.2.fastq.gz"
printf '@d\nGGCC\n+\nIIII\n@e\nGC\n+\nII\n@f\nAC\n+\nII\n' | gzip > "$meta_temp_dir/input.3.fastq.gz"
printf '@g\nAAAAAA\n+\n!!!!!!\n@h\nTT\n+\n!!\n' | gzip > "$meta_temp_dir/input.4.fastq.gz"

# Reverse the first pair to expose accidental sorting.
"$meta_executable" --input "$meta_temp_dir/input.2.fastq.gz;$meta_temp_dir/input.1.fastq.gz" -o "$meta_temp_dir/two.tsv"
"$meta_executable" --input "$meta_temp_dir/input.1.fastq.gz;$meta_temp_dir/input.2.fastq.gz;$meta_temp_dir/input.3.fastq.gz;$meta_temp_dir/input.4.fastq.gz" --out_file "$meta_temp_dir/four.tsv"

for table in two four; do
  awk -F '\t' -v table="$table" '
    function near(value, expected) { return value - expected < 0.051 && expected - value < 0.051 }
    NR == 1 {
      for (i = 1; i <= NF; i++) { if ($i in col) exit 1; col[$i] = i }
      n[1]=2; n[2]=1; n[3]=3; n[4]=2
      total[1]=10; total[2]=8; total[3]=8; total[4]=8
      lo[1]=4; lo[2]=8; lo[3]=2; lo[4]=2
      hi[1]=6; hi[2]=8; hi[3]=4; hi[4]=6
      gc[1]=60; gc[2]=0; gc[3]=87.5; gc[4]=0
      split("file format type num_seqs sum_len min_len max_len avg_len GC(%) Q20(%) Q30(%) Q1 Q2 Q3 N50 N50_num AvgQual", fields, " ")
      for (i in fields) if (!(fields[i] in col)) exit 1
      next
    }
    {
      row=NR-1; id=(table == "two" ? 3-row : row)
      if (id < 1 || id > 4 || (table == "two" && row > 2)) exit 1
      if ($(col["file"]) !~ ("/input\\." id "\\.fastq\\.gz$")) exit 1
      if ($(col["format"]) != "FASTQ" || $(col["type"]) != "DNA") exit 1
      if ($(col["num_seqs"])+0 != n[id] || $(col["sum_len"])+0 != total[id]) exit 1
      if ($(col["min_len"])+0 != lo[id] || $(col["max_len"])+0 != hi[id]) exit 1
      if (!near($(col["avg_len"])+0, total[id]/n[id]) || !near($(col["GC(%)"])+0, gc[id])) exit 1
      quality=(id == 1 || id == 3 ? 100 : 0)
      if (!near($(col["Q20(%)"])+0,quality) || !near($(col["Q30(%)"])+0,quality)) exit 1
      if (quality && $(col["AvgQual"])+0 <= 30) exit 1
      if (!quality && $(col["AvgQual"])+0 != 0) exit 1
      for (i=12; i<=17; i++) if ($(col[fields[i]]) !~ /^[0-9]+([.][0-9]+)?$/) exit 1
    }
    END { if (NR != (table == "two" ? 3 : 5)) exit 1 }
  ' "$meta_temp_dir/$table.tsv"
done
print_test_summary "seqkit tests"
