#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

barcodes="$meta_temp_dir/barcodes.tsv"

# Keep pipes at line ends: the Nextflow runner strips a line up to a leading "|".
samtools view "$par_consensus" |
  awk -F'\t' -v OFS='\t' '
      {
        for (i = 12; i <= NF; i++) {
          if (substr($i, 1, 5) == "CB:Z:") { print $1, substr($i, 6); break }
        }
      }
    ' > "$barcodes"

# The awk scripts here must stay POSIX: awk in the samtools image is BusyBox.
samtools view -h "$par_input" |
  awk -F'\t' -v OFS='\t' -v map="$barcodes" '
      BEGIN {
        while ((getline line < map) > 0) {
          split(line, kv, "\t")
          bc[kv[1]] = kv[2]
        }
      }
      /^@/ { print; next }
      {
        if (!($1 in bc)) {
          print "restore_cell_barcode: no cell barcode for read " $1 " in the consensus BAM" > "/dev/stderr"
          exit 1
        }
        b = bc[$1]
        $1 = $1 "_CB:Z:" b
        print $0 "\tCB:Z:" b
      }
    ' |
  samtools view -b ${meta_cpus:+-@ "$meta_cpus"} -o "$par_output" -
