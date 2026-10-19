#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

# Keep pipes at line ends: the Nextflow runner strips a line up to a leading "|".
# The awk script must stay POSIX: awk in the samtools image is BusyBox.
samtools view -h "$par_input" |
  awk -F'\t' -v OFS='\t' '
      /^@/ { print; next }
      {
        cut = 0; rest = $1
        while ((k = index(rest, "_CB:Z:")) > 0) { cut += k; rest = substr(rest, k + 1) }
        bc = substr($1, cut + 6)
        if (cut == 0 || bc == "") {
          print "tag_cell_barcode: read name " $1 " has no _CB:Z:<barcode> suffix" > "/dev/stderr"
          exit 1
        }
        for (i = 12; i <= NF; i++) {
          if (substr($i, 1, 5) == "CB:Z:") {
            print "tag_cell_barcode: read " $1 " already has a CB tag" > "/dev/stderr"
            exit 1
          }
        }
        print $0 "\tCB:Z:" bc
      }
    ' |
  samtools view -b ${meta_cpus:+-@ "$meta_cpus"} -o "$par_output" -
