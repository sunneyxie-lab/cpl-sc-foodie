#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

tags="$meta_temp_dir/tags.tsv"

# Mates of one pair can have different cD, so the key is the read name plus 1 (0x40) or 2 (0x80).
# Keep pipes at line ends: the Nextflow runner strips a line up to a leading "|".
samtools view "$par_consensus" |
  awk -F'\t' -v OFS='\t' '
      {
        bc = ""; cd = ""
        for (i = 12; i <= NF; i++) {
          if (substr($i, 1, 5) == "CB:Z:") bc = substr($i, 6)
          else if (substr($i, 1, 5) == "cD:i:") cd = substr($i, 6)
        }
        if (bc == "" || cd == "") {
          print "restore_consensus_tags: consensus read " $1 " lacks a CB or cD tag" > "/dev/stderr"
          exit 1
        }
        print $1, (int($2 / 64) % 2 ? 1 : 2), bc, cd
      }
    ' > "$tags"

# The awk scripts here must stay POSIX: awk in the samtools image is BusyBox.
samtools view -h "$par_input" |
  awk -F'\t' -v OFS='\t' -v map="$tags" '
      BEGIN {
        while ((getline line < map) > 0) {
          split(line, kv, "\t")
          bc[kv[1], kv[2]] = kv[3]
          cd[kv[1], kv[2]] = kv[4]
        }
      }
      /^@/ { print; next }
      {
        k = $1 SUBSEP (int($2 / 64) % 2 ? 1 : 2)
        if (!(k in bc)) {
          print "restore_consensus_tags: no consensus mate for read " $1 " (flag " $2 ")" > "/dev/stderr"
          exit 1
        }
        b = bc[k]
        $1 = $1 "_CB:Z:" b
        print $0 "\tCB:Z:" b "\tcD:i:" cd[k]
      }
    ' |
  samtools view -b ${meta_cpus:+-@ "$meta_cpus"} -o "$par_output" -
