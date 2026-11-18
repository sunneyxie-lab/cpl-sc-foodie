#!/bin/bash

## VIASH START
## VIASH END

set -eo pipefail

IFS=';' read -r -a inputs <<< "$par_input"
if [[ ${#inputs[@]} -ne 2 ]]; then
  echo 'Expected exactly two input FASTQs in R1, R2 order.' >&2
  exit 1
fi

for argument in fastqc_html fastqc_zip; do
  variable="par_$argument"
  template="${!variable}"
  if [[ "$template" != *'*'* || "${template#*\*}" == *'*'* || "$template" == *';'* ]]; then
    echo "--$argument must be one output template containing exactly one '*' and no semicolons." >&2
    exit 1
  fi
done

out=$(mktemp -d "$meta_temp_dir/trim_galore.XXXXXX")
cores=1
if [[ ${meta_cpus:-1} -ge 6 ]]; then
  cores=$(( meta_cpus - 4 ))
fi
preset=$(python3 "$meta_resources_dir/select_legacy_adapter.py" "${inputs[0]}")
trim_galore --version
cmd_args=(
  trim_galore --paired -q 20 --stringency 2 -e 0.1 --no_poly_g
  "$preset" --fastqc --basename trimmed --output_dir "$out"
  --cores "$cores" "${inputs[@]}"
)
printf 'Running:'
printf ' %q' "${cmd_args[@]}"
printf '\n'
"${cmd_args[@]}"

# Check all task-local products before relocating any declared output.
for mate in 1 2; do
  gzip -t "$out/trimmed_val_${mate}.fq.gz"
  test -s "$out/trimmed_val_${mate}_fastqc.html"
  test -s "$out/trimmed_val_${mate}_fastqc.zip"
done
output_paths=("$par_output_read1" "$par_output_read2")
for index in 0 1; do
  output_paths+=("${par_fastqc_html/\*/$index}" "${par_fastqc_zip/\*/$index}")
done
for path in "${output_paths[@]}"; do
  mkdir -p -- "$(dirname -- "$path")"
done
mv "$out/trimmed_val_1.fq.gz" "$par_output_read1"
mv "$out/trimmed_val_2.fq.gz" "$par_output_read2"
for index in 0 1; do
  mate=$(( index + 1 ))
  mv "$out/trimmed_val_${mate}_fastqc.html" "${par_fastqc_html/\*/$index}"
  mv "$out/trimmed_val_${mate}_fastqc.zip" "${par_fastqc_zip/\*/$index}"
done
