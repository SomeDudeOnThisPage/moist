#!/bin/bash

epsilon="1e-3"
input="./surfaces/"
output="./volumes/"
csv="trace.csv"

while [[ $# -gt 0 ]]; do
  case $1 in
    -input)
      input="$2"
      shift 2
      ;;
    -output)
      output="$2"
      shift 2
      ;;
    -csv)
      csv="$2"
      shift 2
      ;;
    -epsilon)
      epsilon="$2"
      shift 2
      ;;
    *)
      echo "Unknown parameter: $1"
      exit 1
      ;;
  esac
done


if [[ -z "$input" || ! -d "$output" ]]; then
  echo "Error: '$input' or '$output' is not a directory"
  exit 1
fi

# probably quiet output of fTetWild
for file in "$input"/*; do
  if [[ -f "$file" ]]; then
    filename=$(basename "$file")
    echo "Running fTetWild for $file with eps $epsilon, saving to ${output%/}/${filename%.*}.msh"

    filesize_input_bytes=$(stat -c%s "$file")
    filesize_input_mb=$(awk "BEGIN {printf \"%.5f\", $filesize_input_bytes/1024/1024}")

    # track mem usage (max resident size) and time
    /usr/bin/time -f "%M %e" -o time_stats.txt \
        fTetWild --input "$file" --epsr "$epsilon" \
                --output "${output%/}/${filename%.*}.msh" \
                --no-binary --no-color --max-its "25" \
        >fTetWild.log 2>fTetWild_errors.log

    time_output=$(<time_stats.txt)
    rm time_stats.txt

    ru_maxrss_kb=$(echo "$time_output" | awk '{print $1}')
    runtime_s=$(echo "$time_output" | awk '{print $2}')
    ru_maxrss_mb=$(awk "BEGIN {printf \"%.5f\", $ru_maxrss_kb/1024}")

    # check output filesize
    msh_file="${output%/}/${filename%.*}.msh"
    filesize_bytes=$(stat -c%s "$msh_file")
    filesize_mb=$(awk "BEGIN {printf \"%.5f\", $filesize_bytes/1024/1024}")

    if [[ ! -f "$csv" ]]; then
        echo "$csv"
        echo "testcase;ru_maxrss_mb;runtime_s;input_filesize_mb;output_filesize_mb" >> "$csv"
    fi
    echo "$filename;$ru_maxrss_mb;$runtime_s;$filesize_input_mb;$filesize_mb" >> "$csv"
  fi
done

# we don't need all the supplementary fTetWild files like surface, tracked surface, and csv values.
find "${output%/}" -type f ! -name '*.msh' -delete
