#!/bin/bash

epsilon="1e-3"
input="./surfaces/"
output="./volumes/"

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


for file in "$input"/*; do
  if [[ -f "$file" ]]; then
    filename=$(basename "$file")
    echo "Running fTetWild for $file with eps $epsilon, saving to ${output%/}/${filename%.*}.msh"

    fTetWild --input "$file" --epsr "$epsilon" --output "${output%/}/${filename%.*}.msh" --no-binary --no-color --max-its "25"
  fi
done

# we don't need all the supplementary fTetWild files like surface, tracked surface, and csv values.
find "${output%/}" -type f ! -name '*.msh' -delete
