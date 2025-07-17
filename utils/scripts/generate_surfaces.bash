#!/bin/bash

# Default values
input="./data/{}.tif"
output="./surfaces/{}-{}.off"
first_layer=0
sample_bits=16
isovalue=0.999

num_layers=""
slices=""

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
    -first-layer)
      first_layer="$2"
      shift 2
      ;;
    -num-layers)
      num_layers="$2"
      shift 2
      ;;
    -slices)
      slices="$2"
      shift 2
      ;;
    -sample-bits)
      sample_bits="$2"
      shift 2
      ;;
    -isovalue)
      isovalue="$2"
      shift 2
      ;;
    *)
      echo "Unknown parameter: $1"
      exit 1
      ;;
  esac
done

if [[ -z "$num_layers" ]]; then
  echo "-num-layers is required"
  exit 1
fi

if [[ -z "$slices" ]]; then
  echo "-slices is required"
  exit 1
fi

if (( num_layers % slices != 0 )); then
  echo "-num-layers must be divisible by -slices"
  exit 1
fi

range=$(( num_layers / slices ))

for (( i=0; i < slices; i++ )); do


  if (( i != 0 )); then
    i_first=$(( first_layer + i * range - i ))
  else
    i_first=$(( first_layer + i * range ))
  fi

  echo "Generating slice $(( i_first )): $current_first -> $((i_first + range))"

  # Just assume the program is run in the project source directory....
  ./build/bin/moist-extract \
    --input "$input" \
    --output "$output" \
    --first "$i_first" \
    --amount "$range" \
    --isovalue "$isovalue" \
    --directional-offset "-$first_layer" \
    --invert
done
