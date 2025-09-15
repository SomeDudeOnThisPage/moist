#!/usr/bin/env bash

CSV_FILE="$1"
PROGRAM="./build/bin/moist-insert" # assume we run in root dir

tail -n +2 "$CSV_FILE" | while IFS=',' read -r a b epsilon extent output_a output_b output_metrics output_metrics_tc_name hmin hmax
do
    echo "Running case: $output_metrics_tc_name"

    ./build/bin/moist-insert \
        -a "$a" \
        -b "$b" \
        --epsilon "$epsilon" \
        --extent "$extent" \
        --output-a "$output_a" \
        --output-b "$output_b" \
        --output-metrics "$output_metrics" \
        --output-metrics-tc-name "$output_metrics_tc_name" \
        --hmin "$hmin" \
        --hmax "$hmax" \
        --grid-factor 1
done
