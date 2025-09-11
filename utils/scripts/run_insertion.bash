#!/usr/bin/env bash

#!/usr/bin/env bash

CSV_FILE="$1"
PROGRAM="./build/bin/moist-insert" # assume we run in root dir
GRID_FACTORS=(0.5 1.0 1.5 2.0 2.5)

if [[ ! -f "$CSV_FILE" ]]; then
    echo "CSV file not found: $CSV_FILE"
    exit 1
fi

# Skip header line, read CSV
tail -n +2 "$CSV_FILE" | while IFS=',' read -r a b epsilon extent output_a output_b output_metrics output_metrics_tc_name hmin hmax
do
    for FACTOR in "${GRID_FACTORS[@]}"; do
        echo "Running with grid_factor=$FACTOR"

        "$PROGRAM" \
            -a "$a" \
            -b "$b" \
            --epsilon "$epsilon" \
            --extent "$extent" \
            --output-a "$output_a" \
            --output-b "$output_b" \
            --output-metrics "$output_metrics" \
            --output-metrics-tc-name "${output_metrics_tc_name}" \
            --hmin "$hmin" \
            --hmax "$hmax" \
            --grid-factor "$FACTOR"

        echo "-----------------------------------------"
    done
done
