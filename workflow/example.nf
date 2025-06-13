#!/usr/bin/env nextflow

params.datapath = null
params.slicesize = 50
params.outdir = 'results'

// "/mnt/c/Users/robin/Desktop/Master-Thesis/test-rocks/p4_res23_mu24/merge{:05}_s.tif"

process moist_extract {
    input:
        val _DATA_PATH
        path "b.txt"
    output: path "a.txt"

    script:
    """
    ../build/bin/moist-extract --input $_DATA_PATH -- output "--output", "./test/target/generated-surface/{}-{}.off"
    """
}

workflow {

}
