#!/usr/bin/env nextflow

// Enable DSL2 syntax
nextflow.enable.dsl=2

// Define the first process to create a file
process CreateFileA {

    output:
    path "fileA.txt"

    script:
    """
    which sbatch
    touch fileA.txt
    """
}

// Define the second process to use the output from the first process
process CreateFileB {

    input:
    path fileA

    output:
    path "fileB.txt"

    script:
    """
    sleep 10; cp ${fileA} fileB.txt
    """
}

// Define the workflow execution
workflow {

    // Call the first process to create fileA
    fileA_out = CreateFileA()

    CreateFileB(fileA_out)

    // Use the output of the first process as input for the second process
    //CreateFileB('/home/ld32/data/nextflow/fileA.txt')
}
