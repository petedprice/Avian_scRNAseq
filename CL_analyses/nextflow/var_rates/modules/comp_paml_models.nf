process comp_paml_models {

    cpus { 1 * task.attempt }
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'ignore' }
    maxRetries 6
    memory { 4.GB * task.attempt }

    label 'R'

    input:
    tuple file(paml_outputs), val(sc)
    
    output:
    file("${sc}_model_comparison.txt")
 
    script:
    """
    #!/bin/bash
    Rscript ${baseDir}/scripts/comp_paml_models.R . 2 ${sc}_model_comparison.txt

    """
}
