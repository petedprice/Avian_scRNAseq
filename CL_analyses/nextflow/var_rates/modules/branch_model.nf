process branch_model {

    cpus { 1 * task.attempt }
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'ignore' }
    maxRetries 6
    memory { 4.GB * task.attempt }

    tag {'branch_models' + '_' + og }

    publishDir 'pergene_branch_model_summaries', mode: 'copy', overwrite: true, pattern: '*_paml_branch.txt'

    label 'R'

    input:
    tuple file(phy), val(og), file('paml_branch_trees')
    
    output:
    tuple file("*.txt"), env(sc)
 
    script:
    """
    #!/bin/bash
    cp paml_branch_trees/* .
    cp ${baseDir}/data/PAML_CTLs/branch.ctl .

    sed -i 's/ALN/$phy/g' branch.ctl
    sed -i 's/TREE/tree_paml.txt/g' branch.ctl
    out="\$(basename /$phy .phy)"

    sed -i 's/OUT/paml_branch.txt/g' branch.ctl

    ${baseDir}/software/paml4.8/bin/codeml branch.ctl
    mv paml_branch.txt \${out}_paml_branch.txt

    sc=\${out#"${og}_codon.nogaps_"}


    """
}

