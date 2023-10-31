process mod0_paml {

    cpus { 1 * task.attempt }
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'ignore' }
    maxRetries 6
    memory { 4.GB * task.attempt }

    publishDir 'check_trees', mode: 'copy', overwrite: true, pattern: '*tree.txt'

    input:
    tuple file("${og}_codon.nogaps.phy"), val(og)

    output:
    tuple file("${og}_swamp_analysis"), val(og)

    script:
    """
    #!/bin/bash

    cp ${baseDir}/data/PAML_CTLs/mod0.ctl .
    cp ${baseDir}/data/paml_tree.txt .

    sed -i 's/ALN/${og}_codon.nogaps.phy/g' mod0.ctl
    sed -i 's/TREE/paml_tree.txt/g' mod0.ctl
    sed -i 's/OUT/${og}_mod0.txt/g' mod0.ctl

    ${baseDir}/software/paml-4.10.7/bin/codeml mod0.ctl


    mkdir ${og}_swamp_analysis
    cp 2NG* ${og}_mod0.txt ${og}_codon.nogaps.phy rst rst1 rub lnf 4fold.nuc ${og}_swamp_analysis
    

    """
}
