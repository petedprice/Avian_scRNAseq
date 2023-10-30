process m1avsm2a {

    cpus { 1 * task.attempt }
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'ignore' }
    maxRetries 6
    memory { 4.GB * task.attempt }

    label 'python'

    input:
    tuple file(phy), val(og)
    
    output:
    tuple file("*_paml_mod1a2a.txt"), env(sc)
 
    script:
    """
    #!/bin/bash

    cp ${baseDir}/data/PAML_CTLs/mod1a2a.ctl .
    cp ${baseDir}/data/paml_tree.txt .


    sed -i 's/ALN/$phy/g' mod1a2a.ctl
    sed -i 's/TREE/paml_tree.txt/g' mod1a2a.ctl
    out="\$(basename /$phy .phy)"

    sed -i 's/OUT/paml_mod1a2a.txt/g' mod1a2a.ctl

    ${baseDir}/software/paml-4.10.7/bin/codeml mod1a2a.ctl
    mv paml_mod1a2a.txt \${out}_paml_mod1a2a.txt

    sc=\${out#"${og}_codon.nogaps_"}


    """
}

