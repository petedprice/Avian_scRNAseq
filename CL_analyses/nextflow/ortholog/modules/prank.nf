process prank {
    cpus { 1 * task.attempt }
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'ignore' }
    maxRetries 6
    memory { 4.GB * task.attempt }


    //label 'prank'
    conda 'prank'

    input:
    tuple file("${og}_fastas"), val(og)

    output:


    script:
    """
    #!/bin/bash
    conda --version
    cp ${og}_fastas/* .
    ulimit -c unlimited
    prank -d=pro_species.fa -t=$params.tree -o=output_file -F -showxml
    prank -d=cds_species.fa -t=$params.tree -o=output_cds_codon -F -showxml
    prank -d=cds_species.fa -t=$params.tree -codon -F -o=output_cds_codon
    prank -d=cds_species.fa -o=output_translated -translate -F

    """
}
