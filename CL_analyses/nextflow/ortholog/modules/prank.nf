process prank {
    cpus = 8
    memory = '64 GB'
    time = '2h'

    label 'prank'

    input:
    tuple file("${og}_fastas"), val(og)

    output:


    script:
    """
    #!/bin/bash

    cp ${og}_fastas/* .
    prank -d=pro_species.fa -t=$params.tree -o=output_file -F -showxml
    prank -d=cds_species.fa -t=$params.tree -o=output_cds_codon -F -showxml
    prank -d=cds_species.fa -t=$params.tree -codon -F -o=output_cds_codon
    prank -d=cds_species.fa -o=output_translated -translate -F

    """
}
