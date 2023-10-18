process ortho_cds {
    cpus = 4
    memory = '64 GB'
    time = '8h'


    input:
    tuple file(proteins), file('Orthofinder_Results')

    output:
    tuple file('*_folder'), file('species_order.txt')


    script:
    """
    #!/bin/bash

    man grep | cat

    cat *cds_longest.fna > comp_longest.fna
    grep -F -f Orthofinder_Results/Orthogroups/Orthogroups_SingleCopyOrthologues.txt Orthofinder_Results/Orthogroups/Orthogroups.tsv > single_copy_orthogroups.tsv
    head -1 Orthofinder_Results/Orthogroups/Orthogroups.tsv | cut -f2-> species_order.txt
    for c in \$(cat Orthofinder_Results/Orthogroups/Orthogroups_SingleCopyOrthologues.txt | head)
    do
    grep \$c single_copy_orthogroups.tsv | cut -f2- | tr '\\t' '\\n' | dos2unix > tmp_orthos_name.txt
    grep --no-group-separator -A 1 -F -f tmp_orthos_name.txt comp_longest.fna > \${c}_cds.fa
    mkdir \${c}_folder

    mv \${c}_cds.fa \${c}_folder
    cp Orthofinder_Results/Orthogroup_Sequences/\${c}.fa \${c}_folder/\${c}_pro.fa


    done
    


    """
}
