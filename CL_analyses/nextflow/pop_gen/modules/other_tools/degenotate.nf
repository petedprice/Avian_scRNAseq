process degenotate {

    label 'degen'
    errorStrategy 'retry'

    cpus { 4 * task.attempt }
    memory { 32.GB * task.attempt }

    publishDir 'degen', mode: 'copy', overwrite: true, pattern: 'tmpplaceholder_*'

    input:

    output:
    tuple file("zfold.bed"), file("ffold.bed")


    script:
    """
    #!/bin/bash
    python3 ${projectDir}/software/degenotate/degenotate.py -a ${params.gff} -d " " -g ${params.fasta} -o degen

    awk '\$5 == 0' degen/degeneracy-all-sites.bed | cut -f1,2,3 | uniq > zfold.bed
    awk '\$5 == 4' degen/degeneracy-all-sites.bed | cut -f1,2,3 | uniq > ffold.bed

    """
}

