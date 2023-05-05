process split_bam {

    input:
    tuple val(species), val(sample), file("contig.txt"), file("${sample}_${species}.bam")
    output:
    tuple val(species), val(sample), env(contig), file("subset.bam")

    script:
    """
    #!/bin/bash
    echo
    contig=\$(cat contig.txt)
    touch subset.bam
    """
}
