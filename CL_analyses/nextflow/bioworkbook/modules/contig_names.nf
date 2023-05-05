process contig_names {
    input:
    tuple val(species), val(sample), file("${species}_cellranger_reference"), file("${sample}_${species}.bam")

    output:
    tuple val(species), val(sample), file("*.txt"), file("${sample}_${species}.bam")
    script:
    """
    #!/bin/bash
    mkdir contigs
    for contig in \$(cat ${species}_cellranger_reference)
    do
    echo \$contig > \${contig}_${species}.txt
    done
    """
}
