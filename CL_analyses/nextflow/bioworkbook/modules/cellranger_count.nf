process cellranger_count {

    input:
    tuple val(species), file("${species}_cellranger_reference"), val(sample)

    output:
    //tuple val(species), val(sample), file("${species}_${sample}.bam"), file("${species}_cellranger_reference")
    tuple val(species), val(sample), file("${species}_cellranger_reference"), file("${sample}_${species}.bam")

    script:
    """
    #!/bin/bash
    touch ${sample}_${species}.bam
    """

}
