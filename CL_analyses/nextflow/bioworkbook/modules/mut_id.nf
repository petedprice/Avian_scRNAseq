process mut_id {

    input:
    tuple val(species), val(sample), val(contig), file("subset.bam")

    output:

    script:
    """
    #!/bin/bash
    mv subset.bam ${contig}.bam
    """
}
