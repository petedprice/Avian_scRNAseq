process split_bam {
    //conda './envs/mut_id.yaml'


    input:
    tuple val(species), val(sample), file("sample/outs/possorted_genome_bam.bam"), file("sample/outs/possorted_genome_bam.bam.bai"), file("contig.txt")
    output:
    tuple val(species), val(sample), env(contig), file("subset.bam"), file("subset.bam.bai")

    script:
    """
    #!/bin/bash
    echo 'dog'
    contig=\$(cat contig.txt)
    samtools view -bh ${sample}_${species}.bam /$contig > subset.bam
    samtools index subset.bam
    """
}
