process cellranger_count {

    input:
    tuple val(species), file("${species}_cellranger_reference"), val(sample)

    output:
    tuple val(species), val(sample), file("sample/outs/possorted_genome_bam.bam"), file("sample/outs/possorted_genome_bam.bam.bai")

    script:
    """
    #!/bin/bash
    $params.cellranger count \
	--id=$sample \
   	--fastqs=${params.read_dir}/$sample \
   	--transcriptome=${species}_cellranger_reference
    """

}