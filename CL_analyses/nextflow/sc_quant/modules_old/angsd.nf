process angsd {
    cache 'lenient'

    tag {'angsd_' + species + '_' + sample }

    input:
    tuple val(species), val(sample), file("sample/outs/possorted_genome_bam.bam"), file("sample/outs/possorted_genome_bam.bam.bai")??
	//need to include ref value by merging with metadata 
    
    output:

    script:
    """
    #!/bin/bash
    ls *.bam > bam.list

    ## FILTER BAM FILES ?

    ## CREATE SFS
    ${projectDir}/software/angsd/angsd \
	-bam bam.list \
	-anc $ref.fa \
	-doSaf 1 \
	-out smallFolded \
	-GL 2

    ## CALCULATE TAJIMAS D ETC 

    """

}
