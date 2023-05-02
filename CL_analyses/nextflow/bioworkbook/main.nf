species_ch=channel
    .fromPath(params.metadata)
    .splitCsv()
    .map {row -> tuple(row[1], row[2])}
    .unique()
    .view()



process cellranger_mkref {

    input:
    tuple val(species), val(ref) from species_ch 

    output:
    //publishDir "cellranger_reference"
    tuple file("${species}_cellranger_reference"), val(species) into ref_made


    script:
    """
    #!/bin/bash
    #$params.cellranger mkref --genome=output_genome \
	--fasta=$params.fasta_dir/${ref}.fasta \
	--genes=$params.gtf_dir/${ref}.gtf
    
    touch ${species}_cellranger_reference

	
    """
}


samples = channel
    .fromPath(params.metadata)
    .splitCsv()
    .map {row -> tuple(row[0], row[1])}
    .view()


process cellranger_count {

    input: 
    tuple val(species), file("${species}_cellranger_reference"), val(sample)  from ref_made.combine(samples, by:1)

    output:

    script: 
    """
    #!/bin/bash
    #$params.cellranger count --id=$sample \
	--fastqs=$params.read_dir \ 
	--transcriptome=${species}_cellranger_reference 
    
    """
}




