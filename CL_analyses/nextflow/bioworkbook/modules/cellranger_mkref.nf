process cellranger_mkref {
    queue = "ressexcon.q"
    cpus = 16
    memory = '8 GB'
    time = '4h'
    clusterOptions = { '-P ressexcon' }


    input:
    tuple val(species), val(ref)
    output:
    tuple val(species), file("${species}_cellranger_reference")
    script:
    """
    #!/bin/bash
    #touch ${species}_cellranger_reference
    #echo contig1 > ${species}_cellranger_reference
    #echo contig2 >> ${species}_cellranger_reference
    $params.cellranger mkref \
	--genome=${species}_cellranger_reference \
	--fasta=${params.fasta_dir}/${ref}.fna \
	--genes=${params.gtf_dir}/${ref}.gtf
    """
}
