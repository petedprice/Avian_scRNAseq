nextflow.enable.dsl=2

process cellranger_mkref {
    input:
    tuple val(species), val(ref)
    output:
    tuple val(species), file("${species}_cellranger_reference")
    script:
    """
    #!/bin/bash
    touch ${species}_cellranger_reference
    echo contig1 > ${species}_cellranger_reference
    echo contig2 >> ${species}_cellranger_reference
    """
}


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




workflow {
    species_ch=Channel
	.fromPath(params.metadata)
	.splitCsv()
	.map {row -> tuple(row[1], row[2])}
	.unique()
   
    samples_ch = Channel
	.fromPath(params.metadata)
	.splitCsv()
	.map {row -> tuple(row[1], row[0])}

    ref_made=cellranger_mkref(species_ch)

    counted=cellranger_count(ref_made.combine(samples_ch, by: 0))

    cns=contig_names(counted)
	.transpose()

    splitted=split_bam(cns).view()
    mut_ided = mut_id(splitted)

}




