process idsfs {

    label 'angsd'
    errorStrategy 'retry'

    cpus { 16 * task.attempt }
    errorStrategy 'retry'
    memory { 64.GB * task.attempt }


    input:
    tuple val(contig),  file(bam), val(samples)

    output:
    tuple val(contig), file("${contig}_idsfs_out")

    script:
    """
    #!/bin/bash
    echo \$(ls -1 *bam) > bam_list

    NSAMPLES=\$(cat bam_list | wc -l)
    MIN_SAMPLES=\$(((NSAMPLES*70)/100))

    angsd -dosaf 1 \
	-GL 2 \
	-ref ${params.ref} -anc ${params.ref} \
 	-b bam_list -out 1dsfs \
  	-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 \
        -trim 0 -baq 1 \
        -minMapQ 20 -minQ 20 -minInd \$MIN_SAMPLES \
        -setMinDepth 5 -setMaxDepth 2000 -doCounts 1
   
    mkdir ${contig}_idsfs_out
    mv 1dsfs* ${contig}_idsfs_out

    """
}


