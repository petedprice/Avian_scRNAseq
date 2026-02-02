process thetas {

    label 'angsd'
    errorStrategy 'retry'

    cpus { 16 * task.attempt }
    errorStrategy 'retry'
    memory { 64.GB * task.attempt }


    input:
    file(1dsfs)

    output:

    script:
    """
    #!/bin/bash
    # fold sfs
    angsd realSFS $SFS/${1}.saf.idx -P 16 -fold 1 > ${1}.folded.sfs
    angsd realSFS saf2theta $SFS/${1}.saf.idx -outname ${1} -sfs ${1}.folded.sfs -fold 1
    #Estimate thetas for every Chromosome/scaffold
    
    angsd thetaStat do_stat ${1}.thetas.idx
    #Do a sliding window analysis
    angsd thetaStat do_stat ${1}.thetas.idx -win 100000 -step 10000  -outnames ${1}_theta.thetasWindow.gz


    """
}

