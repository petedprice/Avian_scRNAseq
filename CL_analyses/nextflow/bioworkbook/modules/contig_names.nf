process contig_names {
    queue = "ressexcon.q"
    cpus = 16
    memory = '8 GB'
    time = '4h'
    clusterOptions = { '-P ressexcon' }


    input:
    tuple val(species), val(ref)
    output:
    tuple val(species), file("*.txt")

    script:
    """
    #!/bin/bash

    mkdir contigs
    cat ${params.fasta_dir}/${ref}.fna | grep '>' | cut -d ' ' -f1 | cut -c2- | grep NC* > contigs.txt
    for contig in \$(grep NC contigs.txt)
    do
    echo \$contig > \${contig}_${species}.txt
    done
    rm contigs.txt
    """
}
