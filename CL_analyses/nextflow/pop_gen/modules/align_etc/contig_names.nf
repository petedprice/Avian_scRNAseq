process contig_names {
    cpus = 1
    memory = '4 GB'
    time = '4h'


    input:

    output:
    file("*.txt")

    script:
    """
    #!/bin/bash
    mkdir contigs
    cat ${params.fasta} | grep ">" | cut -f1 -d " " | sed 's/^>//' | grep "NC" | uniq > contigs.txt
    for contig in \$(cat contigs.txt)
    do
    echo \$contig > \${contig}.txt
    done
    rm contigs.txt
    """
}
