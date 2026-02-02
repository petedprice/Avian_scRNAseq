process SnpEff_ncbi {

    errorStrategy 'retry'

    label 'snpf'
    cpus { 1 * task.attempt }
    memory { 16.GB * task.attempt }


    publishDir 'filtered_vcfs_stringent', mode: 'copy', overwrite: true, pattern: '*recode.vcf'

    input:
    tuple val(genbank), val(contig)

    output:

    script:
    """
    #!/bin/bash

    DIR="data/${genbank}"
    GENE_FILE="\$DIR/genes.gbk"

    mkdir -p "\$DIR" >/dev/null 2>&1 || true

    # Download GenBank file
    wget "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=${genbank}&rettype=gbwithparts&retmode=text" > \$GENE_FILE

    echo "${genbank}.genome : ${genbank}" >> snpEff.config

    java -Xmx${(task.memory.toMega() * 0.9).toInteger()}m \
         -jar /opt/custflow/epi2meuser/conda/share/snpeff-5.1-4/snpEff.jar \
         build -genbank -v $genbank



    """
}

