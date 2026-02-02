process SnpEff_builddb {

    errorStrategy 'retry'

    label 'snpf'
    cpus { 1 * task.attempt }
    memory { 16.GB * task.attempt }


    publishDir 'filtered_vcfs_stringent', mode: 'copy', overwrite: true, pattern: '*recode.vcf'

    input:

    output:

    script:
    """
    #!/bin/bash

    echo "species.genome : species" > snpEff.config 

    mkdir data
    mkdir data/species
    mkdir data/genomes

    cat ${params.assembly_report} | cut -f 7,1 | grep -v "^#" > chrom_names.txt

    perl ${projectDir}/software/SnpEff/scripts_build/ncbi/fix_gtf.pl ${params.gtf}.gz chrom_names.txt > data/species/genes.gtf

    cat ${params.fasta} | perl ${projectDir}/software/SnpEff/scripts_build/ncbi/fix_fasta.pl chrom_names.txt > data/species/sequences.fa
    cp data/species/sequences.fa data/genomes/species.fa
    cat ${params.protein} | perl ${projectDir}/software/SnpEff/scripts_build/ncbi/fix_fasta_protein_cds.pl protein_id.map.txt > data/species/protein.fa
    #cat ${params.mrna} | perl ${projectDir}/software/SnpEff/scripts_build/ncbi/fix_fasta_protein_cds.pl protein_id.map.txt > data/species/cds.fa
    echo "Processing RNA FASTA files"
    cat ${params.mrna} | perl -pe 's/^>(\\S+).*/>\$1/' > data/species/cds.fa


    java -Xmx${(task.memory.toMega() * 0.9).toInteger()}m \
         -jar /opt/custflow/epi2meuser/conda/share/snpeff-5.1-4/snpEff.jar \
         build -gtf22 -v species

    #snpEff build -gtf22 -v species


    """
}

