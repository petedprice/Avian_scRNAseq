process get_refs {
    cpus = 1
    memory = '4 GB'
    time = '4h'


    input:
    tuple val(species), val(ref), val(mt_contig)

    output:
    tuple val(species), val(mt_contig),  file("${species}.gtf"), file("${species}.fa")
    
    script:
    """
    #!/bin/bash
    echo $species
    echo $ref
    rsync -t -v rsync://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt ./
    wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq_historical.txt ./
    grep -E ${ref} assembly_summary_refseq_historical.txt | cut -f 20 > ftp_links.txt
   grep -E ${ref} assembly_summary_refseq.txt | cut -f 20 >> ftp_links.txt
    awk 'BEGIN{FS=OFS="/";filesuffix="genomic.gtf.gz"}{ftpdir=\$0;asm=\$10;file=asm"_"filesuffix;print "rsync -t -v "ftpdir,file" ./"}' ftp_links.txt | sed 's/https/rsync/g' > download_gtf_files.sh
    awk 'BEGIN{FS=OFS="/";filesuffix="genomic.fna.gz"}{ftpdir=\$0;asm=\$10;file=asm"_"filesuffix;print "rsync -t -v "ftpdir,file" ./"}' ftp_links.txt | sed 's/https/rsync/g' > download_genome_files.sh

    source download_gtf_files.sh
    source download_genome_files.sh

    gunzip *gz
    mv *genomic.gtf ${species}.gtf
    mv *genomic.fna ${species}.fa
   


    """
}
