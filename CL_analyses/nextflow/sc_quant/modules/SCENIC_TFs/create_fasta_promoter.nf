process create_fasta_promoter {
    cpus = 8
    memory = '64 GB'
    time = '4h'

    label 'agat'

    input:
    tuple val(species), val(mt_contig),  file("${species}.gtf"), file("${species}.fa")

    output:
    tuple val(species), file("${species}_5000up_1000down_TSS.fasta")    

    script:
    """
    #!/bin/bash
    agat_sp_extract_sequences.pl \
	--gff ${species}.gtf \
	--fasta ${species}.fa \
	-t gene --up "5000" --down "1000" \
	--full -o tmp_${species}_5000up_1000down_TSS.fasta

    sed 's/>gene-/>/g' tmp_${species}_5000up_1000down_TSS.fasta > ${species}_5000up_1000down_TSS.fasta

#    agat_sp_extract_sequences.pl \
#	--gff ${species}*gtf \
#	--fasta ${species}*fa \
#	-t gene --up "500" --down "100" \
#	--full -o tmp_${species}_500up_100down_TSS.fasta

#   sed 's/>gene-/>/g' tmp_${species}_500up_100down_TSS.fasta > ${species}_500up_100down_TSS.fasta




    """
}
