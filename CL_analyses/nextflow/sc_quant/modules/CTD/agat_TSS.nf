process agat_TSS {
    cpus = 1
    memory = '4 GB'
    time = '4h'
 
    label 'agat'

    input:

    output:

    
    script:
    """
    #!/bin/bash

    for ref in Duck Pheasant Guineafowl 
    do 
    agat_sp_extract_sequences.pl --gff \${ref}*gff --fasta \${ref}*fna \
	-t gene --up "5000" --down "1000" --full -o tmp_\${ref}_5000up_1000down_TSS.fasta

    sed 's/>gene-/>/g' tmp_\${ref}_5000up_1000down_TSS.fasta > \${ref}_5000up_1000down_TSS.fasta

    agat_sp_extract_sequences.pl --gff \${ref}*gff --fasta \${ref}*fna \
	-t gene --up "500" --down "100" --full -o tmp_\${ref}_500up_100down_TSS.fasta
    sed 's/>gene-/>/g' tmp_\${ref}_500up_100down_TSS.fasta > \${ref}_500up_100down_TSS.fasta

    done
    """
}
