process seurat_filter {

    label 'seurat'

    cpus { 8 * task.attempt }
    errorStrategy 'retry'
    maxRetries 6
    memory { 96.GB * task.attempt }

    //publishDir 'mut_compiled', mode: 'copy', overwrite: true, pattern: '*snp_summarise.txt.gz'

    input: 
    tuple val(sample), val(sex), val(stage), val(species), val(ref)

    output:
    tuple val(sample), val(sex), val(stage), val(species), val(ref), file("filtered_seurat.RData")

    script:
    """
    #!/bin/bash
    mt_contig=\$(cat ${params.gff_dir}/${ref}.gff | grep Name=MT | grep region | grep 1 | cut -f1)
    cat ${params.gtf_dir}/${ref}.gtf | grep \$mt_contig | grep gene_id > ss.gtf
    cut -f9 ss.gtf | cut -d ";" -f1 | cut -d " " -f2 | uniq | sed 's/"//g' > mt_genes.txt



    cat ${params.metadata} | grep $species | grep ,$sex, | grep $stage | grep $sample > metadata_ss.csv
    Rscript ${projectDir}/Rscripts/seurat/subsetted_funcs/1.filtering.R \
	 ${projectDir} \
	 ${params.cellranger_data} \
	 . \
	metadata_ss.csv \
	mt_genes.txt
    mv outdata/filtered_seurat.RData .


    """
}


