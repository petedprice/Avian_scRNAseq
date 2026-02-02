process seurat_filter {

    label 'seurat'

    cpus { 8 * task.attempt }
    errorStrategy 'retry'
    maxRetries 6
    memory { 96.GB * task.attempt }

    publishDir 'cellranger_htmls', mode: 'copy', overwrite: true, pattern: '*html'

    input: 
    tuple val(species), val(sample), val(stage), val(sex), file("${species}_mt_genes.txt"), file("${sample}_crdata")

    output:
    tuple val(species), val(sample), val(stage), val(sex), file("filtered_seurat.RData"), file("${sample}_crdata"), file("${sample}_initial_QC_plots"), file("${sample}_web_summary.html")

    
    script:
    """
    #!/bin/bash

    cat ${params.metadata} | grep '${sample},' > metadata_ss.csv   

    Rscript ${projectDir}/Rscripts/seurat/seurat_filter.R \
	 ${projectDir} \
	 ${sample}_crdata/outs/filtered_feature_bc_matrix \
	 . \
	metadata_ss.csv \
	${species}_mt_genes.txt \
	200 \
	0.2 \
	0
    

    mv plots ${sample}_initial_QC_plots
    mv outdata/filtered_seurat.RData .
    cp ${sample}_crdata/outs/web_summary.html ${sample}_web_summary.html



    """
}


