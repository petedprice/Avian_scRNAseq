process seurat_markers {

    label 'seurat'

    cpus { 26 }
    errorStrategy 'retry'
    maxRetries 6
    memory { 340.GB }
    //publishDir 'mut_compiled', mode: 'copy', overwrite: true, pattern: '*snp_summarise.txt.gz'

    input: 
    tuple val(samples), val(sex), val(stage), val(species), val(ref), file("integrated_seurat.RData")

    output:
    tuple val(samples), file("marker_seurat.RData"), val(sex), val(stage), val(species)

    script:
    """
    #!/bin/bash
    echo dog
    Rscript ${projectDir}/Rscripts/seurat/3.cell_type_ID.R \
	integrated_seurat.RData \
	. \
	${params.celltype_markers} \
	$sex \
	$species

    mv outdata/marker_seurat.RData .	


    """
}


