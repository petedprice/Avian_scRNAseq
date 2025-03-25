nextflow.enable.dsl=2

include { get_refs } from './modules/ncbi/get_refs.nf'
include { get_mt_genes } from './modules/ncbi/get_mt_genes.nf'

include { cellranger_mkref } from './modules/cellranger/cellranger_mkref.nf'
include { cellranger_count } from './modules/cellranger/cellranger_count.nf'

include { seurat_filter } from './modules/seurat/seurat_filter.nf'
include { seurat_doublet } from './modules/seurat/seurat_doublet.nf'
include { seurat_SCT_integrate } from './modules/seurat/seurat_SCT_integrate.nf'
include { seurat_SCT_integrate as seurat_SCT_integrate2 } from './modules/seurat/seurat_SCT_integrate.nf'
include { seurat_SCT_integrate as seurat_SCT_integrate3 } from './modules/seurat/seurat_SCT_integrate.nf'
include { seurat_SCT_integrate as seurat_SCT_integrate4 } from './modules/seurat/seurat_SCT_integrate.nf'

workflow {
    //Channels species name and reference name
    species_ch=Channel
	.fromPath(params.metadata)
	.splitCsv()
	.map {row -> tuple(row[1], row[4], row[5])}
	.unique()
	.view()


    //Channels sample name and species name

    //Get's list of contig names, in:species, ref ;out:species, contig.txt
    seqs=get_refs(species_ch)

   mt_genes=get_mt_genes(seqs)

    //Channels sample name and species name
    samples_ch = Channel
	.fromPath(params.metadata)
	.splitCsv()
	.map {row -> tuple(row[1], row[0])}

    //make cellranger index
    ref_made=cellranger_mkref(seqs)

    //Quantifies cellranger data, in:species, species_cr_ref, sample; out: species, sample, bam, bam_index
    counted=cellranger_count(ref_made.combine(samples_ch, by: 0))
   
    //seurat pipeline
    sex_stage_ch = Channel
		.fromPath(params.metadata)
        	.splitCsv()
        	.map {row -> tuple(row[1], row[0], row[2], row[3])}
		.combine(mt_genes, by: 0)
		.combine(counted, by: [0,1])

    //READ IN CELLRANGER DATA AND FILTER SEURAT OBJECTS
    seurat_filtered=seurat_filter(sex_stage_ch)

    

    //REMOVE DOUBLETS 
    seurat_doubleted=seurat_doublet(seurat_filtered)


    //INTEGRATE individually
    seurat_integrated=seurat_SCT_integrate(seurat_doubleted.groupTuple(by: [0,1,2]))


    //Integrate within species across sexes for each stage
    seurat_integrated2=seurat_SCT_integrate2(
	    seurat_doubleted.groupTuple(by: [0, 1])
		.map { id -> tuple(id[0], id[1], "MIXED", id[3])} 
    )

    //Integrate within species across stages for each sex
    seurat_integrated3=seurat_SCT_integrate3(
            seurat_doubleted.groupTuple(by: [0, 2])
	       	.map { id -> tuple(id[0], "MIXED", id[2], id[3])}
    )

    // Integrate all samples for a species 
    seurat_integrated4=seurat_SCT_integrate4(
            seurat_doubleted.groupTuple(by: [0])
	       	.map { id -> tuple(id[0], "MIXED",  "MIXED", id[3])}
    )


}



