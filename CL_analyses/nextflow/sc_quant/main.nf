nextflow.enable.dsl=2

include { get_refs } from './modules/ncbi/get_refs.nf'
include { get_mt_genes } from './modules/ncbi/get_mt_genes.nf'

include { cellranger_mkref } from './modules/cellranger/cellranger_mkref.nf'
include { cellranger_count } from './modules/cellranger/cellranger_count.nf'

include { seurat_filter } from './modules/seurat/seurat_filter.nf'
include { seurat_doublet } from './modules/seurat/seurat_doublet.nf'
include { seurat_SCT_integrate } from './modules/seurat/seurat_SCT_integrate.nf'
include { seurat_count_matrix } from './modules/seurat/seurat_count_matrix.nf'

//SCENIC
include { orthodb1 } from './modules/SCENIC_orthos/orthodb1.nf'
include { orthodb2 } from './modules/SCENIC_orthos/orthodb2.nf'
include { ortho_gff } from './modules/SCENIC_orthos/ortho_gff.nf'
include { chicken_seqs } from './modules/SCENIC_orthos/chicken_seqs.nf'

include { motif_modify } from './modules/SCENIC_TFs/motif_modify.nf'
include { TF_download } from './modules/SCENIC_TFs/TF_download.nf'
include { create_fasta_promoter } from './modules/SCENIC_TFs/create_fasta_promoter.nf'

include { create_cistarget_motif_databases } from './modules/SCENIC_pipeline/create_cistarget_motif_databases.nf'
include { combine_partials } from './modules/SCENIC_pipeline/combine_partials.nf'
include { grn } from './modules/SCENIC_pipeline/grn.nf'
include { ctx } from './modules/SCENIC_pipeline/ctx.nf'
include { aucell } from './modules/SCENIC_pipeline/aucell.nf'

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
/*
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

    //INTEGRATE
    seurat_integrated=seurat_SCT_integrate(seurat_doubleted.groupTuple(by: [0,1,2]))
*/


    if (params.SCENIC == "run"){   
	    // SCENIC WORKFLOW //

	    //Download chicken reference genome and individually combine with each species reference
    	gffs=chicken_seqs() //Download chicken ref
		.combine(
			seqs.map { id -> (id[0,2])}) //seqs are the sequences from each of the other species 
			.map { id -> (id[1,0,2])} //reordering the outputs 

	    //Get the species taxanomic ID from the metadata file 
    	species_ids= Channel
        	.fromPath(params.metadata)
        	.splitCsv()
        	.map {row -> tuple(row[1], row[6])}
		.unique()

	    //Download orthoDB information for orthologs 
    	orthodbs1=orthodb1()

	    //
    	orthodbs2=orthodb2(species_ids.combine(orthodbs1))

	    orthodb_gff=orthodbs2.combine(gffs, by:0)

	    ortho_gffed=ortho_gff(orthodb_gff)
   
	    TFs=TF_download()
    	modified_motifs=motif_modify(TFs.combine(ortho_gffed))
   
	updownstreamed=create_fasta_promoter(seqs)

	ctds1=create_cistarget_motif_databases(
		TFs.combine(updownstreamed).combine(channel.of(1..250))
		)
    	ctds1_grouped=ctds1.groupTuple(by: 0)
  
    	ctds2=combine_partials(ctds1_grouped) 
/*
    	seurat_counts=seurat_count_matrix(seurat_integrated.combine(modified_motifs))   
    	grned=grn(seurat_counts)
    	ctxed=ctx(grn.combine(ctds2, by: 0))
    	aucelled=aucell(ctxed)
*/
	}
   
}



