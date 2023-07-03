nextflow.enable.dsl=2

include { cellranger_mkref } from './modules/cellranger_mkref.nf'
include { cellranger_count } from './modules/cellranger_count.nf'
include { contig_names } from './modules/contig_names.nf'
include { split_bam } from './modules/split_bam.nf'
include { mut_id } from './modules/mut_id.nf'
include { mut_compiled } from './modules/mut_compiled.nf'
include { nocellranger_split_bam } from './modules/nocellranger_split_bam.nf'
include { create_seq_dict } from './modules/CreateSequenceDictionary.nf'
include { R_mut_filtering } from './modules/R_mut_filtering.nf'
include { seurat_filter } from './modules/seurat/seurat_filter.nf'
include { seurat_SCT_integrate } from './modules/seurat/seurat_SCT_integrate.nf'
include { seurat_markers } from './modules/seurat/seurat_markers.nf'
include { seurat_mutation } from './modules/seurat/seurat_mutation.nf'
include { seurat_doublet } from './modules/seurat/seurat_doublet.nf'

params.seurat_etc="leave"

workflow {
    //Channels species name and reference name
    species_ch=Channel
	.fromPath(params.metadata)
	.splitCsv()
	.map {row -> tuple(row[1], row[2])}
	.unique()

    //Channels sample name and species name
    samples_ch = Channel
	.fromPath(params.metadata)
	.splitCsv()
	.map {row -> tuple(row[1], row[0])}


//    create_seq_dict(species_ch)



    //Get's list of contig names, in:species, ref ;out:species, contig.txt
    cns=contig_names(species_ch)
        .transpose()


    if(params.run_cellranger == "TRUE"){
    	//Makes cellranger ref for each species, in: species, ref; out:species, species_cr_ref 
    	ref_made=cellranger_mkref(species_ch)

    	//Quantifies cellranger data, in:species, species_cr_ref, sample; out: species, sample, bam, bam_index
    	counted=cellranger_count(ref_made.combine(samples_ch, by: 0))
   
        //split bam into individual contigs, in:species, sample, bam, bam_index, contig.txt
        splitted=split_bam(counted.combine(cns, by: 0))
   
    } else {
	//if cellranger data pre made then run as so 
	//splitted=nocellranger_split_bam(samples_ch)
        splitted=nocellranger_split_bam(samples_ch.combine(cns, by:0))
    }
 
	        
    //Mutation rate analysis
    mut_ided = mut_id(splitted.combine(species_ch, by:0))
    Rfiltered = R_mut_filtering(mut_ided)

    mut_comp = mut_compiled(Rfiltered.groupTuple())


    if(params.seurat_etc == "run"){
	sex_stage_ch = Channel
        	.fromPath(params.metadata)
        	.splitCsv()
        	.map {row -> tuple(row[0], row[3], row[4], row[1], row[2])}
	seurat_filtered=seurat_filter(sex_stage_ch)
	seurat_doubleted=seurat_doublet(seurat_filtered).groupTuple(by: [1,2,3,4]).view()

	seurat_integrated=seurat_SCT_integrate(seurat_doubleted)
        seurat_marked=seurat_markers(seurat_integrated)

	sex_stage_mut = mut_comp
		.combine(sex_stage_ch, by:0)
		.groupTuple(by: [2,3,4,5])
		 
	marker_mut = sex_stage_mut
		.combine(seurat_marked, by:[2,3,4])
		.map{it[3,0,1,2,4,7]}.view()
        seurat_mutated=seurat_mutation(marker_mut)
    }


 
}




