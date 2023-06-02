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

workflow {
    //Channels species name and reference name
    species_ch=Channel
	.fromPath(params.metadata)
	.splitCsv()
	.map {row -> tuple(row[1], row[2])}
	.unique().view()

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
        splitted=split_bam(counted.combine(cns, by: 0)).view()
   
    } else {
	//if cellranger data pre made then run as so 
	//splitted=nocellranger_split_bam(samples_ch)
        splitted=nocellranger_split_bam(samples_ch.combine(cns, by:0))
    }
 
	        
    //Mutation rate analysis
    mut_ided = mut_id(splitted.combine(species_ch, by:0))
    Rfiltered = R_mut_filtering(mut_ided)

    mut_comp = mut_compiled(Rfiltered.groupTuple())

 
}




