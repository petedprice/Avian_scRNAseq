nextflow.enable.dsl=2

include { cellranger_mkref } from './modules/cellranger_mkref.nf'
include { cellranger_count } from './modules/cellranger_count.nf'
include { contig_names } from './modules/contig_names.nf'
include { split_bam } from './modules/split_bam.nf'
include { mut_id } from './modules/mut_id.nf'

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

    //Makes cellranger ref for each species, in: species, ref; out:species, species_cr_ref 
    ref_made=cellranger_mkref(species_ch)

    //Quantifies cellranger data, in:species, species_cr_ref, sample; out: species, sample, bam, bam_index
    counted=cellranger_count(ref_made.combine(samples_ch, by: 0))

<<<<<<< HEAD
    //Get's list of contig names, in:species, ref ;out:species, contig.txt
    cns=contig_names(species_ch)

    //cns=contig_names(counted)
    //	.transpose()

=======
    cns=contig_names(counted)
	.transpose()
>>>>>>> develop
    splitted=split_bam(cns).view()
    mut_ided = mut_id(splitted)
 
    //cns2=contig_names(species_ch).transpose()
    //splitted=split_bam(cns2.combine(counted, by: 0)

}




