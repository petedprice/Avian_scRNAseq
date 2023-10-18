nextflow.enable.dsl=2

include { get_refs } from './modules/get_refs.nf'
include { longest_isoform } from './modules/longest_isoform.nf'
include { orthofinder } from './modules/orthofinder.nf'
include { ortho_cds } from './modules/ortho_cds.nf'
include { prank } from './modules/prank.nf'
include { format_fastas } from './modules/format_fastas.nf'

workflow {
    //Channels species name and reference name
    species_ch=Channel
	.fromPath(params.metadata)
	.splitCsv()
	.map {row -> tuple(row[0], row[1])}
	.unique()

    //Channels sample name and species name

    //Get's list of contig names, in:species, ref ;out:species, contig.txt
    seqs=get_refs(species_ch)
    longest=longest_isoform(seqs).collect()
    orthofound=orthofinder(longest)
    of_cds=ortho_cds(orthofound)
	.transpose()
    species_formatted=format_fastas(of_cds)
    pranked=prank(species_formatted)

    // we are now going to have to allign the orthogroups against each other
    // allign both the cds and protein sequence 
    //pranked=prank(orthogroups)



 
}




