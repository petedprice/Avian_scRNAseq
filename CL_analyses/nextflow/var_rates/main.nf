nextflow.enable.dsl=2

include { get_refs } from './modules/get_refs.nf'
include { longest_isoform } from './modules/longest_isoform.nf'
include { orthofinder } from './modules/orthofinder.nf'
include { ortho_cds } from './modules/ortho_cds.nf'
include { prank_allign } from './modules/prank_allign.nf'
include { format_fastas } from './modules/format_fastas.nf'
include { prank_phy } from './modules/prank_phy.nf'
include { remove_gaps } from './modules/remove_gaps.nf'
include { mod0_paml } from './modules/mod0_paml.nf'
include { swamp } from './modules/swamp.nf'



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
    pranked_alligned=prank_allign(species_formatted)
    No_gaps=remove_gaps(pranked_alligned)
    phyed=prank_phy(No_gaps)


    pre_swamp_z_model=mod0_paml(phyed)

    swamp1=swamp(pre_swamp_z_model)


    // we are now going to have to allign the orthogroups against each other
    // allign both the cds and protein sequence 
    //pranked=prank(orthogroups)



 
}




