nextflow.enable.dsl=2

include { get_refs } from './modules/get_refs.nf'
include { longest_isoform } from './modules/longest_isoform.nf'
include { orthofinder } from './modules/orthofinder.nf'
include { ortho_cds } from './modules/ortho_cds.nf'

include { prank_allign } from './modules/prank_allign.nf'
include { format_fastas } from './modules/format_fastas.nf'
include { prank_phy } from './modules/prank_phy.nf'

include { remove_gaps } from './modules/remove_gaps.nf'
include { remove_Ns } from './modules/remove_Ns.nf'

include { mod0_paml } from './modules/mod0_paml.nf'
include { m1avsm2a } from './modules/m1avsm2a.nf'
include { comp_paml_models } from './modules/comp_paml_models.nf'

include { swamp_test } from './modules/swamp_test.nf'
include { swamp_final } from './modules/swamp_final.nf'




workflow {
    //Channels species name and reference name
    if (params.swamp == "checks"){
	    swamp_test_parameters=Channel
        	.fromPath(params.swamp_test_params)
        	.splitCsv()
        	.map {row -> tuple(row[0], row[1], row[2], row[3], row[4], row[5], row[6], row[7])}
        	.unique()
		.view()
    } else if(params.swamp == "yes"){
	    swamp_final_parameters=Channel
                .fromPath(params.swamp_final_params)
                .splitCsv()
                .map {row -> tuple(row[0], row[1], row[2], row[3], row[4], row[5], row[6], row[7])}
                .unique()
                .view()

    }


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

    if (params.swamp !="no"){
	    mod0ed=mod0_paml(phyed)
    }


    if(params.swamp == "checks"){
	   swamped=swamp_test(mod0ed.combine(swamp_test_parameters))
		.transpose()
	   M1aM2aed=m1avsm2a(swamped)
                .groupTuple(by: 1)

	   top_mod=comp_paml_models(M1aM2aed)

    } else if(params.swamp == "yes"){
	   swamped=swamp_final(mod0ed.combine(swamp_final_parameters))
	   NoNs=remove_Ns(swamped)
           M1aM2aed=m1avsm2a(NoNs)
		.groupTuple(by: 1)
           top_mod=comp_paml_models(M1aM2aed)

    } else if(params.swamp == "no"){
           M1aM2aed=m1avsm2a(phyed)
                .groupTuple(by: 1)
           top_mod=comp_paml_models(M1aM2aed)	     

    }
 
}



