nextflow.enable.dsl=2

include { cellranger_mkref } from './modules/cellranger_mkref.nf'
include { cellranger_count } from './modules/cellranger_count.nf'
include { seurat_filter } from './modules/seurat/seurat_filter.nf'
include { seurat_doublet } from './modules/seurat/seurat_doublet.nf'

include { seurat_SCT_integrate } from './modules/seurat/seurat_SCT_integrate.nf'
include { seurat_SCT_integrate_noX } from './modules/seurat/seurat_SCT_integrate_noX.nf'

include { seurat_get_clusters } from './modules/seurat/seurat_get_clusters.nf'
include { seurat_get_barcodes } from './modules/seurat/seurat_get_barcodes.nf'

include { get_samp } from './modules/get_samp.nf'
include { sc_var_call } from './modules/sc_var_call.nf'
include { R_var_filtering } from './modules/R_var_filtering.nf'
include { clus_sc_var_call } from './modules/clus_sc_var_call.nf'
include { clus_R_var_filtering } from './modules/clus_R_var_filtering.nf'


include { contig_names } from './modules/contig_names.nf'
include { split_bam } from './modules/split_bam.nf'
include { split_bam_cluster } from './modules/split_bam_cluster.nf'


include { ws_var_call_cleanup } from './modules/ws_varcall/ws_var_call_cleanup.nf'
include { ws_var_call_splitNCR } from './modules/ws_varcall/ws_var_call_splitNCR.nf'
include { ws_var_call_HC } from './modules/ws_varcall/ws_var_call_HC.nf'
include { ws_var_call_VF_stringent } from './modules/ws_varcall/ws_var_call_VF_stringent.nf'
include { ws_var_call_VF_relaxed } from './modules/ws_varcall/ws_var_call_VF_relaxed.nf'


workflow {
    //Channels species name and reference name
    species_ch=Channel.fromPath(params.metadata)
	.splitCsv()
	.map {row -> tuple(row[1], row[2])}
	.unique()

    param_test=Channel.fromPath(params.filter_params)
	.splitCsv(header: true)
        .map {row -> tuple(row.nfs, row.mtr, row.gu)}
        .unique()
	

    //Channels sample name and species name
    samples_ch = Channel
	.fromPath(params.metadata)
	.splitCsv()
	.map {row -> tuple(row[1], row[0])}

    species_mt_ch=Channel.fromPath(params.metadata)
        .splitCsv()
        .map {row -> tuple(row[1], row[2], row[5])}
        .unique()

    if (params.alevin_quant == 'FALSE'){
	    ref_made=cellranger_mkref(species_ch).view()
	    counted=cellranger_count(ref_made.combine(samples_ch, by: 0))
    } else if (params.alevin_quant == 'TRUE'){
    	    ref_made_alevin=alevin_fry_index(species_ch)
	    counted=alevin_fry_quant(ref_made_alevin.combine(samples_ch, by: 0).view())
    }


    seurat_filtered=seurat_filter(counted.combine(species_mt_ch, by: 0)
	.combine(param_test))


    seurat_doubleted=seurat_doublet(seurat_filtered)
	.groupTuple(by: [1,2,3])

    seurat_integrated=seurat_SCT_integrate(seurat_doubleted)



    if (params.sc_var_call == 'TRUE'){
	cleaned=ws_var_call_cleanup(counted)
	splitNCRed=ws_var_call_splitNCR(cleaned.combine(species_ch, by: 0))
	ws_VCed=ws_var_call_HC(splitNCRed)
	VFs=ws_var_call_VF_stringent(ws_VCed)
	VFr=ws_var_call_VF_relaxed(ws_VCed)
    }

}




