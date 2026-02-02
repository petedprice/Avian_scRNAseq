nextflow.enable.dsl=2

include { samtools_sort_bam_index } from './modules/align_etc/samtools_sort_bam_index.nf'
include { split_bam_contig } from './modules/align_etc/split_bam_contig.nf'

include { contig_names } from './modules/align_etc/contig_names.nf'

include { AddOrReplaceReadGroups } from './modules/gatk/AddOrReplaceReadGroups.nf'
include { MarkDuplicates } from './modules/gatk/MarkDuplicates.nf'
include { SplitNCigarReads } from './modules/gatk/SplitNCigarReads.nf'
include { HaplotypeCaller } from './modules/gatk/HaplotypeCaller.nf'
include { GenomicsDBImport } from './modules/gatk/GenomicsDBImport.nf'
include { GenotypeGVCFs } from './modules/gatk/GenotypeGVCFs.nf'

include { vcftools_filt } from './modules/vcftools/vcftools_filt.nf'
include { vcftools_TajD } from './modules/vcftools/vcftools_TajD.nf'

include { idsfs } from './modules/angsd/idsfs.nf'

include { pixy } from './modules/other_tools/pixy.nf'

workflow {
    bam_ch=Channel
	.fromPath(params.metadata)
	.splitCsv()
	.map {row -> row[0]}

     cns=contig_names().flatten()

     bam_sorted=samtools_sort_bam_index(bam_ch)
     RGs=AddOrReplaceReadGroups(bam_sorted)
     dups=MarkDuplicates(RGs)
     splitNCR=SplitNCigarReads(dups)

     contig_split=split_bam_contig(splitNCR.combine(cns))

     all_bams=contig_split
	.groupTuple(by: 0)

     //idsfs=idsfs(all_bams)


     HCed=HaplotypeCaller(contig_split)
     GDBI=GenomicsDBImport(HCed.groupTuple(by: 0))
     GGVCFS=GenotypeGVCFs(GDBI)
     VT_filt=vcftools_filt(GGVCFS)
//     TajD=vcftools_TajD(VT_filt)

     pixied=pixy(VT_filt)
	
    
}
