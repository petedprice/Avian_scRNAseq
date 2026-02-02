nextflow.enable.dsl=2

include { fastqc } from './modules/QC/fastqc.nf'
include { fastqc as fastqc2 } from './modules/QC/fastqc.nf'
include { multiqc } from './modules/QC/multiqc.nf'
include { trimmomatic } from './modules/QC/trimmomatic.nf'
include { trimmed_multiqc } from './modules/QC/trimmed_multiqc.nf'
include { paired_multiqc } from './modules/QC/paired_multiqc.nf'

include { bwa_index } from './modules/align_etc/bwa_index.nf'
include { bwa } from './modules/align_etc/bwa.nf'
include { samtools_sort_bam_index } from './modules/align_etc/samtools_sort_bam_index.nf'

include { qualimap } from './modules/QC/qualimap.nf'
include { bam_multiqc } from './modules/QC/bam_multiqc.nf'

include { split_bam_contig } from './modules/align_etc/split_bam_contig.nf'

include { contig_names } from './modules/align_etc/contig_names.nf'

include { CreateSequenceDictionary } from './modules/gatk/CreateSequenceDictionary.nf'
include { AddOrReplaceReadGroups } from './modules/gatk/AddOrReplaceReadGroups.nf'
include { MarkDuplicates } from './modules/gatk/MarkDuplicates.nf'
include { SplitNCigarReads } from './modules/gatk/SplitNCigarReads.nf'
include { HaplotypeCaller } from './modules/gatk/HaplotypeCaller.nf'
include { GenomicsDBImport } from './modules/gatk/GenomicsDBImport.nf'
include { GenotypeGVCFs } from './modules/gatk/GenotypeGVCFs.nf'
include { VariantFiltration } from './modules/gatk/VariantFiltration.nf'

include { vcftools_filt } from './modules/bvcftools/vcftools_filt.nf'
include { bcftools_bed } from './modules/bvcftools/bcftools_bed.nf'

include { pixy } from './modules/other_tools/pixy.nf'

include { SnpEff_builddb } from './modules/SnpEff/SnpEff_builddb.nf'
include { SnpEff_ncbi } from './modules/SnpEff/SnpEff_ncbi.nf'
include { degenotate } from './modules/other_tools/degenotate.nf'


workflow {
     sample_ch=Channel
        .fromPath(params.metadata)
        .splitCsv()
        .map {row -> tuple(row[0], row[1])}

     //initial qc of all reads
     fastqc1=fastqc(sample_ch)
     multiqc=multiqc(fastqc1.collect())

     //trim those reads down
     trimmed=trimmomatic(sample_ch)

     //Second qc of all reads
     fastqc2=fastqc2(trimmed)
     trimmed_multiqced=trimmed_multiqc(fastqc2.collect())
     prepost=fastqc1.combine(fastqc2, by:0)

     paired_multiqced=paired_multiqc(prepost)



     ref_index=bwa_index()

     aligned=bwa(trimmed)
     bam_sorted=samtools_sort_bam_index(aligned)


     //bam_qced=qualimap(bam_sorted)
     //bam_multiqc=bam_multiqc(bam_qced.collect())


     cns=contig_names().flatten()
     CreateSequenceDictionary()
     RGs=AddOrReplaceReadGroups(bam_sorted)
     dups=MarkDuplicates(RGs)

     contig_split=split_bam_contig(dups.combine(cns))

     all_bams=contig_split
	.groupTuple(by: 0)

     HCed=HaplotypeCaller(contig_split)
     GDBI=GenomicsDBImport(HCed.groupTuple(by: 0))
     GGVCFS=GenotypeGVCFs(GDBI)

     degen=degenotate().flatten().view()
     VFed=VariantFiltration(GGVCFS)
     VT_filt=vcftools_filt(VFed)
     FOLD_INTERSECT=bcftools_bed(VT_filt.combine(degen))


     gb_ch=Channel
        .fromPath(params.genbank)
        .splitCsv()
        .map {row -> tuple(row[0], row[1])}
  


     SnpEff_builddb()
     pixied=pixy(FOLD_INTERSECT)
    
}
