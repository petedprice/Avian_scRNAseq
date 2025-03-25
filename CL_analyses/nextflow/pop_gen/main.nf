nextflow.enable.dsl=2

include { fasterq_dump } from './modules/ncbi/fasterq_dump.nf'

include { fastqc } from './modules/QC/fastqc.nf'
include { fastqc as fastqc2 } from './modules/QC/fastqc.nf'

include { multiqc } from './modules/QC/multiqc.nf'
include { trimmed_multiqc } from './modules/QC/trimmed_multiqc.nf'

include { trimmomatic } from './modules/QC/trimmomatic.nf'
include { bwa_index } from './modules/align_etc/bwa_index.nf'
include { bwa } from './modules/align_etc/bwa.nf'
include { samtools_sort_bam_index } from './modules/align_etc/samtools_sort_bam_index.nf'

include { qualimap } from './modules/QC/qualimap.nf'
include { bam_multiqc } from './modules/QC/bam_multiqc.nf'

workflow {
    reads_ch=Channel
	.fromPath(params.metadata)
	.splitCsv()
	.map {row -> row[0]}


     reads=fasterq_dump(reads_ch)

     //initial qc of all reads
     fastqc1=fastqc(reads)
     multiqc=multiqc(fastqc1.collect())

     //trim those reads down
     trimmed=trimmomatic(reads)

     //Second qc of all reads 
     fastqc2=fastqc2(trimmed)
     trimmed_multiqced=trimmed_multiqc(fastqc2.collect())

     ref_index=bwa_index()

     aligned=bwa(trimmed)
     bam_sorted=samtools_sort_bam_index(aligned)


     bam_qced=qualimap(bam_sorted)
     bam_multiqc=bam_multiqc(bam_qced.collect())

/*
     bam_files=aligned.collect()
     pop_gen=angsd(bam_files)
  */
}
