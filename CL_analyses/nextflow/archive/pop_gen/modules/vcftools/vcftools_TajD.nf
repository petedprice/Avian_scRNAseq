process vcftools_TajD {

    label 'vcftools'
    errorStrategy 'retry'

    cpus { 4 * task.attempt }
    memory { 32.GB * task.attempt }

    publishDir 'TajD', mode: 'copy', overwrite: true, pattern: 'Taj_*'

    input:
    tuple val(contig), file("${contig}_filt.recode.vcf")

    output:
    tuple val(contig), file("Taj_${contig}")

    script:
    """
    #!/bin/bash
    
    cat ${params.gtf} | awk '{if (\$3 == "gene") print \$0;}' | grep ${contig} | cut -f1,4,5 > tmp.bed
     
    while read -r chrom start end; do
	start_num=\$((start + 0))
	end_num=\$((end + 0))
	window=\$((end - start))
	echo "Running Tajima's D for region: \${chrom}:\${start_num}-\${end_num}"

	vcftools --vcf "\${chrom}_filt.recode.vcf" \
        	--chr "\${chrom}" \
        	--from-bp "\${start_num}" \
        	--to-bp "\${end_num}" \
        	--TajimaD "\${window}" \
        	--out "region_\${chrom}_\${start_num}_\${end_num}"

    done < tmp.bed

    mkdir Taj_${contig}
    mv *.Tajima.D Taj_${contig}

    """
}


