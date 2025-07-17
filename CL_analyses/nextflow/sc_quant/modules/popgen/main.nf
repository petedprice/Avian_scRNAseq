nextflow.enable.dsl=2

workflow {
    //Channels species name and reference name
    species_ch=Channel
	.fromPath(params.metadata)
	.splitCsv()
	.map {row -> tuple(row[1], row[4], row[5])}
	.unique()
	.view()

   
}



