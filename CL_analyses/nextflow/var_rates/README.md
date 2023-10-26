A pipeline for analysing sequence evolution rates across a number of user specified species from NCBI. 

PIPELINE: 
get_refs: pull references down from NCBI 

longest: get longest isoforms for each gene for each species 

orthofinder: find one-to-one orthologs of proteins using orthofinder 

ortho_cds: Get the CDS orthologs 

format_fastas: format fastas so for each gene, the gene name is the species name 

prank_allign: allign orthogroups 

remove_gaps: remove gaps in allignment and stop codon, remove allignments below a certain length and if not multiple of three or if the gaps are in multiples of three

prank_phy: convert fa allignment to .phy for paml input 

mod0_paml: run 0 model for masking by SWAMP




MAKE THE TWO MODEL CTL FILES 
MAKE ALL CTL FILES CHANGEABLE BY THE USER
CHANGE ORTHOFINDER OUTPUTS TO N0
ASK ALISON IF WE NEED THIS REMOVE N SECTION THAT SHE HAS IN HER ORIGINAL SCRIPTS

PROCESSES TO PRODUCE:
- SWAMP, INPUT: INPUTFILES AND SET OF NUMBERS FOR FILTERING 
- save in data, SWAMP_DEFAULT_PARAMETERS.txt
- MODEL, INPUT: PHY FILE AND CTL
- Print alignments in prior process 
- SWAMP AND MODELS

PROCESS BUT LESS PRIORITY:
- CHECK IF CDS FOR EACH GENOME, IF THERE THEN FINE
	[]IF NO CDS FILE, DOWNLOAD GTF AND GENOME FASTA TO CREATE ORFS FOR EACH GENE. 
	[]OR MANUALLY ASK FOR GTF AND GENOME ORF


```
SWAMP = NO, YES or MULTI

if SWAMP %in% NO/YES(	
	IF SWAMP ==YES (
		- Give default parameters, can override with user parameters DATA/SWAMP_DEFAULT_PARAMETERS.txt or USER_PARAMETERS.txt
		- Run SWAMP twice with pair of parameters, user input. 
	)
	- Run the two models for each og (single process) 
	-Get loglikelihoods and report p values
	IF SWAMP == YES (
		-Print all collective top hits to PAML_SUMMARY/MASKEDXXXX_paml_out.out
		-Print all alignments (PHY) to Allignments/masked_${og}XXXXX.phy
	)
	IF SWAMP == NO (
		-Print all collective top hits to PAML_SUMMARY/UNMASKED_paml_out.out
	)
)


if SWAMP == MULTI (
	- Require User input of a table of paremter combinations. 
	- Apply a factorial design of masking strategies 
	- Run SWAMP twice with pairs of parameters, user input. 
	- Run the two models for each og (single process) 
	- Get loglikelihoods and report p values
	- publishDir all PAML outputs and collective top orthogroup hits 
	- Print all collective top hits to PAML_SUMMARY_MULTI/MASKEDXXXX.out
	- Print all alignments (PHY) to MASKED_ALLIGNMENTS_MULTI/masked_${og}XXXXX.phy
)


```

