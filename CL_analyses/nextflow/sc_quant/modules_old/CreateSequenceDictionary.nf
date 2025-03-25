process create_seq_dict {

    input:
    tuple val(species), val(ref)

    script:
    """
    #!/bin/bash
    picard CreateSequenceDictionary R=${params.fasta_dir}/${ref}.fna O=${params.fasta_dir}/${ref}.dict

    """

}
