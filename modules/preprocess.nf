#!/bin/bash nextflow
params.outdir = 'results'

process COMBINE_SEQS {
    container 'edgano/base:latest'
    tag "Combine Seqs on $id"
    publishDir params.outdir

    input:
    tuple id, path(seqs)
    tuple id2, path(refs)

    output:
    file("${id}.fa") 

    script:
    """
    # CREATE A FASTA FILE CONTAINING ALL SEQUENCES (SEQS + REFS)
    t_coffee -other_pg seq_reformat -in ${refs} -output fasta_seq -out refs.tmp.fa
    t_coffee -other_pg seq_reformat -in ${seqs} -output fasta_seq -out seqs.tmp.fa
    cat refs.tmp.fa > completeSeqs.fa
    cat seqs.tmp.fa >> completeSeqs.fa
    # SHUFFLE ORDER OF SEQUENCES
    t_coffee -other_pg seq_reformat -in completeSeqs.fa -output fasta_seq -out ${id}.fa -action +reorder random
    """
}

process GENERATE_DYNAMIC_CONFIG {
    container 'edgano/base:latest'
    tag "Config 4 Dynamic"

    input:
    val (masterAln)
    val (masterSize)
    val (slaveAln)
    val (slaveSize)

    output:
      path "${masterAln}.${masterSize}_${slaveAln}.${slaveSize}.config", emit: configFile
      tuple val(masterAln), val(masterSize), val(slaveAln), val(slaveSize), emit: configValues

    script:
    """
    echo '${masterAln} ${masterSize}' > ${masterAln}.${masterSize}_${slaveAln}.${slaveSize}.config
    echo '${slaveAln} ${slaveSize}' >> ${masterAln}.${masterSize}_${slaveAln}.${slaveSize}.config
    """
}

process PRECOMPUTE_BLAST {
    container 'edgano/tcoffee:pdb'
    tag "BLAST on $id"
    publishDir "${params.blastOutdir}", mode: 'copy', overwrite: true, pattern: '*tmp.gz'
    label 'process_medium'

    input:
    tuple val(id), path(seqs)

    output:
    val "${id}", emit: id
    path "*tmp.gz", emit: blast

    script:
    """
    t_coffee -other_pg seq_reformat -in ${seqs} -action +compress +db ${params.database_path} +seq2blast
    """
}

process INTRAMOL_MATRIX_GENERATION {
    container 'edgano/tcoffee:protocols'
    tag "Intramolecular matrix on $id"
    publishDir "${params.outdir}/matrices", mode: 'copy', overwrite: true

    input:
    tuple val(id), file(fasta), file(template)

    output:
    tuple val (id), path("${id}.matrices"), emit: id_Matrix

    script:
    """
    export THREED_TREE_MODE=${params.threedTreeMode}
    export THREED_TREE_NO_WEIGHTS=${params.threedTreeNoWeights}
    export THREED_TREE_MODE_EXP=${params.threedTreeModeExp}

    t_coffee -other_pg seq_reformat -in ${fasta} -in2 ${template} -action +tree replicates ${params.replicatesNum} \
        +evaluate3D ${params.evaluate3DVal} +tree2bs first +print_replicates -output dm > ${id}.matrices 
    """
}

process SELECTED_PAIRS_OF_COLUMNS_MATRIX {     //TODO
    container 'edgano/tcoffee:protocols'
    tag "Selected Columns on $id"
    publishDir "${params.outdir}/matrices", mode: 'copy', overwrite: true, pattern: '*tmp.gz'

    input:
    tuple val(id), path(fasta)
    tuple val(id), path (pair)
    set val(id), file(fasta), file(template)
    val (numReplicates)

    output:
    tuple val (id), path("${id}.matrices"), emit: id_Matrix

    script:
    """
    export THREED_TREE_MODE=${params.threedTreeMode}
    export THREED_TREE_NO_WEIGHTS=${params.threedTreeNoWeights}
    export THREED_TREE_MODE_EXP=${params.threedTreeModeExp}

    t_coffee -other_pg seq_reformat -in ${fasta} -in2 ${template} -action +columns4tree ${pair} +tree replicates ${numReplicates} \
        +evaluate3D ${params.evaluate3DVal} +tree2bs first +print_replicates -output dm > ${pair}.matrices
    """
}

process LIBRARY_GENERATION {
    container 'edgano/tcoffee:protocols'
    tag "on $id"
    publishDir "${params.outdir}/${id}/libraries", mode: 'copy', overwrite: true

    input:
    set val(id), file(fasta), file(template)

    output:
    tuple val (id), path("${id}_sap.lib"), emit: sap_lib
    tuple val (id), path("${id}_tmalign.lib"), emit: tmalign_lib
    tuple val (id), path("${id}_mustang.lib"), emit: mustang_lib

    script:
    """
    t_coffee -seq ${fasta} -template_file ${template} -method sap_pair -out_lib ${id}_sap.lib
    t_coffee -seq ${fasta} -template_file ${template} -method TMalign_pair -out_lib ${id}_tmalign.lib
    t_coffee -seq ${fasta} -template_file ${template} -method mustang_pair -out_lib ${id}_mustang.lib
    """
}

process ALN_2_PHYLIP {
    container 'edgano/tcoffee:protocols'
    tag "on $id"
    publishDir "${params.outdir}/phylip", mode: 'copy', overwrite: true

    input:
    tuple val(id), path(fasta)

    output:
    tuple val (id), path("${id}.ph"), emit: id_phylip

    script:
    """
    t_coffee -other_pg seq_reformat -in ${fasta} -output phylip_aln > ${id}.ph
    """
}

