#!/bin/bash nextflow
//params.outdir = 'results'

include set_templates_path from './functions.nf'
path_templates = set_templates_path()

process REG_ALIGNER {
    container 'edgano/tcoffee:pdb'
    tag "$align_method - $tree_method - $bucket_size on $id"
    publishDir "${params.outdir}/alignments", pattern: '*.aln'

    input:
    tuple val(id), val(tree_method), path(seqs), path(guide_tree)
    each align_method
    each bucket_size

    output:
    val align_method, emit: alignMethod
    val tree_method, emit: treeMethod
    val bucket_size, emit: bucketSize
    tuple val (id), path ("${id}.reg_${bucket_size}.${align_method}.with.${tree_method}.tree.aln"), emit: alignmentFile
    path "${id}.homoplasy", emit: homoplasyFile
    path ".command.trace", emit: metricFile

    script:
    template "${path_templates}/regressive_align/reg_${align_method}.sh"   
}

process PROG_ALIGNER {
    container 'edgano/tcoffee:pdb'
    tag "$align_method - $tree_method on $id"
    publishDir "${params.outdir}/alignments", pattern: '*.aln'

    input:
    tuple val(id), val(tree_method), path(seqs), path(guide_tree)
    each align_method

    output:
    val align_method, emit: alignMethod
    val tree_method, emit: treeMethod
    tuple val (id), path ("${id}.prog.*.tree.aln"), emit: alignmentFile
    path ".command.trace", emit: metricFile

    script:
    template "${path_templates}/progressive_align/prog_${align_method}.sh"
}

process SLAVE_ALIGNER {
    container 'edgano/tcoffee:pdb'
    tag "$align_method - $tree_method - $slave_method - $bucket_size on $id"
    publishDir "${params.outdir}/alignments", pattern: '*.aln'

    input:
    tuple val(id), val(tree_method), path(seqs), path(guide_tree)
    each align_method
    each bucket_size
    each slave_method   

    output:
    val align_method, emit: alignMethod
    val "${tree_method}_${slave_method}", emit: treeMethod
    val bucket_size, emit: bucketSize
    tuple val (id), path ("${id}.slave_${bucket_size}.${align_method}.with.${tree_method}_${slave_method}.tree.aln"), emit: alignmentFile
    path "${id}.homoplasy", emit: homoplasyFile
    path ".command.trace", emit: metricFile

    script:
    template "${path_templates}/slave_align/slave_${align_method}.sh"
}

process DYNAMIC_ALIGNER {
    container 'edgano/tcoffee:pdb'
    tag "$align_method - $tree_method on $id"
    publishDir "${params.outdir}/alignments", pattern: '*.aln'
    label 'process_medium'

    input:
    tuple val(id), val(tree_method), path(seqs), path(guide_tree)
    val align_method
    each bucket_size
    each dynamicX
    path (dynamicConfig)
    tuple val(masterAln), val(masterSize), val(slaveAln), val(slaveSize)
    val dynamicValues

    output:
    val "${dynamicValues}_${params.db}", emit: alignMethod
    val tree_method, emit: treeMethod
    val "${bucket_size}_${dynamicX}", emit: bucketSize
    tuple val (id), path("*.with.${tree_method}.tree.aln"), emit: alignmentFile 
    path "${id}.homoplasy", emit: homoplasyFile
    path ".command.trace", emit: metricFile

    script:
    // template can be default or using a config file
    template "${path_templates}/dynamic_align/dynamic_${align_method}.sh"
}

process POOL_ALIGNER {
    container 'edgano/tcoffee:pdb'
    tag "$align_method - $tree_method - $bucket_size on $id"
    publishDir "${params.outdir}/alignments", pattern: '*.aln'

    input:
    tuple val(id), val(tree_method), path(seqs), path(guide_tree)
    each align_method
    each bucket_size

    output:
    val align_method, emit: alignMethod
    val tree_method, emit: treeMethod
    val bucket_size, emit: bucketSize
    tuple val (id), path ("${id}.pool_${bucket_size}.${align_method}.with.${tree_method}.tree.aln"), emit: alignmentFile
    path "${id}.homoplasy", emit: homoplasyFile
    path ".command.trace", emit: metricFile

    script:
    template "${path_templates}/pool_align/pool_${align_method}.sh"
}

process TCOFFEE_ALIGNER{
    container '44d0bbaf546c'//'edgano/tcoffee:protocols'
    tag "$tcoffee_mode  on $id"
    publishDir "${params.outdir}/alignments", pattern: '*.aln'

    input:
    tuple val(id), path(seqs)
    each tcoffee_mode
    set val(id), file(fasta), file(template)
    //path(lib)
    val(fakeId)             //to ensure the process before when precompute Blast
    //tuple val(id), path (template)
    //tuple val(id), path (pdbFile)

    output:
    val tcoffee_mode, emit: tcMode
    tuple val (id), path ("${id}.tcoffee.*.aln"), emit: alignmentFile
    path ".command.trace", emit: metricFile

    script:    
    if( tcoffee_mode == 'default' )
        """
        t_coffee -seq $seqs -outfile ${id}.tcoffee.${tcoffee_mode}.aln ${params.params4default}
        """ 
    else if( tcoffee_mode == 'quickaln' )
        """
        t_coffee -seq $seqs -mode quickaln -outfile ${id}.tcoffee.${tcoffee_mode}.aln ${params.params4quickalign}
        """    
    else if( tcoffee_mode == 'mcoffee' )
        """
        t_coffee -seq $seqs -mode mcoffee -outfile ${id}.tcoffee.${tcoffee_mode}.aln ${params.params4mcoffee}
        """ 
    else if( tcoffee_mode == 'accurate' )
    // TODO not working
    //t_coffee -in=data_93c5fbb0.in -mode=accurate -blast=LOCAL -pdb_db=/db/pdb/derived_data_format/blast/latest/pdb_seqres.fa -protein_db=/db/ncbi/201511/blast/db/nr.fa 
        """
        t_coffee -seq $seqs -mode accurate -blast LOCAL -pdb_db ${params.database_path} -outfile ${id}.tcoffee.${tcoffee_mode}.aln
        """            
    else if( tcoffee_mode == 'fmcoffee' )
        """
        t_coffee -seq $seqs -mode fmcoffee -outfile ${id}.tcoffee.${tcoffee_mode}.aln ${params.paramsfmcoffee}
        """
    else if( tcoffee_mode == 'psicoffee' )
        """
        t_coffee -seq $seqs -mode psicoffee -blast_server LOCAL -protein_db ${params.database_path} -cache ${params.blastOutdir} -outfile ${id}.tcoffee.${tcoffee_mode}.aln ${params.params4psicoffee}
        """
    else if( tcoffee_mode == 'expresso' )
        """
        t_coffee -seq $seqs -mode expresso -pdb_type d -blast LOCAL -pdb_db ${params.database_path} -outfile ${id}.tcoffee.${tcoffee_mode}.aln  ${params.params4expresso}
        """ 
    else if( tcoffee_mode == 'procoffee' )
        """
        t_coffee -seq $seqs -mode procoffee -outfile ${id}.tcoffee.${tcoffee_mode}.aln ${params.params4procoffee}
        """        
    else if( tcoffee_mode == '3dcoffee' )
        """
        t_coffee -seq $seqs -method sap_pair -template_file ${templates} -outfile ${id}.tcoffee.${tcoffee_mode}.aln ${params.params43dcoffee}
        """
    else if( tcoffee_mode == 'trmsd' )
    // TODO -> check to output
        """
        t_coffee -seq $seqs -method mustang_pair -template_file ${templates} -outfile ${id}.tcoffee.${tcoffee_mode}.aln ${params.params4trmsd}

        t_coffee -other_pg trmsd -aln ${id}.tcoffee.${tcoffee_mode}.aln -template_file ${templates}
        """
    else if( tcoffee_mode == 'rcoffee' )
        """
        t_coffee -seq $seqs -mode rcoffee -outfile ${id}.tcoffee.${tcoffee_mode}.aln ${params.params4rcoffee}
        """
    else if( tcoffee_mode == 'rcoffee_consan' )
        """
        t_coffee -seq $seqs -mode rcoffee_consan -outfile ${id}.tcoffee.${tcoffee_mode}.aln ${params.params4rcoffee_cosan}
        """
    else if( (tcoffee_mode == '3d_align') || (tcoffee_mode == '3dM_align') )     
    //3d align Leila CLANS                                -out ${id}_3dcoffee.fa
        """
        t_coffee -seq $seqs -lib ${lib} -output fasta_aln -outfile ${id}.tcoffee.${tcoffee_mode}.aln
        """    
    else
        error "Invalid alignment mode: ${tcoffee_mode}"
}