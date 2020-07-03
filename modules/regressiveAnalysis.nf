#!/bin/bash nextflow

include { REG_CLUSTALO ; REG_FAMSA ; REG_MAFFT_FFTNS1 }     from './regressiveAlignment.nf'   
include EVAL_ALIGNMENT    from './evaluateAlignment.nf'  
include EASEL_INFO        from './evaluateAlignment.nf'  
include HOMOPLASY         from './evaluateAlignment.nf'  
include METRICS           from './evaluateAlignment.nf'  

workflow ALIGN_REG {
    take:
        seqs_and_trees
        align_method
        tree_method
        bucket_size

    main:
        align_method = align_method.toString().replace("[", "").replace("]", "")
        println align_method    

        if (align_method == "CLUSTALO"){
            REG_CLUSTALO(seqs_and_trees, bucket_size)

            reg_alignment_method = REG_CLUSTALO.out.alignMethod
            reg_tree_method = REG_CLUSTALO.out.treeMethod
            reg_alignmentFile = REG_CLUSTALO.out.alignmentFile
            reg_bucket = REG_CLUSTALO.out.bucketSize
            reg_homoplasy = REG_CLUSTALO.out.homoplasyFile
            reg_metrics = REG_CLUSTALO.out.metricFile
            flavour = "regressive"
        }
        if (align_method == "FAMSA"){
            REG_FAMSA(seqs_and_trees, bucket_size)

            reg_alignment_method = REG_FAMSA.out.alignMethod
            reg_tree_method = REG_FAMSA.out.treeMethod
            reg_alignmentFile = REG_FAMSA.out.alignmentFile
            reg_bucket = REG_FAMSA.out.bucketSize
            reg_homoplasy = REG_FAMSA.out.homoplasyFile
            reg_metrics = REG_FAMSA.out.metricFile
            flavour = "regressive"
        }
        if (align_method == "MAFFT-FFTNS1"){
            REG_MAFFT_FFTNS1(seqs_and_trees, bucket_size)

            reg_alignment_method = REG_MAFFT_FFTNS1.out.alignMethod
            reg_tree_method = REG_MAFFT_FFTNS1.out.treeMethod
            reg_alignmentFile = REG_MAFFT_FFTNS1.out.alignmentFile
            reg_bucket = REG_MAFFT_FFTNS1.out.bucketSize
            reg_homoplasy = REG_MAFFT_FFTNS1.out.homoplasyFile
            reg_metrics = REG_MAFFT_FFTNS1.out.metricFile
            flavour = "regressive"
        }
    emmit:  
        //reg_alignment_method
        //reg_tree_method
        //reg_alignmentFile
        //reg_bucket
        //reg_homoplasy 
        //reg_metrics
        //flavour 
        reg_alignmentFile = REG_CLUSTALO.out.alignmentFile
}

workflow EVAL_REG {
    take:
        refs_ch
        //aligmentFile
        
        //flavour
        //align_method
        //tree_method
        //bucket_size
    main:
    refs_ch.voew()
        /*refs_ch
            .cross (aligmentFile)
            .map { it -> [ it[1][0], it[1][1], it[0][1] ] }
            .view()
            .set { alignment_and_ref }
        /*
        flavour = flavour.toString()
        align_method = align_method.toString()
        tree_method = tree_method.toString()
        println flavour    
        println align_method    
        println tree_method    

        EVAL_ALIGNMENT (flavour, alignment_and_ref, align_method, tree_method, bucket_size)
        EVAL_ALIGNMENT.out.tcScore
                        .map{ it ->  "${it[0]};${it[1]};${it[2]};${it[3]};${it[4]};${it[5].text}" }
                        .collectFile(name: "${workflow.runName}.${flavour}.tcScore.csv", newLine: true, storeDir:"${params.outdir}/CSV/${workflow.runName}/")
        EVAL_ALIGNMENT.out.spScore
                        .map{ it ->  "${it[0]};${it[1]};${it[2]};${it[3]};${it[4]};${it[5].text}" }
                        .collectFile(name: "${workflow.runName}.${flavour}.spScore.csv", newLine: true, storeDir:"${params.outdir}/CSV/${workflow.runName}/")
        EVAL_ALIGNMENT.out.colScore
                        .map{ it ->  "${it[0]};${it[1]};${it[2]};${it[3]};${it[4]};${it[5].text}" }
                        .collectFile(name: "${workflow.runName}.${flavour}.colScore.csv", newLine: true, storeDir:"${params.outdir}/CSV/${workflow.runName}/")
    emmit:  
        alignment_and_ref    */
     
}

workflow REG_ANALYSIS {
  take:
    seqs_and_trees
    refs_ch
    align_method
    tree_method
    bucket_size
    
  main: 
    //for(align in align_method){
        ALIGN_REG(seqs_and_trees,align_method,tree_method,bucket_size)
        
        //ALIGN_REG.out.reg_alignmentFile.view()
        
        //EVAL_REG(refs_ch)
        //ALIGN_REG.out.reg_alignmentFile, refs_ch, ALIGN_REG.out.flavour, ALIGN_REG.out.reg_alignmentM, ALIGN_REG.out.reg_tree_method, ALIGN_REG.out.reg_bucket)

    //}
}
    