#!/bin/bash nextflow

include { REG_CLUSTALO ; REG_FAMSA ; REG_MAFFT_FFTNS1 ; REG_MAFFT_GINSI ; REG_MAFFT_SPARSECORE ; REG_MAFFT ; REG_MSAPROBS ; REG_MUSCLE ; REG_PROBCONS ; REG_PSICOFFEE ; REG_UPP}     from './modules_regressiveAlignment.nf'   
workflow ALIGNMENT {
    take:
        seqs_and_trees
        align_method
        tree_method
        bucket_size

    main:
        //align_method = align_method.toString().replace("[", "").replace("]", "")
        //println align_method    

        if (align_method == "CLUSTALO"){
            REG_CLUSTALO(seqs_and_trees, bucket_size)

            alignment_method = REG_CLUSTALO.out.alignMethod
            tree_method = REG_CLUSTALO.out.treeMethod
            alignmentFile = REG_CLUSTALO.out.alignmentFile
            bucket = REG_CLUSTALO.out.bucketSize
            homoplasy = REG_CLUSTALO.out.homoplasyFile
            metrics = REG_CLUSTALO.out.metricFile
        }
        else if (align_method == "FAMSA"){
            REG_FAMSA(seqs_and_trees, bucket_size)

            alignment_method = REG_FAMSA.out.alignMethod
            tree_method = REG_FAMSA.out.treeMethod
            alignmentFile = REG_FAMSA.out.alignmentFile
            bucket = REG_FAMSA.out.bucketSize
            homoplasy = REG_FAMSA.out.homoplasyFile
            metrics = REG_FAMSA.out.metricFile
        }
        else if (align_method == "MAFFT-FFTNS1"){
            REG_MAFFT_FFTNS1(seqs_and_trees, bucket_size)

            alignment_method = REG_MAFFT_FFTNS1.out.alignMethod
            tree_method = REG_MAFFT_FFTNS1.out.treeMethod
            alignmentFile = REG_MAFFT_FFTNS1.out.alignmentFile
            bucket = REG_MAFFT_FFTNS1.out.bucketSize
            homoplasy = REG_MAFFT_FFTNS1.out.homoplasyFile
            metrics = REG_MAFFT_FFTNS1.out.metricFile
        }
        else if (align_method == "MAFFT-GINSI"){
            REG_MAFFT_GINSI(seqs_and_trees, bucket_size)

            alignment_method = REG_MAFFT_GINSI.out.alignMethod
            tree_method = REG_MAFFT_GINSI.out.treeMethod
            alignmentFile = REG_MAFFT_GINSI.out.alignmentFile
            bucket = REG_MAFFT_GINSI.out.bucketSize
            homoplasy = REG_MAFFT_GINSI.out.homoplasyFile
            metrics = REG_MAFFT_GINSI.out.metricFile
        }
        else if (align_method == "MAFFT-SPARSECORE"){
            REG_MAFFT_SPARSECORE(seqs_and_trees, bucket_size)

            alignment_method = REG_MAFFT_SPARSECORE.out.alignMethod
            tree_method = REG_MAFFT_SPARSECORE.out.treeMethod
            alignmentFile = REG_MAFFT_SPARSECORE.out.alignmentFile
            bucket = REG_MAFFT_SPARSECORE.out.bucketSize
            homoplasy = REG_MAFFT_SPARSECORE.out.homoplasyFile
            metrics = REG_MAFFT_SPARSECORE.out.metricFile
        }
        else if (align_method == "MAFFT"){
            REG_MAFFT(seqs_and_trees, bucket_size)

            alignment_method = REG_MAFFT.out.alignMethod
            tree_method = REG_MAFFT.out.treeMethod
            alignmentFile = REG_MAFFT.out.alignmentFile
            bucket = REG_MAFFT.out.bucketSize
            homoplasy = REG_MAFFT.out.homoplasyFile
            metrics = REG_MAFFT.out.metricFile
        }
        else if (align_method == "MSAPROBS"){
            REG_MSAPROBS(seqs_and_trees, bucket_size)

            alignment_method = REG_MSAPROBS.out.alignMethod
            tree_method = REG_MSAPROBS.out.treeMethod
            alignmentFile = REG_MSAPROBS.out.alignmentFile
            bucket = REG_MSAPROBS.out.bucketSize
            homoplasy = REG_MSAPROBS.out.homoplasyFile
            metrics = REG_MSAPROBS.out.metricFile
        }
        else if (align_method == "MUSCLE"){
            REG_MUSCLE(seqs_and_trees, bucket_size)

            alignment_method = REG_MUSCLE.out.alignMethod
            tree_method = REG_MUSCLE.out.treeMethod
            alignmentFile = REG_MUSCLE.out.alignmentFile
            bucket = REG_MUSCLE.out.bucketSize
            homoplasy = REG_MUSCLE.out.homoplasyFile
            metrics = REG_MUSCLE.out.metricFile
        }
        else if (align_method == "PROBCONS"){
            REG_PROBCONS(seqs_and_trees, bucket_size)

            alignment_method = REG_PROBCONS.out.alignMethod
            tree_method = REG_PROBCONS.out.treeMethod
            alignmentFile = REG_PROBCONS.out.alignmentFile
            bucket = REG_PROBCONS.out.bucketSize
            homoplasy = REG_PROBCONS.out.homoplasyFile
            metrics = REG_PROBCONS.out.metricFile
        }     
        else if (align_method == "PSICOFFEE"){
            REG_PSICOFFEE(seqs_and_trees, bucket_size)

            alignment_method = REG_PSICOFFEE.out.alignMethod
            tree_method = REG_PSICOFFEE.out.treeMethod
            alignmentFile = REG_PSICOFFEE.out.alignmentFile
            bucket = REG_PSICOFFEE.out.bucketSize
            homoplasy = REG_PSICOFFEE.out.homoplasyFile
            metrics = REG_PSICOFFEE.out.metricFile
        }    
        else if (align_method == "UPP"){
            REG_UPP(seqs_and_trees, bucket_size)

            alignment_method = REG_UPP.out.alignMethod
            tree_method = REG_UPP.out.treeMethod
            alignmentFile = REG_UPP.out.alignmentFile
            bucket = REG_UPP.out.bucketSize
            homoplasy = REG_UPP.out.homoplasyFile
            metrics = REG_UPP.out.metricFile
        } 
        //else{
        //    error "Invalid alignment mode: ${alignment_method}"
        //}

        flavour = "regressive"

    emit:  
        alignment_method
        tree_method
        alignmentFile
        bucket
        homoplasy 
        metrics
        flavour 
}

include { EVALUATION ; METRICS ; HOMOPLASY ; EASEL}         from './analysis_evaluation.nf'  

workflow REG_ANALYSIS {
  take:
    seqs_and_trees
    refs_ch
    align_method
    tree_method
    bucket_size
    
  main: 
    //println align_method 

    for(align in align_method){
        //println align 
        ALIGNMENT(seqs_and_trees,align,tree_method,bucket_size)

        if (params.evaluate){       
            EVALUATION(refs_ch, 
                ALIGNMENT.out.alignmentFile, 
                ALIGNMENT.out.flavour, 
                align,
                ALIGNMENT.out.tree_method, 
                ALIGNMENT.out.bucket)
        }
        if (params.metrics){
            METRICS(ALIGNMENT.out.alignmentFile, 
                ALIGNMENT.out.flavour, 
                align,
                ALIGNMENT.out.tree_method, 
                ALIGNMENT.out.bucket,
                ALIGNMENT.out.metrics)
        }
        if (params.homoplasy){  
            HOMOPLASY( ALIGNMENT.out.alignmentFile, 
                ALIGNMENT.out.flavour, 
                align,
                ALIGNMENT.out.tree_method, 
                ALIGNMENT.out.bucket,
                ALIGNMENT.out.homoplasy)
        }
        if (params.easel){ 
            EASEL( ALIGNMENT.out.alignmentFile, 
                ALIGNMENT.out.flavour, 
                align,
                ALIGNMENT.out.tree_method, 
                ALIGNMENT.out.bucket) 
        }
    }
}