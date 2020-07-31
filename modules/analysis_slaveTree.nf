#!/bin/bash nextflow

include { SLAVE_CLUSTALO ; SLAVE_FAMSA ; 
            SLAVE_MAFFT_FFTNS1 ; SLAVE_MAFFT_GINSI ; 
            SLAVE_MAFFT_SPARSECORE ; SLAVE_MAFFT ; 
            SLAVE_MSAPROBS ; SLAVE_MUSCLE ; 
            SLAVE_PROBCONS ; SLAVE_TCOFFEE ; SLAVE_UPP}     from './modules_slaveTreeAlignment.nf'   

workflow ALIGNMENT {
    take:
        seqs_and_trees
        align_method
        tree_method
        slave_method
        bucket_size

    main:
        //align_method = align_method.toString().replace("[", "").replace("]", "")
        //println align_method    

        if (align_method == "CLUSTALO"){
            SLAVE_CLUSTALO(seqs_and_trees, bucket_size,slave_method)

            alignment_method = SLAVE_CLUSTALO.out.alignMethod
            tree_method = SLAVE_CLUSTALO.out.treeMethod
            alignmentFile = SLAVE_CLUSTALO.out.alignmentFile
            bucket = SLAVE_CLUSTALO.out.bucketSize
            homoplasy = SLAVE_CLUSTALO.out.homoplasyFile
            metrics = SLAVE_CLUSTALO.out.metricFile
        }
        else if (align_method == "FAMSA"){
            SLAVE_FAMSA(seqs_and_trees, bucket_size,slave_method)

            alignment_method = SLAVE_FAMSA.out.alignMethod
            tree_method = SLAVE_FAMSA.out.treeMethod
            alignmentFile = SLAVE_FAMSA.out.alignmentFile
            bucket = SLAVE_FAMSA.out.bucketSize
            homoplasy = SLAVE_FAMSA.out.homoplasyFile
            metrics = SLAVE_FAMSA.out.metricFile
        }
        else if (align_method == "MAFFT-FFTNS1"){
            SLAVE_MAFFT_FFTNS1(seqs_and_trees, bucket_size,slave_method)

            alignment_method = SLAVE_MAFFT_FFTNS1.out.alignMethod
            tree_method = SLAVE_MAFFT_FFTNS1.out.treeMethod
            alignmentFile = SLAVE_MAFFT_FFTNS1.out.alignmentFile
            bucket = SLAVE_MAFFT_FFTNS1.out.bucketSize
            homoplasy = SLAVE_MAFFT_FFTNS1.out.homoplasyFile
            metrics = SLAVE_MAFFT_FFTNS1.out.metricFile
        }
        else if (align_method == "MAFFT-GINSI"){
            SLAVE_MAFFT_GINSI(seqs_and_trees, bucket_size,slave_method)

            alignment_method = SLAVE_MAFFT_GINSI.out.alignMethod
            tree_method = SLAVE_MAFFT_GINSI.out.treeMethod
            alignmentFile = SLAVE_MAFFT_GINSI.out.alignmentFile
            bucket = SLAVE_MAFFT_GINSI.out.bucketSize
            homoplasy = SLAVE_MAFFT_GINSI.out.homoplasyFile
            metrics = SLAVE_MAFFT_GINSI.out.metricFile
        }
        else if (align_method == "MAFFT-SPARSECORE"){
            SLAVE_MAFFT_SPARSECORE(seqs_and_trees, bucket_size,slave_method)

            alignment_method = SLAVE_MAFFT_SPARSECORE.out.alignMethod
            tree_method = SLAVE_MAFFT_SPARSECORE.out.treeMethod
            alignmentFile = SLAVE_MAFFT_SPARSECORE.out.alignmentFile
            bucket = SLAVE_MAFFT_SPARSECORE.out.bucketSize
            homoplasy = SLAVE_MAFFT_SPARSECORE.out.homoplasyFile
            metrics = SLAVE_MAFFT_SPARSECORE.out.metricFile
        }
        else if (align_method == "MAFFT"){
            SLAVE_MAFFT(seqs_and_trees, bucket_size,slave_method)

            alignment_method = SLAVE_MAFFT.out.alignMethod
            tree_method = SLAVE_MAFFT.out.treeMethod
            alignmentFile = SLAVE_MAFFT.out.alignmentFile
            bucket = SLAVE_MAFFT.out.bucketSize
            homoplasy = SLAVE_MAFFT.out.homoplasyFile
            metrics = SLAVE_MAFFT.out.metricFile
        }
        else if (align_method == "MSAPROBS"){
            SLAVE_MSAPROBS(seqs_and_trees, bucket_size,slave_method)

            alignment_method = SLAVE_MSAPROBS.out.alignMethod
            tree_method = SLAVE_MSAPROBS.out.treeMethod
            alignmentFile = SLAVE_MSAPROBS.out.alignmentFile
            bucket = SLAVE_MSAPROBS.out.bucketSize
            homoplasy = SLAVE_MSAPROBS.out.homoplasyFile
            metrics = SLAVE_MSAPROBS.out.metricFile
        }
        else if (align_method == "MUSCLE"){
            SLAVE_MUSCLE(seqs_and_trees, bucket_size,slave_method)

            alignment_method = SLAVE_MUSCLE.out.alignMethod
            tree_method = SLAVE_MUSCLE.out.treeMethod
            alignmentFile = SLAVE_MUSCLE.out.alignmentFile
            bucket = SLAVE_MUSCLE.out.bucketSize
            homoplasy = SLAVE_MUSCLE.out.homoplasyFile
            metrics = SLAVE_MUSCLE.out.metricFile
        }
        else if (align_method == "PROBCONS"){
            SLAVE_PROBCONS(seqs_and_trees, bucket_size,slave_method)

            alignment_method = SLAVE_PROBCONS.out.alignMethod
            tree_method = SLAVE_PROBCONS.out.treeMethod
            alignmentFile = SLAVE_PROBCONS.out.alignmentFile
            bucket = SLAVE_PROBCONS.out.bucketSize
            homoplasy = SLAVE_PROBCONS.out.homoplasyFile
            metrics = SLAVE_PROBCONS.out.metricFile
        }     
        else if (align_method == "TCOFFEE"){
            SLAVE_TCOFFEE(seqs_and_trees, bucket_size,slave_method)

            alignment_method = SLAVE_TCOFFEE.out.alignMethod
            tree_method = SLAVE_TCOFFEE.out.treeMethod
            alignmentFile = SLAVE_TCOFFEE.out.alignmentFile
            bucket = SLAVE_TCOFFEE.out.bucketSize
            homoplasy = SLAVE_TCOFFEE.out.homoplasyFile
            metrics = SLAVE_TCOFFEE.out.metricFile
        }    
        else if (align_method == "UPP"){
            SLAVE_UPP(seqs_and_trees, bucket_size,slave_method)

            alignment_method = SLAVE_UPP.out.alignMethod
            tree_method = SLAVE_UPP.out.treeMethod
            alignmentFile = SLAVE_UPP.out.alignmentFile
            bucket = SLAVE_UPP.out.bucketSize
            homoplasy = SLAVE_UPP.out.homoplasyFile
            metrics = SLAVE_UPP.out.metricFile
        } 
        //else{
        //    error "Invalid alignment mode: ${alignment_method}"
        //}

        flavour = "slave"

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

workflow SLAVE_ANALYSIS {
  take:
    seqs_and_trees
    refs_ch
    align_method
    tree_method
    bucket_size
    slave_method
    
  main: 
    //println align_method 

    for(align in align_method){
        //println align 
        for (tree in slave_method){
            ALIGNMENT(seqs_and_trees,align,tree_method,bucket_size,tree)

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
}