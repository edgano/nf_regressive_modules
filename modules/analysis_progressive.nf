#!/bin/bash nextflow

include { PROG_CLUSTALO ; PROG_FAMSA ; PROG_MAFFT_FFTNS1 ; PROG_MAFFT_GINSI ; PROG_MAFFT_SPARSECORE ; PROG_MAFFT ; PROG_MSAPROBS ; PROG_MUSCLE ; PROG_PROBCONS ; PROG_UPP}     from './modules_progressiveAlignment.nf' 
include { TCOFFEE_DEFAULT ; TCOFFEE_QUICKALN ; TCOFFEE_MCOFFEE ; TCOFFEE_ACCURATE ; TCOFFEE_FMCOFFEE ; TCOFFEE_PSICOFFEE ; TCOFFEE_EXPRESSO ; TCOFFEE_PROCOFFEE ; TCOFFEE_3DCOFFEE ; TCOFFEE_TRMSD ; TCOFFEE_RCOFFEE ; TCOFFEE_RCOFFEE_CONSAN ; TCOFFEE_3DALIGN ; TCOFFEE_3DMALIGN }     from './modules_tcoffeeModes.nf' 

workflow ALIGNMENT {
    take:
        seqs_and_trees
        align_method
        tree_method
        bucket_size

    main:
        align_method = align_method.toString().replace("[", "").replace("]", "")
        //println align_method    

        if (align_method == "CLUSTALO"){
            PROG_CLUSTALO(seqs_and_trees)

            alignment_method = PROG_CLUSTALO.out.alignMethod
            tree_method = PROG_CLUSTALO.out.treeMethod
            alignmentFile = PROG_CLUSTALO.out.alignmentFile
            metrics = PROG_CLUSTALO.out.metricFile
        }
        if (align_method == "FAMSA"){
            PROG_FAMSA(seqs_and_trees)

            alignment_method = PROG_FAMSA.out.alignMethod
            tree_method = PROG_FAMSA.out.treeMethod
            alignmentFile = PROG_FAMSA.out.alignmentFile
            metrics = PROG_FAMSA.out.metricFile
        }
        if (align_method == "MAFFT-FFTNS1"){
            PROG_MAFFT_FFTNS1(seqs_and_trees)

            alignment_method = PROG_MAFFT_FFTNS1.out.alignMethod
            tree_method = PROG_MAFFT_FFTNS1.out.treeMethod
            alignmentFile = PROG_MAFFT_FFTNS1.out.alignmentFile
            metrics = PROG_MAFFT_FFTNS1.out.metricFile
        }
        if (align_method == "MAFFT-GINSI"){
            PROG_MAFFT_GINSI(seqs_and_trees)

            alignment_method = PROG_MAFFT_GINSI.out.alignMethod
            tree_method = PROG_MAFFT_GINSI.out.treeMethod
            alignmentFile = PROG_MAFFT_GINSI.out.alignmentFile
            metrics = PROG_MAFFT_GINSI.out.metricFile
        }
        if (align_method == "MAFFT-SPARSECORE"){
            PROG_MAFFT_SPARSECORE(seqs_and_trees)

            alignment_method = PROG_MAFFT_SPARSECORE.out.alignMethod
            tree_method = PROG_MAFFT_SPARSECORE.out.treeMethod
            alignmentFile = PROG_MAFFT_SPARSECORE.out.alignmentFile
            metrics = PROG_MAFFT_SPARSECORE.out.metricFile
        }
        if (align_method == "MAFFT"){
            PROG_MAFFT(seqs_and_trees)

            alignment_method = PROG_MAFFT.out.alignMethod
            tree_method = PROG_MAFFT.out.treeMethod
            alignmentFile = PROG_MAFFT.out.alignmentFile
            metrics = PROG_MAFFT.out.metricFile
        }
        if (align_method == "MSAPROBS"){
            PROG_MSAPROBS(seqs_and_trees)

            alignment_method = PROG_MSAPROBS.out.alignMethod
            tree_method = PROG_MSAPROBS.out.treeMethod
            alignmentFile = PROG_MSAPROBS.out.alignmentFile
            metrics = PROG_MSAPROBS.out.metricFile
        }
        if (align_method == "MUSCLE"){
            PROG_MUSCLE(seqs_and_trees)

            alignment_method = PROG_MUSCLE.out.alignMethod
            tree_method = PROG_MUSCLE.out.treeMethod
            alignmentFile = PROG_MUSCLE.out.alignmentFile
            metrics = PROG_MUSCLE.out.metricFile
        }
        if (align_method == "PROBCONS"){
            PROG_PROBCONS(seqs_and_trees)

            alignment_method = PROG_PROBCONS.out.alignMethod
            tree_method = PROG_PROBCONS.out.treeMethod
            alignmentFile = PROG_PROBCONS.out.alignmentFile
            metrics = PROG_PROBCONS.out.metricFile
        }     
        /*if (align_method == "PSICOFFEE"){
            TCOFFEE_PSICOFFEE(seqs_and_trees)       //tuple val(id), file(seqs), file(template), file(library)

            alignment_method = TCOFFEE_PSICOFFEE.out.alignMethod
            tree_method = TCOFFEE_PSICOFFEE.out.treeMethod
            alignmentFile = TCOFFEE_PSICOFFEE.out.alignmentFile
            metrics = TCOFFEE_PSICOFFEE.out.metricFile
        } */
        if (align_method == "UPP"){
            PROG_UPP(seqs_and_trees)       //tuple val(id), file(seqs), file(template), file(library)

            alignment_method = PROG_UPP.out.alignMethod
            tree_method = PROG_UPP.out.treeMethod
            alignmentFile = PROG_UPP.out.alignmentFile
            metrics = PROG_UPP.out.metricFile
        }    
        //else{
        //    error "Invalid alignment method: ${align_method}" 
        //}       
    flavour = "progressive"
        
    emit:  
        alignment_method
        tree_method
        alignmentFile
        metrics
        flavour 
}

include { EVALUATION ; METRICS ; GAPS ; EASEL}         from './analysis_evaluation.nf'  
workflow PROG_ANALYSIS {
  take:
    seqs_and_trees
    refs_ch
    align_method
    tree_method
    bucket_size
    
   main: 
    //align_method = align_method.toString().replace("[", "").replace("]", "")
    //println align_method 

    for(align in align_method){
        println align
        ALIGNMENT(seqs_and_trees,align,tree_method,bucket_size)

        if (params.evaluate){       
            EVALUATION(refs_ch, 
                ALIGNMENT.out.alignmentFile, 
                ALIGNMENT.out.flavour, 
                align,
                ALIGNMENT.out.tree_method, 
                "NA")
        }
        if (params.metrics){
            METRICS(ALIGNMENT.out.alignmentFile, 
                ALIGNMENT.out.flavour, 
                align,
                ALIGNMENT.out.tree_method, 
                "NA",
                ALIGNMENT.out.metrics)
        }
        if (params.gapCount){  
            GAPS(ALIGNMENT.out.alignmentFile, 
                ALIGNMENT.out.flavour, 
                align,
                ALIGNMENT.out.tree_method, 
                "NA")
        }
        if (params.easel){ 
            EASEL( ALIGNMENT.out.alignmentFile, 
                ALIGNMENT.out.flavour, 
                align,
                ALIGNMENT.out.tree_method, 
                "NA") 
        }
    }
}