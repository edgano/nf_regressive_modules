#!/bin/bash nextflow

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

        if (align_method == "TCOFFEE"){
            TCOFFEE(seqs_and_trees)

            alignment_method = TCOFFEE.out.alignMethod
            tree_method = TCOFFEE.out.treeMethod
            alignmentFile = TCOFFEE.out.alignmentFile
            metrics = TCOFFEE.out.metricFile
        }
        if (align_method == "TCOFFEE_QUICKALN"){
            TCOFFEE_QUICKALN(seqs_and_trees)

            alignment_method = TCOFFEE_QUICKALN.out.alignMethod
            tree_method = TCOFFEE_QUICKALN.out.treeMethod
            alignmentFile = TCOFFEE_QUICKALN.out.alignmentFile
            metrics = TCOFFEE_QUICKALN.out.metricFile
        }
        if (align_method == "TCOFFEE_MCOFFEE"){
            TCOFFEE_MCOFFEE(seqs_and_trees)

            alignment_method = TCOFFEE_MCOFFEE.out.alignMethod
            tree_method = TCOFFEE_MCOFFEE.out.treeMethod
            alignmentFile = TCOFFEE_MCOFFEE.out.alignmentFile
            metrics = TCOFFEE_MCOFFEE.out.metricFile
        }
        if (align_method == "TCOFFEE_ACCURATE"){
            TCOFFEE_ACCURATE(seqs_and_trees)

            alignment_method = TCOFFEE_ACCURATE.out.alignMethod
            tree_method = TCOFFEE_ACCURATE.out.treeMethod
            alignmentFile = TCOFFEE_ACCURATE.out.alignmentFile
            metrics = TCOFFEE_ACCURATE.out.metricFile
        }
        if (align_method == "TCOFFEE_FMCOFFEE"){
            TCOFFEE_FMCOFFEE(seqs_and_trees)

            alignment_method = TCOFFEE_FMCOFFEE.out.alignMethod
            tree_method = TCOFFEE_FMCOFFEE.out.treeMethod
            alignmentFile = TCOFFEE_FMCOFFEE.out.alignmentFile
            metrics = TCOFFEE_FMCOFFEE.out.metricFile
        }
        if (align_method == "TCOFFEE_PSICOFFEE"){
            TCOFFEE_PSICOFFEE(seqs_and_trees)

            alignment_method = TCOFFEE_PSICOFFEE.out.alignMethod
            tree_method = TCOFFEE_PSICOFFEE.out.treeMethod
            alignmentFile = TCOFFEE_PSICOFFEE.out.alignmentFile
            metrics = TCOFFEE_PSICOFFEE.out.metricFile
        }
        if (align_method == "TCOFFEE_EXPRESSO"){
            TCOFFEE_EXPRESSO(seqs_and_trees)

            alignment_method = TCOFFEE_EXPRESSO.out.alignMethod
            tree_method = TCOFFEE_EXPRESSO.out.treeMethod
            alignmentFile = TCOFFEE_EXPRESSO.out.alignmentFile
            metrics = TCOFFEE_EXPRESSO.out.metricFile
        }
        if (align_method == "TCOFFEE_PROCOFFEE"){
            TCOFFEE_PROCOFFEE(seqs_and_trees)

            alignment_method = TCOFFEE_PROCOFFEE.out.alignMethod
            tree_method = TCOFFEE_PROCOFFEE.out.treeMethod
            alignmentFile = TCOFFEE_PROCOFFEE.out.alignmentFile
            metrics = TCOFFEE_PROCOFFEE.out.metricFile
        }
        if (align_method == "TCOFFEE_3DCOFFEE"){
            TCOFFEE_3DCOFFEE(seqs_and_trees)

            alignment_method = TCOFFEE_3DCOFFEE.out.alignMethod
            tree_method = TCOFFEE_3DCOFFEE.out.treeMethod
            alignmentFile = TCOFFEE_3DCOFFEE.out.alignmentFile
            metrics = TCOFFEE_3DCOFFEE.out.metricFile
        }     
        if (align_method == "TCOFFEE_TRMSD"){
            TCOFFEE_TRMSD(seqs_and_trees)       //tuple val(id), file(seqs), file(template), file(library)

            alignment_method = TCOFFEE_TRMSD.out.alignMethod
            tree_method = TCOFFEE_TRMSD.out.treeMethod
            alignmentFile = TCOFFEE_TRMSD.out.alignmentFile
            metrics = TCOFFEE_TRMSD.out.metricFile
        }  
        if (align_method == "TCOFFEE_RCOFFEE"){
            TCOFFEE_RCOFFEE(seqs_and_trees)

            alignment_method = TCOFFEE_RCOFFEE.out.alignMethod
            tree_method = TCOFFEE_RCOFFEE.out.treeMethod
            alignmentFile = TCOFFEE_RCOFFEE.out.alignmentFile
            metrics = TCOFFEE_RCOFFEE.out.metricFile
        }
        if (align_method == "TCOFFEE_RCOFFEE_CONSAN"){
            TCOFFEE_RCOFFEE_CONSAN(seqs_and_trees)

            alignment_method = TCOFFEE_RCOFFEE_CONSAN.out.alignMethod
            tree_method = TCOFFEE_RCOFFEE_CONSAN.out.treeMethod
            alignmentFile = TCOFFEE_RCOFFEE_CONSAN.out.alignmentFile
            metrics = TCOFFEE_RCOFFEE_CONSAN.out.metricFile
        }
        if (align_method == "TCOFFEE_3DALIGN"){
            TCOFFEE_3DALIGN(seqs_and_trees)

            alignment_method = TCOFFEE_3DALIGN.out.alignMethod
            tree_method = TCOFFEE_3DALIGN.out.treeMethod
            alignmentFile = TCOFFEE_3DALIGN.out.alignmentFile
            metrics = TCOFFEE_3DALIGN.out.metricFile
        }     
        if (align_method == "TCOFFEE_3DMALIGN"){
            TCOFFEE_3DMALIGN(seqs_and_trees)       //tuple val(id), file(seqs), file(template), file(library)

            alignment_method = TCOFFEE_3DMALIGN.out.alignMethod
            tree_method = TCOFFEE_3DMALIGN.out.treeMethod
            alignmentFile = TCOFFEE_3DMALIGN.out.alignmentFile
            metrics = TCOFFEE_3DMALIGN.out.metricFile
        }    

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
    for(align in align_method){
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