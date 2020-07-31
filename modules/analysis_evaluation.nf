#!/bin/bash nextflow

include {   EVAL_ALIGNMENT as EVAL_CLUSTALO ; 
            EVAL_ALIGNMENT as EVAL_FAMSA ; 
            EVAL_ALIGNMENT as EVAL_FFTNS1 ;
            EVAL_ALIGNMENT as EVAL_MAFFT_GINSI ; 
            EVAL_ALIGNMENT as EVAL_MAFFT_SPARSECORE ; 
            EVAL_ALIGNMENT as EVAL_MAFFT ;
            EVAL_ALIGNMENT as EVAL_MSAPROBS ; 
            EVAL_ALIGNMENT as EVAL_MUSCLE ; 
            EVAL_ALIGNMENT as EVAL_PROBCONS ;
            EVAL_ALIGNMENT as EVAL_PSICOFFEE ;
            EVAL_ALIGNMENT as EVAL_UPP} from './modules_evaluateAlignment.nf'  

workflow EVALUATION {
    take:
        refs_ch
        aligmentFile
        flavour
        align_method
        tree_method
        bucket_size

    main:
        refs_ch
            .cross (aligmentFile)
            .map { it -> [ it[1][0], it[1][1], it[0][1] ] }
            .set { alignment_and_ref }  

        if (align_method == "CLUSTALO"){
            EVAL_CLUSTALO(flavour, alignment_and_ref, align_method, tree_method, bucket_size)

            tcScore = EVAL_CLUSTALO.out.tcScore
            spScore = EVAL_CLUSTALO.out.spScore
            colScore = EVAL_CLUSTALO.out.colScore
        } 
        if (align_method == "FAMSA"){
            EVAL_FAMSA(flavour, alignment_and_ref, align_method, tree_method, bucket_size)

            tcScore = EVAL_FAMSA.out.tcScore
            spScore = EVAL_FAMSA.out.spScore
            colScore = EVAL_FAMSA.out.colScore
        } 
        if (align_method == "MAFFT-FFTNS1"){
            EVAL_FFTNS1(flavour, alignment_and_ref, align_method, tree_method, bucket_size)

            tcScore = EVAL_FFTNS1.out.tcScore
            spScore = EVAL_FFTNS1.out.spScore
            colScore = EVAL_FFTNS1.out.colScore
        }
        if (align_method == "MAFFT-GINSI"){
            EVAL_MAFFT_GINSI(flavour, alignment_and_ref, align_method, tree_method, bucket_size)

            tcScore = EVAL_MAFFT_GINSI.out.tcScore
            spScore = EVAL_MAFFT_GINSI.out.spScore
            colScore = EVAL_MAFFT_GINSI.out.colScore
        }
        if (align_method == "MAFFT-SPARSECORE"){
            EVAL_MAFFT_SPARSECORE(flavour, alignment_and_ref, align_method, tree_method, bucket_size)

            tcScore = EVAL_MAFFT_SPARSECORE.out.tcScore
            spScore = EVAL_MAFFT_SPARSECORE.out.spScore
            colScore = EVAL_MAFFT_SPARSECORE.out.colScore            
        }
        if (align_method == "MAFFT"){
            EVAL_MAFFT(flavour, alignment_and_ref, align_method, tree_method, bucket_size)

            tcScore = EVAL_MAFFT.out.tcScore
            spScore = EVAL_MAFFT.out.spScore
            colScore = EVAL_MAFFT.out.colScore            
        }
        if (align_method == "MSAPROBS"){
            EVAL_MSAPROBS(flavour, alignment_and_ref, align_method, tree_method, bucket_size)

            tcScore = EVAL_MSAPROBS.out.tcScore
            spScore = EVAL_MSAPROBS.out.spScore
            colScore = EVAL_MSAPROBS.out.colScore   
        }
        if (align_method == "MUSCLE"){
            EVAL_MUSCLE(flavour, alignment_and_ref, align_method, tree_method, bucket_size)

            tcScore = EVAL_MUSCLE.out.tcScore
            spScore = EVAL_MUSCLE.out.spScore
            colScore = EVAL_MUSCLE.out.colScore   
        }
        if (align_method == "PROBCONS"){
            EVAL_PROBCONS(flavour, alignment_and_ref, align_method, tree_method, bucket_size)

            tcScore = EVAL_PROBCONS.out.tcScore
            spScore = EVAL_PROBCONS.out.spScore
            colScore = EVAL_PROBCONS.out.colScore   
        }
        /*if (align_method == "PSICOFFEE"){
            EVAL_PSICOFFEE(flavour, alignment_and_ref, align_method, tree_method, bucket_size)

            tcScore = EVAL_PSICOFFEE.out.tcScore
            spScore = EVAL_PSICOFFEE.out.spScore
            colScore = EVAL_PSICOFFEE.out.colScore   
        }*/
        if (align_method == "UPP"){
            EVAL_UPP(flavour, alignment_and_ref, align_method, tree_method, bucket_size)

            tcScore = EVAL_UPP.out.tcScore
            spScore = EVAL_UPP.out.spScore
            colScore = EVAL_UPP.out.colScore   
        }
        //else{
        //    error "Invalid alignment method: ${align_method}" 
        //}    
    emit:  
        tcScore 
        spScore  
        colScore
}

include {   METRICS as METRIC_CLUSTALO ; 
            METRICS as METRIC_FAMSA ; 
            METRICS as METRIC_FFTNS1 ;
            METRICS as METRIC_MAFFT_GINSI ; 
            METRICS as METRIC_MAFFT_SPARSECORE ; 
            METRICS as METRIC_MAFFT ;
            METRICS as METRIC_MSAPROBS ; 
            METRICS as METRIC_MUSCLE ; 
            METRICS as METRIC_PROBCONS ;
            METRICS as METRIC_PSICOFFEE } from './modules_evaluateAlignment.nf'  

workflow METRICS {
    take:
        aligmentFile
        flavour
        align_method
        tree_method
        bucket_size
        metricFile

    main:
        if (align_method == "CLUSTALO"){
            METRIC_CLUSTALO(  flavour, 
                        aligmentFile, 
                        align_method, 
                        tree_method, 
                        bucket_size, 
                        metricFile)

            metricFiles = METRIC_CLUSTALO.out.metricFiles
        } 
        if (align_method == "FAMSA"){
            METRIC_FAMSA(  flavour, 
                        aligmentFile, 
                        align_method, 
                        tree_method, 
                        bucket_size, 
                        metricFile)
            metricFiles = METRIC_FAMSA.out.metricFiles

        } 
        if (align_method == "MAFFT-FFTNS1"){
            METRIC_FFTNS1(  flavour, 
                        aligmentFile, 
                        align_method, 
                        tree_method, 
                        bucket_size, 
                        metricFile) 
            
            metricFiles = METRIC_FFTNS1.out.metricFiles
        }
        if (align_method == "MAFFT-GINSI"){
            METRIC_GINSI(  flavour, 
                        aligmentFile, 
                        align_method, 
                        tree_method, 
                        bucket_size, 
                        metricFile) 
            
            metricFiles = METRIC_GINSI.out.metricFiles
        }        
        if (align_method == "MAFFT-SPARSECORE"){
            METRIC_SPARSECORE(  flavour, 
                        aligmentFile, 
                        align_method, 
                        tree_method, 
                        bucket_size, 
                        metricFile) 
            
            metricFiles = METRIC_SPARSECORE.out.metricFiles
        }
        if (align_method == "MAFFT"){
            METRIC_MAFFT(  flavour, 
                        aligmentFile, 
                        align_method, 
                        tree_method, 
                        bucket_size, 
                        metricFile) 
            
            metricFiles = METRIC_MAFFT.out.metricFiles
        }
        if (align_method == "MSAPROBS"){
            METRIC_MSAPROBS(  flavour, 
                        aligmentFile, 
                        align_method, 
                        tree_method, 
                        bucket_size, 
                        metricFile) 
            
            metricFiles = METRIC_MSAPROBS.out.metricFiles
        }
        if (align_method == "MUSCLE"){           
            METRIC_MUSCLE(  flavour, 
                        aligmentFile, 
                        align_method, 
                        tree_method, 
                        bucket_size, 
                        metricFile) 
            
            metricFiles = METRIC_MUSCLE.out.metricFiles
        }        
        if (align_method == "PROBCONS"){
            METRIC_PROBCONS(  flavour, 
                        aligmentFile, 
                        align_method, 
                        tree_method, 
                        bucket_size, 
                        metricFile) 
            
            metricFiles = METRIC_PROBCONS.out.metricFiles
        }        
        if (align_method == "PSICOFFEE"){
            METRIC_PSICOFFEE(  flavour, 
                        aligmentFile, 
                        align_method, 
                        tree_method, 
                        bucket_size, 
                        metricFile) 
            
            metricFiles = METRIC_PSICOFFEE.out.metricFiles
        }  

    emit:
        metricFiles
}

include {   HOMOPLASY as HOMO_CLUSTALO ; 
            HOMOPLASY as HOMO_FAMSA ; 
            HOMOPLASY as HOMO_MAFFT_FFTNS1 ;
            HOMOPLASY as HOMO_MAFFT_GINSI ; 
            HOMOPLASY as HOMO_MAFFT_SPARSECORE ; 
            HOMOPLASY as HOMO_MAFFT ;
            HOMOPLASY as HOMO_MSAPROBS ; 
            HOMOPLASY as HOMO_MUSCLE ; 
            HOMOPLASY as HOMO_PROBCONS ;
            HOMOPLASY as HOMO_PSICOFFEE } from './modules_evaluateAlignment.nf'  
workflow HOMOPLASY {
    take:
        aligmentFile
        flavour
        align_method
        tree_method
        bucket_size
        homoplasyFile

    main:
        if (align_method == "CLUSTALO"){
            HOMO_CLUSTALO(  flavour, 
                        aligmentFile, 
                        align_method, 
                        tree_method, 
                        bucket_size, 
                        homoplasyFile)

            homoFiles = HOMO_CLUSTALO.out.homoFiles
        } 
        if (align_method == "FAMSA"){
            HOMO_FAMSA(  flavour, 
                        aligmentFile, 
                        align_method, 
                        tree_method, 
                        bucket_size, 
                        homoplasyFile)
            
            homoFiles = HOMO_FAMSA.out.homoFiles
        } 
        if (align_method == "MAFFT-FFTNS1"){
            HOMO_MAFFT_FFTNS1(  flavour, 
                        aligmentFile, 
                        align_method, 
                        tree_method, 
                        bucket_size, 
                        homoplasyFile)  

            homoFiles = HOMO_MAFFT_FFTNS1.out.homoFiles
        }
        if (align_method == "MAFFT-GINSI"){
            HOMO_MAFFT_GINSI(  flavour, 
                        aligmentFile, 
                        align_method, 
                        tree_method, 
                        bucket_size, 
                        homoplasyFile) 

            homoFiles = HOMO_MAFFT_GINSI.out.homoFiles
        }
        if (align_method == "MAFFT-SPARSECORE"){
            HOMO_MAFFT_SPARSECORE(  flavour, 
                        aligmentFile, 
                        align_method, 
                        tree_method, 
                        bucket_size, 
                        homoplasyFile) 

            homoFiles = HOMO_MAFFT_SPARSECORE.out.homoFiles
        }
        if (align_method == "MAFFT"){
            HOMO_MAFFT(  flavour, 
                        aligmentFile, 
                        align_method, 
                        tree_method, 
                        bucket_size, 
                        homoplasyFile) 

            homoFiles = HOMO_MAFFT.out.homoFiles
        }            
        if (align_method == "MSAPROBS"){   
            HOMO_MSAPROBS(  flavour, 
                        aligmentFile, 
                        align_method, 
                        tree_method, 
                        bucket_size, 
                        homoplasyFile) 

            homoFiles = HOMO_MSAPROBS.out.homoFiles
        }          
        if (align_method == "MUSCLE"){     
            HOMO_MUSCLE(  flavour, 
                        aligmentFile, 
                        align_method, 
                        tree_method, 
                        bucket_size, 
                        homoplasyFile) 

            homoFiles = HOMO_MUSCLE.out.homoFiles
        }            
        if (align_method == "PROBCONS"){
            HOMO_PROBCONS(  flavour, 
                        aligmentFile, 
                        align_method, 
                        tree_method, 
                        bucket_size, 
                        homoplasyFile) 

            homoFiles = HOMO_PROBCONS.out.homoFiles
        } 
        if (align_method == "PSICOFFEE"){
            HOMO_PSICOFFEE(  flavour, 
                        aligmentFile, 
                        align_method, 
                        tree_method, 
                        bucket_size, 
                        homoplasyFile) 

            homoFiles = HOMO_PSICOFFEE.out.homoFiles
        } 

    emit:
        homoFiles
}

include {   EASEL_INFO as EASEL_CLUSTALO ; 
            EASEL_INFO as EASEL_FAMSA ; 
            EASEL_INFO as EASEL_MAFFT_FFTNS1 ;
            EASEL_INFO as EASEL_MAFFT_GINSI ; 
            EASEL_INFO as EASEL_MAFFT_SPARSECORE ; 
            EASEL_INFO as EASEL_MAFFT ;
            EASEL_INFO as EASEL_MSAPROBS ; 
            EASEL_INFO as EASEL_MUSCLE ; 
            EASEL_INFO as EASEL_PROBCONS ;
            EASEL_INFO as EASEL_PSICOFFEE } from './modules_evaluateAlignment.nf'   
workflow EASEL {
    take:
        aligmentFile
        flavour
        align_method
        tree_method
        bucket_size

    main:
        if (align_method == "CLUSTALO"){
            EASEL_CLUSTALO(  flavour, 
                        aligmentFile, 
                        align_method, 
                        tree_method, 
                        bucket_size)

            easelFiles = EASEL_CLUSTALO.out.easelFiles
        } 
        if (align_method == "FAMSA"){
            EASEL_FAMSA(  flavour, 
                        aligmentFile, 
                        align_method, 
                        tree_method, 
                        bucket_size)

            easelFiles = EASEL_FAMSA.out.easelFiles
        } 
        if (align_method == "MAFFT-FFTNS1"){
            EASEL_MAFFT_FFTNS1(  flavour, 
                        aligmentFile, 
                        align_method, 
                        tree_method, 
                        bucket_size)

            easelFiles = EASEL_MAFFT_FFTNS1.out.easelFiles
        }
        if (align_method == "MAFFT-GINSI"){
            EASEL_MAFFT_GINSI(  flavour, 
                        aligmentFile, 
                        align_method, 
                        tree_method, 
                        bucket_size)

            easelFiles = EASEL_MAFFT_GINSI.out.easelFiles
        }
        if (align_method == "MAFFT-SPARSECORE"){
            EASEL_MAFFT_SPARSECORE(  flavour, 
                        aligmentFile, 
                        align_method, 
                        tree_method, 
                        bucket_size)

            easelFiles = EASEL_MAFFT_SPARSECORE.out.easelFiles
        }
        if (align_method == "MAFFT"){
            EASEL_MAFFT(  flavour, 
                        aligmentFile, 
                        align_method, 
                        tree_method, 
                        bucket_size)

            easelFiles = EASEL_MAFFT.out.easelFiles
        }            
        if (align_method == "MSAPROBS"){   
            EASEL_MSAPROBS(  flavour, 
                        aligmentFile, 
                        align_method, 
                        tree_method, 
                        bucket_size)

            easelFiles = EASEL_MSAPROBS.out.easelFiles
        }          
        if (align_method == "MUSCLE"){     
            EASEL_MUSCLE(  flavour, 
                        aligmentFile, 
                        align_method, 
                        tree_method, 
                        bucket_size)

            easelFiles = EASEL_MUSCLE.out.easelFiles
        }            
        if (align_method == "PROBCONS"){
            EASEL_PROBCONS(  flavour, 
                        aligmentFile, 
                        align_method, 
                        tree_method, 
                        bucket_size)

            easelFiles = EASEL_PROBCONS.out.easelFiles
        } 
        if (align_method == "PSICOFFEE"){
            EASEL_PSICOFFEE(  flavour, 
                        aligmentFile, 
                        align_method, 
                        tree_method, 
                        bucket_size)

            easelFiles = EASEL_PSICOFFEE.out.easelFiles
        }
}

include {   GAPS_PROGRESSIVE as GAPS_CLUSTALO ; 
            GAPS_PROGRESSIVE as GAPS_FAMSA ; 
            GAPS_PROGRESSIVE as GAPS_MAFFT_FFTNS1 ;
            GAPS_PROGRESSIVE as GAPS_MAFFT_GINSI ; 
            GAPS_PROGRESSIVE as GAPS_MAFFT_SPARSECORE ; 
            GAPS_PROGRESSIVE as GAPS_MAFFT ;
            GAPS_PROGRESSIVE as GAPS_MSAPROBS ; 
            GAPS_PROGRESSIVE as GAPS_MUSCLE ; 
            GAPS_PROGRESSIVE as GAPS_PROBCONS ;
            GAPS_PROGRESSIVE as GAPS_PSICOFFEE } from './modules_evaluateAlignment.nf'   
workflow GAPS {
    take:
        aligmentFile
        flavour
        align_method
        tree_method
        bucket_size
    
    main:
        if (align_method == "CLUSTALO"){
            GAPS_CLUSTALO(  flavour, 
                        aligmentFile, 
                        align_method, 
                        tree_method, 
                        bucket_size)

            gapFiles = GAPS_CLUSTALO.out.gapFiles
        } 
        if (align_method == "FAMSA"){
            GAPS_FAMSA(  flavour, 
                        aligmentFile, 
                        align_method, 
                        tree_method, 
                        bucket_size)
                
            gapFiles = GAPS_FAMSA.out.gapFiles        
        } 
        if (align_method == "MAFFT-FFTNS1"){
            GAPS_MAFFT_FFTNS1(  flavour, 
                        aligmentFile, 
                        align_method, 
                        tree_method, 
                        bucket_size)
            gapFiles = GAPS_MAFFT_FFTNS1.out.gapFiles       
        }
        if (align_method == "MAFFT-GINSI"){
             GAPS_MAFFT_GINSI(  flavour, 
                        aligmentFile, 
                        align_method, 
                        tree_method, 
                        bucket_size)
            gapFiles = GAPS_MAFFT_GINSI.out.gapFiles              
        } 
        if (align_method == "MAFFT-SPARSECORE"){
            GAPS_MAFFT_SPARSECORE(  flavour, 
                        aligmentFile, 
                        align_method, 
                        tree_method, 
                        bucket_size)
            gapFiles = GAPS_MAFFT_SPARSECORE.out.gapFiles  
        }
        if (align_method == "MAFFT"){
            GAPS_MAFFT(  flavour, 
                        aligmentFile, 
                        align_method, 
                        tree_method, 
                        bucket_size)
            gapFiles = GAPS_MAFFT.out.gapFiles  
        }            
        if (align_method == "MSAPROBS"){   
            GAPS_MSAPROBS(  flavour, 
                        aligmentFile, 
                        align_method, 
                        tree_method, 
                        bucket_size)
            gapFiles = GAPS_MSAPROBS.out.gapFiles  
        }          
        if (align_method == "MUSCLE"){     
            GAPS_MUSCLE(  flavour, 
                        aligmentFile, 
                        align_method, 
                        tree_method, 
                        bucket_size)
            gapFiles = GAPS_MUSCLE.out.gapFiles  
        }            
        if (align_method == "PROBCONS"){
            GAPS_PROBCONS(  flavour, 
                        aligmentFile, 
                        align_method, 
                        tree_method, 
                        bucket_size)
            gapFiles = GAPS_PROBCONS.out.gapFiles  
        } 
        if (align_method == "PSICOFFEE"){
            GAPS_PSICOFFEE(  flavour, 
                        aligmentFile, 
                        align_method, 
                        tree_method, 
                        bucket_size)
            gapFiles = GAPS_PSICOFFEE.out.gapFiles  
        }  

    emit:
        gapFiles
}

include {   IRMSD as IRMSD_CLUSTALO ; 
            IRMSD as IRMSD_FAMSA ; 
            IRMSD as IRMSD_MAFFT_FFTNS1 ;
            IRMSD as IRMSD_MAFFT_GINSI ; 
            IRMSD as IRMSD_MAFFT_SPARSECORE ; 
            IRMSD as IRMSD_MAFFT ;
            IRMSD as IRMSD_MSAPROBS ; 
            IRMSD as IRMSD_MUSCLE ; 
            IRMSD as IRMSD_PROBCONS ;
            IRMSD as IRMSD_PSICOFFEE } from './modules_evaluateAlignment.nf'   
workflow IRMSD {
    take:
        aligmentFile
        flavour
        align_method
        tree_method
        bucket_size
    
    main:
        if (align_method == "CLUSTALO"){
            IRMSD_CLUSTALO(  flavour, 
                        aligmentFile, 
                        align_method, 
                        tree_method, 
                        bucket_size)

            irmsd = IRMSD_CLUSTALO.out.irmsd
        } 
        if (align_method == "FAMSA"){
            IRMSD_FAMSA(  flavour, 
                        aligmentFile, 
                        align_method, 
                        tree_method, 
                        bucket_size)
                
            irmsd = IRMSD_FAMSA.out.irmsd        
        } 
        if (align_method == "MAFFT-FFTNS1"){
            IRMSD_MAFFT_FFTNS1(  flavour, 
                        aligmentFile, 
                        align_method, 
                        tree_method, 
                        bucket_size)
            irmsd = IRMSD_MAFFT_FFTNS1.out.irmsd       
        }
        if (align_method == "MAFFT-GINSI"){
             IRMSD_MAFFT_GINSI(  flavour, 
                        aligmentFile, 
                        align_method, 
                        tree_method, 
                        bucket_size)
            irmsd = IRMSD_MAFFT_GINSI.out.irmsd              
        } 
        if (align_method == "MAFFT-SPARSECORE"){
            IRMSD_MAFFT_SPARSECORE(  flavour, 
                        aligmentFile, 
                        align_method, 
                        tree_method, 
                        bucket_size)
            irmsd = IRMSD_MAFFT_SPARSECORE.out.irmsd  
        }
        if (align_method == "MAFFT"){
            IRMSD_MAFFT(  flavour, 
                        aligmentFile, 
                        align_method, 
                        tree_method, 
                        bucket_size)
            irmsd = IRMSD_MAFFT.out.irmsd  
        }            
        if (align_method == "MSAPROBS"){   
            IRMSD_MSAPROBS(  flavour, 
                        aligmentFile, 
                        align_method, 
                        tree_method, 
                        bucket_size)
            irmsd = IRMSD_MSAPROBS.out.irmsd  
        }          
        if (align_method == "MUSCLE"){     
            IRMSD_MUSCLE(  flavour, 
                        aligmentFile, 
                        align_method, 
                        tree_method, 
                        bucket_size)
            irmsd = IRMSD_MUSCLE.out.irmsd  
        }            
        if (align_method == "PROBCONS"){
            IRMSD_PROBCONS(  flavour, 
                        aligmentFile, 
                        align_method, 
                        tree_method, 
                        bucket_size)
            irmsd = IRMSD_PROBCONS.out.irmsd  
        } 
        if (align_method == "PSICOFFEE"){
            IRMSD_PSICOFFEE(  flavour, 
                        aligmentFile, 
                        align_method, 
                        tree_method, 
                        bucket_size)
            irmsd = IRMSD_PSICOFFEE.out.irmsd  
        }  

    emit:
        irmsd
}