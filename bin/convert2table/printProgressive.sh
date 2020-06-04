#!/bin/bash

#################
## FAMILY NAME ##
#################

# >> HOMFAM
declare -a all=(seatoxin hip     scorptoxin      cyt3    rnasemam        bowman  toxin   ghf11   TNF     sti     Stap_Strp_toxin profilin        ricin   ghf22   ChtBD   ins     trfl    slectin phoslip ltn     il8     az      kringle cryst   DEATH   cah     mmp     rub     ghf10   tgfb    sodcu   KAS     DMRL_synthase   tms     GEL     kunitz  Sulfotransfer   mofe    Ald_Xan_dh_2    ghf5    phc     aadh    annexin serpin  cytb    asp     oxidored_q6     hpr     hormone_rec     hr      tim     glob    ace     cys     ghf1    sodfe   peroxidase      uce     flav    HMG_box OTCace  msb     icd     proteasome      cyclo   LIM     HLH     ldh     subt    int     lyase_1 gpdh    egf     blm     gluts   myb_DNA-binding tRNA-synt_2b    biotin_lipoyl   hom     ghf13   aldosered       hla     Rhodanese       PDZ     blmb    rhv     p450    adh     aat     rrm     Acetyltransf    sdr     zf-CCHH rvp)

# declare -a over1000=(ltn     il8     az      kringle cryst   DEATH   cah     mmp     rub     ghf10   tgfb    sodcu   KAS     DMRL_synthase   tms     GEL     kunitz  Sulfotransfer   mofe    Ald_Xan_dh_2    ghf5    phc     aadh    annexin serpin  cytb    asp     oxidored_q6     hpr     hormone_rec     hr      tim     glob    ace     cys     ghf1    sodfe   peroxidase      uce     flav    HMG_box OTCace  msb     icd     proteasome      cyclo   LIM     HLH     ldh     subt    int     lyase_1 gpdh    egf     blm )

# declare -a top20=(gluts   myb_DNA-binding tRNA-synt_2b    biotin_lipoyl   hom     ghf13   aldosered       hla     Rhodanese       PDZ     blmb    rhv     p450    adh     aat     rrm     Acetyltransf    sdr     zf-CCHH rvp)

# >> extHOMFAM
# declare -a all=(Cyclotide Colicin adeno_fiber PSI_PsaE  lta B DNA_PPF hip bowman  leukocidin  Colicin_Pyocin  ghf33 UcrQ  HTH_AraC  TMV_coat  E2_C  scorptoxin  mthina  mthinb  UCR_TM  FN1 cyt3  CBD_3 cks bv  Huristasin  pot planttoxin  Binary_toxA toxin ALBUMIN NTP_transf_2_2  ricin NO_synthase MCR_alpha_N ghf14 ghf6  rnase P Bcl-2 sti cytprime  MIF rnasemam  ins ghf11 EF1BD TNF ghf22 Propep_M14  ngf ets ENTH  profilin  ChtBD Stap_Strp_toxin laminin_G trfl  phoslip mmob  lamb  CBM_20  RuBisCO_small Cys-protease-3C il8 kazal Bacterial_PQQ cystatin  LuxS  pilin BIR DHHA2 TRANSTHYRETIN THAUMATIN intb  tbpc  Cu_bind_like  DISIN SAP pec_lyase tryp_alpha_amyl plantltp  cbm12 slectin PTB kringle rub az  GATA  aldose  ltn DEATH parv  cbp lipocalin fabp  tgfb  recombinase Invasin cryst svmp  ANK ADF fbpase  Fork_head phc ghf10 Bac_rhodopsin S_100 Band_41_M Band_41_N mmp tln YgbB  C1  hemopexin ppase igC1  MHC_II_C  igV igcon igvar-h igvar-l sodcu TIL porin DMRL_synthase Band_41_NME igI cah ATP-gua_Ptrans  Cu_amine_oxid hr  B_lectin  NiFeSe_Hases  RRF Fapy_DNA_glyco  RGS SRF-TF  hexokinase  mofe  oxidored_q6 alpha-amylase_C fibrinogen_C  Ribosomal_L7Ae  ghf7  rnr tms dhfr  acid_phosphat gpr sugbp hpr glob  WW  Pro_CA  ndk chromo  GEL Lipase_3  kunitz  Filamin Ribosomal_L1  KAS glnasn  lipase  MHC_II_alpha_NC Sulfotransfer Adenylsucc_synt Sm  RNase_HII ligase-CoA  hormone_rec phs pgk dsrm  kex RING  cat HMG_box DUF170  annexin sh2 serpin  asp dutpase Glyco_hydro_2 cpa cytb  Lum_binding DHOdehase ptpase  ghf5  bromodomain ACPS  mle TIG Bac_DNA_binding fer4  DHH PEP-utilizers Ribosomal_L6_D  tim EF_TS UBQ LMWPc MCR_alpha sh3 Asp_Glu_race_D  histone NADHdh_2  MM_CoA_mutase Orn_DAP_Arg_deC flav  DAHP_synth_1  tRNA_bind ferritin  CH  cys MCR_alpha_NC  ARM HLH csp Methyltransf_1  citrate_synt  igps  galpha  cat3  int ung FAD-oxidase_C peroxidase  SBP_bacterial_1 cyclin  uce LIM CUB cyt5  cytc  Cu_nir  Glyco_hydro_18  Glyco_hydro_18_D1 Glyco_hydro_18_D2 ghf18 DNA_pol3_beta LRR inositol_P  proteasome  UPF0076 msb PK  MoCF_biosynth sodfe aadh  KISc  ace hom EFTU_C  CPSase_L_chain  COX3  rnh ghf1  EGF_Lam icd cyclo myb_DNA-binding fkbp  grs LDLa  Cys_Met_Meta_PP pnp TPR dpo1  PDZ Toprim  GAF MHC_II_N  PH  HGTP_anticodon  dCMP_cyt_deam transketolase_C hla Mur_ligase_C  chorismate_bind RuBisCO_large_N ccH lacI  PAS GAL4  adk fer2  lyase_1 OTCace  MHC_II_beta_NC  tRNA-synt_2b  rhv Extradiol_dioxy AP_endonuc_2  xia transket_pyr  HMA egf blm Ald_Xan_dh_2  ghf34 AhpC-TSA  Nitroreductase  biotin_lipoyl subt  S4  bac_luciferase  Peptidase_M24 serbact sermam  cox tRNA-synt_1 FAD_binding_4 CPS EFTU_2  gpdh  gluts Ald_Xan_dh_1  NTP_transferase C2  ldh PGAM  prt Rhodanese GATase_2  alpha-amylase oat GATase  reductases  asprs ech aabp  rep FAD-oxidase_NC  thiored cdh alpha-amylase_NC  transketolase_PC  Haloperoxidase  zf-CCHH blmb  actin COX2  aldosered RuBisCO_large Epimerase adh apbact  PA_phosphatase  fn3 Cyclodex_gly_tran hexapep molybdopterin AAA ATP-synt  rrm pdc Phage_integrase Acetyltransf  RuBisCO_large_NC  CBS DEAD  aldedh  p450  aat helicase_C  gtp mdd rvp sdr HATPase_c TyrKc helicase_NC kinase  response_reg  ABC_tranCyclotide Colicin adeno_fiber PSI_PsaE  lta B DNA_PPF hip bowman  leukocidin  Colicin_Pyocin  ghf33 UcrQ  HTH_AraC  TMV_coat  E2_C  scorptoxin  mthina  mthinb  UCR_TM  FN1 cyt3  CBD_3 cks bv  Huristasin  pot planttoxin  Binary_toxA toxin ALBUMIN NTP_transf_2_2  ricin NO_synthase MCR_alpha_N ghf14 ghf6  rnase P Bcl-2 sti cytprime  MIF rnasemam  ins ghf11 EF1BD TNF ghf22 Propep_M14  ngf ets ENTH  profilin  ChtBD Stap_Strp_toxin laminin_G trfl  phoslip mmob  lamb  CBM_20  RuBisCO_small Cys-protease-3C il8 kazal Bacterial_PQQ cystatin  LuxS  pilin BIR DHHA2 TRANSTHYRETIN THAUMATIN intb  tbpc  Cu_bind_like  DISIN SAP pec_lyase tryp_alpha_amyl plantltp  cbm12 slectin PTB kringle rub az  GATA  aldose  ltn DEATH parv  cbp lipocalin fabp  tgfb  recombinase Invasin cryst svmp  ANK ADF fbpase  Fork_head phc ghf10 Bac_rhodopsin S_100 Band_41_M Band_41_N mmp tln YgbB  C1  hemopexin ppase igC1  MHC_II_C  igV igcon igvar-h igvar-l sodcu TIL porin DMRL_synthase Band_41_NME igI cah ATP-gua_Ptrans  Cu_amine_oxid hr  B_lectin  NiFeSe_Hases  RRF Fapy_DNA_glyco  RGS SRF-TF  hexokinase  mofe  oxidored_q6 alpha-amylase_C fibrinogen_C  Ribosomal_L7Ae  ghf7  rnr tms dhfr  acid_phosphat gpr sugbp hpr glob  WW  Pro_CA  ndk chromo  GEL Lipase_3  kunitz  Filamin Ribosomal_L1  KAS glnasn  lipase  MHC_II_alpha_NC Sulfotransfer Adenylsucc_synt Sm  RNase_HII ligase-CoA  hormone_rec phs pgk dsrm  kex RING  cat HMG_box DUF170  annexin sh2 serpin  asp dutpase Glyco_hydro_2 cpa cytb  Lum_binding DHOdehase ptpase  ghf5  bromodomain ACPS  mle TIG Bac_DNA_binding fer4  DHH PEP-utilizers Ribosomal_L6_D  tim EF_TS UBQ LMWPc MCR_alpha sh3 Asp_Glu_race_D  histone NADHdh_2  MM_CoA_mutase Orn_DAP_Arg_deC flav  DAHP_synth_1  tRNA_bind ferritin  CH  cys MCR_alpha_NC  ARM HLH csp Methyltransf_1  citrate_synt  igps  galpha  cat3  int ung FAD-oxidase_C peroxidase  SBP_bacterial_1 cyclin  uce LIM CUB cyt5  cytc  Cu_nir  Glyco_hydro_18  Glyco_hydro_18_D1 Glyco_hydro_18_D2 ghf18 DNA_pol3_beta LRR inositol_P  proteasome  UPF0076 msb PK  MoCF_biosynth sodfe aadh  KISc  ace hom EFTU_C  CPSase_L_chain  COX3  rnh ghf1  EGF_Lam icd cyclo myb_DNA-binding fkbp  grs LDLa  Cys_Met_Meta_PP pnp TPR dpo1  PDZ Toprim  GAF MHC_II_N  PH  HGTP_anticodon  dCMP_cyt_deam transketolase_C hla Mur_ligase_C  chorismate_bind RuBisCO_large_N ccH lacI  PAS GAL4  adk fer2  lyase_1 OTCace  MHC_II_beta_NC  tRNA-synt_2b  rhv Extradiol_dioxy AP_endonuc_2  xia transket_pyr  HMA egf blm Ald_Xan_dh_2  ghf34 AhpC-TSA  Nitroreductase  biotin_lipoyl subt  S4  bac_luciferase  Peptidase_M24 serbact sermam  cox tRNA-synt_1 FAD_binding_4 CPS EFTU_2  gpdh  gluts Ald_Xan_dh_1  NTP_transferase C2  ldh PGAM  prt Rhodanese GATase_2  alpha-amylase oat GATase  reductases  asprs ech aabp  rep FAD-oxidase_NC  thiored cdh alpha-amylase_NC  transketolase_PC  Haloperoxidase  zf-CCHH blmb  actin COX2  aldosered RuBisCO_large Epimerase adh apbact  PA_phosphatase  fn3 Cyclodex_gly_tran hexapep molybdopterin AAA ATP-synt  rrm pdc Phage_integrase Acetyltransf  RuBisCO_large_NC  CBS DEAD  aldedh  p450  aat helicase_C  gtp mdd rvp sdr HATPase_c TyrKc helicase_NC kinase  response_reg  ABC_tran)

################
##  ALIGNERS  ##
################
declare -a aligner=(CLUSTALO MAFFT-FFTNS1 FAMSA)

################
##    TREES   ##
################   
declare -a tree=(CLUSTALO FAMSA MAFFT_PARTTREE)
#             (codnd dpparttreednd1 dpparttreednd2 dpparttreednd2size fastaparttreednd fftns1dnd fftns1dndmem fftns2dnd fftns2dndmem mafftdnd parttreednd0 parttreednd1 parttreednd2 parttreednd2size)

###############
##   Nseq    ##
###############
declare -a bucket=(NA) #(1000 3000 5000 10000 20000)
#               (NA)        -> for PROG
#               (1000 3000 5000)

##############
## Prog/REG ##
##############
declare -a flavour="progressive" #"prog_align"  #"reg_align"  

printf "\t\t########################\n"
printf "\t\t######### TC ###########\n"
printf "\t\t########################\n"

#####
###   ALIG |  CO  |  CO  |  CO  |  CO  |
###   TREE | mBed | mBed |  PT  |  PT  |
###   NSeq | 1000 | 2000 | 1000 | 2000 |
#####

## print 1st line -> ALIGN
printf "ALIG;"
for x in ${aligner[@]} 
do
  for y in ${tree[@]}
  do
    for z in ${bucket[@]}
    do
      printf ${x}";"
    done
  done
done
printf "\n"

#print 2nd line -> TREE
printf "TREE;"
for x in ${aligner[@]}
do
  for y in ${tree[@]}
  do
    for z in ${bucket[@]}
    do
      printf ${y}";"
    done
  done
done
printf "\n"

##print 3r line -> nSeq
printf "nSeq;"
for x in ${aligner[@]}
do
  for y in ${tree[@]}
  do
    for z in ${bucket[@]}
    do
      printf ${z}";"
    done
  done
done
printf "\n\n"

for family in ${all[@]} ## loop for all families
do
  printf "${family};"
	
	for align_method in ${aligner[@]} ## loop each alignment
    do
        for tree_method in ${tree[@]} ## loop each tree
        do
          for nSeq in ${bucket[@]}  ## loop all the buckets
          do
            cat ../results/score/tc/${family}.${flavour}.${nSeq}.${align_method}.with.${tree_method}.tree.tc | tr '' ';'| tr -d "[:space:]"
            printf ";" 
          done
       	done
    done
    printf "\n"
done

printf "\n"
printf "\t\t########################\n"
printf "\t\t## Homoplasy -- EASEL ##\n"
printf "\t\t########################\n"
printf "\n"
## print 1st line -> ALIGN
printf "ALIG;"
for x in ${aligner[@]} 
do
  for y in ${tree[@]}
  do
    for z in ${bucket[@]}
    do
      for i in {1..4}
      do
        printf ${x}";"
      done
    done
  done
done
printf "\n"

#print 2nd line -> TREE
printf "TREE;"
for x in ${aligner[@]}
do
  for y in ${tree[@]}
  do
    for z in ${bucket[@]}
    do
      for i in {1..4}
      do
        printf ${y}";"
      done
    done
  done
done
printf "\n"

##print 3r line -> nSeq
printf "nSeq;"
for x in ${aligner[@]}
do
  for y in ${tree[@]}
  do
    for z in ${bucket[@]}
    do
      for i in {1..4}
      do
        printf ${z}";"
      done
    done
  done
done
printf "\n"
##print 4r line -> homoplasy
printf "nSeq;"
for x in ${aligner[@]}
do
  for y in ${tree[@]}
  do
    for z in ${bucket[@]}
    do
      printf "avgID;alnLen;totGap;numSeq;"
    done
  done
done

printf "\n\n"

for family in ${all[@]} ## loop for all families
do
  printf "${family};"
  
	for align_method in ${aligner[@]} ## loop each alignment
    do
        for tree_method in ${tree[@]} ## loop each tree
        do
          for nSeq in ${bucket[@]}  ## loop all the buckets
          do
            cat ../results/easel/${family}.${flavour}.${nSeq}.${align_method}.with.${tree_method}.tree.avgId | tr '' ';'| tr -d "[:space:]"
            printf ";" 
            cat ../results/gaps/${family}.${flavour}.${nSeq}.${align_method}.with.${tree_method}.tree.alnLen | tr '' ';'| tr -d "[:space:]"
            printf ";" 
            cat ../results/gaps/${family}.${flavour}.${nSeq}.${align_method}.with.${tree_method}.tree.totGap | tr '' ';'| tr -d "[:space:]"
            printf ";" 
            cat ../results/gaps/${family}.${flavour}.${nSeq}.${align_method}.with.${tree_method}.tree.numSeq | tr '' ';'| tr -d "[:space:]"
            printf ";" 
          done
       	done
    done
    printf "\n"
done