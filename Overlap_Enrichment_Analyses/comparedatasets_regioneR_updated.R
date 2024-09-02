#Try to use the new software that Arnar sent
library(dplyr)
library(sqldf)
library(regioneR)  #browseVignettes("regioneR") #For infos

#Set seed for reproducible results
set.seed(123)
#ATH: Always specify mc.set.seed=FALSE when doing permutation, otherwise the function overrides the seed?

############################################# LOAD DATA ######################################################

#Load GRanges objects of the datasets of interest. (made with thescript: C:/Users/sebma/Desktop/Rscripts/makeGRanges_outofdatasets.R)
    #Whole genome assembly (NOT WGseq): Could be used with RegioneR if no resampling.
    genomecharr <- readRDS("/Users/sebma/Desktop/GRanges_Objects/genome_noNW_nomit_GR.rds")

  #Methylation  (ATH: These do not have a peak around them. Just an Irange of 2bp)
    #significant CpGs by morph (Morph, MorphxTime, MorphxSex)
    glm27_signifmorph_noNW_nomit_GR <- readRDS("/Users/sebma/Desktop/GRanges_Objects/glm27_signifmorphall_noNW_nomit_GR.rds")
    #non-significant CpGs:
    glm27_nonsignifmorphall_noNW_nomit_GR <- readRDS("/Users/sebma/Desktop/GRanges_Objects/glm27_nonsignifmorphall_noNW_nomit_GR.rds")
    #whole dataset but without the mitochondrial scaffold
    glm27_allPos_noNW_noMit_GR <- readRDS("/Users/sebma/Desktop/GRanges_Objects/glm27_allPos_noNW_noMit_GR.rds")
    #All signifCpGs by Time (Time, TIme xMorph, TimexSex) for negative control
    glm27_signiftime_noNW_nomit_GR <- readRDS("/Users/sebma/Desktop/GRanges_Objects/glm27_signiftimeall_noNW_nomit_GR.rds")
    
    
  #RADseq data  
    #significant RADseq SNPs (outside 2sigmas). 20kb peak.
    RAD_2sigmas_20kpeak_noNW_nomit_GR <- readRDS("/Users/sebma/Desktop/GRanges_Objects/RAD_2sigmas_20kpeak_noNW_nomit_GR.rds")
    #nonsignif RADseq positions: With 20kpeak (Trimmed)
    RADSNPs_nonsignif_outof2sigma_20kpeak_noNW_nomit_GR_trim <- readRDS("/Users/sebma/Desktop/GRanges_Objects/RADSNPs_nonsignif_outof2sigma_20kpeak_noNW_nomit_GR_trim.rds")
    #All RADseq positions on NC_ with 20kpeak (Trimmed) 
    RADallpos_20kpeak_noNW_nomit_GR_trim <- readRDS("/Users/sebma/Desktop/GRanges_Objects/RADallpos_20kpeak_noNW_nomit_GR_trim.rds")
    
  #WGseq data (100kb peaks)
    #Significant regions outside 2 sigmas (trimmed)
    WGSNPs_2sigmas_NConly_nomit_100k_GR_trim <- readRDS("/Users/sebma/Desktop/GRanges_Objects/WG_2sigma_NConly_nomit_100k_GR_trim.rds") 
    #Non-significant regions (Trimmed)
    WG_nonsignificant_outof2sigma_NConly_nomit_100k_GR_trim <- readRDS("/Users/sebma/Desktop/GRanges_Objects/WG_nonsignificant_outof2sigma_NConly_nomit_100k_GR_trim.rds")
    #All WGseq regions (Trimmed)
    WGpeaks_NConly_nomit_allwindows_GR_trim <- readRDS("/Users/sebma/Desktop/GRanges_Objects/WGpeaks_NConly_nomit_allwindows_GR_trim.rds")
       
  #Transcriptome data (Expression) ATH: Check how I made these.
    #ATH: 146 DE genes are not in the list of genes in the expr_fullref.tsv
    #Need to figure out why
    #With peaks 20kb (10kb upstream of start and 10kb downstream of end)
    #Signif DE
    DEgenesFullInfo_NConly_uniq_20kbpeak_GR <- readRDS("/Users/sebma/Desktop/GRanges_Objects/DEgenesFullInfo_NConly_uniq_20kbpeak_GR.rds")
    #DE genes but without the 146 that are missing from the whole table
    DEgenes_NConly_uniq_20kbpeak_withoutmissinggenes_GR <- readRDS("/Users/sebma/Desktop/GRanges_Objects/DEgenes_NConly_uniq_20kbpeak_withoutmissinggenes_GR.rds")
    
    #All dataset (Trimmed)
    allgenesFullInfo_NConly_uniq_20kbpeak_GR_trim <- readRDS("/Users/sebma/Desktop/GRanges_Objects/allgenesFullInfo_NConly_uniq_20kbpeak_GR_trim.rds")
    #Non signif genes 
    #nonDEgenes_NConly_uniq_20kbpeak_GR <- readRDS("/Users/sebma/Desktop/GRanges_Objects/nonDEgenes_NConly_uniq_20kbpeak_GR_trim.rds")
    

  
################### Overlaps between whole datasets, randomize over whole genome ######################### 
    ############### OLD STUFF
    
    #RRBS vs RAD
    RRBSvsRAD <- overlapPermTest(A=glm27_allPos_noNW_noMit_GR, B=RADallpos_20kpeak_noNW_nomit_GR_trim, ntimes=50, genome=genomecharr, mc.set.seed=FALSE)
    plot(RRBSvsRAD)
      
    #RRBS vs WG (not sure if interesting: WGseq has the whole genome more or less. Also should re-generate the WGseq GR files)
    RRBSvsWG <- overlapPermTest(A=glm27_allPos_noNW_noMit_GR, B=WGpeaks_NConly_nomit_allwindows_GR_trim, ntimes=50, genome=genomecharr, mc.set.seed=FALSE)
    plot(RRBSvsWG)
    
    #RRBS vs Transcriptome genes
    RRBSvsTrans <- overlapPermTest(A=glm27_allPos_noNW_noMit_GR, B=allgenesFullInfo_NConly_uniq_20kbpeak_GR_trim, ntimes=50, genome=genomecharr, mc.set.seed=FALSE)
    plot(RRBSvsTrans)
    
    #RAD vs WG
    RADvsWG <- overlapPermTest(A=RADallpos_20kpeak_noNW_nomit_GR_trim, B=WGpeaks_NConly_nomit_allwindows_GR_trim, ntimes=50, genome=genomecharr, mc.set.seed=FALSE)
    plot(RADvsWG)
    
    #RAD vs Transcriptome genes
    RADvsTrans <- overlapPermTest(A=RADallpos_20kpeak_noNW_nomit_GR_trim, B=allgenesFullInfo_NConly_uniq_20kbpeak_GR_trim, ntimes=50, genome=genomecharr, mc.set.seed=FALSE)
    plot(RADvsTrans)
    
    #WG vs Transcriptome genes
    WGvsTrans <- overlapPermTest(A=WGpeaks_NConly_nomit_allwindows_GR_trim, B=allgenesFullInfo_NConly_uniq_20kbpeak_GR_trim, ntimes=50, genome=genomecharr, mc.set.seed=FALSE)
    plot(WGvsTrans)
    
    
   
##################################### RESAMPLE REGIONS ########################################

    #signif CpGs vs signif RAD
    ptsignifCpG_signifRAD <- permTest(A=glm27_signifmorph_noNW_nomit_GR, ntimes=1000, randomize.function=resampleRegions, universe=glm27_allPos_noNW_noMit_GR,
                                      evaluate.function=numOverlaps, B=RAD_2sigmas_20kpeak_noNW_nomit_GR, verbose=FALSE, mc.set.seed=FALSE)  
    plot(ptsignifCpG_signifRAD)
    
    #Signif CpGs vs nonsignif RAD
    ptsignifCpG_nonsignifRAD <- permTest(A=glm27_signifmorph_noNW_nomit_GR, ntimes=5000, randomize.function=resampleRegions, universe=glm27_allPos_noNW_noMit_GR,
                                         evaluate.function=numOverlaps, B=RADSNPs_nonsignif_outof2sigma_20kpeak_noNW_nomit_GR_trim, verbose=FALSE, mc.set.seed=FALSE)  
    plot(ptsignifCpG_nonsignifRAD)
    
    #Try signif CpGs vs the whole RAD dataset
    ptsignifCpG_allRAD <- permTest(A=glm27_signifmorph_noNW_nomit_GR, ntimes=1000, randomize.function=resampleRegions, universe=glm27_allPos_noNW_noMit_GR,
                                      evaluate.function=numOverlaps, B=RADallpos_20kpeak_noNW_nomit_GR_trim, verbose=FALSE, mc.set.seed=FALSE)  
    plot(ptsignifCpG_allRAD)
    
    #Non sgnif CpGs vs the whole RAD dataset
    ptnonsignifCpG_allRAD <- permTest(A=glm27_nonsignifmorphall_noNW_nomit_GR, ntimes=50, randomize.function=resampleRegions, universe=glm27_allPos_noNW_noMit_GR,
                                   evaluate.function=numOverlaps, B=RADallpos_20kpeak_noNW_nomit_GR_trim, verbose=FALSE, mc.set.seed=FALSE)  
    plot(ptnonsignifCpG_allRAD)
    
    #Just for fun, what happens with signifCpGstime ==> Same as with signifCpGsMorph.
    ptsignifCpGTime_allRAD <- permTest(A=glm27_signiftime_noNW_nomit_GR, ntimes=50, randomize.function=resampleRegions, universe=glm27_allPos_noNW_noMit_GR,
                                   evaluate.function=numOverlaps, B=RADallpos_20kpeak_noNW_nomit_GR_trim, verbose=FALSE, mc.set.seed=FALSE)  
    plot(ptsignifCpGTime_allRAD)
    
    #NEED TO DO SIGNIF RAD vs SIGNIF CPG WITH UNIVERSE = ALL RAD ########################
    ptsignifRAD_signifCpG <- permTest(A=RAD_2sigmas_20kpeak_noNW_nomit_GR, ntimes=1000, randomize.function=resampleRegions, universe=RADallpos_20kpeak_noNW_nomit_GR_trim,
                                      evaluate.function=numOverlaps, B=glm27_signifmorph_noNW_nomit_GR, verbose=FALSE, mc.set.seed=FALSE)  
    plot(ptsignifRAD_signifCpG)
    # ==> Non significant
    
    
##################   
    #Signif CpGs vs WG
    ptsignifCpG_signifWG <- permTest(A=glm27_signifmorph_noNW_nomit_GR, ntimes=50, randomize.function=resampleRegions, universe=glm27_allPos_noNW_noMit_GR,
                                      evaluate.function=numOverlaps, B=WGSNPs_2sigmas_NConly_nomit_100k_GR_trim, verbose=FALSE, mc.set.seed=FALSE)  
    plot(ptsignifCpG_signifWG)
    
    #Signif CpGs vs non-signif WG
    ptsignifCpG_nonsignifWG <- permTest(A=glm27_signifmorph_noNW_nomit_GR, ntimes=1000, randomize.function=resampleRegions, universe=glm27_allPos_noNW_noMit_GR,
                                     evaluate.function=numOverlaps, B=WG_nonsignificant_outof2sigma_NConly_nomit_100k_GR_trim, verbose=FALSE, mc.set.seed=FALSE)  
    plot(ptsignifCpG_nonsignifWG)
    
    #Negative control with Signif Timepoints
    ptsignifCpGTime_signifWG <- permTest(A=glm27_signiftime_noNW_nomit_GR, ntimes=1000, randomize.function=resampleRegions, universe=glm27_allPos_noNW_noMit_GR,
                                     evaluate.function=numOverlaps, B=WGSNPs_2sigmas_NConly_nomit_100k_GR_trim, verbose=FALSE, mc.set.seed=FALSE)  
    plot(ptsignifCpGTime_signifWG)
    
    
##################    
    #Signif CpGs vs DE genes
    ptsignifCpG_DE <- permTest(A=glm27_signifmorph_noNW_nomit_GR, ntimes=50, randomize.function=resampleRegions, universe=glm27_allPos_noNW_noMit_GR,
                                     evaluate.function=numOverlaps, B=DEgenesFullInfo_NConly_uniq_20kbpeak_GR, verbose=FALSE)  
    plot(ptsignifCpG_DE)
    
    #Non-signif CpGs vs DE genes
    ptnonsignifCpG_DE <- permTest(A=glm27_nonsignifmorphall_noNW_nomit_GR, ntimes=50, randomize.function=resampleRegions, universe=glm27_allPos_noNW_noMit_GR,
                               evaluate.function=numOverlaps, B=DEgenesFullInfo_NConly_uniq_20kbpeak_GR, verbose=FALSE)  
    plot(ptnonsignifCpG_DE)
    
    #signif CpGs vs nonDE genes
    ptsignifCpG_nonDE <- permTest(A=glm27_signifmorph_noNW_nomit_GR, ntimes=50, randomize.function=resampleRegions, universe=glm27_allPos_noNW_noMit_GR,
                               evaluate.function=numOverlaps, B=allgenesFullInfo_NConly_uniq_20kbpeak_GR_trim, verbose=FALSE)  
    plot(ptsignifCpG_nonDE)
  
    
##################
    #Signif RAD vs DE genes
    ptsignifRAD_DE <- permTest(A=RAD_2sigmas_20kpeak_noNW_nomit_GR, ntimes=1000, randomize.function=resampleRegions, universe=RADallpos_20kpeak_noNW_nomit_GR_trim,
                                     evaluate.function=numOverlaps, B=DEgenesFullInfo_NConly_uniq_20kbpeak_GR, verbose=FALSE)  
    plot(ptsignifRAD_DE)
    
    #Nonsignif RAD vs DE
    ptNonsignifRAD_DE <- permTest(A=RADSNPs_nonsignif_outof2sigma_20kpeak_noNW_nomit_GR_trim, ntimes=1000, randomize.function=resampleRegions, universe=RADallpos_20kpeak_noNW_nomit_GR_trim,
                               evaluate.function=numOverlaps, B=DEgenesFullInfo_NConly_uniq_20kbpeak_GR, verbose=FALSE)  
    plot(ptNonsignifRAD_DE)
    
    #Basically it is non significant. 
   
    #Try DE vs RAD even though it is not the best
    ptDE_signifRAD <- permTest(A=DEgenesFullInfo_NConly_uniq_20kbpeak_GR, ntimes=1000, randomize.function=resampleRegions, universe=allgenesFullInfo_NConly_uniq_20kbpeak_GR_trim,
                               evaluate.function=numOverlaps, B=RAD_2sigmas_20kpeak_noNW_nomit_GR, verbose=FALSE)  
    plot(ptDE_signifRAD)
   
##################
    #DE genes vs Signif WG
    ptsignifDE_WG <- permTest(A=DEgenesFullInfo_NConly_uniq_20kbpeak_GR, ntimes=1000, randomize.function=resampleRegions, universe=allgenesFullInfo_NConly_uniq_20kbpeak_GR_trim,
                              evaluate.function=numOverlaps, B=WGSNPs_2sigmas_NConly_nomit_100k_GR_trim, verbose=FALSE)  
    plot(ptsignifDE_WG)
    
    #DE genes vs Signif WG but without the 146 genes that are missing from the whole list
    ptsignifDE2_WG <- permTest(A=DEgenes_NConly_uniq_20kbpeak_withoutmissinggenes_GR, ntimes=1000, randomize.function=resampleRegions, universe=allgenesFullInfo_NConly_uniq_20kbpeak_GR_trim,
                              evaluate.function=numOverlaps, B=WGSNPs_2sigmas_NConly_nomit_100k_GR_trim, verbose=FALSE)  
    plot(ptsignifDE2_WG)
    
    #DE genes vs non-signif WG
    ptsignifDE_nonsignifWG <- permTest(A=DEgenesFullInfo_NConly_uniq_20kbpeak_GR, ntimes=50, randomize.function=resampleRegions, universe=allgenesFullInfo_NConly_uniq_20kbpeak_GR_trim,
                              evaluate.function=numOverlaps, B=WG_nonsignificant_outof2sigma_NConly_nomit_100k_GR_trim, verbose=FALSE)  
    plot(ptsignifDE_nonsignifWG)
    
    #DE genes vs whole WG data
    ptsignifDE_WGwhole <- permTest(A=DEgenesFullInfo_NConly_uniq_20kbpeak_GR, ntimes=1000, randomize.function=resampleRegions, universe=allgenesFullInfo_NConly_uniq_20kbpeak_GR_trim,
                              evaluate.function=numOverlaps, B=WGpeaks_NConly_nomit_allwindows_GR_trim, verbose=FALSE)  
    plot(ptsignifDE_WGwhole)
    
    #That is great result. 
    #I should do it in the same way I did the RAD vs DE though. But I can't think straight into which one is the best way to do stuff
    ptsignifWG_DE <- permTest(A=WGSNPs_2sigmas_NConly_nomit_100k_GR_trim, ntimes=1000, randomize.function=resampleRegions, universe=WGpeaks_NConly_nomit_allwindows_GR_trim,
                              evaluate.function=numOverlaps, B=DEgenesFullInfo_NConly_uniq_20kbpeak_GR, verbose=FALSE)  
    plot(ptsignifWG_DE)
    
    #Try WG vs the whole genes dataset
    ptsignifWG_allgenes <- permTest(A=WGSNPs_2sigmas_NConly_nomit_100k_GR_trim, ntimes=1000, randomize.function=resampleRegions, universe=WGpeaks_NConly_nomit_allwindows_GR_trim,
                              evaluate.function=numOverlaps, B=allgenesFullInfo_NConly_uniq_20kbpeak_GR_trim, verbose=FALSE)  
    plot(ptsignifWG_allgenes)
    #Not good
    
    
#########################    
    
    #Signif RAD vs Signif WG
    ptsignifRAD_signifWG_2sig <- permTest(A=RAD_2sigmas_20kpeak_noNW_nomit_GR, ntimes=1000, randomize.function=resampleRegions, universe=RADallpos_20kpeak_noNW_nomit_GR_trim,
                                     evaluate.function=numOverlaps, B=WGSNPs_2sigmas_NConly_nomit_100k_GR_trim, verbose=FALSE)  
    plot(ptsignifRAD_signifWG_2sig)
    #Same thing wooohooo. So I just use these (2sigma data) now instead. 
    
    
    #Signif RAD vs nonsignifWG
    ptsignifRAD_nonsignifWG <- permTest(A=RAD_2sigmas_20kpeak_noNW_nomit_GR, ntimes=50, randomize.function=resampleRegions, universe=RADallpos_20kpeak_noNW_nomit_GR_trim,
                                     evaluate.function=numOverlaps, B=WG_nonsignificant_outof2sigma_NConly_nomit_100k_GR_trim, verbose=FALSE)  
    plot(ptsignifRAD_nonsignifWG)
    
    #THIS IS WHERE I STOPPED 
    "red" "red""red""red""red""red""red""red""red""red""red""red""red""red""red""red""red""red""red""red"

    
##################
    #DE genes and RADSnps
    ptDE_signifRAD <- permTest(A=DEgenesFullInfo_NConly_uniq_20kbpeak_GR, ntimes=1000, randomize.function=resampleRegions, universe=allgenesFullInfo_NConly_uniq_20kbpeak_GR_trim,
                               evaluate.function=numOverlaps, B=RAD_2sigmas_20kpeak_noNW_nomit_GR, verbose=FALSE)  
    plot(ptDE_signifRAD)
    
    #DE genes and RADSnps, but without the missing genes 
    ptDE2_signifRAD <- permTest(A=DEgenes_NConly_uniq_20kbpeak_withoutmissinggenes_GR, ntimes=1000, randomize.function=resampleRegions, universe=allgenesFullInfo_NConly_uniq_20kbpeak_GR_trim,
                               evaluate.function=numOverlaps, B=RAD_2sigmas_20kpeak_noNW_nomit_GR, verbose=FALSE)  
    plot(ptDE2_signifRAD)
    
    #DE genes and non-signif RADsnps
    ptDE_nonsignifRAD <- permTest(A=DEgenesFullInfo_NConly_uniq_20kbpeak_GR, ntimes=50, randomize.function=resampleRegions, universe=allgenesFullInfo_NConly_uniq_20kbpeak_GR_trim,
                               evaluate.function=numOverlaps, B=RADallpos_20kpeak_noNW_GR_trim, verbose=FALSE)  
    plot(ptDE_nonsignifRAD)
    #Weird. I guess I can't take this into account..
    
    #DE genes and the whole RAD dataset
    ptDE_allRAD <- permTest(A=DEgenesFullInfo_NConly_uniq_20kbpeak_GR, ntimes=500, randomize.function=resampleRegions, universe=allgenesFullInfo_NConly_uniq_20kbpeak_GR_trim,
                               evaluate.function=numOverlaps, B=RADallpos_20kpeak_noNW_GR_trim, verbose=FALSE)  
    plot(ptDE_allRAD)
    
    
    
    #DE genes and RADSnps but with 100kbp RAD so I can compare it to WG
    ptDE_signifRAD100kb <- permTest(A=DEgenesFullInfo_NConly_uniq_20kbpeak_GR, ntimes=500, randomize.function=resampleRegions, universe=allgenesFullInfo_NConly_uniq_20kbpeak_GR_trim,
                               evaluate.function=numOverlaps, B=RAD0p2_100kpeak_GR_trim, verbose=FALSE)  
    plot(ptDE_signifRAD100kb)
    
    #DE genes and nonsignif RAD (or the whole RAD dataset) but with 100kbp
    ptDE_nonsignifRAD100kb <- permTest(A=DEgenesFullInfo_NConly_uniq_20kbpeak_GR, ntimes=500, randomize.function=resampleRegions, universe=allgenesFullInfo_NConly_uniq_20kbpeak_GR_trim,
                                    evaluate.function=numOverlaps, B=RADpeaks100kb_nonsignif_GR_trim, verbose=FALSE)  
    plot(ptDE_nonsignifRAD100kb)
    
    #DE genes and the whole RAD dataset with 100kb
    
    
#################
    #DE genes and WGpeaks
    ptDE_signifWG <- permTest(A=DEgenesFullInfo_NConly_uniq_20kbpeak_GR, ntimes=50, randomize.function=resampleRegions, universe=allgenesFullInfo_NConly_uniq_20kbpeak_GR_trim,
                               evaluate.function=numOverlaps, B=WGpeaks_NConly_Signif0p2_GR_trim, verbose=FALSE)  
    plot(ptDE_signifWG)
    
    #DE genes and nonsignif WGpeaks
    ptDE_nonsignifWG <- permTest(A=DEgenesFullInfo_NConly_uniq_20kbpeak_GR, ntimes=50, randomize.function=resampleRegions, universe=allgenesFullInfo_NConly_uniq_20kbpeak_GR_trim,
                              evaluate.function=numOverlaps, B=WGpeaks_NConly_nonsignif_GR_trim, verbose=FALSE)  
    plot(ptDE_nonsignifWG)
    
    #WGpeaks and DE genes
    ptsignifWG_DE <- permTest(A=WGpeaks_NConly_Signif0p2_GR_trim, ntimes=50, randomize.function=resampleRegions, universe=WGpeaks_NConly_allwindows_GR_trim,
                               evaluate.function=numOverlaps, B=DEgenesFullInfo_NConly_uniq_20kbpeak_GR, verbose=FALSE)  
    plot(ptsignifWG_DE)
    
    #DE genes and non signif WGpeaks
    ptDE_nonsignifWG <- permTest(A=DEgenesFullInfo_NConly_uniq_20kbpeak_GR, ntimes=50, randomize.function=resampleRegions, universe=allgenesFullInfo_NConly_uniq_20kbpeak_GR_trim,
                              evaluate.function=numOverlaps, B=WGpeaks_NConly_nonsignif_GR_trim, verbose=FALSE)  
    plot(ptDE_nonsignifWG)
    
    #DE genes and the whole WG dataset?
    ptDE_AllWG <- permTest(A=DEgenesFullInfo_NConly_uniq_20kbpeak_GR, ntimes=50, randomize.function=resampleRegions, universe=allgenesFullInfo_NConly_uniq_20kbpeak_GR_trim,
                              evaluate.function=numOverlaps, B=WGpeaks_NConly_allwindows_GR_trim, verbose=FALSE)  
    plot(ptDE_AllWG)
        

######### Evaluate whether there is congruence in terms of morph specific peaks in significant datasets    

#Load morph specific data:
    #RAD data (ATH OLD DATA, BIAS BECAUSE OF SNPs WITH SAME FST)
    # ddRAD_SNPs_LBonly_2sigmas_NConly_nomit_20k_GR <- readRDS("/Users/sebma/Desktop/GRanges_Objects/ddRAD_SNPs_LBonly_2sigmas_NConly_nomit_20k_GR.rds")
    # ddRAD_SNPs_SBonly_2sigmas_NConly_nomit_20k_GR <- readRDS("/Users/sebma/Desktop/GRanges_Objects/ddRAD_SNPs_SBonly_2sigmas_NConly_nomit_20k_GR.rds")
    # ddRAD_SNPs_PLonly_2sigmas_NConly_nomit_20k_GR <- readRDS("/Users/sebma/Desktop/GRanges_Objects/ddRAD_SNPs_PLonly_2sigmas_NConly_nomit_20k_GR.rds")
    #NEW RAD 
    ddRAD_SNPs_LBonly_2sigmas_NConly_nomit_20k_GR <- readRDS("/Users/sebma/Desktop/GRanges_Objects/ddRAD_SNPs_LBonly_NOBIAS_2sigmas_NConly_nomit_20k_GR.rds")
    ddRAD_SNPs_SBonly_2sigmas_NConly_nomit_20k_GR <- readRDS("/Users/sebma/Desktop/GRanges_Objects/ddRAD_SNPs_SBonly_NOBIAS_2sigmas_NConly_nomit_20k_GR.rds")
    ddRAD_SNPs_PLonly_2sigmas_NConly_nomit_20k_GR <- readRDS("/Users/sebma/Desktop/GRanges_Objects/ddRAD_SNPs_PLonly_NOBIAS_2sigmas_NConly_nomit_20k_GR.rds")
    
    #WG data
    WGSNPs_LBonly_2sigmas_NConly_nomit_100k_GR_trim <- readRDS("/Users/sebma/Desktop/GRanges_Objects/WGSNPs_LBonly_2sigmas_NConly_nomit_100k_GR_trim.rds")
    WGSNPs_SBonly_2sigmas_NConly_nomit_100k_GR_trim <- readRDS("/Users/sebma/Desktop/GRanges_Objects/WGSNPs_SBonly_2sigmas_NConly_nomit_100k_GR_trim.rds")
    WGSNPs_PLonly_2sigmas_NConly_nomit_100k_GR_trim <- readRDS("/Users/sebma/Desktop/GRanges_Objects/WGSNPs_PLonly_2sigmas_NConly_nomit_100k_GR_trim.rds")
    
    #DE data
    DEgenesFullInfo_NConly_uniq_LBonly_20k_GR <- readRDS("/Users/sebma/Desktop/GRanges_Objects/DEgenesFullInfo_NConly_uniq_LBonly_20k_GR.rds")
    DEgenesFullInfo_NConly_uniq_SBonly_20k_GR <- readRDS("/Users/sebma/Desktop/GRanges_Objects/DEgenesFullInfo_NConly_uniq_SBonly_20k_GR.rds")
    DEgenesFullInfo_NConly_uniq_PLonly_20k_GR <- readRDS("/Users/sebma/Desktop/GRanges_Objects/DEgenesFullInfo_NConly_uniq_PLonly_20k_GR.rds")

    DEwithoutmissinggenes_LBonly_20k_GR <- readRDS("/Users/sebma/Desktop/GRanges_Objects/DEgenesFullInfo_NConly_uniq_LBonly_nomissinggene_20k_GR.rds")
    DEwithoutmissinggenes_SBonly_20k_GR <- readRDS("/Users/sebma/Desktop/GRanges_Objects/DEgenesFullInfo_NConly_uniq_SBonly_nomissinggene_20k_GR.rds")
    DEwithoutmissinggenes_PLonly_20k_GR <- readRDS("/Users/sebma/Desktop/GRanges_Objects/DEgenesFullInfo_NConly_uniq_PLonly_nomissinggene_20k_GR.rds")
    
    
#Load significant datasets for universe
    #significant RADseq SNPs (outside 2sigmas). 20kb peak.
    RAD_2sigmas_20kpeak_noNW_nomit_GR <- readRDS("/Users/sebma/Desktop/GRanges_Objects/RAD_2sigmas_20kpeak_noNW_nomit_GR.rds")
    #All DE genes. 20kb peak
    DEgenesFullInfo_NConly_uniq_20kbpeak_GR <- readRDS("/Users/sebma/Desktop/GRanges_Objects/DEgenesFullInfo_NConly_uniq_20kbpeak_GR.rds")
    #All DE genes but without the 146 that are missing from the whole table
    DEgenes_NConly_uniq_20kbpeak_withoutmissinggenes_GR <- readRDS("/Users/sebma/Desktop/GRanges_Objects/DEgenes_NConly_uniq_20kbpeak_withoutmissinggenes_GR.rds")
    
    
#Perform comparison analyses
    RADLB_WGLB <- permTest(A=ddRAD_SNPs_LBonly_2sigmas_NConly_nomit_20k_GR, ntimes=1000, randomize.function=resampleRegions, universe=RAD_2sigmas_20kpeak_noNW_nomit_GR,
                                 evaluate.function=numOverlaps, B=WGSNPs_LBonly_2sigmas_NConly_nomit_100k_GR_trim, verbose=FALSE)  
    
    RADSB_WGSB <- permTest(A=ddRAD_SNPs_SBonly_2sigmas_NConly_nomit_20k_GR, ntimes=1000, randomize.function=resampleRegions, universe=RAD_2sigmas_20kpeak_noNW_nomit_GR,
                           evaluate.function=numOverlaps, B=WGSNPs_SBonly_2sigmas_NConly_nomit_100k_GR_trim, verbose=FALSE)  
    
    RADPL_WGPL <- permTest(A=ddRAD_SNPs_PLonly_2sigmas_NConly_nomit_20k_GR, ntimes=1000, randomize.function=resampleRegions, universe=RAD_2sigmas_20kpeak_noNW_nomit_GR,
                           evaluate.function=numOverlaps, B=WGSNPs_PLonly_2sigmas_NConly_nomit_100k_GR_trim, verbose=FALSE)  
    
    
    plot(RADLB_WGLB)
    plot(RADSB_WGSB)
    plot(RADPL_WGPL)
    
    DELB_WGLB <- permTest(A=DEgenesFullInfo_NConly_uniq_LBonly_20k_GR, ntimes=1000, randomize.function=resampleRegions, universe=DEgenesFullInfo_NConly_uniq_20kbpeak_GR,
                          evaluate.function=numOverlaps, B=WGSNPs_LBonly_2sigmas_NConly_nomit_100k_GR_trim, verbose=FALSE)  
    
    DELB_WGLB2 <- permTest(A=DEwithoutmissinggenes_LBonly_20k_GR, ntimes=1000, randomize.function=resampleRegions, universe=DEgenes_NConly_uniq_20kbpeak_withoutmissinggenes_GR,
                          evaluate.function=numOverlaps, B=WGSNPs_LBonly_2sigmas_NConly_nomit_100k_GR_trim, verbose=FALSE)  
    
    DESB_WGSB <- permTest(A=DEgenesFullInfo_NConly_uniq_SBonly_20k_GR, ntimes=1000, randomize.function=resampleRegions, universe=DEgenesFullInfo_NConly_uniq_20kbpeak_GR,
                          evaluate.function=numOverlaps, B=WGSNPs_SBonly_2sigmas_NConly_nomit_100k_GR_trim, verbose=FALSE)  
    
    DESB_WGSB2 <- permTest(A=DEwithoutmissinggenes_SBonly_20k_GR, ntimes=1000, randomize.function=resampleRegions, universe=DEgenes_NConly_uniq_20kbpeak_withoutmissinggenes_GR,
                          evaluate.function=numOverlaps, B=WGSNPs_SBonly_2sigmas_NConly_nomit_100k_GR_trim, verbose=FALSE) 
    
    DEPL_WGPL <- permTest(A=DEgenesFullInfo_NConly_uniq_PLonly_20k_GR, ntimes=1000, randomize.function=resampleRegions, universe=DEgenesFullInfo_NConly_uniq_20kbpeak_GR,
                          evaluate.function=numOverlaps, B=WGSNPs_PLonly_2sigmas_NConly_nomit_100k_GR_trim, verbose=FALSE)  
    
    DEPL_WGPL2 <- permTest(A=DEwithoutmissinggenes_PLonly_20k_GR, ntimes=1000, randomize.function=resampleRegions, universe=DEgenes_NConly_uniq_20kbpeak_withoutmissinggenes_GR,
                          evaluate.function=numOverlaps, B=WGSNPs_PLonly_2sigmas_NConly_nomit_100k_GR_trim, verbose=FALSE)
    
    plot(DELB_WGLB)
    plot(DELB_WGLB2)
    plot(DESB_WGSB)
    plot(DESB_WGSB2)
    plot(DEPL_WGPL)
    plot(DEPL_WGPL2)
    

############ Calculate NumOverlaps between the different datasets
#LOAD DATA
#whole RRBS dataset but without the mitochondrial scaffold
glm27_allPos_noNW_noMit_GR <- readRDS("/Users/sebma/Desktop/GRanges_Objects/glm27_allPos_noNW_noMit_GR.rds")
#All RADseq positions on NC_ with 20kpeak (Trimmed) 
RADallpos_20kpeak_noNW_nomit_GR_trim <- readRDS("/Users/sebma/Desktop/GRanges_Objects/RADallpos_20kpeak_noNW_nomit_GR_trim.rds")
#All transcriptome genes on NC_ (I think there is no mitochondrial data in transcriptome)
allgenesFullInfo_NConly_uniq_20kbpeak_GR_trim <- readRDS("/Users/sebma/Desktop/GRanges_Objects/allgenesFullInfo_NConly_uniq_20kbpeak_GR_trim.rds")
#WGseq data, just to have an idea. ATH: I can't have an idea because there are reps..
WGpeaks_NConly_allwindows_GR_trim <- readRDS("/Users/sebma/Desktop/GRanges_Objects/WGpeaks_NConly_allwindows_GR_trim.rds")

numOverlaps(glm27_allPos_noNW_noMit_GR, RADallpos_20kpeak_noNW_nomit_GR_trim) #2210 overlaps, 
# out of 8000 that's actually still pretty good.. So I should find the same results.
numOverlaps(RADallpos_20kpeak_noNW_nomit_GR_trim, allgenesFullInfo_NConly_uniq_20kbpeak_GR_trim) #8540 overlaps
numOverlaps(glm27_allPos_noNW_noMit_GR,allgenesFullInfo_NConly_uniq_20kbpeak_GR_trim)
numOverlaps(glm27_allPos_noNW_noMit_GR,WGpeaks_NConly_allwindows_GR_trim)

#I feel like this should already be enough to judge whether the significant ones are enriched.
#Check again what it says with the significant ones.
#signif WG
WGSNPs_2sigmas_NConly_nomit_100k_GR_trim <- readRDS("/Users/sebma/Desktop/GRanges_Objects/WG_2sigma_NConly_nomit_GR_trim.rds") 
#Signif RRBS.
glm27_signifmorphall_noNW_nomit_GR <- readRDS("/Users/sebma/Desktop/GRanges_Objects/glm27_signifmorphall_noNW_nomit_GR.rds") 
#signif RAD
RAD_2sigmas_20kpeak_noNW_nomit_GR <- readRDS("/Users/sebma/Desktop/GRanges_Objects/RAD_2sigmas_20kpeak_noNW_nomit_GR.rds")
#signif Transcriptome
DEgenesFullInfo_NConly_uniq_20kbpeak_GR <- readRDS("/Users/sebma/Desktop/GRanges_Objects/DEgenesFullInfo_NConly_uniq_20kbpeak_GR.rds")


#Signif CpGs vs Signif WG
ptsignifCpG_signifWG <- permTest(A=glm27_signifmorphall_noNW_nomit_GR, ntimes=50, randomize.function=resampleRegions, universe=glm27_allPos_noNW_noMit_GR,
                                 evaluate.function=numOverlaps, B=WGSNPs_2sigmas_NConly_nomit_100k_GR_trim, verbose=FALSE)  
plot(ptsignifCpG_signifWG)

#Signif CpGs vs Signif RADseq 
ptsignifCpG_signifRAD <- permTest(A=glm27_signifmorphall_noNW_nomit_GR, ntimes=50, randomize.function=resampleRegions, universe=glm27_allPos_noNW_noMit_GR,
                                 evaluate.function=numOverlaps, B=RAD_2sigmas_20kpeak_noNW_nomit_GR, verbose=FALSE)  
plot(ptsignifCpG_signifRAD)

#Nonsignif CpGs vs Signif RADseq
ptnonsignifCpG_signifRAD <- permTest(A=glm27_nonsignifmorphall_noNW_nomit_GR, ntimes=50, randomize.function=resampleRegions, universe=glm27_allPos_noNW_noMit_GR,
                                  evaluate.function=numOverlaps, B=RAD_2sigmas_20kpeak_noNW_nomit_GR, verbose=FALSE)  
plot(ptnonsignifCpG_signifRAD)

#Need to check the opposite: Signif CpGs vs non signif RAD












