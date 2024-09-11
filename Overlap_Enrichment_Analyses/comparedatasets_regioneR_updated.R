library(dplyr)
library(sqldf)
library(regioneR)  #browseVignettes("regioneR") #For infos

#Note: Setting seed for reproducible results does not seem to work with RegioneR resampling of regions? 
#It is not a big deal when using as big as 1000 permutations (even 50 permutations already gives consistent results).
#But because of this, while the result of the analyses will be the same, do not be surprised if the histogram bars on the graphs change 
#slightly from one run of the script to the other.

############################################# LOAD DATA ######################################################
#Load GRanges objects of the datasets of interest located in Superimposition-Paper/GRanges_Objects/ 
#These were made with thescript: Superimposition-Paper/GRanges_Objects/makeGRanges_outofdatasets.R

#Notes:
#"no_NW" or "NC_only" in file name indicates that the dataset does not comprise data from unplaced scaffolds
#"nomit" or "noMit" in file name indicates that the dataset does not comprise data from the mitochondrial scaffold
#"Trimmed" or "_trim" in file name indicates that the GRanges object were trimmed to remove sequences outside the boundaries of a scaffold.
    #For instance, if a SNP is located 500bp before the end of a scaffold, then if I set-up a peak that goes from 10kb 
    #downstream to 10kb upstream of this base, the GRanges object will actually exceed the end of the scaffold by 9500bp.
    #Because I did not know whether this could confuse the permutation software or not, I decided to trim all sequences outside the scaffolds.
    #Thus, GRanges objects do not all have exactly the same size. (Some are more equal than others #AnimalFarm).
#Because we will be using the resampleRegions randomizing function, we do not actually need the GRanges objects coresponding to the non-significant parts of the data:
    #Indeed, you will compare positions of the significant parts of the data with resampling done on the ENTIRETY of the data (comprising significant 
    #and non-significant regions). Still, I have included to this script and to Superimposition-Paper/GRanges_Objects/ GRanges objects corresponding
    #to the non-significant parts of the data. Feel free to use them as negative controls or if you want to test other things.

  #Methylation  (These do not have a peak around them. Just an Irange of 2bp)
    #significant CpGs by morph (Morph, MorphxTime, MorphxSex)
    glm27_signifmorph_noNW_nomit_GR <- readRDS("~/glm27_signifmorphall_noNW_nomit_GR.rds")
    #non-significant CpGs:
    glm27_nonsignifmorphall_noNW_nomit_GR <- readRDS("~/glm27_nonsignifmorphall_noNW_nomit_GR.rds")
    #whole RRBS dataset (all CpGs)
    glm27_allPos_noNW_noMit_GR <- readRDS("~/glm27_allPos_noNW_noMit_GR.rds")
    #All signifCpGs by Time (Time, Time xMorph, TimexSex) for negative control
    glm27_signiftime_noNW_nomit_GR <- readRDS("~/glm27_signiftimeall_noNW_nomit_GR.rds")

                #MORPH SPECIFIC methylation data
                #Does not exist (not possible to establish congruence)   
    
  #RADseq data (20kb peak: from 10kb downstream to 10kb upstream the SNP)
    #significant RADseq SNPs (outside 2sigmas)
    RAD_2sigmas_20kpeak_noNW_nomit_GR <- readRDS("~/RAD_2sigmas_20kpeak_noNW_nomit_GR.rds")
    #nonsignif RADseq positions (Trimmed)
    RADSNPs_nonsignif_outof2sigma_20kpeak_noNW_nomit_GR_trim <- readRDS("~/RADSNPs_nonsignif_outof2sigma_20kpeak_noNW_nomit_GR_trim.rds")
    #All RADseq positions on NC_ with 20kpeak (Trimmed) 
    RADallpos_20kpeak_noNW_nomit_GR_trim <- readRDS("~/RADallpos_20kpeak_noNW_nomit_GR_trim.rds")

                #MORPH SPECIFIC ddRAD significant data
                ddRAD_SNPs_LBonly_2sigmas_NConly_nomit_20k_GR <- readRDS("~/ddRAD_SNPs_LBonly_NOBIAS_2sigmas_NConly_nomit_20k_GR.rds")
                ddRAD_SNPs_SBonly_2sigmas_NConly_nomit_20k_GR <- readRDS("~/ddRAD_SNPs_SBonly_NOBIAS_2sigmas_NConly_nomit_20k_GR.rds")
                ddRAD_SNPs_PLonly_2sigmas_NConly_nomit_20k_GR <- readRDS("~/ddRAD_SNPs_PLonly_NOBIAS_2sigmas_NConly_nomit_20k_GR.rds")

  #WGseq data (100kb peaks, as identified by Alexander)
    #Significant regions outside 2 sigmas (trimmed)
    WGSNPs_2sigmas_NConly_nomit_100k_GR_trim <- readRDS("~/WG_2sigma_NConly_nomit_100k_GR_trim.rds") 
    #Non-significant regions (Trimmed)
    WG_nonsignificant_outof2sigma_NConly_nomit_100k_GR_trim <- readRDS("~/WG_nonsignificant_outof2sigma_NConly_nomit_100k_GR_trim.rds")
    #All WGseq regions (Trimmed)
    WGpeaks_NConly_nomit_allwindows_GR_trim <- readRDS("~/WGpeaks_NConly_nomit_allwindows_GR_trim.rds")

                #MORPH SPECIFIC WGS significant data
                WGSNPs_LBonly_2sigmas_NConly_nomit_100k_GR_trim <- readRDS("~/WGSNPs_LBonly_2sigmas_NConly_nomit_100k_GR_trim.rds")
                WGSNPs_SBonly_2sigmas_NConly_nomit_100k_GR_trim <- readRDS("~/WGSNPs_SBonly_2sigmas_NConly_nomit_100k_GR_trim.rds")
                WGSNPs_PLonly_2sigmas_NConly_nomit_100k_GR_trim <- readRDS("~/WGSNPs_PLonly_2sigmas_NConly_nomit_100k_GR_trim.rds")

  #Transcriptome data (20kb peaks)
    #List of all genes in the transcriptome, made from the realigned transcriptome data to the charr genome by Alexander
    allgenesFullInfo_NConly_uniq_20kbpeak_GR_trim <- readRDS("~/allgenesFullInfo_NConly_uniq_20kbpeak_GR_trim.rds")
    #List of the DE genes, as identified by Alexander
        #DEgenesFullInfo_NConly_uniq_20kbpeak_GR <- readRDS("~/DEgenesFullInfo_NConly_uniq_20kbpeak_GR.rds")
        #HOWEVER, 146 of these DE genes are not in the full gene list. Why? I do not know. Did Alexander use different raw files to generate both?
        #Is there an issue with the file I used to get the ID of the genes? (i.e., is the file I used incomplete?)
        #Because of this issue, I have also generated a list of DE genes without the 146 genes which are not in the whole gene list
    #List of the DE genes, as identified by Alexander, minus the 146 weird genes
    DEgenes_NConly_uniq_20kbpeak_withoutmissinggenes_GR <- readRDS("~/DEgenes_NConly_uniq_20kbpeak_withoutmissinggenes_GR.rds")
    #I personnally prefer using this file. It removes a bit of data, but it is better for comparisons between datasets in my opinion.
    #Feel free to use the file above if you prefer. It might be completely valid too. I just don't like that I can't figure out why there is a difference.

                #MORPH SPECIFIC transcriptome data
                #Again, two options, I prefer the datasets without the missing genes.
                #But I leave these too in case you are curious.
                #DEgenesFullInfo_NConly_uniq_LBonly_20k_GR <- readRDS("~/DEgenesFullInfo_NConly_uniq_LBonly_20k_GR.rds")
                #DEgenesFullInfo_NConly_uniq_SBonly_20k_GR <- readRDS("~/DEgenesFullInfo_NConly_uniq_SBonly_20k_GR.rds")
                #DEgenesFullInfo_NConly_uniq_PLonly_20k_GR <- readRDS("~/DEgenesFullInfo_NConly_uniq_PLonly_20k_GR.rds")
                
                DEwithoutmissinggenes_LBonly_20k_GR <- readRDS("~/DEgenesFullInfo_NConly_uniq_LBonly_nomissinggene_20k_GR.rds")
                DEwithoutmissinggenes_SBonly_20k_GR <- readRDS("~/DEgenesFullInfo_NConly_uniq_SBonly_nomissinggene_20k_GR.rds")
                DEwithoutmissinggenes_PLonly_20k_GR <- readRDS("~/DEgenesFullInfo_NConly_uniq_PLonly_nomissinggene_20k_GR.rds")


#BONUS
#Whole genome assembly (NOT WGseq, the whole assembly). This can be used with RegioneR if you don't want to use the resampleRegions function.
#You can see an example of this at the end of this script. I do not think it should be used, as the resampleRegions function is,
#to me, the most accurate. But anyway, if you want to use the assembly as a base for other functions, here it is as a GRanges object.
   #genomecharr <- readRDS("~/genome_noNW_nomit_GR.rds")


####################################### OVERLAP ANALYSES USING THE RESAMPLEREGIONS FUNCTION ######################################################
#To perform these analyses, you will always need three datasets:
#Two datasets that you want to compare the overlap of, for instance we could call them "SIGNIF_EGG" and "SIGNIF_POTATOES"
#A third dataset called universe with which the software will perform permutation analyses (the number of permutations is signified by "ntimes=")
#In our example, we could call this file "UNIVERSE_EGG". It contains positions for all the eggs in our egg basket, and not just the significant eggs.
#The function then works like this 
#result <- permTest(A=SIGNIF_EGG, ntimes=1000, randomize.function=resampleRegions, universe=UNIVERSE_EGG,
#                   evaluate.function=numOverlaps, B=SIGNIF_POTATOES, verbose=FALSE, mc.set.seed=FALSE)
#It will count the number of overlaps between "SIGNIF_EGG" and "SIGNIF_POTATOES". 
#It will create 1000 subsets of "UNIVERSE_EGG" which will all have the same size as "SIGNIF_EGG".
#It will count the number of overlaps between these subsets and "SIGNIF POTATOES".
#At the end, number of overlaps between the subsets of universe and "SIGNIF_POTATOES" will be compared to the number of overlaps between
#"SIGNIF_EGG" and "SIGNIF_POTATOES".
#Results can be directly visualized with plot(result)


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

#Load morph specific data: ATH: I should move this UP.
   
    
    
    
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
    


#Are these not above?
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

######################################################
################### Overlaps between whole datasets, randomize over whole genome ######################### 
    ############### OLD STUFF, Leaving here just for example
    
    #RRBS vs RAD
    RRBSvsRAD <- overlapPermTest(A=glm27_allPos_noNW_noMit_GR, B=RADallpos_20kpeak_noNW_nomit_GR_trim, ntimes=50, genome=genomecharr, mc.set.seed=FALSE)
    plot(RRBSvsRAD)

############ Calculate NumOverlaps between the different datasets
numOverlaps(glm27_allPos_noNW_noMit_GR, RADallpos_20kpeak_noNW_nomit_GR_trim) #2210 overlaps, 







