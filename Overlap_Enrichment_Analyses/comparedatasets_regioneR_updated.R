library(dplyr)
library(sqldf)
library(regioneR)

#Note: Setting seed for reproducible results does not seem to work with RegioneR resampling of regions? 
#It is not a big deal when using as big as 1000 permutations (even 50 permutations already gives consistent results).
#But because of this, while the result of the analyses will be the same, do not be surprised if the histogram bars on the graphs change 
#slightly from one run of the script to the other.

############################################# I. LOAD DATA ######################################################
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
                #Note: Sometimes, two morphs have the same Fst value for the same SNP. In this case, the SNP appears in both files. 
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


####################################### II. OVERLAP ANALYSES USING THE RESAMPLEREGIONS FUNCTION ######################################################
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

#To simplify and clarify this script, I have only kept here the analyses leading to a Figure or Supplementary Figure in the paper.
#But feel free to try other comparisons yourself, try other controls, switch things around, etc...

#Figure 4 (Tests of enrichment between DE genes and regions of genetic differentiation (WGS), within each of the morphs)
    #A) LB
        DELB_WGLB <- permTest(A=DEwithoutmissinggenes_LBonly_20k_GR, ntimes=1000, randomize.function=resampleRegions, universe=DEgenes_NConly_uniq_20kbpeak_withoutmissinggenes_GR,
                          evaluate.function=numOverlaps, B=WGSNPs_LBonly_2sigmas_NConly_nomit_100k_GR_trim, verbose=FALSE)
        plot(DELB_WGLB)

    #B) SB
        DESB_WGSB <- permTest(A=DEwithoutmissinggenes_SBonly_20k_GR, ntimes=1000, randomize.function=resampleRegions, universe=DEgenes_NConly_uniq_20kbpeak_withoutmissinggenes_GR,
                          evaluate.function=numOverlaps, B=WGSNPs_SBonly_2sigmas_NConly_nomit_100k_GR_trim, verbose=FALSE)
        plot(DESB_WGSB)

    #C) PL
        DEPL_WGPL <- permTest(A=DEwithoutmissinggenes_PLonly_20k_GR, ntimes=1000, randomize.function=resampleRegions, universe=DEgenes_NConly_uniq_20kbpeak_withoutmissinggenes_GR,
                          evaluate.function=numOverlaps, B=WGSNPs_PLonly_2sigmas_NConly_nomit_100k_GR_trim, verbose=FALSE)
        plot(DEPL_WGPL)

#Supplementary Figure 6 (overlap and morph congruence WGS VS RAD):
    #A) Regions of genetic differentiation between morphs identified through ddRAD-seq and WGS overlap more than expected by chance.
        ptsignifRAD_signifWG_2sig <- permTest(A=RAD_2sigmas_20kpeak_noNW_nomit_GR, ntimes=1000, randomize.function=resampleRegions, universe=RADallpos_20kpeak_noNW_nomit_GR_trim,
                                         evaluate.function=numOverlaps, B=WGSNPs_2sigmas_NConly_nomit_100k_GR_trim, verbose=FALSE)  
        plot(ptsignifRAD_signifWG_2sig)

    #B-D) Moreover, congruence was observed between regions of genetic differentiation identified by WGS and SNPs from ddRADseq for all three sympatric morphs
        RADLB_WGLB <- permTest(A=ddRAD_SNPs_LBonly_2sigmas_NConly_nomit_20k_GR, ntimes=1000, randomize.function=resampleRegions, universe=RAD_2sigmas_20kpeak_noNW_nomit_GR,
                                     evaluate.function=numOverlaps, B=WGSNPs_LBonly_2sigmas_NConly_nomit_100k_GR_trim, verbose=FALSE)  
        plot(RADLB_WGLB)
        
        RADSB_WGSB <- permTest(A=ddRAD_SNPs_SBonly_2sigmas_NConly_nomit_20k_GR, ntimes=1000, randomize.function=resampleRegions, universe=RAD_2sigmas_20kpeak_noNW_nomit_GR,
                               evaluate.function=numOverlaps, B=WGSNPs_SBonly_2sigmas_NConly_nomit_100k_GR_trim, verbose=FALSE)  
        plot(RADSB_WGSB)
    
        RADPL_WGPL <- permTest(A=ddRAD_SNPs_PLonly_2sigmas_NConly_nomit_20k_GR, ntimes=1000, randomize.function=resampleRegions, universe=RAD_2sigmas_20kpeak_noNW_nomit_GR,
                               evaluate.function=numOverlaps, B=WGSNPs_PLonly_2sigmas_NConly_nomit_100k_GR_trim, verbose=FALSE)
        plot(RADPL_WGPL)

#Supplementary Figure 7
    #A) Regions of genetic differentiation (WGS) and DE genes overlap more than expected by chance. (resampling on all WG regions) 
        ptsignifWG_DE <- permTest(A=WGSNPs_2sigmas_NConly_nomit_100k_GR_trim, ntimes=1000, randomize.function=resampleRegions, universe=WGpeaks_NConly_nomit_allwindows_GR_trim,
                              evaluate.function=numOverlaps, B=DEgenes_NConly_uniq_20kbpeak_withoutmissinggenes_GR, verbose=FALSE)  
        plot(ptsignifWG_DE)

                ATH: signif WGS VS all genes also gives enrichment (but lower Z-score):
                ptsignifWG_allgenes <- permTest(A=WGSNPs_2sigmas_NConly_nomit_100k_GR_trim, ntimes=1000, randomize.function=resampleRegions, universe=WGpeaks_NConly_nomit_allwindows_GR_trim,
                              evaluate.function=numOverlaps, B=allgenesFullInfo_NConly_uniq_20kbpeak_GR_trim, verbose=FALSE)  
                plot(ptsignifWG_allgenes)
                

    #B) Regions of genetic differentiation (ddRAD) and DE genes tend towards overlap enrichment.
        ptsignifRAD_DE <- permTest(A=RAD_2sigmas_20kpeak_noNW_nomit_GR, ntimes=1000, randomize.function=resampleRegions, universe=RADallpos_20kpeak_noNW_nomit_GR_trim,
                                     evaluate.function=numOverlaps, B=DEgenes_NConly_uniq_20kbpeak_withoutmissinggenes_GR, verbose=FALSE)  
        plot(ptsignifRAD_DE)

#Supplementary Figure 8: 
#Tests for overlap of CpGs differently methylated by morphs (tested with linear model in R, Matlosz et al. 2022) and 
#regions of genetic differentiation (from WGS data).
        ptsignifCpG_signifWG <- permTest(A=glm27_signifmorph_noNW_nomit_GR, ntimes=1000, randomize.function=resampleRegions, universe=glm27_allPos_noNW_noMit_GR,
                                      evaluate.function=numOverlaps, B=WGSNPs_2sigmas_NConly_nomit_100k_GR_trim, verbose=FALSE, mc.set.seed=FALSE)  
        plot(ptsignifCpG_signifWG)

#Supplementary Figure XX10: 
#CpGs differently methylated by timepoint (“Time”, “Time X Morph” and “Time X Sex”) are not enriched in regions of genetic differentiation (WGS).
        ptsignifCpGTime_signifWG <- permTest(A=glm27_signiftime_noNW_nomit_GR, ntimes=1000, randomize.function=resampleRegions, universe=glm27_allPos_noNW_noMit_GR,
                                     evaluate.function=numOverlaps, B=WGSNPs_2sigmas_NConly_nomit_100k_GR_trim, verbose=FALSE, mc.set.seed=FALSE)  
        plot(ptsignifCpGTime_signifWG)

#Supplementary Figure XX9: 
#Differently expressed genes and differently methylated CpGs do not overlap more or less than expected.
        ptsignifCpG_DE <- permTest(A=glm27_signifmorph_noNW_nomit_GR, ntimes=1000, randomize.function=resampleRegions, universe=glm27_allPos_noNW_noMit_GR,
                                     evaluate.function=numOverlaps, B=DEgenes_NConly_uniq_20kbpeak_withoutmissinggenes_GR, verbose=FALSE)  
        plot(ptsignifCpG_DE)

#Supplemental Figure 9: Differently methylated CpGs are enriched in all regions sequenced by RADseq, regardless of FST value. 
    #A) Differentially methylated CpGs between morphs are enriched close to outlier SNPs from the ddRAD. 
        ptsignifCpG_signifRAD <- permTest(A=glm27_signifmorph_noNW_nomit_GR, ntimes=1000, randomize.function=resampleRegions, universe=glm27_allPos_noNW_noMit_GR,
                                      evaluate.function=numOverlaps, B=RAD_2sigmas_20kpeak_noNW_nomit_GR, verbose=FALSE, mc.set.seed=FALSE)  
        plot(ptsignifCpG_signifRAD)

    #B) Differentially methylated CpGs between morphs are also enriched close to all genomic positions sequenced by ddRAD, regardless of FST value.
        ptsignifCpG_allRAD <- permTest(A=glm27_signifmorph_noNW_nomit_GR, ntimes=1000, randomize.function=resampleRegions, universe=glm27_allPos_noNW_noMit_GR,
                                      evaluate.function=numOverlaps, B=RADallpos_20kpeak_noNW_nomit_GR_trim, verbose=FALSE, mc.set.seed=FALSE)  
        plot(ptsignifCpG_allRAD)

#Supplementary Figure xxxx: 
#ddRAD SNPs are not enriched in regions close to differentially methylated CpGs, compared with the rest of the ddRAD dataset.
        ptsignifRAD_signifCpG <- permTest(A=RAD_2sigmas_20kpeak_noNW_nomit_GR, ntimes=1000, randomize.function=resampleRegions, universe=RADallpos_20kpeak_noNW_nomit_GR_trim,
                                      evaluate.function=numOverlaps, B=glm27_signifmorph_noNW_nomit_GR, verbose=FALSE, mc.set.seed=FALSE)  
        plot(ptsignifRAD_signifCpG)


####################################### III. MISCELLANEOUS / OLD / EXAMPLES ######################################################
# Overlaps between whole datasets, randomize over the whole genome 
# If you do not want to use the resampleRegions function, you could decide to do the comparison with randomisation from the whole genome
# In this case use the assembly as GRanges as "genome"
    #Example: RRBS vs RAD
        RRBSvsRAD <- overlapPermTest(A=glm27_allPos_noNW_noMit_GR, B=RADallpos_20kpeak_noNW_nomit_GR_trim, ntimes=50, genome=genomecharr, mc.set.seed=FALSE)
        plot(RRBSvsRAD)

# If you just want to calculate the number of overlaps between two datasets, use the numOverlaps() function
    #Example: overlaps between all CpGs in RRBS and the whole RAD dataset
        numOverlaps(glm27_allPos_noNW_noMit_GR, RADallpos_20kpeak_noNW_nomit_GR_trim) 

