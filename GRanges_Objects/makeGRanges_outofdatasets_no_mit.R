#### This script takes the dfs of interest and transforms them into Formal class GRanges object
library(genomation)
library(GenomeInfoDb) #For SeqInfo
library(dplyr)
library(GenomicRanges) #For makeGRangesFromDataFrame
library(ggplot2)

#Make the seqinfo (i.e. the name of each scaffold and their lengths)
    #Load the table made/extracted by Lea
    ACchrom <- read.csv("/Users/sebma/Desktop/GRanges_Objects/AC_chr_lengths.csv")
    #Only keep info for the placed scaffolds (NOTE: NC_000861.1 is the mitochondrial genome)
    ACchrom_noNW <- ACchrom %>% filter(grepl('NC_', NCBI))
    #Remove mitochondrial scaffold
    ACchrom_noNW_nomit <- ACchrom_noNW %>% filter(!grepl('NC_000861.1', NCBI))
    #change NC_ into chr
    ACchrom_noNW_nomit$chrom <- gsub("NC_","chr",ACchrom_noNW_nomit$NCBI)
    #Make the Seqinfo object
    charrgenomeseqinfo <- Seqinfo(seqnames = ACchrom_noNW_nomit$chrom, seqlengths = ACchrom_noNW_nomit$Length, genome = "canada_charr")

    
#Make the function that I can use so I don't have to always update it:
    makeGRfunction <- function(dfname) {
      
     newdfGR <- makeGRangesFromDataFrame(dfname,
                               keep.extra.columns=TRUE,
                               ignore.strand=FALSE,
                               seqinfo=charrgenomeseqinfo,
                               seqnames.field=c("seqnames", "seqname",
                                                "chromosome", "chrom",
                                                "chr", "chromosome_name",
                                                "seqid"),
                               start.field="start",
                               end.field=c("end", "stop"),
                               strand.field="strand",
                               starts.in.df.are.0based=FALSE)
      
      return(newdfGR)
      
    }
    
    
####################### REFERENCE GENOME
#Make the whole genome as a Granges object (needed for RegioneR to work) 
    #actually not sure if needed if we do the resample regions but will generate the data anyways
    genomedf <- data.frame("chrom"=ACchrom_noNW_nomit$chrom, "start"=1,"end"=ACchrom_noNW_nomit$Length)
 
    genome_charr_nomit_GR <- makeGRfunction(genomedf)
        # saveRDS(genome_charr_nomit_GR, "/Users/sebma/Desktop/GRanges_Objects/genome_noNW_nomit_GR.rds")
    
    
    
################################## METH Here, I want to combine signif morph with signif morph*time and morph*sex
#Load methylation data
    #the ones that are signif by morph with glm (181)
        signifMorph <- read.csv("/Users/sebma/Desktop/27samples/glm27_signif_Morphs_noNC_nomit.csv")
        
    #The ones that are signif by morph*time (27)   
        signifMorphxTime <- read.csv("/Users/sebma/Desktop/27samples/glm27_signif_MorphsxTime_noNC_nomit.csv")
        
    #The ones signif by MorphxSex (28)
        signifMorphxSex <- read.csv("/Users/sebma/Desktop/27samples/glm27_signif_MorphsxSex_noNC_nomit.csv")
        
    #Combine them, remove duplicates (left with 217)
        signifMorphall <- rbind(signifMorph,signifMorphxTime,signifMorphxSex)
        signifMorphall_nodup <- signifMorphall[!duplicated(signifMorphall[ , "Residue"]),]

        #write.csv(signifMorphall_nodup, "/Users/sebma/Desktop/27samples/glm27_signif_MorphsAll_noNC_nomit.csv", row.names = FALSE)
        
                  #ALREADY DONE PRIOR #Now to compare these to RAD and WG, better to have only the residues in NC_  (195 CpGs) 
                                      #Also remove mitochondrial positions
                                      #signifMorph_NConly_mit <- signifMorphall_nodup %>% filter(grepl('NC_', NW))
                                      #signifMorph_NConly <- signifMorph_NConly_mit %>% filter(!grepl('NC_000861.1', NCBI))
        
        signifMorph_NConly <- signifMorphall_nodup
        
        #Make a recap table with just scaffold and position
        recapsignif <- data.frame(signifMorph_NConly$NCBI,signifMorph_NConly$start)
        colnames(recapsignif) <- c("NCBI_MET","POS_MET")
        recapsignif$ID_MET <- paste(recapsignif$NCBI_MET,recapsignif$POS_MET,sep="_")
    
    #Or the full table    
        #load the full list of 14427 residues
        allCpG <- read.table("/Users/sebma/Desktop/27samples/methmin1_27_noPIno100.csv", header=TRUE, sep=",")
        allCpG_POS <- allCpG[,c(1,2)]
        
        #Format the chr into NC/NW
        allCpG_POS$temp <- gsub("chr","",allCpG_POS$chr)
        allCpG_POS$NW <- ifelse(nchar(allCpG_POS$temp)==11, "NW_", "NC_")
        allCpG_POS$NCBI_nonmet <- paste(allCpG_POS$NW,allCpG_POS$temp, sep="")
        allCpG_POS$ID1 <- paste(allCpG_POS$NCBI_nonmet,allCpG_POS$start,sep="_")
        #Only keep the residues in NC_ and not mit (9599 CpGs)
        allCpG_POS_NConlymit <- allCpG_POS %>% filter(grepl('NC_', NW))
        allCpG_POS_NConly <- allCpG_POS_NConlymit %>% filter(!grepl('NC_000861.1', NCBI_nonmet))

        
        #The start and chrom are already there. Just need to make the end column
        colnames(allCpG_POS_NConly)[1] <- "chrom"
        allCpG_POS_NConly$end <- allCpG_POS_NConly$start+1
        
        allCpG_POS_NConlyGR <- makeGRfunction(allCpG_POS_NConly)
        # saveRDS(allCpG_POS_NConlyGR, "/Users/sebma/Desktop/GRanges_Objects/glm27_allPos_noNW_nomit_GR.rds")
        
        
        
    #OR the ones that are "non-signif". 
        #For that, load the full list of 14427 residues
        allCpG <- read.table("/Users/sebma/Desktop/27samples/methmin1_27_noPIno100.csv", header=TRUE, sep=",")
        allCpG_POS <- allCpG[,c(1,2)]
        
        #Format the chr into NC/NW
        allCpG_POS$temp <- gsub("chr","",allCpG_POS$chr)
        allCpG_POS$NW <- ifelse(nchar(allCpG_POS$temp)==11, "NW_", "NC_")
        allCpG_POS$NCBI_nonmet <- paste(allCpG_POS$NW,allCpG_POS$temp, sep="")
        allCpG_POS$ID1 <- paste(allCpG_POS$NCBI_nonmet,allCpG_POS$start,sep="_")
        
        #Only keep the residues in NC_ and not mit (9599 CpGs)
        allCpG_POS_NConlymit <- allCpG_POS %>% filter(grepl('NC_', NW))
        allCpG_POS_NConly <- allCpG_POS_NConlymit %>% filter(!grepl('NC_000861.1', NCBI_nonmet))
        
        #make a vector of the positions in the methylated data
        diffmeth <- recapsignif$ID_MET 
        #Substract them from the full residue list. to have non methylation data (9382 residues). #SHOULD I DO THAT?
        nondiffmethCpG <- allCpG_POS_NConly %>% filter(!grepl((paste(diffmeth,collapse="|")), ID1))

                
        
# Transform methylation data into GRanges. Need to turn NC back in chr, make start and end
    #First the signif CpGs on NC_ only
        recapsignif$chrom <- gsub("NC_","chr",recapsignif$NCBI_MET)
        recapsignif$start <- recapsignif$POS_MET        
        recapsignif$end <- recapsignif$start+1

        recapsignifGR <- makeGRfunction(recapsignif)
        #saveRDS(recapsignifGR, "/Users/sebma/Desktop/GRanges_Objects/glm27_signifmorphall_noNW_nomit_GR.rds")
    
    #Second the non-signif CpGs on NC_ only  
        colnames(nondiffmethCpG)[1] <- "chrom"
        nondiffmethCpG$end <- nondiffmethCpG$start+1
        
        nondiffmethCpG_GR <- makeGRfunction(nondiffmethCpG)
        saveRDS(nondiffmethCpG_GR, "/Users/sebma/Desktop/GRanges_Objects/glm27_nonsignifmorphall_noNW_nomit_GR.rds")
        


    #### NEGATIVE CONTROL: SIGNIF METHYLATION BY TIME AND TIMExSEX
        #Load methylation data
        #the ones that are signif by time with glm (545)
        signifTime <- read.csv("/Users/sebma/Desktop/27samples/glm27_signif_Time_noNC_nomit.csv")
        
        #The ones that are signif by morph*time (27)   #This one is called bis just because it has an extra column that allows for rbind. But data is the same
        signifMorphxTime <- read.csv("/Users/sebma/Desktop/27samples/glm27_signif_MorphsxTime_noNC_nomit_bis.csv")
        
        #The ones signif by time*sex (11)
        signifTimexSex <- read.csv("/Users/sebma/Desktop/27samples/glm27_signif_TimexSex_noNC_nomit.csv")
        
        #Combine them, remove duplicates (left with 569)
        signifTimeall <- rbind(signifTime,signifMorphxTime,signifTimexSex)
        signifTimeall_nodup <- signifTimeall[!duplicated(signifTimeall[ , "Residue"]),]
        
        signifTime_NConly <- signifTimeall_nodup
        
        #Make a recap table with just scaffold and position
        recapsignif <- data.frame(signifTime_NConly$NCBI,signifTime_NConly$start)
        colnames(recapsignif) <- c("NCBI_MET","POS_MET")
        recapsignif$ID_MET <- paste(recapsignif$NCBI_MET,recapsignif$POS_MET,sep="_")
        
    # Transform methylation data into GRanges. Need to turn NC back in chr, make start and end
        #First the signif CpGs on NC_ only
        recapsignif$chrom <- gsub("NC_","chr",recapsignif$NCBI_MET)
        recapsignif$start <- recapsignif$POS_MET        
        recapsignif$end <- recapsignif$start+1
        
        recapsignifGR <- makeGRfunction(recapsignif)
        #saveRDS(recapsignifGR, "/Users/sebma/Desktop/GRanges_Objects/glm27_signiftimeall_noNW_nomit_GR.rds")
        
        
                
               
################################## RAD-seq
    #Load the RADseq SNPs
        newRADlist <- read.csv("/Users/sebma/Desktop/SNP_SLH-main/SNP_SLH-main/methylation/vcftools_output_merged.tsv", sep="\t")
        #Keep only NC_ and no mit
        newRADlist_NConlymit <- newRADlist %>% filter(grepl('NC_', NCBI)) #8480 SNPs
        newRADlist_NConly <- newRADlist_NConlymit %>% filter(!grepl('NC_000861.1', NCBI)) #8472 SNPs
        
     #Need to find SNPs which are lying outside 2sigmas of the distribution.
                  #The thing is, I don't know how to select them, because the PL will always outweight the others.
                  #I guess let's just take the highest of the three each time.
                  #Then do the distribution with the high numbers and see what lies above and below 2sigmas.
                  
                  newRADlist_NConly$highest <- pmax(newRADlist_NConly$SB,newRADlist_NConly$LB,newRADlist_NConly$PL)
                  #hist(newRADlist_NConly$highest)
                  
                  #Find what two sigmas of the distribution represents.
                  max <- mean(newRADlist_NConly$highest) + 2*sd(newRADlist_NConly$highest)
                  ddRAD_SNPs <- filter(newRADlist_NConly, highest > max)
                  
                  ####ATH! This is not a normal distrib. I think it is folded normal.
                  #### But technically it shouldn't change anything. I think I can still keep two sigmas
                  ddRAD_SNPs$chrom <- gsub("NC_","chr",ddRAD_SNPs$NCBI) 
                  #Make the peaks. Here 10 000bp each direction so 20000bp peaks
                  ddRAD_SNPs$start <- ddRAD_SNPs$POS - 10000
                  #Make sure there is no negative numbers for the SNPs that are close to the beginning of a contig
                  for(i in 1:nrow(ddRAD_SNPs)){
                    if(ddRAD_SNPs$start[i] <= 0){
                      ddRAD_SNPs$start[i] = 1
                    }
                  }
                  #Make end of peak:
                  ddRAD_SNPs$end <- ddRAD_SNPs$POS + 10000
                  #write.csv(ddRAD_SNPs, "/Users/sebma/Desktop/GRanges_Objects/RAD_2sigmas_20kpeak_noNW_nomit.csv", row.names = FALSE)
                  
                  RADSNPs_2sigma_20kpeak_GR <- makeGRfunction(ddRAD_SNPs)
                  saveRDS(RADSNPs_2sigma_20kpeak_GR, "/Users/sebma/Desktop/GRanges_Objects/RAD_2sigmas_20kpeak_noNW_nomit_GR.rds")
        
                  
                  #Plot the distribution, just for info:
                  ggplot(newRADlist_NConly, aes(x=highest)) + geom_histogram(bins = 30) + theme_bw() + labs(x="Fst", y="Frequency", title = "Fst distribution RADseq positions") + 
                    geom_vline(aes(xintercept = max), color= "red") + geom_text(aes(label = "mean + 2σ", x=0.52, y= 1600), color="red", size=3) + theme(plot.title = element_text(hjust = 0.5))
 
#########################                  
              #Separate the significant dataset to make one dataset per morph ATH: The way I was doing it before is biased
                  #By the SNPs with same Fst value in all morphs
                  # ddRAD_SNPs$HighMorph <- if_else(ddRAD_SNPs$LB == ddRAD_SNPs$highest, "LB","")
                  # ddRAD_SNPs$HighMorph <- if_else(ddRAD_SNPs$SB == ddRAD_SNPs$highest, "SB",ddRAD_SNPs$HighMorph)
                  # ddRAD_SNPs$HighMorph <- if_else(ddRAD_SNPs$PL == ddRAD_SNPs$highest, "PL",ddRAD_SNPs$HighMorph)
                  
                  #Sometimes, two morphs have the same Fst value for the same SNP. Thus I do not know whether to consider 
                  #the SNP from one morph or another. Because of this, I consider it as significant for both morphs.
                  ddRAD_SNPs_LB <- ddRAD_SNPs %>% filter(ddRAD_SNPs$highest == ddRAD_SNPs$LB)
                  ddRAD_SNPs_SB <- ddRAD_SNPs %>% filter(ddRAD_SNPs$highest == ddRAD_SNPs$SB)
                  ddRAD_SNPs_PL <- ddRAD_SNPs %>% filter(ddRAD_SNPs$highest == ddRAD_SNPs$PL)
                  
                  #Make GRange objects (no need to trim here)
                  ddRAD_SNPs_LBonly_2sigmas_NConly_nomit_20k_GR <- makeGRfunction(ddRAD_SNPs_LB)
                  saveRDS(ddRAD_SNPs_LBonly_2sigmas_NConly_nomit_20k_GR, "/Users/sebma/Desktop/GRanges_Objects/ddRAD_SNPs_LBonly_NOBIAS_2sigmas_NConly_nomit_20k_GR.rds")
                  
                  ddRAD_SNPs_SBonly_2sigmas_NConly_nomit_20k_GR <- makeGRfunction(ddRAD_SNPs_SB)
                  saveRDS(ddRAD_SNPs_SBonly_2sigmas_NConly_nomit_20k_GR, "/Users/sebma/Desktop/GRanges_Objects/ddRAD_SNPs_SBonly_NOBIAS_2sigmas_NConly_nomit_20k_GR.rds")
                  
                  ddRAD_SNPs_PLonly_2sigmas_NConly_nomit_20k_GR <- makeGRfunction(ddRAD_SNPs_PL)
                  saveRDS(ddRAD_SNPs_PLonly_2sigmas_NConly_nomit_20k_GR, "/Users/sebma/Desktop/GRanges_Objects/ddRAD_SNPs_PLonly_NOBIAS_2sigmas_NConly_nomit_20k_GR.rds")
                  
                  #Then I also need to make the inverse datasets.. ATH I DO NOT NEED TO DO THIS
                  # ddRAD_SNPs_PLSB <- ddRAD_SNPs %>% filter(!grepl('LB', HighMorph))
                  # ddRAD_SNPs_PLLB <- ddRAD_SNPs %>% filter(!grepl('SB', HighMorph))
                  # ddRAD_SNPs_SBLB <- ddRAD_SNPs %>% filter(!grepl('PL', HighMorph))
                  # 
                  # ddRAD_SNPs_PLSB_2sigmas_NConly_nomit_100k_GR <- makeGRfunction(ddRAD_SNPs_PLSB)
                  # saveRDS(ddRAD_SNPs_PLSB_2sigmas_NConly_nomit_100k_GR, "/Users/sebma/Desktop/GRanges_Objects/ddRAD_SNPs_PLSB_2sigmas_NConly_nomit_100k_GR.rds")
                  # 
                  # ddRAD_SNPs_PLLB_2sigmas_NConly_nomit_100k_GR <- makeGRfunction(ddRAD_SNPs_PLLB)
                  # saveRDS(ddRAD_SNPs_PLLB_2sigmas_NConly_nomit_100k_GR, "/Users/sebma/Desktop/GRanges_Objects/ddRAD_SNPs_PLLB_2sigmas_NConly_nomit_100k_GR.rds")
                  # 
                  # ddRAD_SNPs_SBLB_2sigmas_NConly_nomit_100k_GR <- makeGRfunction(ddRAD_SNPs_SBLB)
                  # saveRDS(ddRAD_SNPs_SBLB_2sigmas_NConly_nomit_100k_GR, "/Users/sebma/Desktop/GRanges_Objects/ddRAD_SNPs_SBLB_2sigmas_NConly_nomit_100k_GR.rds")
                  # 
                 
######################## 
      #Keep only non-significant SNPs:  with Fst < 0.2 in all morphs. (6800 positions? That doesn't add up) 
                  ddRAD_nonSNPs <- filter(newRADlist_NConly, highest < max)
                  ddRAD_nonSNPs$chrom <- gsub("NC_","chr",ddRAD_nonSNPs$NCBI) 
                  #Make the peaks. Here 10 000bp each direction so 20000bp peaks
                  ddRAD_nonSNPs$start <- ddRAD_nonSNPs$POS - 10000
                  #Make sure there is no negative numbers for the SNPs that are close to the beginning of a contig
                  for(i in 1:nrow(ddRAD_nonSNPs)){
                    if(ddRAD_nonSNPs$start[i] <= 0){
                      ddRAD_nonSNPs$start[i] = 1
                    }
                  }
                  #Make end of peak:
                  ddRAD_nonSNPs$end <- ddRAD_nonSNPs$POS + 10000
                  RADSNPs_nonsignif_outof2sigma_20kpeak_GR <- makeGRfunction(ddRAD_nonSNPs)
                  #Need to trim
                  RADSNPs_nonsignif_outof2sigma_20kpeak_GR_trim <- trim(RADSNPs_nonsignif_outof2sigma_20kpeak_GR, use.names=TRUE)
                  saveRDS(RADSNPs_nonsignif_outof2sigma_20kpeak_GR_trim, "/Users/sebma/Desktop/GRanges_Objects/RADSNPs_nonsignif_outof2sigma_20kpeak_noNW_nomit_GR_trim.rds")

        
      #keep all "SNPs"
        newRADlist <- read.csv("/Users/sebma/Desktop/SNP_SLH-main/SNP_SLH-main/methylation/vcftools_output_merged.tsv", sep="\t")
        #Keep only NC_ and no mit
        newRADlist_NConlymit <- newRADlist %>% filter(grepl('NC_', NCBI))
        newRADlist_NConly <- newRADlist_NConlymit %>% filter(!grepl('NC_000861.1', NCBI))
        
        newRADlist_NConly$chrom <- gsub("NC_","chr",newRADlist_NConly$NCBI) 
        #Make the peaks. Here 10 000bp each direction so 20000bp peaks
        newRADlist_NConly$start <- newRADlist_NConly$POS - 10000
        #Make sure there is no negative numbers for the SNPs that are close to the beginning of a contig
        for(i in 1:nrow(newRADlist_NConly)){
          if(newRADlist_NConly$start[i] <= 0){
            newRADlist_NConly$start[i] = 1
          }
        }
        #Make end of peak:
        newRADlist_NConly$end <- newRADlist_NConly$POS + 10000
        
        #ATH, this produces a couple of out of range SNPs because the end goes further than the end of the scaffold
        newRADlist_NConly_GR <- makeGRfunction(newRADlist_NConly)
        newRADlist_NConly_GR_trim <- trim(newRADlist_NConly_GR, use.names=TRUE)
        saveRDS(newRADlist_NConly_GR_trim, "/Users/sebma/Desktop/GRanges_Objects/RADallpos_20kpeak_noNW_nomit_GR_trim.rds")
        
        # newRADlist100kb_NConly_GR <- makeGRfunction(newRADlist_NConly)
        # newRADlist100kb_NConly_GR_trim <- trim(newRADlist100kb_NConly_GR, use.names=TRUE)
        # saveRDS(newRADlist100kb_NConly_GR, "/Users/sebma/Desktop/GRanges_Objects/RADallpos_100kpeak_noNW_GR_trim.rds")

        
                
########################### WG-seq (there is no sequences on Mitochondria scaffold here so its the same)
    #Note, this is the WG data I am using, but I am still a little bit unsure as to where it came from
      #Load the data:
        WGpeaks <- read.table("/Users/sebma/Desktop/RAD_WG_Meth/fst_window_WG.tsv", sep=",", header = TRUE)
      #Only keep the data on NC_ and no mit
        WGpeaks_NConlymit <- WGpeaks %>% filter(grepl('NC_', CHROM)) #(151540 obs)
        WGpeaks_NConly <- WGpeaks_NConlymit %>% filter(!grepl('NC_000861.1', CHROM)) #(151540 obs)
        
    ########ATH I realize there are some negative Fst values. 
        #From what I can see online, these negative Fst should be considered 0 I think.
        #And they show that there are more variations within than between pop? I mean we are at 2 samples per morph so it 
        #would make sense that it is not great
        
        #Put negative values at 0.
        WGpeaks_NConly[WGpeaks_NConly < 0] <- 0
        
        
        #Need to find SNPs which are lying outside 2sigmas of the distribution.
        WGpeaks_NConly$highest <- pmax(WGpeaks_NConly$fst_LB,WGpeaks_NConly$fst_SB,WGpeaks_NConly$fst_PL)
        # hist(WGpeaks_NConly$highest)
        
        # write.csv(WGpeaks_NConly, "/Users/sebma/Desktop/thirdpaperfigures/tables/listofWGregions_NConly.csv", row.names = F)
        
        #Find what two sigmas of the distribution represents.
        maxWG <- mean(WGpeaks_NConly$highest) + 2*sd(WGpeaks_NConly$highest)
        WG_SNPs <- filter(WGpeaks_NConly, highest > maxWG)
        
        WG_SNPs$chrom <- gsub("NC_","chr",WG_SNPs$CHROM)
        colnames(WG_SNPs)[1] <- "NCBI" #Otherwise there are two columns chr and its ambiguous
        #These are actually 100kb windows that slide every 10kb.
        #I believe POS is the middle of each window. Because the first window on each LG always is POS=50000.
        #That makes way bigger regions than the RADseq
        WG_SNPs$start <- WG_SNPs$POS - 49999
        WG_SNPs$end <- WG_SNPs$POS + 50000
        
        #write.csv(WG_SNPs, "/Users/sebma/Desktop/27samples/WGSNPs_2sigmas_NConly_nomit_100k.csv", row.names = FALSE)
        
        WGSNPs_2sigmas_NConly_nomit_100k_GR <- makeGRfunction(WG_SNPs)
        #Trim the windows out of range
        WGSNPs_2sigmas_NConly_nomit_100k_GR_trim <- trim(WGSNPs_2sigmas_NConly_nomit_100k_GR, use.names=TRUE)
        # saveRDS(WGSNPs_2sigmas_NConly_nomit_100k_GR_trim, "/Users/sebma/Desktop/GRanges_Objects/WG_2sigma_NConly_nomit_100k_GR_trim.rds")
       
        #Plot the distribution, just for info:
        ggplot(WGpeaks_NConly, aes(x=highest)) + geom_histogram(bins = 30) + theme_bw() + labs(x="Fst", y="Frequency", title = "Fst distribution WGseq regions") + 
          geom_vline(aes(xintercept = maxWG), color= "red") + geom_text(aes(label = "mean + 2σ", x=0.30, y= 30000), color="red", size=3) + theme(plot.title = element_text(hjust = 0.5))
        
##############        
    #Separate the significant dataset and create one dataset per morph
        WG_SNPs$HighMorph <- if_else(WG_SNPs$fst_LB == WG_SNPs$highest, "LB","")
        WG_SNPs$HighMorph <- if_else(WG_SNPs$fst_SB == WG_SNPs$highest, "SB",WG_SNPs$HighMorph)
        WG_SNPs$HighMorph <- if_else(WG_SNPs$fst_PL == WG_SNPs$highest, "PL",WG_SNPs$HighMorph)
        
        WG_SNPs_LB <- WG_SNPs %>% filter(grepl('LB', HighMorph))
        WG_SNPs_SB <- WG_SNPs %>% filter(grepl('SB', HighMorph))
        WG_SNPs_PL <- WG_SNPs %>% filter(grepl('PL', HighMorph))
        
        #Make GRange objects and trim 
        WGSNPs_LBonly_2sigmas_NConly_nomit_100k_GR <- makeGRfunction(WG_SNPs_LB)
        WGSNPs_LBonly_2sigmas_NConly_nomit_100k_GR_trim <- trim(WGSNPs_LBonly_2sigmas_NConly_nomit_100k_GR, use.names=TRUE)
        saveRDS(WGSNPs_LBonly_2sigmas_NConly_nomit_100k_GR_trim, "/Users/sebma/Desktop/GRanges_Objects/WGSNPs_LBonly_2sigmas_NConly_nomit_100k_GR_trim.rds")
        
        WGSNPs_SBonly_2sigmas_NConly_nomit_100k_GR <- makeGRfunction(WG_SNPs_SB)
        WGSNPs_SBonly_2sigmas_NConly_nomit_100k_GR_trim <- trim(WGSNPs_SBonly_2sigmas_NConly_nomit_100k_GR, use.names = TRUE)
        saveRDS(WGSNPs_SBonly_2sigmas_NConly_nomit_100k_GR_trim, "/Users/sebma/Desktop/GRanges_Objects/WGSNPs_SBonly_2sigmas_NConly_nomit_100k_GR_trim.rds")
        
        WGSNPs_PLonly_2sigmas_NConly_nomit_100k_GR <- makeGRfunction(WG_SNPs_PL)
        WGSNPs_PLonly_2sigmas_NConly_nomit_100k_GR_trim <- trim(WGSNPs_PLonly_2sigmas_NConly_nomit_100k_GR, use.names = TRUE)
        saveRDS(WGSNPs_PLonly_2sigmas_NConly_nomit_100k_GR_trim, "/Users/sebma/Desktop/GRanges_Objects/WGSNPs_PLonly_2sigmas_NConly_nomit_100k_GR_trim.rds")
        
##############        
        
                
    #Keep only nonsignificant windows 
        WG_nonSNPs <- filter(WGpeaks_NConly, highest < maxWG)
        WG_nonSNPs$chrom <- gsub("NC_","chr",WG_nonSNPs$CHROM)
        colnames(WG_nonSNPs)[1] <- "NCBI" #Otherwise there are two columns chr and its ambiguous
        #These are actually 100kb windows that slide every 10kb.
        #I believe POS is the middle of each window. Because the first window on each LG always is POS=50000.
        #That makes way bigger regions than the RADseq
        WG_nonSNPs$start <- WG_nonSNPs$POS - 49999
        WG_nonSNPs$end <- WG_nonSNPs$POS + 50000
        
        WGSNPs_nonsignificant_outof2sigmas_NConly_nomit_100k_GR <- makeGRfunction(WG_nonSNPs)
        #Trim the windows out of range
        WGSNPs_nonsignificant_outof2sigmas_NConly_nomit_100k_GR_trim <- trim(WGSNPs_nonsignificant_outof2sigmas_NConly_nomit_100k_GR, use.names=TRUE)
        saveRDS(WGSNPs_nonsignificant_outof2sigmas_NConly_nomit_100k_GR_trim, "/Users/sebma/Desktop/GRanges_Objects/WG_nonsignificant_outof2sigma_NConly_nomit_100k_GR_trim.rds")
        
        
     #Keep all windows on NC_
        WGpeaks_NConly <- WGpeaks %>% filter(grepl('NC_', CHROM)) #(151540 obs)
        
        WGpeaks_NConly$chrom <- gsub("NC_","chr",WGpeaks_NConly$CHROM)
        colnames(WGpeaks_NConly)[1] <- "NCBI" #Otherwise there are two columns chr and its ambiguous
        
        WGpeaks_NConly$start <- WGpeaks_NConly$POS - 49999
        WGpeaks_NConly$end <- WGpeaks_NConly$POS + 50000
        
        WGpeaks_NConly_allwindows_GR <- makeGRfunction(WGpeaks_NConly)
        #Trim the windows out of range
        WGpeaks_NConly_nomit_allwindows_GR_trim <- trim(WGpeaks_NConly_allwindows_GR, use.names=TRUE)
        saveRDS(WGpeaks_NConly_nomit_allwindows_GR_trim, "/Users/sebma/Desktop/GRanges_Objects/WGpeaks_NConly_nomit_allwindows_GR_trim.rds")
        
        
########################### DE genes
        #Use C:/Users/sebma/Desktop/Rscripts/getGeneInfo_Transc.R to get the Scaffold and position of the genes
        
      # 1) All genes in transcriptome data
        allgenesFullInfo <- read.csv("/Users/sebma/Desktop/Alexander stuff/allgeneswithPosition.csv")
        #Only keep the data on NC_
        allgenesFullInfo_NConly <- allgenesFullInfo %>% filter(grepl('NC_', chromosome_RefSeq))
        #Only keep uniq geneID (basically only keep one of the transcripts if there is more than one)
        allgenesFullInfo_NConly_uniq <- allgenesFullInfo_NConly[!duplicated(allgenesFullInfo_NConly$geneID),]
        #Take the middle of the gene? the start?
        #Some genes are huge..
        #OK what about in that case I do the peak around the CpG?
        #And I just make the GRange for the genes whatever size they are supposed to be, at least for a specific transcript.
        allgenesFullInfo_NConly_uniq$chrom <- gsub("NC_","chr",allgenesFullInfo_NConly_uniq$chromosome_RefSeq) 
        
        #allgenesFullInfo_NConly_uniq_GR <- WGpeaks_NC_Signif0p2(allgenesFullInfo_NConly_uniq)
        #saveRDS(allgenesFullInfo_NConly_uniq_GR, "/Users/sebma/Desktop/GRanges_Objects/allgenesFullInfo_NConly_uniq_GR.rds")

      # 2) All genes but this time with a 20kb peak around them:
        allgenesFullInfo_NConly_uniq$start <- allgenesFullInfo_NConly_uniq$start - 10000
        allgenesFullInfo_NConly_uniq$stop <- allgenesFullInfo_NConly_uniq$stop + 10000
        #Make sure there is no negative numbers for the genes that are close to the beginning of a contig
        for(i in 1:nrow(allgenesFullInfo_NConly_uniq)){
          if(allgenesFullInfo_NConly_uniq$start[i] <= 0){
            allgenesFullInfo_NConly_uniq$start[i] = 1
          }
        }
        
        allgenesFullInfo_NConly_uniq_20kbpeak_GR <- makeGRfunction(allgenesFullInfo_NConly_uniq)
        #Trim.
        allgenesFullInfo_NConly_uniq_20kbpeak_GR_trim <- trim(allgenesFullInfo_NConly_uniq_20kbpeak_GR, use.names=TRUE)
        #saveRDS(allgenesFullInfo_NConly_uniq_20kbpeak_GR_trim, "/Users/sebma/Desktop/GRanges_Objects/allgenesFullInfo_NConly_uniq_20kbpeak_GR_trim.rds")
        
        
        
        #3) Only DE genes
        DEgenesFullInfo <- read.csv("/Users/sebma/Desktop/Alexander stuff/DEgeneswithPosition.csv")
        #Only keep the data on NC_
        DEgenesFullInfo_NConly <- DEgenesFullInfo %>% filter(grepl('NC_', chromosome_RefSeq))
                  #There is nothing on mitochondria so no need to remove these
                  #DEgenesFullInfo_NConly_nomit <- DEgenesFullInfo_NConly %>% filter(!grepl('NC_000861.1', chromosome_RefSeq))
        #Only keep uniq geneID (need to ask what is the best way of doing things.)
        DEgenesFullInfo_NConly_uniq <- DEgenesFullInfo_NConly[!duplicated(DEgenesFullInfo_NConly$gene_name),]
        #Replace NC_ with chr
        DEgenesFullInfo_NConly_uniq$chrom <- gsub("NC_","chr",DEgenesFullInfo_NConly_uniq$chromosome_RefSeq) 
        
        #DEgenesFullInfo_NConly_uniq_GR <- makeGRfunction(DEgenesFullInfo_NConly_uniq)
        #saveRDS(DEgenesFullInfo_NConly_uniq_GR, "/Users/sebma/Desktop/GRanges_Objects/DEgenesFullInfo_NConly_uniq_GR.rds")
        
        
        #4) Only DE genes but this time with a 20kb peak around them:
              DEgenesFullInfo_NConly_uniq$start <- DEgenesFullInfo_NConly_uniq$start - 10000
              DEgenesFullInfo_NConly_uniq$stop <- DEgenesFullInfo_NConly_uniq$stop + 10000
        #Make sure there is no negative numbers for the genes that are close to the beginning of a contig
        for(i in 1:nrow(DEgenesFullInfo_NConly_uniq)){
          if(DEgenesFullInfo_NConly_uniq$start[i] <= 0){
            DEgenesFullInfo_NConly_uniq$start[i] = 1
          }
        }
  
        #write.csv(DEgenesFullInfo_NConly_uniq, "/Users/sebma/Desktop/GRanges_Objects/DEgenesFullInfo_NConly_uniq_20kbpeak.csv", row.names = FALSE)
        
        DEgenesFullInfo_NConly_uniq_20kbpeak_GR <- makeGRfunction(DEgenesFullInfo_NConly_uniq)
  #No need to trim this time, as there is no warning of out of bounds
  #DEgenesFullInfo_NConly_uniq_20kbpeak_GR_trim <- trim(DEgenesFullInfo_NConly_uniq_20kbpeak_GR, use.names=TRUE)
  
  saveRDS(DEgenesFullInfo_NConly_uniq_20kbpeak_GR, "/Users/sebma/Desktop/GRanges_Objects/DEgenesFullInfo_NConly_uniq_20kbpeak_GR.rds")

####################                  
  #Separate the significant dataset to make one dataset per morph.
  #Here, I have the choice to take either the uniq (which is the highest) or uniq + all. 
  #I will only take uniq, as I have taken the highest SNP for RAD and WG, even when there was more than 1morph.
  #I could do uniq+all in other datasets too and use uniq+all instead of uniq.
  DEgenesFullInfo_NConly_uniq_LBonly <- DEgenesFullInfo_NConly_uniq %>% filter(grepl('LB', direction_unique))
  DEgenesFullInfo_NConly_uniq_SBonly <- DEgenesFullInfo_NConly_uniq %>% filter(grepl('SB', direction_unique))
  DEgenesFullInfo_NConly_uniq_PLonly <- DEgenesFullInfo_NConly_uniq %>% filter(grepl('PL', direction_unique))
  
  #Make GRange objects (no need to trim here)
  DEgenesFullInfo_NConly_uniq_LBonly_20k_GR <- makeGRfunction(DEgenesFullInfo_NConly_uniq_LBonly)
  saveRDS(DEgenesFullInfo_NConly_uniq_LBonly_20k_GR, "/Users/sebma/Desktop/GRanges_Objects/DEgenesFullInfo_NConly_uniq_LBonly_20k_GR.rds")
  
  DEgenesFullInfo_NConly_uniq_SBonly_20k_GR <- makeGRfunction(DEgenesFullInfo_NConly_uniq_SBonly)
  saveRDS(DEgenesFullInfo_NConly_uniq_SBonly_20k_GR, "/Users/sebma/Desktop/GRanges_Objects/DEgenesFullInfo_NConly_uniq_SBonly_20k_GR.rds")
  
  DEgenesFullInfo_NConly_uniq_PLonly_20k_GR <- makeGRfunction(DEgenesFullInfo_NConly_uniq_PLonly)
  saveRDS(DEgenesFullInfo_NConly_uniq_PLonly_20k_GR, "/Users/sebma/Desktop/GRanges_Objects/DEgenesFullInfo_NConly_uniq_PLonly_20k_GR.rds")
  

  
  
  #Then I also need to make the inverse datasets.. WAIT I DO NOT NEED TO DO THIS
  
  
  ########################   
  
  
    
        #5) Non significant genes, with a 20kbp around them
            allgenesFullInfo <- read.csv("/Users/sebma/Desktop/Alexander stuff/allgeneswithPosition.csv")
            #Only keep the data on NC_
            allgenesFullInfo_NConly <- allgenesFullInfo %>% filter(grepl('NC_', chromosome_RefSeq))
            #Only keep uniq geneID (need to ask what is the best way of doing things.) #23627 genes
            allgenesFullInfo_NConly_uniq <- allgenesFullInfo_NConly[!duplicated(allgenesFullInfo_NConly$geneID),]
        
            DEgenesFullInfo <- read.csv("/Users/sebma/Desktop/Alexander stuff/DEgeneswithPosition.csv")
            #Only keep the data on NC_
            DEgenesFullInfo_NConly <- DEgenesFullInfo %>% filter(grepl('NC_', chromosome_RefSeq))
            #Only keep uniq geneID (need to ask what is the best way of doing things.) #1141 genes
            DEgenesFullInfo_NConly_uniq <- DEgenesFullInfo_NConly[!duplicated(DEgenesFullInfo_NConly$gene_name),]
       
            
            #make a vector of the geneID in the DE dataframe
            DEgenesid <- DEgenesFullInfo_NConly_uniq$gene_id
            #Substract them from the full gene list. #22632 obs... Something is not quite right. 
            #146 DE genes are missing from the whole gene table
            nonDEgenes <- allgenesFullInfo_NConly_uniq %>% filter(!grepl((paste(DEgenesid,collapse="|")), geneID))
            
            test<- allgenesFullInfo_NConly_uniq %>% filter(grepl((paste(DEgenesid,collapse="|")), geneID))
            testgeneiD <- test$geneID
            missinggenes <- DEgenesFullInfo_NConly_uniq %>% filter(!grepl((paste(testgeneiD,collapse="|")), gene_id))
  ###############################
      #Remove these from the DE list (done on 26/02/24)
            DEwithoutmissinggenes <- DEgenesFullInfo_NConly_uniq %>% filter(grepl((paste(testgeneiD,collapse="|")), gene_id))
            DEwithoutmissinggenes$chrom <- gsub("NC_","chr",DEwithoutmissinggenes$chromosome_RefSeq)
            
            DEwithoutmissinggenes$start <- DEwithoutmissinggenes$start - 10000
            DEwithoutmissinggenes$stop <- DEwithoutmissinggenes$stop + 10000
            #Make sure there is no negative numbers for the genes that are close to the beginning of a contig
            for(i in 1:nrow(DEwithoutmissinggenes)){
              if(DEwithoutmissinggenes$start[i] <= 0){
                DEwithoutmissinggenes$start[i] = 1
              }
            }
            
            DEgenes_NConly_uniq_20kbpeak_withoutmissinggenes_GR <- makeGRfunction(DEwithoutmissinggenes)
            #trim
            #DEgenes_NConly_uniq_20kbpeak_withoutmissinggenes_GR_trim <- trim(DEgenes_NConly_uniq_20kbpeak_withoutmissinggenes_GR, use.names=TRUE)
            saveRDS(DEgenes_NConly_uniq_20kbpeak_withoutmissinggenes_GR, "/Users/sebma/Desktop/GRanges_Objects/DEgenes_NConly_uniq_20kbpeak_withoutmissinggenes_GR.rds")
            
            #Separate the significant dataset to make one dataset per morph.
            DEwithoutmissinggenes_LBonly <- DEwithoutmissinggenes %>% filter(grepl('LB', direction_unique))
            DEwithoutmissinggenes_SBonly <- DEwithoutmissinggenes %>% filter(grepl('SB', direction_unique))
            DEwithoutmissinggenes_PLonly <- DEwithoutmissinggenes %>% filter(grepl('PL', direction_unique))
            
            #Make GRange objects (no need to trim here)
            DEwithoutmissinggenes_LBonly_20k_GR <- makeGRfunction(DEwithoutmissinggenes_LBonly)
            saveRDS(DEwithoutmissinggenes_LBonly_20k_GR, "/Users/sebma/Desktop/GRanges_Objects/DEgenesFullInfo_NConly_uniq_LBonly_nomissinggene_20k_GR.rds")
            
            DEwithoutmissinggenes_SBonly_20k_GR <- makeGRfunction(DEwithoutmissinggenes_SBonly)
            saveRDS(DEwithoutmissinggenes_SBonly_20k_GR, "/Users/sebma/Desktop/GRanges_Objects/DEgenesFullInfo_NConly_uniq_SBonly_nomissinggene_20k_GR.rds")
            
            DEwithoutmissinggenes_PLonly_20k_GR <- makeGRfunction(DEwithoutmissinggenes_PLonly)
            saveRDS(DEwithoutmissinggenes_PLonly_20k_GR, "/Users/sebma/Desktop/GRanges_Objects/DEgenesFullInfo_NConly_uniq_PLonly_nomissinggene_20k_GR.rds")
            
            
            
            
  ###############################          
            #Now I have the list of the 146 genes that are not in the whole table and I need to figure out why
            #If I don't figure out why, I need to figure out what to do about it: Should I just remove them from the data
            #
            
            
            nonDEgenes$chrom <- gsub("NC_","chr",nonDEgenes$chromosome_RefSeq)
            
            nonDEgenes$start <- nonDEgenes$start - 10000
            nonDEgenes$stop <- nonDEgenes$stop + 10000
            #Make sure there is no negative numbers for the genes that are close to the beginning of a contig
            for(i in 1:nrow(nonDEgenes)){
              if(nonDEgenes$start[i] <= 0){
                nonDEgenes$start[i] = 1
              }
            }
            
            nonDEgenes_NConly_uniq_20kbpeak_GR <- makeGRfunction(nonDEgenes)
            #trim
            nonDEgenes_NConly_uniq_20kbpeak_GR_trim <- trim(nonDEgenes_NConly_uniq_20kbpeak_GR, use.names=TRUE)
            saveRDS(nonDEgenes_NConly_uniq_20kbpeak_GR, "/Users/sebma/Desktop/GRanges_Objects/nonDEgenes_NConly_uniq_20kbpeak_GR_trim.rds")
            