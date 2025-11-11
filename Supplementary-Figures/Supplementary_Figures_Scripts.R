#Script to make supplementary figures with the right colour code
library(dplyr)
library(ggplot2)

#######################

#Supplementary Figure 2
#Load ddRAD data
newRADlist <- read.csv("/Users/sebma/Desktop/SNP_SLH-main/SNP_SLH-main/methylation/vcftools_output_merged.tsv", sep="\t")
#Keep only NC_ and no mit
newRADlist_NConlymit <- newRADlist %>% filter(grepl('NC_', NCBI)) #8480 SNPs
newRADlist_NConly <- newRADlist_NConlymit %>% filter(!grepl('NC_000861.1', NCBI))
#Find out the highest Fst Value for a window out of all three morphs
newRADlist_NConly$highest <- pmax(newRADlist_NConly$SB,newRADlist_NConly$LB,newRADlist_NConly$PL)

hist_ddRAD_LB<-ggplot(newRADlist_NConly, aes(x=LB)) + geom_histogram(color="#74add1", fill="#74add1") + 
  theme_bw() + coord_cartesian((xlim = c(0, 0.9)),(ylim=c(0,5000)))
hist_ddRAD_LB

hist_ddRAD_SB<-ggplot(newRADlist_NConly, aes(x=SB)) + geom_histogram(color="#313695", fill="#313695") + theme_bw() + 
  coord_cartesian((xlim = c(0, 0.9)),(ylim=c(0,5000)))
hist_ddRAD_SB

hist_ddRAD_PL<-ggplot(newRADlist_NConly, aes(x=PL)) + geom_histogram(color="#a50026", fill="#a50026") + theme_bw() +
  coord_cartesian((xlim = c(0, 0.9)),(ylim=c(0,5000)))
hist_ddRAD_PL

#Plot distribution of highest values, add intercept line at mean+2sigmas:
hist_ddRAD_highest<-ggplot(newRADlist_NConly, aes(x=highest)) + geom_histogram() + theme_bw() + 
  coord_cartesian((xlim = c(0, 0.9)),(ylim=c(0,5000))) + geom_vline(xintercept=0.415, colour="red")
hist_ddRAD_highest

#Export all plots in .pdf and put them together in Inkscape
pdf("/Users/sebma/Desktop/hist_ddRAD_LB.pdf",10,10)
print(hist_ddRAD_LB)
dev.off() 

pdf("/Users/sebma/Desktop/hist_ddRAD_SB.pdf",10,10)
print(hist_ddRAD_SB)
dev.off() 

pdf("/Users/sebma/Desktop/hist_ddRAD_PL.pdf",10,10)
print(hist_ddRAD_PL)
dev.off() 

pdf("/Users/sebma/Desktop/hist_ddRAD_highest.pdf",10,10)
print(hist_ddRAD_highest)
dev.off() 

##########################

#Supplementary Figure 3:
#Load WGS data
WGpeaks <- read.table("/Users/sebma/Desktop/RAD_WG_Meth/fst_window_WG.tsv", sep=",", header = TRUE)
#Only keep the data on NC_ and no mit
WGpeaks_NConlymit <- WGpeaks %>% filter(grepl('NC_', CHROM)) #(151540 obs)
WGpeaks_NConly <- WGpeaks_NConlymit %>% filter(!grepl('NC_000861.1', CHROM)) #(151540 obs)
########ATH I realize there are some negative Fst values. 
#Put negative values at 0.
WGpeaks_NConly[WGpeaks_NConly < 0] <- 0
#Find out the highest Fst Value for a window out of all three morphs
WGpeaks_NConly$highest <- pmax(WGpeaks_NConly$fst_LB,WGpeaks_NConly$fst_SB,WGpeaks_NConly$fst_PL)

#Make histogram with ggplot2
hist_WGS_LB<-ggplot(WGpeaks_NConly, aes(x=fst_LB)) + geom_histogram(color="#74add1", fill="#74add1") + 
  theme_bw() + coord_cartesian(xlim = c(0, 0.9))
hist_WGS_LB

hist_WGS_SB<-ggplot(WGpeaks_NConly, aes(x=fst_SB)) + geom_histogram(color="#313695", fill="#313695") + theme_bw() + 
  coord_cartesian(xlim = c(0, 0.9))
hist_WGS_SB

hist_WGS_PL<-ggplot(WGpeaks_NConly, aes(x=fst_PL)) + geom_histogram(color="#a50026", fill="#a50026") + theme_bw() +
  coord_cartesian(xlim = c(0, 0.9))
hist_WGS_PL

#Plot distribution of highest values, add intercept line at mean+2sigmas:
hist_WGS_highest<-ggplot(WGpeaks_NConly, aes(x=highest)) + geom_histogram() + theme_bw() + 
  coord_cartesian((xlim = c(0, 0.9)),(ylim=c(0,100000))) + geom_vline(xintercept=0.199, colour="red")
hist_WGS_highest

#Export all plots in .pdf and put them together in Inkscape
pdf("/Users/sebma/Desktop/hist_WGS_LB.pdf",10,10)
print(hist_WGS_LB)
dev.off() 

pdf("/Users/sebma/Desktop/hist_WGS_SB.pdf",10,10)
print(hist_WGS_SB)
dev.off() 

pdf("/Users/sebma/Desktop/hist_WGS_PL.pdf",10,10)
print(hist_WGS_PL)
dev.off() 

pdf("/Users/sebma/Desktop/hist_WGS_highest.pdf",10,10)
print(hist_WGS_highest)
dev.off() 

###############################################################################################

#Supplementary Figure 4
#ddRAD
    #Load data
    newRADlist <- read.csv("/Users/sebma/Desktop/SNP_SLH-main/SNP_SLH-main/methylation/vcftools_output_merged.tsv", sep="\t")
    #Keep only NC_ and no mit
    newRADlist_NConlymit <- newRADlist %>% filter(grepl('NC_', NCBI)) #8480 SNPs
    newRADlist_NConly <- newRADlist_NConlymit %>% filter(!grepl('NC_000861.1', NCBI))

    #Number of individual SNPs that pass thresholds (>0.2 or >0.5)
    ddRADoutliers_LB_0p2 <- newRADlist_NConly %>% filter(newRADlist_NConly$LB > 0.2)
    ddRADoutliers_SB_0p2 <- newRADlist_NConly %>% filter(newRADlist_NConly$SB > 0.2)
    ddRADoutliers_PL_0p2 <- newRADlist_NConly %>% filter(newRADlist_NConly$PL > 0.2)
    
    ddRADoutliers_LB_0p5 <- newRADlist_NConly %>% filter(newRADlist_NConly$LB > 0.5)
    ddRADoutliers_SB_0p5 <- newRADlist_NConly %>% filter(newRADlist_NConly$SB > 0.5)
    ddRADoutliers_PL_0p5 <- newRADlist_NConly %>% filter(newRADlist_NConly$PL > 0.5)
    
    #Count them
    num_RAD_LB_0p2 <- nrow(ddRADoutliers_LB_0p2)
    num_RAD_SB_0p2 <- nrow(ddRADoutliers_SB_0p2)
    num_RAD_PL_0p2 <- nrow(ddRADoutliers_PL_0p2)
    
    num_RAD_LB_0p5 <- nrow(ddRADoutliers_LB_0p5)
    num_RAD_SB_0p5 <- nrow(ddRADoutliers_SB_0p5)
    num_RAD_PL_0p5 <- nrow(ddRADoutliers_PL_0p5)
    
    #Make df out of them for plotting purposes
    numberpassingthreshold_RAD <- c(num_RAD_LB_0p2,num_RAD_SB_0p2,num_RAD_PL_0p2,num_RAD_LB_0p5,num_RAD_SB_0p5,num_RAD_PL_0p5)
    threshold <- c(">0.2",">0.2",">0.2",">0.5",">0.5",">0.5")
    morph <- c("LB","SB","PL","LB","SB","PL")
    
    ddRAD_Table <- data.frame(numberpassingthreshold_RAD,threshold,morph)

    #Make barplot:
    color <- c("#74add1","#313695","#a50026")
    names(color) <- c("LB","SB","PL")
    
    RAD_supp4 <- ggplot(ddRAD_Table, aes(x=threshold, y=numberpassingthreshold_RAD, fill=morph)) +
    geom_bar(stat = "identity", position = "dodge") + theme_bw() + scale_fill_manual(values = color)
    RAD_supp4

    #Save .pdf
    pdf("/Users/sebma/Desktop/ddRAD_suppFig4.pdf",10,10)
    print(RAD_supp4)
    dev.off()
    
#WG    
    #Load WGS data
    WGpeaks <- read.table("/Users/sebma/Desktop/RAD_WG_Meth/fst_window_WG.tsv", sep=",", header = TRUE)
    #Only keep the data on NC_ and no mit
    WGpeaks_NConlymit <- WGpeaks %>% filter(grepl('NC_', CHROM)) #(151540 obs)
    WGpeaks_NConly <- WGpeaks_NConlymit %>% filter(!grepl('NC_000861.1', CHROM)) #(151540 obs)
    #Put negative values at 0.
    WGpeaks_NConly[WGpeaks_NConly < 0] <- 0

    #Number of individual SNPs that pass thresholds (>0.2 or >0.5)
    WGoutliers_LB_0p2 <- WGpeaks_NConly %>% filter(WGpeaks_NConly$fst_LB > 0.2)
    WGoutliers_SB_0p2 <- WGpeaks_NConly %>% filter(WGpeaks_NConly$fst_SB > 0.2)
    WGoutliers_PL_0p2 <- WGpeaks_NConly %>% filter(WGpeaks_NConly$fst_PL > 0.2)
    
    WGoutliers_LB_0p5 <- WGpeaks_NConly %>% filter(WGpeaks_NConly$fst_LB > 0.5)
    WGoutliers_SB_0p5 <- WGpeaks_NConly %>% filter(WGpeaks_NConly$fst_SB > 0.5)
    WGoutliers_PL_0p5 <- WGpeaks_NConly %>% filter(WGpeaks_NConly$fst_PL > 0.5)
    
    #Count them
    num_WG_LB_0p2 <- nrow(WGoutliers_LB_0p2)
    num_WG_SB_0p2 <- nrow(WGoutliers_SB_0p2)
    num_WG_PL_0p2 <- nrow(WGoutliers_PL_0p2)
    
    num_WG_LB_0p5 <- nrow(WGoutliers_LB_0p5)
    num_WG_SB_0p5 <- nrow(WGoutliers_SB_0p5)
    num_WG_PL_0p5 <- nrow(WGoutliers_PL_0p5)
    
    numberpassingthreshold_WG <- c(num_WG_LB_0p2,num_WG_SB_0p2,num_WG_PL_0p2,num_WG_LB_0p5,num_WG_SB_0p5,num_WG_PL_0p5)
    threshold <- c(">0.2",">0.2",">0.2",">0.5",">0.5",">0.5")
    morph <- c("LB","SB","PL","LB","SB","PL")
    
    WG_Table <- data.frame(numberpassingthreshold_WG,threshold,morph)
    
    #Make barplot:
    color <- c("#74add1","#313695","#a50026")
    names(color) <- c("LB","SB","PL")
    
    WG_supp4 <- ggplot(WG_Table, aes(x=threshold, y=numberpassingthreshold_WG, fill=morph)) +
      geom_bar(stat = "identity", position = "dodge") + theme_bw() + scale_fill_manual(values = color)
    WG_supp4
    
    #Save .pdf
    pdf("/Users/sebma/Desktop/WG_suppFig4.pdf",10,10)
    print(WG_supp4)
    dev.off()

    #Then put the Figure together in Inkscape

###################
    
#Supplementary Figure 5:
#For now I have no idea how it was done, but it was something similar to the script TopWindows2.R

