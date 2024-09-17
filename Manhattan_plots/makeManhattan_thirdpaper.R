library(dplyr)
library(sqldf)
library(ggplot2)

############### LOAD INFO FOR EACH LINKAGE GROUP (SIZE, NAME, ETC...) ###################
    ACchrom <- read.csv("~/AC_chr_lengths.csv")
#Remove unplaced scaffolds
    ACchrom_noNW <- ACchrom %>% filter(grepl('NC_', NCBI))
#Remove mitochondrial scaffold
    ACchrom_noNW_nomit <- ACchrom_noNW %>% filter(!grepl('NC_000861.1', NCBI))
    names(ACchrom_noNW_nomit) <- c("NAME","NCBI","LENGTH","START","LG","ABLENGTH")
#Only keep name, start and length, add end
    chromosome_info <- sqldf('SELECT LG, min(START), ABLENGTH
                             FROM ACchrom_noNW_nomit
                             GROUP BY LG')
    names(chromosome_info) <- c("LG","START","LENGTH")
    chromosome_info <- chromosome_info[order(chromosome_info$START),]
    chromosome_info$END <- chromosome_info$START + chromosome_info$LENGTH
    chromosome_info$MIDDLE <- (chromosome_info$START+chromosome_info$END)/2
    starting <- chromosome_info$START
    ending <- chromosome_info$END
#Make a pattern that will be used for the alternate coloring on the plot
    chromosome_info$pattern <- rep(c(TRUE,FALSE), length.out=nrow(chromosome_info))


########################## LOAD RADseq DATA #################################
    newRADlist <- read.csv("~/vcftools_output_merged.tsv", sep="\t")
#Keep only NC_ (placed scaffolds)
    newRADlist_NConlymit <- newRADlist %>% filter(grepl('NC_', NCBI)) #8480 SNPs
#Remove mitochondrial scaffold
    newRADlist_NConly <- newRADlist_NConlymit %>% filter(!grepl('NC_000861.1', NCBI))
#Identify the highest fst value out of the three morphs
    newRADlist_NConly$highest <- pmax(newRADlist_NConly$SB,newRADlist_NConly$LB,newRADlist_NConly$PL)

#ATH: Here, I realized that for some SNPs, 1 morph has an Fst of 0 while the other 2 have the same value.
#Should these SNPs be taken into account? I remember us talking about it but can't remember what we ended up going with
#I think I removed them from the graph, but kept them in the general analysis 
#(see script *Superimposition-Paper/GRanges_Objects/makeGRanges_outofdatasets_nomit.R**)
#To remove them:
    #Identify rows where this happens
        newRADlist_NConly$row <- row.names(newRADlist_NConly)
        SBLB <-  filter(newRADlist_NConly, SB == LB, SB >0, LB>0)
        SBPL <- filter(newRADlist_NConly, SB == PL, SB >0, PL>0)
        LBPL <- filter(newRADlist_NConly, LB == PL, LB >0, PL>0)
        badrows <- as.numeric(c(SBLB$row,SBPL$row,LBPL$row))
    #Remove these rows from the table
        newRADlist_NConly_good <- newRADlist_NConly[-badrowsnodup,]

#Identify what two sigmas of the distribution represents
#This will be used to put the threshold on the plot
    max <- mean(newRADlist_NConly$highest) + 2*sd(newRADlist_NConly$highest) #Max for RADseq = 0.41497
#If you do not want these weird SNP mentioned before to be taken into account when calculating 2 sigmas, use the other table instead:
#    max <- mean(newRADlist_NConly_good$highest) + 2*sd(newRADlist_NConly_good$highest) #Max for RADseq = 0.41387
#I believe that I still went with a threshold of 0.41497 (calculated with all SNPs)

#Now for plotting purposes, I do not want to see these SNPs
#Separate table per morph, then rbind it.
    RAD_SB <- data.frame(newRADlist_NConly_good[,c(1:4)])
    RAD_LB <- data.frame(newRADlist_NConly_good[,c(1:3,5)])
    RAD_PL <- data.frame(newRADlist_NConly_good[,c(1:3,6)])

#Prepare colname for rbind
    colnames(RAD_SB)[4] <- "Fst"
    colnames(RAD_LB)[4] <- "Fst"
    colnames(RAD_PL)[4] <- "Fst"

#Add the name of the morph to prepare for the color
    RAD_SB$Comparison <- "SB vs LB/PL"
    RAD_LB$Comparison <- "LB vs SB/PL"
    RAD_PL$Comparison <- "PL vs SB/LB"

#Bind everything together
    RAD_cleanTable <- rbind(RAD_SB,RAD_LB,RAD_PL)

#If Fst less than max, change the Color to "Non-significant"
    RAD_cleanTable$Comparison <- ifelse(RAD_cleanTable$Fst < max, "Non-significant",RAD_cleanTable$Comparison)
#Relevel so I have non-significant at the end
    RAD_cleanTable$Comparison <- factor(RAD_cleanTable$Comparison)
    RAD_cleanTable$Comparison <- relevel(RAD_cleanTable$Comparison,"Non-significant")
    RAD_cleanTable$Comparison <- relevel(RAD_cleanTable$Comparison,"PL vs SB/LB")
    RAD_cleanTable$Comparison <- relevel(RAD_cleanTable$Comparison,"SB vs LB/PL")
    RAD_cleanTable$Comparison <- relevel(RAD_cleanTable$Comparison,"LB vs SB/PL")
















#I need to calculate the absolute Position of each dot. Because here it is just the basic position.
mergeData <- function(datatable){
  datatable$LG <- ACchrom_noNW_nomit$LG[match(datatable$NCBI,ACchrom_noNW_nomit$NCBI)]
  datatable$START <- ACchrom_noNW_nomit$START[match(datatable$NCBI,ACchrom_noNW_nomit$NCBI)]
  datatable$ABPOS <- datatable$START + datatable$POS
  return(datatable)
}

RAD_cleanTable2 <- mergeData(RAD_cleanTable)

#Now I have all columns of interest: NCBI, POS, Fst and Color, with Fst and POS numeric.
#Let's try to plot.

color <- c("#00BA38FF","#619CFFFF","#F8766DFF","black")
names(color) <- c("LB vs SB/PL","SB vs LB/PL","PL vs SB/LB","Non-significant")

LGplot <- ggplot()+ 
  geom_rect(aes(NULL,NULL,xmin=starting,xmax=ending,fill=chromosome_info$pattern),
            ymin=-Inf,ymax=Inf,data=chromosome_info) +
  scale_fill_manual(values = alpha(c("white","gray76"), .5)) +
  scale_x_continuous(label=chromosome_info$LG, breaks = chromosome_info$MIDDLE) +
  theme_classic()+ xlab("Linkage group") + ylab("Fst") + ggtitle("Positions of RADseq SNPs on charr genome") +
  theme(axis.text.x = element_text(angle = 90, size = 8, face = "bold", vjust = 0.5),
        axis.title.x = element_text(size=16,
                                    margin = margin(t=15,r=0,b=10,l=0)),
        axis.title.y = element_text(size=16,
                                    margin = margin(t=0,r=10,b=0,l=10)),
        panel.border = element_rect(colour = "black",fill = "NA", size = 1),
        legend.background = element_blank(),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.box.background = element_rect(color = "black")) +
  geom_point(aes(x=RAD_cleanTable2$ABPOS, y=RAD_cleanTable2$Fst, color=RAD_cleanTable2$Comparison), size=1) +
  ylim(0,(max(RAD_cleanTable2$Fst)+0.01))  + scale_colour_manual(values = color,name="Comparison") +
  guides(fill=FALSE) + theme(plot.title = element_text(hjust = 0.5))#+ geom_hline(yintercept=max)


LGplot

#Save as very long pdf that is easily readable.
pdf("/Users/sebma/Desktop/Manhattan_RADseq_2sigmas_NOWEIRDSNPS_280224.pdf", 20,10)
LGplot
dev.off()

################  WGseq DATA - WGseq DATA - WGseq DATA - WGseq DATA - WGseq DATA ############################
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

###27/02/24
###Test, remove all SNPs that have a 0 value, because it means the other 2 are biased
WGdf_filtered <- filter(WGpeaks_NConly, fst_SB == fst_LB, fst_SB >0, fst_LB>0)
WGdf_filtered2 <- filter(WGpeaks_NConly, fst_SB == fst_PL, fst_SB >0, fst_PL>0)
WGdf_filtered3 <- filter(WGpeaks_NConly, fst_LB == fst_PL, fst_LB >0, fst_PL>0)

############# Just doing some counting for barplots 26/02/24
sum(WGpeaks_NConly$fst_SB > 0.2) # 2971
sum(WGpeaks_NConly$fst_LB > 0.2) # 2769
sum(WGpeaks_NConly$fst_PL > 0.2) # 3432

sum(WGpeaks_NConly$fst_SB > 0.5) # 19
sum(WGpeaks_NConly$fst_LB > 0.5) # 46
sum(WGpeaks_NConly$fst_PL > 0.5) # 40


#Now doing some barplots
test <- read.csv("/Users/sebma/Desktop/thirdpaperfigures/tables/numberSNPs_aboveThreshold.csv", sep=",")
library(ggplot2)
color <- c("#00BA38FF","#619CFFFF","#F8766DFF")
names(color) <- c("LB","SB","PL")


#Use labs(x= "", y = "Relative expression") for ARL16
#And labs(x = "Developmental stage (ts)", y = "") for Rhoguanin
  graphddRAD <- ggplot(test, aes(x=test$Threshold, y=test$Number.ddSNPs, fill=test$Morph)) + 
    geom_bar(stat="identity", position = position_dodge()) +
    labs(title="ddRAD",x="Threshold",y="Number of SNPs passing threshold", fill="Morph") +
    theme_bw() + theme(plot.title = element_text(hjust=0.5)) + scale_fill_manual(values=color)

  graphddRAD
  
  graphWG <- ggplot(test, aes(x=test$Threshold, y=test$Number.WG.SNP, fill=test$Morph)) + 
    geom_bar(stat="identity", position = position_dodge()) +
    labs(title="WG",x="Threshold",y="Number of windows passing threshold", fill="Morph") +
    theme_bw() + theme(plot.title = element_text(hjust=0.5)) + scale_fill_manual(values=color)
  graphWG
  
########################

#Need to find SNPs which are lying outside 2sigmas of the distribution.
WGpeaks_NConly$highest <- pmax(WGpeaks_NConly$fst_LB,WGpeaks_NConly$fst_SB,WGpeaks_NConly$fst_PL)
maxWG <- mean(WGpeaks_NConly$highest) + 2*sd(WGpeaks_NConly$highest) # max = 0.1989775727

#Separate table per morph, then rbind it.
WG_SB <- data.frame(WGpeaks_NConly[,c(1:2,4)])
WG_LB <- data.frame(WGpeaks_NConly[,c(1:3)])
WG_PL <- data.frame(WGpeaks_NConly[,c(1:2,5)])

#Prepare colname for rbind
colnames(WG_SB)[3] <- "Fst"
colnames(WG_LB)[3] <- "Fst"
colnames(WG_PL)[3] <- "Fst"

#Add the name of the morph to prepare for the color
WG_SB$Comparison <- "SB vs LB/PL"
WG_LB$Comparison <- "LB vs SB/PL"
WG_PL$Comparison <- "PL vs SB/LB"

#Bind everything together
WG_cleanTable <- rbind(WG_SB,WG_LB,WG_PL)

#If Fst less than max, change the Color to "Non-significant"
WG_cleanTable$Comparison <- ifelse(WG_cleanTable$Fst < maxWG, "Non-significant",WG_cleanTable$Comparison)
#Relevel so I have non-significant at the end
WG_cleanTable$Comparison <- factor(WG_cleanTable$Comparison)
WG_cleanTable$Comparison <- relevel(WG_cleanTable$Comparison,"Non-significant")
WG_cleanTable$Comparison <- relevel(WG_cleanTable$Comparison,"PL vs SB/LB")
WG_cleanTable$Comparison <- relevel(WG_cleanTable$Comparison,"SB vs LB/PL")
WG_cleanTable$Comparison <- relevel(WG_cleanTable$Comparison,"LB vs SB/PL")

#Change colnames so it's the same as for RAD
colnames(WG_cleanTable)[1] <- "NCBI"

#I need to calculate the absolute Position of each dot. Because here it is just the basic position.
mergeData <- function(datatable){
  datatable$LG <- ACchrom_noNW_nomit$LG[match(datatable$NCBI,ACchrom_noNW_nomit$NCBI)]
  datatable$START <- ACchrom_noNW_nomit$START[match(datatable$NCBI,ACchrom_noNW_nomit$NCBI)]
  datatable$ABPOS <- datatable$START + datatable$POS
  return(datatable)
}

WG_cleanTable2 <- mergeData(WG_cleanTable)

#Now I have all columns of interest: NCBI, POS, Fst and Color, with Fst and POS numeric.
#Let's try to plot.

color <- c("#00BA38FF","#619CFFFF","#F8766DFF","black")
names(color) <- c("LB vs SB/PL","SB vs LB/PL","PL vs SB/LB","Non-significant")

LGplot_WG <- ggplot()+ 
  geom_rect(aes(NULL,NULL,xmin=starting,xmax=ending,fill=chromosome_info$pattern),
            ymin=-Inf,ymax=Inf,data=chromosome_info) +
  scale_fill_manual(values = alpha(c("white","gray76"), .5)) +
  scale_x_continuous(label=chromosome_info$LG, breaks = chromosome_info$MIDDLE) +
  theme_classic()+ xlab("Linkage group") + ylab("Fst") + ggtitle("Positions of WGseq 'SNPs' on charr genome") +
  theme(axis.text.x = element_text(angle = 90, size = 8, face = "bold", vjust = 0.5),
        axis.title.x = element_text(size=16,
                                    margin = margin(t=15,r=0,b=10,l=0)),
        axis.title.y = element_text(size=16,
                                    margin = margin(t=0,r=10,b=0,l=10)),
        panel.border = element_rect(colour = "black",fill = "NA", size = 1),
        legend.background = element_blank(),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.box.background = element_rect(color = "black")) +
  geom_point(aes(x=WG_cleanTable2$ABPOS, y=WG_cleanTable2$Fst, color=WG_cleanTable2$Comparison), size=1) +
  ylim(0,(max(WG_cleanTable2$Fst)+0.01))  + scale_colour_manual(values = color,name="Comparison") +
  guides(fill=FALSE) + theme(plot.title = element_text(hjust = 0.5))#+ geom_hline(yintercept=max)


LGplot_WG


#Save as very long pdf that is easily readable. Obviously, the WGseq Manhattan is very crowded because we have dots
#Every 50kb on the whole genome.
pdf("/Users/sebma/Desktop/Manhattan_WGseq_2sigmas_201023.pdf", 20,10)
LGplot_WG
dev.off()


### OLD STUFF ###
######## Just doing some counting for barplots 26/02/24
sum(newRADlist_NConly$SB > 0.2) # 723
