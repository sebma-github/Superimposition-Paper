#### PACKAGES ####
#The plan is to make a Manhattan plot that shows the significant SNPs/regions
#For both RADseq and WGseq
#But this time, only about 1MB around the top 5
library(dplyr)
library(sqldf)
library(ggplot2)

#Identify the top 5 highest SNPs first:

        #For RADseq
          newRADlist <- read.csv("/Users/sebma/Desktop/SNP_SLH-main/SNP_SLH-main/methylation/vcftools_output_merged.tsv", sep="\t")
          #Keep only NC_ and no mit
          newRADlist_NConlymit <- newRADlist %>% filter(grepl('NC_', NCBI)) #8480 SNPs
          newRADlist_NConly <- newRADlist_NConlymit %>% filter(!grepl('NC_000861.1', NCBI))
          newRADlist_NConly$highest <- pmax(newRADlist_NConly$SB,newRADlist_NConly$LB,newRADlist_NConly$PL)
          #Only keep the 5 highest rows
          highestSNPs <- newRADlist_NConly %>% arrange(desc(highest)) %>% slice(1:5)                                     # Top N highest values by group
          #This works, but you can see that 3 of the 5 highest SNPs are in the same region
          #Keep the highest SNPs until we get the top 5 regions
          highestSNPs <- newRADlist_NConly %>% arrange(desc(highest)) %>% slice(1:10)
          highestSNPstop5regions <- highestSNPs[c(1,3,4,6,10),]
          #That will give me the position of the 5 regions I want to look at
          #They asked me to look at this on a 1MB basis
          highestSNPstop5regions$NCBIHighest <- highestSNPstop5regions$NCBI
          highestSNPstop5regions$startWindow <- highestSNPstop5regions$POS - 500000
          highestSNPstop5regions$endWindow <- highestSNPstop5regions$POS + 500000
          #Remove all columns that have the same name
          highestSNPstop5regionsCLEAN <- highestSNPstop5regions[,c(8:10)]
          #Now that I have the coordinates, I need to keep all positions in these regions
          All5Regions <- sqldf('SELECT NCBI, POS, ID, SB, LB, PL, highest
                         FROM newRADlist_NConly, highestSNPstop5regionsCLEAN
                         WHERE NCBI = NCBIHighest AND POS > startWindow AND POS < endWindow')
          #Put the name of Morph for colouring
          #Separate table per morph, then rbind it.
          RAD_SB <- data.frame(All5Regions[,c(1:4)])
          RAD_LB <- data.frame(All5Regions[,c(1:3,5)])
          RAD_PL <- data.frame(All5Regions[,c(1:3,6)])
          
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
          #If Fst less than 2sigma threshold (0.41497) for RADseq, change the Color by "Non-significant"
          RAD_cleanTable$Comparison <- ifelse(RAD_cleanTable$Fst < 0.41497, "Non-significant",RAD_cleanTable$Comparison)
          #Relevel so I have non-significant at the end
          RAD_cleanTable$Comparison <- factor(RAD_cleanTable$Comparison)
          RAD_cleanTable$Comparison <- relevel(RAD_cleanTable$Comparison,"Non-significant")
          RAD_cleanTable$Comparison <- relevel(RAD_cleanTable$Comparison,"PL vs SB/LB")
          RAD_cleanTable$Comparison <- relevel(RAD_cleanTable$Comparison,"SB vs LB/PL")
          RAD_cleanTable$Comparison <- relevel(RAD_cleanTable$Comparison,"LB vs SB/PL")
          
          #Separate it in regions:
          Region1 <- RAD_cleanTable %>% filter(grepl('NC_036875.1', NCBI))
          Region2 <- RAD_cleanTable %>% filter(grepl('NC_036850.1', NCBI))
          Region3 <- RAD_cleanTable %>% filter(grepl('NC_036841.1', NCBI))
          Region4 <- RAD_cleanTable %>% filter(grepl('NC_036859.1', NCBI))
          Region5 <- RAD_cleanTable %>% filter(grepl('NC_036843.1', NCBI))
          
          #Now for the plotting
          #I don't need to calculate the absolute position because only 1 LG per graph
          color <- c("#00BA38FF","#619CFFFF","#F8766DFF","black")
          names(color) <- c("LB vs SB/PL","SB vs LB/PL","PL vs SB/LB","Non-significant")
          
          graphFct <- function(datatable,LG,startwindow,endwindow,title) {
            ggplot()+ theme_classic()+ xlab(LG) + ylab("Fst") + ggtitle(title) +
            geom_point(aes(x=datatable$POS, y=datatable$Fst, color=datatable$Comparison), size=1) + 
            theme(axis.text.x = element_text(angle = 0, size = 8, vjust=0.6)) + xlim(startwindow,endwindow) +
            ylim(0,0.9) + scale_colour_manual(values = color,name="Comparison") + theme(plot.title = element_text(hjust = 0.5))
          }
          
          graph1 <- graphFct(Region1,"NC_036875.1",29708000,30708000,"Region1")
          graph2 <- graphFct(Region2,"NC_036850.1",20926835,21926835,"Region2")
          graph3 <- graphFct(Region3,"NC_036841.1",9230820,10230820,"Region3")
          graph4 <- graphFct(Region4,"NC_036859.1",30994438,31994438,"Region4")
          graph5 <- graphFct(Region5,"NC_036843.1",21567680,22567680,"Region5")
          
          graph1
          graph2
          graph3
          graph4
          graph5
          
          #For Region1 and Region2, need to make a graph with a smaller window:
          #For Region 1: 100bp is enough
            graph1bis <- ggplot()+ theme_classic()+ xlab("NC_036875.1") + ylab("Fst") + ggtitle("Region 1 zoomed in") +
            geom_point(aes(x=Region1$POS, y=Region1$Fst, color=Region1$Comparison), size=1) + xlim(30207950,30208050) +
            ylim(0,0.9) + theme(axis.text.x = element_text(angle = 0, size = 8, vjust=0.6)) + 
            scale_colour_manual(values = color,name="Comparison") + theme(plot.title = element_text(hjust = 0.5))
          graph1bis
          
          graph2bis <- ggplot()+ theme_classic()+ xlab("NC_036850.1") + ylab("Fst") + ggtitle("Region 2 zoomed in") +
            geom_point(aes(x=Region2$POS, y=Region2$Fst, color=Region2$Comparison), size=1) + xlim(21426785,21426885) +
            ylim(0,0.9) + theme(axis.text.x = element_text(angle = 0, size = 8, vjust=0.6)) +
            scale_colour_manual(values = color,name="Comparison") + theme(plot.title = element_text(hjust = 0.5))
          graph2bis
          #For Region1 and Region2, need to make a graph with a smaller window: lets say 50kb before and 50kb after
          
          #ATH: ISSUE WITH THE ddRAD SNPs: Some SNPs only have values for 2 morphs: Thus impossible to know which morph separates
          #Check how many for each pair of morphs:
          LBSB <- filter(newRADlist_NConly, LB == SB, LB >0, SB>0)
          LBPL <- filter(newRADlist_NConly, LB == PL, LB >0, PL>0)
          SBPL <- filter(newRADlist_NConly, SB == PL, SB >0, PL>0)
          
pdf("/Users/sebma/Desktop/Top5Regions.pdf", 5,5)
          graph1          
          graph1bis
          graph2
          graph2bis
          graph3
          graph4
          graph5
          dev.off()

pdf("/Users/sebma/Desktop/Region2_zoomedin.pdf", 5,5)   
graph2bis
dev.off()

##############################
        #For WGseq
        WGpeaks <- read.table("/Users/sebma/Desktop/RAD_WG_Meth/fst_window_WG.tsv", sep=",", header = TRUE)
        #Only keep the data on NC_ and no mit
        WGpeaks_NConlymit <- WGpeaks %>% filter(grepl('NC_', CHROM)) #(151540 obs)
        WGpeaks_NConly <- WGpeaks_NConlymit %>% filter(!grepl('NC_000861.1', CHROM)) #(151540 obs)
        ########ATH I realize there are some negative Fst values. 
        #Put negative values at 0.
        WGpeaks_NConly[WGpeaks_NConly < 0] <- 0
        WGpeaks_NConly$highest <- pmax(WGpeaks_NConly$fst_LB,WGpeaks_NConly$fst_SB,WGpeaks_NConly$fst_PL)
        #Only keep the 5 highest rows
        highestSNPsWG <- WGpeaks_NConly %>% arrange(desc(highest)) %>% slice(1:15)                                   
        #Keep the highest SNPs until we get the top 5 regions
        highestSNPsWGtop5regions <- highestSNPsWG[c(1,2,6,9,12),]
        #That will give me the position of the 5 regions I want to look at
        #They asked me to look at this on a 1MB basis
        highestSNPsWGtop5regions$CHROM2 <- highestSNPsWGtop5regions$CHROM
        highestSNPsWGtop5regions$startWindow <- highestSNPsWGtop5regions$POS - 500000
        highestSNPsWGtop5regions$endWindow <- highestSNPsWGtop5regions$POS + 500000
        
        
        #Remove all columns that have the same name
        highestSNPsWGtop5regionsCLEAN <- highestSNPsWGtop5regions[,c(7:9)]
        #Now that I have the coordinates, I need to keep all positions in these regions
        All5RegionsWG <- sqldf('SELECT CHROM, POS, fst_SB, fst_LB, fst_PL, highest
                         FROM WGpeaks_NConly, highestSNPsWGtop5regionsCLEAN
                         WHERE CHROM = CHROM2 AND POS > startWindow AND POS < endWindow')
        #Put the name of Morph for colouring
        #Separate table per morph, then rbind it.
        WG_SB <- data.frame(All5RegionsWG[,c(1:3)])
        WG_LB <- data.frame(All5RegionsWG[,c(1:2,4)])
        WG_PL <- data.frame(All5RegionsWG[,c(1:2,5)])
        
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
        #If Fst less than 2sigma threshold (0.1989775727) for WGseq, change the Color by "Non-significant"
        WG_cleanTable$Comparison <- ifelse(WG_cleanTable$Fst < 0.1989775727, "Non-significant",WG_cleanTable$Comparison)

        #Relevel so I have non-significant at the end
        WG_cleanTable$Comparison <- factor(WG_cleanTable$Comparison)
        WG_cleanTable$Comparison <- relevel(WG_cleanTable$Comparison,"Non-significant")
        WG_cleanTable$Comparison <- relevel(WG_cleanTable$Comparison,"PL vs SB/LB")
        WG_cleanTable$Comparison <- relevel(WG_cleanTable$Comparison,"SB vs LB/PL")
        WG_cleanTable$Comparison <- relevel(WG_cleanTable$Comparison,"LB vs SB/PL")
        
        #Separate it in regions:
        Region1 <- WG_cleanTable %>% filter(grepl('NC_036848.1', CHROM))
        Region2 <- WG_cleanTable %>% filter(grepl('NC_036854.1', CHROM))
        Region3 <- WG_cleanTable %>% filter(grepl('NC_036847.1', CHROM))
        Region4 <- WG_cleanTable %>% filter(grepl('NC_036849.1', CHROM))
        Region5 <- WG_cleanTable %>% filter(grepl('NC_036860.1', CHROM))
        
        #Now for the plotting
        #I don't need to calculate the absolute position because only 1 LG per graph
        color <- c("#00BA38FF","#619CFFFF","#F8766DFF","black")
        names(color) <- c("LB vs SB/PL","SB vs LB/PL","PL vs SB/LB","Non-significant")
        
        graphFct <- function(datatable,LG,startwindow,endwindow,title) {
          ggplot()+ theme_classic()+ xlab(LG) + ylab("Fst") + ggtitle(title) +
            geom_point(aes(x=datatable$POS, y=datatable$Fst, color=datatable$Comparison), size=1) + 
            theme(axis.text.x = element_text(angle = 0, size = 8, vjust=0.6)) + xlim(startwindow,endwindow) +
            ylim(0,0.9) + scale_colour_manual(values = color,name="Comparison") + theme(plot.title = element_text(hjust = 0.5))
        }
        
        graph1 <- graphFct(Region1,"NC_036848.1",54390000,55390000,"Region1")
        graph2 <- graphFct(Region2,"NC_036854.1",48960000,49960000,"Region2")
        graph3 <- graphFct(Region3,"NC_036847.1",11810000,12810000,"Region3")
        graph4 <- graphFct(Region4,"NC_036849.1",14230000,15230000,"Region4")
        graph5 <- graphFct(Region5,"NC_036860.1",19950000,20950000,"Region5")
        
        graph1
        graph2
        graph3
        graph4
        graph5
        
        pdf("/Users/sebma/Desktop/Top5Regions_WG.pdf", 5,5)
        graph1          
        graph2
        graph3
        graph4
        graph5
        dev.off()
        
        
        
        
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
