library(dplyr)
library(sqldf)
library(ggplot2)

#Aim: make Manhattan plots that shows what happens around the top significant SNPs/regions

########################## LOAD RADseq DATA #################################
        newRADlist <- read.csv("~/vcftools_output_merged.tsv", sep="\t")
#Keep only NC_ (placed scaffolds)
        newRADlist_NConlymit <- newRADlist %>% filter(grepl('NC_', NCBI)) #8480 SNPs
#Remove mitochondrial scaffold
        newRADlist_NConly <- newRADlist_NConlymit %>% filter(!grepl('NC_000861.1', NCBI))
#Identify the highest fst value out of the three morphs
        newRADlist_NConly$highest <- pmax(newRADlist_NConly$SB,newRADlist_NConly$LB,newRADlist_NConly$PL)

#Only keep the top n highest Fst SNPs (for instance n=10):
        highestSNPs <- newRADlist_NConly %>% arrange(desc(highest)) %>% slice(1:10)
#Here, out fo the 10 SNPs, only 5 are in different regions, so keep these five to avoid redundancy (optional)
        highestSNPstop5regions <- highestSNPs[c(1,3,4,6,10),]
#Alternatively, find out the SNPs of interest to you through any way you'd like but keep them in this format

#Rename the scaffold column with a different name so it doesn't confuse the sqldf() function down the line
        colnames(highestSNPstop5regions)[1] <- "NCBIHighest"
#Establish the start and end of the window of interest:
#For instance, if you want to plot a 1Mb region:
          highestSNPstop5regions$startWindow <- highestSNPstop5regions$POS - 500000
          highestSNPstop5regions$endWindow <- highestSNPstop5regions$POS + 500000
#Only keep columns with: name of scaffold, start of window and end of window
          highestSNPstop5regionsCLEAN <- highestSNPstop5regions[,c(1,8:9)]

#Use sqldf to mine data regarding the SNPs that fall in these regions
          All5Regions <- sqldf('SELECT NCBI, POS, ID, SB, LB, PL, highest
                         FROM newRADlist_NConly, highestSNPstop5regionsCLEAN
                         WHERE NCBI = NCBIHighest AND POS > startWindow AND POS < endWindow')

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
          
#If Fst less than 2sigma threshold (0.41497) for RADseq, change the Color to "Non-significant"
          RAD_cleanTable$Comparison <- ifelse(RAD_cleanTable$Fst < 0.41497, "Non-significant",RAD_cleanTable$Comparison)
#Relevel so I have non-significant at the end
          RAD_cleanTable$Comparison <- factor(RAD_cleanTable$Comparison)
          RAD_cleanTable$Comparison <- relevel(RAD_cleanTable$Comparison,"Non-significant")
          RAD_cleanTable$Comparison <- relevel(RAD_cleanTable$Comparison,"PL vs SB/LB")
          RAD_cleanTable$Comparison <- relevel(RAD_cleanTable$Comparison,"SB vs LB/PL")
          RAD_cleanTable$Comparison <- relevel(RAD_cleanTable$Comparison,"LB vs SB/PL")
          
#Separate it in regions:
#ATH: If you decide to plot regions that are on the same scaffold, you will need to find other ways to filter the data
          Region1 <- RAD_cleanTable %>% filter(grepl('NC_036875.1', NCBI))
          Region2 <- RAD_cleanTable %>% filter(grepl('NC_036850.1', NCBI))
          Region3 <- RAD_cleanTable %>% filter(grepl('NC_036841.1', NCBI))
          Region4 <- RAD_cleanTable %>% filter(grepl('NC_036859.1', NCBI))
          Region5 <- RAD_cleanTable %>% filter(grepl('NC_036843.1', NCBI))
          
#Now for the plotting
#I don't need to calculate the absolute position because only 1 LG per graph
        #Old color scheme          
                #color <- c("#00BA38FF","#619CFFFF","#F8766DFF","black")
        #New colour scheme:
                color <- c("#74add1","#313695","#a50026","black")
                names(color) <- c("LB vs SB/PL","SB vs LB/PL","PL vs SB/LB","Non-significant")

#Plotting function
#This finds the position of the highest SNP (the middle of the window) and uses it to calculate the xlims
        graphFct <- function(datatable,LG,title) {
            row_with_highest_Fst_SNP <- datatable %>% filter(Fst == max(Fst))
            xlim1 <- row_with_highest_Fst_SNP$POS - 500000
            xlim2 <- row_with_highest_Fst_SNP$POS + 500000
            ggplot()+ theme_classic()+ xlab(LG) + ylab("Fst") + ggtitle(title) +
              geom_point(aes(x=datatable$POS, y=datatable$Fst, color=datatable$Comparison), size=1) + 
              theme(axis.text.x = element_text(angle = 0, size = 8, vjust=0.6)) + xlim(xlim1,xlim2) +
              ylim(0,0.9) + scale_colour_manual(values = color,name="Comparison") + theme(plot.title = element_text(hjust = 0.5))
          }
          
          graph1 <- graphFct(Region1,"NC_036875.1","Region1")
          graph2 <- graphFct(Region2,"NC_036850.1","Region2")
          graph3 <- graphFct(Region3,"NC_036841.1","Region3")
          graph4 <- graphFct(Region4,"NC_036859.1","Region4")
          graph5 <- graphFct(Region5,"NC_036843.1","Region5")
          
          graph1
          graph2
          graph3
          graph4
          graph5

#As you can see, there are multiple peaks in region1, so you can change the function to have a smaller window around the highest SNP:
#For instance, only a 100bp window
        graphFct_100bp <- function(datatable,LG,title) {
                    row_with_highest_Fst_SNP <- datatable %>% filter(Fst == max(Fst))
                    xlim1 <- row_with_highest_Fst_SNP$POS - 50
                    xlim2 <- row_with_highest_Fst_SNP$POS + 50
                    ggplot()+ theme_classic()+ xlab(LG) + ylab("Fst") + ggtitle(title) +
                      geom_point(aes(x=datatable$POS, y=datatable$Fst, color=datatable$Comparison), size=1) + 
                      theme(axis.text.x = element_text(angle = 0, size = 8, vjust=0.6)) + xlim(xlim1,xlim2) +
                      ylim(0,0.9) + scale_colour_manual(values = color,name="Comparison") + theme(plot.title = element_text(hjust = 0.5))
                  }
        graph1_bis <- graphFct_100bp(Region1,"NC_036875.1","Region1")

#Or you could also specify the xlim() directly if it is a one time thing or if the size that makes the most sense changes all the time 
       
#Save your plots in .pdf if needed
        pdf("/Users/sebma/Desktop/Top5Regions.pdf", 5,5)
                  graph1          
                  graph1bis
                  graph2
                  graph3
                  graph4
                  graph5
        dev.off()


########################## LOAD WGseq DATA #################################










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
