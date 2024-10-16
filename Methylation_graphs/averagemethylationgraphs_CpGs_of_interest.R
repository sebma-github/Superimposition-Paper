#Average methylation graphs for CpGs that are close to SNPs
library(stringr)
library(dplyr)
library(ggplot2)

############################### LOAD DATA #################################
#Load the methylation data:
meth.min1df <- read.table("/Users/sebma/Desktop/27samples/methmin1_27_noPIno100.csv", sep=",", header=TRUE)

  ##Use the annotations to create an ID we will use later
    annotation <- meth.min1df[,c(1,2,3,4)]
    annotation$ID <- paste(annotation$chr, annotation$start, sep="_")

  #remove the first columns that I don't want.
    meth.min1df <- meth.min1df[,-c(1,2,3,4)] #144 #81
  #Remove the number of Ts columns (every 3 rows)
    meth.min1df <- meth.min1df[,-seq(0,81,3)] #96 #54

  #Separate the data in two: on one side the number of counted cytosines, on the other the total of bases counted (C + T)
    numCdf <- meth.min1df[,seq(0,54,2)]
    allbasesdf <- meth.min1df[,seq(1,54,2)]
  
  #Calculate the methylation percentage for each position by dividing the number of cytosines by the number of total bases.
    total <- ((numCdf*100)/allbasesdf)

  #Put back the annotation in terms of row name. 
    rownames(total) <- annotation$ID  
    
#Now the "total" df will have the scaffold and position of each SNP as row.name
#and each of the 27 samples as columns (order of samples and sample names can be found a bit below)

#Now, which CpGs you are interested in?
    #For instance, are you interested in CpGs because of their position around SNPs?
    #Here is an example of CpGs in a 30k peak around RAD SNPs with an Fst > 0.2:
    interestPOS <- read.csv("/Users/sebma/Desktop/27samples/Position_CpGvsRADSNP_0p2_peak30k.csv")
        #Again, this is just an example, I encourage you to make your own list of CpGs that are of
        #interest to you. Either with other scripts by comparing with other datasets,
        #Or you could also just note the position of regions of genetic divergence and manually make
        #a list of a range of positions you'd be interested in
    
    #Reformat the IDs so that it matches the row.names in the total df
      #replace "NC_" by "chr"
        interestPOS$NCBI_MET <- str_replace(interestPOS$NCBI_MET, "NC_","chr")
      #paste together scaffold and position
        interestPOS$ID <- paste(interestPOS$NCBI_MET, interestPOS$POS_MET, sep="_")
      #In this example, I had the same CpG registering twice because it appears in two peaks around two SNPs
      # So remove duplicates (this might not be necessary in your case)
        POSofINTEREST <- unique(interestPOS$ID)

    #Filter the data for these CpGs of interest in the total df
    CpGofInterest <- total %>% filter(grepl((paste(POSofINTEREST,collapse="|")), row.names(total)))

    #Transpose the data so that it is in the right way.
    totalt <- as.data.frame(t(CpGofInterest))

    #Add a column for morph and timepoint.
    #Sample names and order can be found at "/Users/sebma/Desktop/27samples/samplerecapdf_samplenames_27.csv"
    #But you can also trust these vectors and reuse them (as long as you work with the 27 samples methylation data, and not the 48 samples methylation data from my first paper)
    totalt$Morph <- c("PL","PL","PL","LB","LB","LB","SB","SB","SB",
                  "PL","PL","LB","LB","SB","PL","SB","SB","LB",
                  "PL","PL","LB","LB","SB","PL","SB","SB","LB")
    totalt$Timepoint <- c("200ts","200ts","200ts","200ts","200ts","200ts","200ts","200ts","200ts",
                      "150ts","150ts","150ts","150ts","150ts","150ts","150ts","150ts","150ts",
                      "50ts","50ts","50ts","50ts","50ts","50ts","50ts","50ts","50ts")


#Now we have all the data we need.
  #(ATH: At this point in the script, I have hard-coded a bunch of things, which is not ideal)
  #Separate the CpGs based on their location in the genome, to make a graph per location
  #For example: let's just do LG4q and LG5 for now
    #LG4:
      LG4q.1_29 <- data.frame(totalt$chr036842.1_62932148,totalt$chr036842.1_62932283,totalt$chr036842.1_62932301,
                              totalt$Morph,totalt$Timepoint)
      colnames(LG4q.1_29) <- c("Meth1","Meth2","Meth3","Morph","Timepoint")
        #Note: These 3 methylation points are three DIFFERENT CpGs, nothing to do with biological replicates.
      #add a mean of these 3 CpGs in the region
        LG4q.1_29$Meth <- (LG4q.1_29$Meth1+LG4q.1_29$Meth2+LG4q.1_29$Meth3)/3
      #add a column with both morph and time data in the same string
        LG4q.1_29$MT <- paste(LG4q.1_29$Morph,LG4q.1_29$Timepoint,sep="")
      
    #LG5: 
      LG5 <- data.frame(totalt$chr036844.1_16136638,totalt$Morph,totalt$Timepoint)
      colnames(LG5) <- c("Meth","Morph","Timepoint")
        #Note: only one CpG here, no mean calculation
      #add a column with both morph and time data in the same string
        LG5$MT <- paste(LG5$Morph,LG5$Timepoint,sep="")
        
#To simplify/automate the process, you could 1) make one graph per CpG, by taking each methylation
#column of totalt and bind it with the column morph and Timepoint.
#You could also 2) extract all columns with a name that start with a specific scaffold string,
#then add the morph and timepoint data and do a big Mean on the whole scaffold.
#3) you just do as I did and hard code the specific CpGs that you want together.
      

############################ PLOTTING #############################
  #old color scheme:
    #colors <- c("LB" = "#00BA38FF", "SB" = "#619CFFFF", "PL" = "#F8766DFF")
  #new color scheme:
    colors <- c("LB" = "#74add1", "SB" = "#313695", "PL" = "#a50026") 

#Make a function to calculate stats:
      calculateStats <- function(df) {
        #calculate Mean and sd for each morph*time combination
          newdf <- data.frame("Mean"= tapply(df$Meth, df$MT, mean, na.rm=TRUE), "Sd"=tapply(df$Meth, df$MT, sd)) 
        #keep the first two characters of row names as morph (i.e. "PL", "SB", etc...)
          newdf$Morph <-  substr(row.names(newdf),1,2)
        #keep the last three to five characters for as time (i.e. "200ts", "150ts"...)
          newdf$Time <- substr(row.names(newdf),3,7)
        
        #Relevel the time so it is in order in the legend when we plot
          newdf$Time <- factor(newdf$Time)
          newdf$Time <- relevel(newdf$Time, "200ts")
          newdf$Time <- relevel(newdf$Time, "150ts")
          newdf$Time <- relevel(newdf$Time, "50ts")
        return(newdf)
      }
      
      LG4q.1_29_stats <- calculateStats(LG4q.1_29)
      LG5_stats <- calculateStats(LG5)

#Make a function for plotting:
      plotting <- function(df, name) {
          
          plot <- ggplot(df, aes(x=Time, y=Mean, color=Morph)) +geom_point(position=position_dodge(width=0.5)) + 
            geom_errorbar(aes(ymin=Mean-Sd, ymax=Mean+Sd), width=.2,position=position_dodge(.5)) +
            labs(x = "Developmental stage (ts)", y = "Methylation %age", color = "Morph", title = name) + 
            scale_color_manual(values = colors) + theme_bw() + theme(plot.title = element_text(hjust = 0.5))
          return(plot)
        }
#Make graphs
LG4q.1_29_plot <- plotting(LG4q.1_29_stats, "CpGs on LG4q.1_29")
LG5_plot <- plotting(LG5_stats, "CpG on LG5")

#Display graphs
LG4q.1_29_plot
LG5_plot



