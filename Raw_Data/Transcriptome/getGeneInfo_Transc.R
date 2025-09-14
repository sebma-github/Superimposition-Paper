## Get the genomic information (Scaffold and position) of genes in the transcriptome data
#Based on gene ID and using the assembly
library(sqldf)

#ATH there is a list of all exons in /data1/RNA_charr/ref/exon_id_full
#That Alexander used to do the htseq counts with htseq.sh by using the "Aligned_filtered" files in /star.
#From there, he used the DEseq script to calculate the DEgenes

#OK so basically, I think the files that Alexander provided us with might be good. 
#BUT, there was no LOC number on the allgenes file.
#So I had to try to match the gene ID to a LOC number and that's where things get dicey because I do not know which
#Dataset he used for that (annotation)
#But actually if I use the files I re-generated with his scripts, I have LOC numbers everywhere.
#So because of that I should be able to compare them.
#However, I will need to use a way to get the start and end positions of these genes..

#For that, I might have to use StortData3 and separate the LOC number from the rest of the V9 column.
#Otherwise the canada_genome_protein has a bunch of LOC as well. I should try to see if there are some that are 
#not in this database
#Question 1:
#Are there geneIDs that are in the LRT results/all genes that are not in the StortData3?
#Are there LOCIDs that are in the LRT results/all genes that are not in the StortData3?
#Which ones, and why. And can I find them in the canada_genome_dbx_exon stuff



#Load the dataframes of interest.
    #All transcripts
        #Alexander version
        allgenes <- read.table("/Users/sebma/Desktop/Alexander stuff/transcriptomeData/expr_fullref.tsv", sep="\t", header=T)
        colnames(allgenes)[1] <- "geneID"
        
        #Regenerated version
        allgenes <- read.table("/Users/sebma/Desktop/Alexander stuff/DESeq_allgenes_LRT_SM_140723.tsv", sep="\t", header=T)
        colnames(allgenes)[1] <- "geneID"

#This is the same dataset
# test <- read.table("/Users/sebma/Desktop/Alexander stuff/WGCNA/fullreference/exprLogNorm_full_m160.csv", sep=",", header=T)

    #Only DE genes
        #Alexander version
        DEgenes <- read.table("/Users/sebma/Desktop/Alexander stuff/transcriptomeData/DESeq_m160_results_LRT.tsv", sep="\t", header=T)
        #First problem I see: There is no geneID here but the full LOC number, and in some cases the real name
        #Actually, shouldn't be a problem because that info is also available in the file "canada_genome_overview_protein.tsv"

        #Regenerated version
        DEgenes <- read.table("/Users/sebma/Desktop/Alexander stuff/DESeq_signif_results_LRT_SM_140723.tsv", sep="\t", header=T)


        
geneLOC <- DEgenes$gene_name
test <- stortData[stortData$gene %in% geneLOC,] 
#As you can see, now we have more genes but it is actually because of duplicates. 
#If you remove the duplicates, then there is less. Meaning that some of the genes
#that are in the DE are not in the annotation table.
geneLOCstort <- stortData$gene
test2 <- DEgenes[!DEgenes$gene_name %in% geneLOCstort,] #356 LOC are not in StortData
test3 <- allgenes[!allgenes$gene_name %in% geneLOCstort,] #9717 out of 46120 LOC are not in stortData

#But are these in StortData3?
test4 <- grepl('LOC112069209',stortData3$V9)
#It is!! And on an NW_ scaffold. 
#What does this tell me? is canada_genome_overview_protein only on placed?
#It would be good to download the datasets from the NCBI genome again and check there as well.
#I might need to get info from StortData3


#I could just remove all the genes that are not in the annotation db for both datasets and go from there.
#But that is kind of weird.


        

#Load the information about the genes: Here I probably have different methods
stortData <- read.table(file='/Users/sebma/Desktop/annotated_charr_genome/canada_genome_overview_protein.tsv', sep = '\t', header=TRUE)

stortData2 <- read.table(file='/Users/sebma/Desktop/annotated_charr_genome/canada_genome_overview_mrna.tsv', sep = '\t', header=TRUE)

stortData3 <- read.table(file='/Users/sebma/Desktop/Alexander stuff/alexanderstuff for DE/canada_exons_dbxref.gff', sep = '\t', header=F)

stortData4 <- read.table(file='/Users/sebma/Desktop/Alexander stuff/alexanderstuff for DE/exon_id_full.gtf', sep = '\t', header=F)

stortData5 <- read.table(file='/Users/sebma/Desktop/GCF_002910315.2/ncbi_dataset/data/GCF_002910315.2/genomic.gff', sep = '\t', header=F)

#This is the latest dl from Salmobase ==> its exactly the same as the one from NCBI
stortData6 <- read.table(file='/Users/sebma/Desktop/annotated_charr_genome/GCF_002910315.2_ASM291031v2_genomic.gff', sep = '\t', header=F)


#Maybe I would have more chances with that.
#Need to split the column x into multiple columns

test <- separate(stortData3, "V9", c("ID","Parent","GeneIDtemp","Key","LOC","Product","TranscriptID"), sep=";", remove=F)
#This would work if the V9 column was always the same but it's not. Sometimes there are additionnal warning
#messages in the column that get split as well and it's a mess.

#So I need to do this better.
library(stringr)
stortData3$geneID <- str_extract(stortData3$V9, "(?<=GeneID:)[^;]*(?=,|$)")
#Works when there is a comma after GeneID:xxx but sometimes it is a semi colon......
sum(is.na(stortData3$geneID)) #23383

#Test with the LOC instead
stortData3$LOCID <- str_extract(stortData3$V9, "(?<=gene=)[^;]*(?=;|$)")
#Works when there is a comma after GeneID:xxx but sometimes it is a semi colon......
sum(is.na(stortData3$LOCID)) #1463 better.

#But OK let's assume that I can identify each LOC to an exon position, then it is kind of awkward to take the 
#Exon positions as gene positions.
#That is why the transcript lists were good right?
#Because it gives me a position for the actual mRNA or gene. 
#Maybe if I take the full genomic.gff file and REMOVE the exons and only keep the CDS/mRNA?
#And then check if I can still find all transcripts there.

stortData5$LOCID <- str_extract(stortData5$V9, "(?<=gene=)[^;]*(?=;|$)")
sum(is.na(stortData5$LOCID)) #31542
#Remove all lines where V3 is exon
stortData5_noexon <- stortData5[!grepl("exon", stortData5$V3), ]
stortData5_noexon_noCDS <- stortData5_noexon[!grepl("CDS", stortData5_noexon$V3), ]

#That would be a good df. Now look into that if I can find the LOCs that were missing before, i.e. LOC112069209
#It is there.
#Now look if I am still missing some genes from the DEseq or allgenes datasets

geneLOC <- DEgenes$gene_name
test <- stortData5_noexon_noCDS[stortData5_noexon_noCDS$LOCID %in% geneLOC,] 
#As you can see, now we have more genes but it is actually because of duplicates. 
#If you remove the duplicates, then there is less. Meaning that some of the genes
#that are in the DE are not in the annotation table.
geneLOCstort <- stortData5_noexon_noCDS$LOCID
geneLOCstort_withexon <- stortData5$LOCID 
test2 <- DEgenes[!DEgenes$gene_name %in% geneLOCstort,] #908 LOC are not in stortData5_noexon_noCDS
test2bis <- DEgenes[!DEgenes$gene_name %in% geneLOCstort_withexon,] #906 are not in stortData5
test3 <- allgenes[!allgenes$geneID %in% geneLOCstort,] #22504 are not in stortData5_noexon_noCDS

#Honestly, even this annotation looks shady as fuck. 
#Look at row 317855 in stortData5 (LOC111979690 is apparently somewhere in V9)
#I might end up having to say: the annotation is weird, we only took genes where the rows in the gff file looked ok.
#But then, it might be best to rerun the htseq.sh with a file that doesn't have these weird columns.

#Also, how many of the LOCIDs are na in this datasets
sum(is.na(stortData5_noexon_noCDS$LOCID)) #31525
#Check what these are
NA_LOC <- stortData5_noexon_noCDS[is.na(stortData5_noexon_noCDS$LOCID),]
library(dplyr)
count(NA_LOC, V3) #(1 D_loop, 27107 cDNA_match, 1 origin of replication, 2rRNA, 4399 regions, 15 tRNA)




#What about the exon_id_full.gtf. Are all DEgenes and allgenes in this one?
stortData4$LOCID <- str_extract(stortData4$V9, "(?<=gene_name )[^;]*(?=;|$)").
sum(is.na(stortData4$LOCID)) #0 no NAs in this dataset

geneLOCstort4 <- stortData4$LOCID

test2 <- DEgenes[!DEgenes$gene_name %in% geneLOCstort,] #0 LOC are not in StortData4
test3 <- allgenes[!allgenes$geneID %in% geneLOCstort4,] #0 are not in StortData4


#That really makes me think that he used the exon_id_full.gtf file to make the htseq counts.
#But how exactly did he make this file is a mystery to me.
#And it has the 

#Should rerun this script only with stortData4














for (i in 1:nrow(stortData3)) {
  if (is.na(stortData3$geneID[i])) 
      {
    stortData3$geneID[i] <- str_extract(stortData3$V9[i], "(?<=GeneID:)[^;]*(?=;|$)")
  }
}
  
sum(is.na(stortData3$geneID)) #1463.... WTF. 

test <- stortData3[is.na(stortData3$geneID),]
factor(stortData3$V3) #There is a bunch of them that is not exons.




is.na(stortData3$geneID[1])



a <- "DP=26;AN=2;DB=1;AC=1;MQ=56;MZ=0;ST=5:10,7:2;CQ=SYNONYMOUS_CODING;GN=NOC2L;PA=1^1:0.720&2^1:0"
str_extract(a, "(?<=GN=)[^;]*(?=;|$)")






#Keep the information for the genes of interest (no need to keep the transcript id I think)
#Note, there is more rows than for allgenes because multiple transcripts in the ref while only 
# 1 geneID in expr_fullref.tsv
allgenesFullInfo <- sqldf('SELECT geneID, transcript_id, chromosome_RefSeq, start, stop
                 FROM stortData, allgenes 
                 WHERE geneID=gene_id')

#ATH, remove duplicates. When I remove the duplicates, I do not end back up with 46120... Which means some of the gene IDs are
#lost when I compare them to the big table. So the canada_genome_overview_protein.tsv file might not be the one.
#I should probably use another one..
allgenesfullinfo_nodup <- allgenesFullInfo[!duplicated(allgenesFullInfo$geneID),]

# write.csv(allgenesFullInfo, "/Users/sebma/Desktop/Alexander stuff/allgeneswithPosition.csv")

#I will need to decide what I do with the different transcripts: i.e. do I count them as one (i think so) or do I separate them
#Separating them would lead to counting the gene more than once? Not sure.

DEgenesFullInfo <- sqldf('SELECT gene_name, gene_id, chromosome_RefSeq, start, stop, direction, direction_unique, log2FoldChange_LB_vs_PL, log2FoldChange_LB_vs_SB, log2FoldChange_PL_vs_SB, pvalue, padj
                 FROM stortData, DEgenes 
                 WHERE gene_name=gene')

DEgenesFullInfo

#write.csv(DEgenesFullInfo, "/Users/sebma/Desktop/Alexander stuff/DEgeneswithPosition.csv")
