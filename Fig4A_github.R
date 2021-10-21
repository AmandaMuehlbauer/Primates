#Model 6 
###This script makes Fig. 4A 

library("hexbin")

library("DESeq2")

library('ggplot2')
library("pheatmap")
library("RColorBrewer")
library("PoiClaClu")
library("dplyr")
library("vsn")
library("AnnotationDbi")
library("org.Hs.eg.db")
source("https://bioconductor.org/biocLite.R")
biocLite("goseq")
library("lattice")
library("hash")
library("biomaRt")
library("gplots")

library("UpSetR")
library("tidyr")

library("compositions")
library(scales)
###I did not change any of the variable names, so primate_genera is kind of a misnomer
##Should be "primate_species"

setwd("/Volumes/GoogleDrive/My\ Drive/Primates_Coculture/Sequencing20180530/ShotgunSequencing/Tables_Combined")
primate_genera <- t(read.table("Merged_Species.txt", sep = "\t", header = TRUE, row.names = 1, check.names = FALSE))

setwd("/Volumes/GoogleDrive/My\ Drive/Primates_Coculture/Primates/MyAnalysis")

primate_meta16s <- read.table("primate_metagenomic_meta.txt", sep = "\t", header = TRUE, row.names ="SampleID")

#deal with metadata formatting first
primate_meta16s$sampname <- rownames(primate_meta16s)
#Get rid of UMGC blank
primate_meta16s <- primate_meta16s[!(primate_meta16s$Microbiome == "BLANK"),]

#deal with sequencing data formatting
newrows <- vapply(strsplit(rownames(primate_genera),"_"), `[`, 1, FUN.VALUE=character(1))
rownames(primate_genera) <- newrows
#Rename the luca samples
rownames(primate_genera)[16] <- "1UNCULTURED"
rownames(primate_genera)[17] <- "65UNCULTURED"
rownames(primate_genera)[18] <- "110UNCULTURED"
rownames(primate_genera)[19] <- "111UNCULTURED"

#Subset the genera so it is between 50-100 using column sums.
gen1 <-primate_genera


#Subset abundance data so that I have only the uncultured samples 
#Get rid of low abundance genera
gen1 <-as.data.frame(gen1)
gen1.1<- gen1[grepl("UNCULTURED",rownames(gen1)),]
#I want to keep genera that are represented in at least half the samples (8)
###Try running for every microbial species detected
cutoff_nsamples <- colSums(gen1.1>0) >9
gen1.11<- gen1.1[,cutoff_nsamples]
#gen1.11 <- gen1.1

#get select genera by abundance per species
#After filtering on number of species with the taxa, this step is superfluous
#gen1.2 <-gen1.11[,colSums(gen1.11)>0.2]
gen1.2 <-gen1.11
gen1.2[gen1.2 ==0] <- 0.0065
gen1.3 <-clr(gen1.2)
gen1.3 <-as.data.frame(gen1.3)
gen1.3$MBID <- rownames(gen1.3)

##Also need to subset the metadata so that it I only have the "uncultured" samples
metagen <- as.data.frame(primate_meta16s)
metagen1 <- metagen[grepl("UNCULTURED", rownames(metagen)),]
#Get rid of  microbiome samples that were not in the colonocyte experiment
metagen1.2 <- metagen1[rownames(metagen1) %in% rownames(gen1.3),]
metagen1.2$MBID <- rownames(metagen1.2)

gen2 <-merge(gen1.3, metagen1.2, by = "MBID", all.y = TRUE)
colnames(gen2)

#These indicies are hardcoded and need to be changed if code is updated!!!!!
gen2[, 2:38][is.na(gen2[, 2:38])] <- 0

gen3 <-gen2[,!grepl("Unassigned", colnames(gen2))]


#MAke a new directory for the results
 systempath <- "/Volumes/GoogleDrive/My\ Drive/Primates_Coculture/Primates/MyAnalysis_v02"
 setwd(systempath)
 dir.create(file.path(systempath, format(Sys.time(), "Fig4Aii_%F")))
# 
 rootpath <- paste0(file.path(systempath, format(Sys.time(), "Fig4Aii_%F")))
 setwd(paste0(file.path(systempath, format(Sys.time(), "Fig4Aii_%F"))))
# 
# 
#setwd("~/Google Drive/Primates_Coculture/Sequencing20180530/ShotgunSequencing/Tables_Combined")

 ######This is from the Final_Analyses_vo4.R script
 #set directory and load data
 setwd("/Volumes/GoogleDrive/My\ Drive/Primates_Coculture/RNAseq_Gopher/R_workspace")
 load("/Volumes/GoogleDrive/My\ Drive/Primates_Coculture/RNAseq_Gopher/R_workspace/subread_DESeq2-data.RDATA")
 
 
 #Filter for protein coding genes
 ########### Only keep protein-coding genes in the read.counts matrix ############## 
 #First retrieve gene ids
 listMarts(host="www.ensembl.org")
 geneIDs <- rownames(data)
 human <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "grch37.ensembl.org")
 
 #ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl") # for current ensembl Db
 ensembl_gene_table <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol', 'gene_biotype'), #output fields
                             filters = c('ensembl_gene_id','biotype'), #filter by
                             values = list(geneIDs,"protein_coding"), #filtering values
                             mart = human)
 dim(ensembl_gene_table)
 
 #Ensure all protein-coding
 all(ensembl_gene_table$gene_biotype=="protein_coding") #should be TRUE
 
 #subset original expression matrix by protein coding gene IDs in ensembl_gene_table
 data <- subset(data,rownames(data) %in% ensembl_gene_table$ensembl_gene_id)
 dim(data) #should have 15,000-20,000 genes (for human data)
 
 meta$Species <- as.character(meta$Species)
 meta$Microbiome <- as.character(meta$Microbiome)
 meta$Species[is.na(meta$Species)] <- "AA"
 meta$Microbiome[is.na(meta$Microbiome)] <- "AA"
 
 meta$Species_H[meta$Species == "Human"] <- "Human"
 meta$Species_H[meta$Microbiome == "AA"] <- "AA"
 meta$Species_H[is.na(meta$Species_H)] <- "BB"
 #meta$Species_H[meta$Species_H == "AA"] <- " "
 
 meta$Species_C[meta$Species == "Chimp"] <- "Chimp"
 meta$Species_C[meta$Microbiome == "AA"] <- "AA"
 meta$Species_C[is.na(meta$Species_C)] <- "BB"
 #meta$Species_C[meta$Species_C == "AA"] <- " "
 
 meta$Species_G[meta$Species == "Gorilla"] <- "Gorilla"
 meta$Species_G[meta$Microbiome == "AA"] <- "AA"
 meta$Species_G[is.na(meta$Species_G)] <- "BB"
 #meta$Species_G[meta$Species_G == "AA"] <- " "
 
 meta$Species_O[meta$Species == "Orangutan"] <- "Orangutan"
 meta$Species_O[meta$Microbiome == "AA"] <- "AA"
 meta$Species_O[is.na(meta$Species_O)] <- "BB"
 #meta$Species_O[meta$Species_O == "AA"] <- " "
 
 meta$Species_All[meta$Species == "AA"] <- "AA"
 meta$Species_All[meta$Species!= "AA"] <- "Primate"
 
 meta$Control.ID <- rep("AA",nrow(meta))
 meta$Control.ID <- factor(meta$Control.ID)
 ControlLevels <- levels(meta$Control.ID)
 
 
 
 meta$Species_H <- factor(meta$Species_H)
 meta$Species_C <- factor(meta$Species_C)
 meta$Species_G <- factor(meta$Species_G)
 meta$Species_O <- factor(meta$Species_O)
 meta$Species_All <- factor(meta$Species_All)

 #meta <- subset(meta, Species != "AA")
 #data <- subset(data, )
 
 #Modifying the metadata so that I can also use Individual as a variable in the model
 meta$Species <- as.character(meta$Species)
 meta$Microbiome <- as.character(meta$Microbiome)
 meta$Individual <- as.character(meta$Individual)
 meta$Individual[meta$Species == "Human" & meta$Microbiome == "1"] <- "Human1"
 meta$Individual[meta$Species == "Human" & meta$Microbiome == "65"] <- "Human65"
 meta$Individual[meta$Species == "Human" & meta$Microbiome == "110"] <- "Human110"
 meta$Individual[meta$Species == "Human" & meta$Microbiome == "111"] <- "Human111"
 meta$Individual[meta$Species == "AA"] <- "AA"
 
 
 meta$Species <- factor(meta$Species)
 meta$ExperimentPlate <- factor(meta$ExperimentPlate)
 meta$Microbiome <-factor(meta$Microbiome)
 
 meta$rnaseqlabel <- rownames(meta)
 
 #remove Gorilla 6 (Schroeder replicate). This is the one that looks like an orangutan. This corresponds to HT2 for the RNA seq data
 meta <- meta[!(meta$Microbiome %in% c("6")),]
 data <- subset(data, select = -CP4.HT2 )
 write.table(data, file = "Genes.txt", sep  = "\t")
#####End FinalAnalyses_vo4.R script segment######
 
#Create metadata table for the DESeq2 model. Need to combine the gen3 table with the meta table
library("gtools")
 
gen4 <- merge(gen3, meta, by = c("Microbiome", "Species", "ExperimentPlate", "Individual"), all.y = TRUE) 
gen4 <- subset(gen4, select = -s__unclassified)
colnames(gen4)

rownames(gen4) <- gen4$rnaseqlabel
gen4[, 6:41][is.na(gen4[, 6:41])] <- 0
gen4 <- gen4[order(row.names(gen4)),]

setwd("/Volumes/GoogleDrive/My Drive/Primates_Coculture/Primates/MyAnalysis_v02/SpeciesResults_2018-10-25")
write.table(gen4, file = "MetadataPathways.txt", sep = "\t")

##Get all OTUs will test
test <- colnames(gen4)[6:41]

 dir.create(file.path(rootpath, "/DESeq2_model"))
 dir.create(file.path(rootpath, "/Plots/"))
 dataDir <- file.path(rootpath, "/DESeq2_model")
# 
 runmany <- function(OTU,stuff,cov) {
   cat("#",OTU,"\n")
   names(cov)[names(cov) == eval(as.symbol("OTU"))] <- 'example'
   ddsFull <- DESeqDataSetFromMatrix(
     countData = round(stuff),
     colData = cov,
     design  = ~ ExperimentPlate + Species_All + example)
   cat("#ran\n")
   keep <- rowSums(counts(ddsFull)) > 20
   ddsFull <- ddsFull[keep,]
   ddsFull <- DESeq(ddsFull)
   res <- results(ddsFull)
   temp <- data.frame(res@rownames, res$'log2FoldChange', res$'lfcSE', res$'pvalue', res$'padj', stringsAsFactors=FALSE)
   temp2 <- temp[!is.na(temp$res.padj),]
   colnames(temp2) <- c("rownames","log2FoldChange","lfcSE","pvalue","padj")
   save(list=c("res", "ddsFull", "temp2"), file=paste(dataDir, '/DESeq2_',OTU,'.RData', sep=''))
 }
 
 sapply(test,runmany,stuff=data,cov=gen4)
###Test analysis
 # dds_test <- DESeqDataSetFromMatrix(countData = data, colData = gen4, design=~ ExperimentPlate + Species_All +  g__Anaerotruncus )
 # dds_test <- DESeq(dds_test)
 # res_Test <- results(dds_test)


##Run the analysis

##load data object



rootpath1 <- "/Volumes/GoogleDrive/My\ Drive/Primates_Coculture/Primates/MyAnalysis_v02/SpeciesResults_2018-10-25/DESeq2_model/"
setwd("/Volumes/GoogleDrive/My\ Drive/Primates_Coculture/Primates/MyAnalysis_v02/SpeciesResults_2018-10-25/DESeq2_model/")
fileList=dir(rootpath1,"*.RData$")

openRData <- function(file)  {
  x <- gsub("DESeq2_","",file)
  example <- gsub(".RData","",x)
  
  load(paste0(rootpath1, file))
  tempe <- data.frame(res@rownames, res$'log2FoldChange', res$'lfcSE', res$'pvalue', res$'padj', stringsAsFactors=FALSE)
  colnames(tempe) <- c("rownames","log2FoldChange","lfcSE","pvalue","padj")
  tempe$OTU <- rep(example,nrow(tempe))
  tempe
}

tempe <- lapply(fileList,openRData)
tempe3 <- do.call(rbind,tempe)

#write.table(temp3,"test_1.txt",sep="\t")



colnames(tempe3) <- c("ENSG","log2FoldChange","lfcSE","pvalue","padj","OTU")

##remove those lines with NA in the pval because they are useless in the next analysis
tempe5 <- tempe3[!is.na(tempe3$pvalue),]

##Add ensg


##Recalculate adjusted p for combined group of gene-taxa pairs
tempe5$newpadj <- p.adjust(tempe5$pvalue,method="BH")


##I think I'll need to do more test corrections but let's just take a look
sig <- tempe5[tempe5$padj < 0.1,]
signew <- tempe5[tempe5$newpadj < 0.1,]


##Get important OTUs
otuCount <- plyr::count(signew, vars = "OTU") 
rownames(otuCount) <- otuCount$OTU
sig2 <- merge(signew,otuCount,by="OTU",all=TRUE) ##12820 rows
length(unique(sig2$OTU)) ## 36
length(unique(sig2$ENSG)) ##231 genes
table(sig2$freq) ##Take each value on the bottom and divide by the top!!!



##Let's get a look at a few examples
##First we'll need the separated gene expression data
library ('tidyr')

#################################################
##Adapted from ExpressionHeatMap_fixed_032317.R##
#################################################
##Read in expression data
# FileFolder='~/piquelab/scratch/Allison/microNEW/expression/DESeq2/out_data_CP3_separate/stats/'
# fileList=dir(FileFolder,'CP3_DEG_stats_.*\\.txt')

##Get sample names so can extract from RData
x <- strsplit(fileList,"__")
z <- sapply(x,function(y) {paste0(y[2])})
samples <- gsub("[.]RData","",z)



DEG <- aux6[aux6$count > 0,]
library("compositions")
##Load microbiome data
primate_genera1 <- primate_genera
primate_genera1[primate_genera1 ==0] <- 0.0065
otu <- t(clr(primate_genera1))
otu <- otu[,1:19]
names <- primate_meta16s
names <- names[names$Treatment=="Uncultured",]


##Focus on OTU stuff that is also DEG with normal individual model
sig3 <- sig2[sig2$ENSG %in% DEG$ENSG,] ##588 -> 393 rows
length(unique(sig3$OTU)) ##30 OTUs
length(unique(sig3$ENSG)) ##127 transcripts

###Want to set a new cutoff -> goal is to get ~25 unique species and 96 genes - > let's use 0.05 

sigt <- sig3[sig3$newpadj<0.05,]

length(unique(sigt$OTU)) ##25
length(unique(sigt$ENSG)) ##80
sig3 <- sigt

#Getting ready to plot things, need to load baseline data (use COL samples)
###Work on this


otu2 <- otu[rownames(otu) %in% unique(sig3$OTU),]
otu3 <- otu2[,grepl("UNCULTURED",colnames(otu2))]

dds_tc <- DESeqDataSetFromMatrix(countData = data, colData = meta, design = ~ ExperimentPlate + Species_All)

dds_tc <- dds_tc[ rowSums(counts(dds_tc)) > 20, ]
dds_tc <- DESeq(dds_tc)

vst1<- varianceStabilizingTransformation(dds_tc, blind = TRUE)
vst2 <- as.data.frame(assay(vst1))
vst2$symbol <- mapIds(org.Hs.eg.db,
                      keys=row.names(vst2),
                      column="SYMBOL",
                      keytype="ENSEMBL",
                      multiVals="first")



rootpath <- paste0(file.path(systempath, format(Sys.time(), "Fig4Aii_%F")))
setwd(paste0(file.path(systempath, format(Sys.time(), "Fig4Aii_%F"))))
dir.create(paste0(rootpath, "/Plots/"))
setwd(paste0(rootpath, "/Plots/"))
spearp<- vector()
spearrho <- vector()
pearp <- vector()
pearr <-vector()


plotshit <- function(transcript,data,sig,otusamp,otubase,num2){
  setup <- data[data$ENSG==eval(as.symbol("transcript")),]
  setup2 <- gather(setup,sample,expr) ##this is not actually a log fold change, it is an expression
  #################################Be careful with this next step
  setup2 <- setup2[grepl("CP4",setup2$sample),]
  setup2$mic <- meta[setup2$sample,"Microbiome"]
  setup2$spec <- meta[setup2$sample, "Species"]
  setup2$match <- paste0(setup2$mic,"UNCULTURED")
  otu <- data.frame(value = otubase, dummy = 1)
  setup2$base <- otu[setup2$match,"value"]
  setup2$base[is.na(setup2$base)] <- 0
  ##make line graph
  setup3 <- setup2[order(setup2$spec),]
  setup3$mic <- factor(setup3$mic,levels=unique(setup3$mic))
  symbol <- setup[setup$ENSG == transcript, "symbol"]
  ##  setup3$base <- as.numeric(setup3$base)
  setup3$expr <- as.numeric(setup3$expr)
  setup4 <- setup3[setup3$mic != "AA",]
  int <- coef(lm(formula=expr ~ base, data = setup4))[[1]][1]
  sl <- coef(lm(formula=expr ~ base, data = setup4))[[2]][1]
  stuff <- cor.test(setup4$base,setup4$expr,method="spearman")
  stuff2 <- cor.test(setup4$base,setup4$expr,method="pearson")
  
  ##collect stats from correlation test
  spearp <<- c(spearp, stuff$p.value)
  spearrho <<- c(spearrho, stuff$estimate)
  pearp <<- c(pearp, stuff2$p.value)
  pearr <<- c(pearr, stuff2$estimate)
  
  
  ##decide where to annotate
  yval <- ifelse(sl > 0, min(setup2$expr),max(setup2$expr))
  xval <- ifelse(sl > 0, max(setup2$expr),min(setup2$expr))
  numname <- num2[num2$OTU==eval(as.symbol("otusamp")),"num"]
  x <- strsplit(otusamp,";")
   y <- sapply(x,function(y) {ifelse(y[2]!="g__",y[2])})
                                     # ifelse(y[5]!="f__",y[5],
                                     #        ifelse(y[4]!="o__",y[4],
                                     #               ifelse(y[3]!="c__",y[3],
                                     #                      ifelse(y[2]!="p__",y[2],"WTF")))))})
  
  
  y <- gsub("\\[","",y)
  y <- gsub("\\]","",y)
  z <- strsplit(y,"__")
  zz <- sapply(z,function(y) {z[[1]][2]})
  otusamp <- zz
 
  label <- x
  setup3$facet <- ifelse(setup3$mic=="AA","C"," ")
  setup3$mic <- as.character(setup3$mic)
  setup3$spec <- as.character(setup3$spec)
  setup3 <- rbind(setup3,c("dummy",mean(setup3$expr),"1000","dummy","dummy", max(setup4$base)/15,"C"))
  setup3$expr <- as.numeric(setup3$expr)
  setup3 <- rbind(setup3,c("dummy",mean(setup3$expr),"1000","dummy","dummy",max(setup4$base)/-15,"C"))
  setup3$expr <- as.numeric(setup3$expr)
  setup3$base <- as.numeric(setup3$base)
  ##make dataframe for line
  dummy2 <- data.frame(X = c("C", " "), intercept = c(0,int), slope=c(0,sl))
  dummy2$facet <- as.character(dummy2$X)
  setup3$sample <- factor(setup3$sample,levels=unique(setup3$sample))
  file1 <- paste0(numname,"_",symbol,"_12132017.pdf")
  pdf(file1)
  p<-  ggplot(setup3, aes(x=base, y=expr, color=spec)) +
    geom_point(size=8) +
    labs(title = paste0(transcript,"-",symbol,"\n","r = ",round(as.numeric(stuff2$estimate),digits = 6),"; p-value = ", round(as.numeric(stuff2$p.value), digits = 6),"; \nrho = ",round(as.numeric(stuff$estimate),digits = 6),"; p-value = ",round(as.numeric(stuff$p.value),digits = 6)),x = paste0("CLR transform\n",label),y="Normalized Gene Expression") +
    facet_grid(. ~ facet,scales = "free_x",space="free_x") +
    geom_abline(data=dummy2,aes(intercept = intercept, slope = slope)) +
    scale_color_manual(values=c("gray", "#F8766D", "white", "#7CAE00", "#00BFC4", "#C77CFF")) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.text.x = element_text(size=18),
          axis.text.y = element_text(size=18),
          axis.title = element_text(size = 20,face="bold"),
          strip.text.x = element_text(size=18),
          panel.border = element_rect(size=2, color="black"),
          axis.ticks = element_line(size=2, color="black"),
          strip.background = element_rect(color = "black", size = 2))
  print(p)
  dev.off()
}

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
cols = gg_color_hue(4)


foreachOTU <- function(um,inddata,genusdata,otutab,num) {
  transcripts <- genusdata[genusdata$OTU==um,"ENSG"]
  table3 <- inddata[transcripts,]
  table3$ENSG <- rownames(table3)
  otutab2 <- otutab[eval(as.symbol("um")),]
  sapply(transcripts,plotshit,data=table3,otusamp=um,otubase=otutab2,num2=num)
}

##Make otu table for self so names aren't so ridiculous
def <- data.frame(unique(sig3$OTU))
colnames(def) <- c("OTU")
def$num <- c(1:nrow(def))

sapply(unique(sig3$OTU),foreachOTU,inddata=vst2,genusdata=sig3,otutab=otu3,num=def)




sigoth <- cbind(spearp, spearrho, pearp, pearr)

# 
 sig3 <- cbind(sig3, sigoth)
 sig3$padjustSpearman <- p.adjust(sig3$spearp, method = "BH")
 sig3$padjustPearson <- p.adjust(sig3$pearp, method = "BH")

sig3.2 <- sig3[sig3$newpadj<0.1,]
length(unique(sig3.2$OTU)) #25
length(unique(sig3.2$ENSG)) #80


sig4 <- merge(sig3.2,def,by="OTU")  

setwd(rootpath)
write.table(sig4,"significantRows_genusModel_07092018.txt",sep="\t",quote=FALSE,row.names=FALSE)

#########Plot everything############
###Heatmaps
length(unique(sig4$OTU))
#25
length(unique(sig4$ENSG))
#80

sensg <- unique(sig4$ENSG)
sotu <- unique(sig4$OTU)

gg1 <- as.character(sensg)
aa1 <- as.character(sotu)

data2 <- vst2[,1:31]

dat1 <- t(data2)
bb2 <-gg1
gen5 <- cbind(gen4, dat1)
gen6 <- gen5[1:19,]
bb3 <- vector()
dog <- data.frame()

taxgen <- vector()
taxp <- vector()
# ##just testing
# aa1 <- head(aa1,2)
# bb2 <- head(bb2)

taxgenlist <-list()
taxplist <- list()
for(i in 1:length(aa1)){
  ##For all the taxa
  taxgen <- vector()
  taxpp <-vector()
  for(j in 1:length(bb2)){
    ##For all genes
    pp <- cor.test(gen6[,aa1[i]], gen6[,bb2[j]], method = "pearson")
    taxgen <- c(taxgen, pp$estimate)
    taxpp <- c(taxpp, pp$p.value)
    
  }
  taxgenlist[[i]]<- taxgen
  taxplist[[i]] <-taxpp
}

taxacor <- do.call(rbind, taxgenlist)
taxapval <- do.call(rbind, taxplist)

bb2 <- as.data.frame(bb2)

bb2$symbol <- mapIds(org.Hs.eg.db,  
                     keys=as.character(bb2$bb2),
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")
bb2
#colnames(taxacor) <- bb2$symbol
colnames(taxacor) <- bb2$symbol
rownames(taxacor) <- aa1

#colnames(taxapval) <-bb2$symbol
colnames(taxapval) <- bb2$symbol
rownames(taxapval) <- aa1

write.table(taxacor, file = "MicrobialSpecies_Correlations.txt", sep="\t")

library("reshape2")


###Now plot the complex heatmap with annotations for which genes are species specific

setwd("/Volumes/GoogleDrive/My Drive/Primates_Coculture/Primates/MyAnalysis_v02/Results_2019-02-04/DE_Model2")
Modcombine1 <- read.table("Modcombine_2.csv", sep = ",")

gs2 <- Modcombine1

gs2$symbol <- mapIds(org.Hs.eg.db,  
                     keys=as.character(rownames(gs2)),
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")

gs3 <- gs2[gs2$symbol %in% colnames(taxacor),]


gs3$Human <- factor(gs3$Human)
gs3$Orangutan <- factor(gs3$Orangutan)
gs3$Gorilla <- factor(gs3$Gorilla)
gs3$Chimp <- factor(gs3$Chimp)

rownames(gs3) <- gs3$symbol
gs4 <- gs3[,1:4]

gs4 <-gs4[order(match(rownames(gs4),colnames(taxacor))),]


library("ComplexHeatmap")
col_ha <- HeatmapAnnotation(df = gs4, 
                            col = list(Chimp = c("1"="#F8766D","0"="gray90"),
                                       Gorilla = c("1" = "#7CAE00","0"="gray90"),
                                       Human = c("1" = "#00BFC4" , "0"="gray90"),
                                       Orangutan =c("1" = "#C77Cff", "0"="gray90")),
                            show_legend = FALSE)
draw(col_ha, 1:10) 
#To draw the barplots on the side for the genera
#I'm going to use gen1.11 because I have removed all of the non UNCULTURED samples and I haven't transformed any of the abundance values
#get averages of each microbial species for each primate species



allchimps <- c("100UNCULTURED", "400UNCULTURED", "300UNCULTURED", "500UNCULTURED")
allgorillas <- c("2UNCULTURED", "4UNCULTURED","14UNCULTURED","26UNCULTURED","40UNCULTURED","19UNCULTURED","24UNCULTURED")
allhumans <- c("1UNCULTURED","65UNCULTURED","110UNCULTURED","111UNCULTURED")
allorangutans <- c("8UNCULTURED","18UNCULTURED","36UNCULTURED","28UNCULTURED")

setwd("/Volumes/GoogleDrive/My\ Drive/Primates_Coculture/Primates/MyAnalysis")
primate_meta16s <- read.table("primate_metagenomic_meta.txt", sep = "\t", header = TRUE, row.names ="SampleID")
setwd("/Volumes/GoogleDrive/My\ Drive/Primates_Coculture/Sequencing20180530/ShotgunSequencing/Tables_Combined")
primate_genera <- t(read.table("Merged_Species.txt", sep = "\t", header = TRUE, row.names = 1, check.names = FALSE))

#deal with metadata formatting first
primate_meta16s$sampname <- rownames(primate_meta16s)
#Get rid of UMGC blank
primate_meta16s <- primate_meta16s[!(primate_meta16s$Microbiome == "BLANK"),]

#deal with sequencing data formatting
newrows <- vapply(strsplit(rownames(primate_genera),"_"), `[`, 1, FUN.VALUE=character(1))
rownames(primate_genera) <- newrows
#Rename the luca samples
rownames(primate_genera)[16] <- "1UNCULTURED"
rownames(primate_genera)[17] <- "65UNCULTURED"
rownames(primate_genera)[18] <- "110UNCULTURED"
rownames(primate_genera)[19] <- "111UNCULTURED"

#Subset the species so it is between 50-100 using column sums.
gen1 <-primate_genera

#Subset abundance data so that I have only the uncultured samples 
#Get rid of low abundance species
gen1 <-as.data.frame(gen1)
gen1.1<- gen1[grepl("UNCULTURED",rownames(gen1)),]
#I want to keep species that are represented in at least half the samples (8)
cutoff_nsamples <- colSums(gen1.1>0) > 9
gen1.11<- gen1.1[,cutoff_nsamples]

microc <-gen1.11[rownames(gen1.11) %in% allchimps,]
microg <- gen1.11[rownames(gen1.11) %in% allgorillas,]
microh <- gen1.11[rownames(gen1.11) %in% allhumans,]
microo <- gen1.11[rownames(gen1.11) %in% allorangutans,]

avgc <- as.data.frame(log2(colMeans(microc)))
avgg <- as.data.frame(log2(colMeans(microg)))
avgh <- as.data.frame(log2(colMeans(microh)))
avgo <- as.data.frame(log2(colMeans(microo)))

microc1 <- microc[,colnames(microc) %in% rownames(taxacor)]
microg1 <- microg[,colnames(microg) %in% rownames(taxacor)]
microh1 <- microh[,colnames(microh) %in% rownames(taxacor)]
microo1 <- microo[,colnames(microo) %in% rownames(taxacor)]

microc1[microc1==0] <- 0.65
microg1[microg1==0] <- 0.65
microh1[microh1==0] <- 0.65
microo1[microo1==0] <- 0.65

microc1 <- log2(microc1)
microg1 <- log2(microg1)
microh1 <- log2(microh1)
microo1 <- log2(microo1)

microg1[microg1== -Inf]<-0
microc1[microc1== -Inf]<-0
microh1[microh1== -Inf]<-0
microo1[microo1== -Inf]<-0

microg1[microg1== Inf]<-0
microc1[microc1== Inf]<-0
microh1[microh1== Inf]<-0
microo1[microo1== Inf]<-0

micall <- cbind(avgc, avgg, avgh, avgo)
tch <-rownames(micall)
#remove infs
micall <-do.call(data.frame,lapply(micall, function(x) replace(x, is.infinite(x),0)))

colnames(micall) <- c("Chimp", "Gorilla", "Human", "Orangutan")
rownames(micall) <- tch
micall1 <- micall[rownames(micall) %in% rownames(taxacor),]
micallc <- as.data.frame(micall1$Chimp)
rownames(micallc) <- rownames(micall1)
colnames(micallc) <- c("Chimp")

micallg <- as.data.frame(micall1$Gorilla)
rownames(micallg) <- rownames(micall1)
colnames(micallg) <- c("Gorilla")

micallh <- as.data.frame(micall1$Human)
rownames(micallh) <- rownames(micall1)
colnames(micallh) <- c("Human")

micallo <- as.data.frame(micall1$Chimp)

rownames(micallo) <- rownames(micall1)
colnames(micallo) <- c("Orangutan")

library("stringr")

rownames(taxacor) <- gsub('s__', '', rownames(taxacor))



rownames(taxapval) <- gsub('s__', '', rownames(taxapval))


row_ha <- rowAnnotation(micall1 = row_anno_boxplot(micall1, 
                                                   axis = TRUE, 
                                                   gp = gpar(fill = c(Chimp = "#F8766D", Gorilla ="#7CAE00",Human = "#00BFC4", Orangutan  = "#C77Cff"))), 
                        width = unit(1.5, "cm"))

microc1 <- microc1[,order(match(colnames(microc1),rownames(taxacor)))]
microg1 <- microg1[,order(match(colnames(microg1),rownames(taxacor)))]
microh1 <- microh1[,order(match(colnames(microh1),rownames(taxacor)))]
microo1 <- microo1[,order(match(colnames(microo1),rownames(taxacor)))]

colnames(microc1) <- rownames(taxacor)
colnames(microg1) <- rownames(taxacor)
colnames(microh1) <- rownames(taxacor)
colnames(microo1) <- rownames(taxacor)
##Make barplot for chimp abundances
row_hc <- HeatmapAnnotation(boxplot = row_anno_boxplot(microc1, outline = FALSE,
                                                       axis = TRUE, 
                                                     ylim = c(-9,5),
                                                       gp = gpar(fill = c(Chimp = "#F8766D"))), 
                            width = unit(2, "cm"),which = "row" )

draw(row_hc, 1:10)

#Make barplot for gorilla abundances
row_hg <- HeatmapAnnotation(boxplot = row_anno_boxplot(microg1, outline = FALSE,
                                                       axis = TRUE, 
                                                       ylim=c(-9,5),
                                                       gp = gpar(fill = c(Gorilla ="#7CAE00"))), 
                            width = unit(2, "cm"),which = "row")


draw(row_hg, 1:10)

#Make barplot for human abundances
row_hh <- HeatmapAnnotation(boxplot = row_anno_boxplot(microh1, outline=FALSE,
                                                       axis = TRUE, 
                                                       ylim=c(-9,5),
                                                       gp = gpar(fill = c(Human = "#00BFC4"))), 
                            width = unit(2, "cm"), which="row")

draw(row_hh, 1:10)

#Make barplot for orangutan abundances
row_ho <- HeatmapAnnotation(boxplot = row_anno_boxplot(microo1, outline = FALSE,
                                                       axis = TRUE, 
                                                       ylim=c(-9,5),
                                                       gp = gpar(fill = c(Orangutan  = "#C77Cff"))), 
                            width = unit(2, "cm"), which="row")

draw(row_ho, 1:10)
systempath <- "/Volumes/GoogleDrive/My\ Drive/Primates_Coculture/Primates/MyAnalysis_v02"
setwd(systempath)
dir.create(file.path(systempath, format(Sys.time(), "Fig4A_%F")))
rootpath <- paste0(file.path(systempath, format(Sys.time(), "Fig4A_%F")))

##Set colors for heatmap correlations
col1 = colorRampPalette(c( "#999999", "#999999","#ffffff","#ef8a62","#ef8a62"))(n = 299)

setwd(rootpath)

pdf("CompHeatmap_sig.pdf")
p <- Heatmap(taxacor, 
             col = col1,
             top_annotation = col_ha, 
             km = 1, 
             row_names_gp = gpar(fontsize = 10, fontface="italic", base_family="Arial"), 
             column_names_gp = gpar(fontsize=3), width = unit(12.1, "cm")) + row_hc + row_hg+  row_hh + row_ho
p <- draw(p, show_heatmap_legend=FALSE)
print(p)
dev.off()

##Draw heatmap without barplots
setwd(rootpath)
pdf("CompHeatmap_nobars.pdf")
p <- Heatmap(taxacor, 
             col=col1,
             top_annotation = col_ha, 
             km = 1, 
             row_names_gp = gpar(fontsize = 10, fontface="italic", base_family="Arial"), 
             column_names_gp = gpar(fontsize=3), width = unit(10, "cm"))
p <- draw(p, show_heatmap_legend=FALSE)
print(p)
dev.off()

#Heatmap with only chimp bars
pdf("CompHeatmap_chimpbars.pdf")
p <- Heatmap(taxacor, 
             col=col1,
             top_annotation = col_ha, 
             km = 1, 
             row_names_gp = gpar(fontsize = 4, fontface="italic", base_family="Arial"), 
             column_names_gp = gpar(fontsize=3), width = unit(12.1, "cm")) + row_hc 
p <- draw(p, show_heatmap_legend=FALSE)
print(p)
dev.off()

#HEatmap with only gorilla bars
pdf("CompHeatmap_gorillabars.pdf")
p <- Heatmap(taxacor, 
             col=col1,
             top_annotation = col_ha, 
             top_annotation_height = unit(1, "cm"), 
             km = 1, 
             row_names_gp = gpar(fontsize = 4, fontface="italic", base_family="Arial"), 
             column_names_gp = gpar(fontsize=3), width = unit(12.1, "cm")) + row_hg
p <- draw(p, show_heatmap_legend=FALSE)
print(p)
dev.off()

#Heatmap with only human bars
pdf("CompHeatmap_humanbars.pdf")
p <- Heatmap(taxacor, 
             col=col1,
             top_annotation = col_ha, 
             top_annotation_height = unit(1, "cm"), 
             km = 1, 
             row_names_gp = gpar(fontsize = 4, fontface="italic", base_family="Arial"), 
             column_names_gp = gpar(fontsize=3), width = unit(12.1, "cm")) + row_hh
p <- draw(p, show_heatmap_legend=FALSE)
print(p)
dev.off()


#Heatmap with only orangutan bars
pdf("CompHeatmap_orangutanbars.pdf")
p <- Heatmap(taxacor, 
             col=col1,
             top_annotation = col_ha, 
             top_annotation_height = unit(1, "cm"), 
             km = 1, 
             row_names_gp = gpar(fontsize = 4, fontface="italic", base_family="Arial"), 
             column_names_gp = gpar(fontsize=3), width = unit(12.1, "cm")) + row_ho
p <- draw(p, show_heatmap_legend=FALSE)
print(p)
dev.off()

q <- Heatmap(taxacor, top_annotation = col_ha, top_annotation_height = unit(1, "cm"), km = 1, row_names_gp = gpar(fontsize = 5.5), column_names_gp = gpar(fontsize=5.5)) + row_ha
draw(q, show_heatmap_legend=FALSE)


pdf("CompHeatmap_justbars.pdf")
p <- Heatmap(taxacor, 
             col=col1,
             top_annotation = col_ha, 
             top_annotation_height = unit(1, "cm"), 
             km = 1, 
             row_names_gp = gpar(fontsize = 4, fontface="italic", base_family="Arial"), 
             column_names_gp = gpar(fontsize=3), width = unit(1, "cm")) + row_hh + row_hc + row_hg + row_ho
print(p)
dev.off()

#setwd("/Volumes/GoogleDrive/My\ Drive/Primates_Coculture/Primates/MyAnalysis_v02/Fig4Aii_2019-02-27")
setwd(rootpath)

pvals <- melt(taxapval)
rhos <- melt(taxacor)


tabs1 <- pvals %>% right_join(rhos, by=c("Var1","Var2"))

colnames(tabs1) <-c("Microbial Species", "Gene", "P-value", "Correlation")
tabs1$Qvalue <- p.adjust(tabs1$`P-value`, method = "BH")

write.table(tabs1,file="SupplementalTable5_CorrSpecies.txt", sep = "\t")

sigprint <- sig4

sigprint$symbol <- mapIds(org.Hs.eg.db,  
                     keys=as.character(sigprint$ENSG),
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")

sigprint$OTU <- str_remove_all(sigprint$OTU, "s__")

write.table(sigprint, file = "MicrobesSpeciesGenesModel.txt", sep = "\t")




