#Model 6 
###This script makes Fig. 4B heatmaps
#Load data first.
library("hexbin")
library("parathyroidSE")
library("DESeq2")
library('lmtest')
library('ggplot2')
library("parathyroidSE")
library("pheatmap")
library("RColorBrewer")
library("PoiClaClu")
library("dplyr")
library("vsn")
library("AnnotationDbi")
library("org.Hs.eg.db")
#source("https://bioconductor.org/biocLite.R")
#biocLite("goseq")
library("lattice")
library("hash")
library("biomaRt")
library("gplots")
library("Cormotif")
library("UpSetR")
library("tidyr")
##install.packages("Hmisc", dependencies=T)
library("Hmisc")
library("compositions")



setwd("/Volumes/GoogleDrive/My\ Drive/Primates_Coculture/Sequencing20180530/ShotgunSequencing/Humann2")

pab1 <- read.table("PrimatePathAbundance.tsv", sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
pcv1 <- read.table("PrimatePathCoverage.tsv", sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)



#deal with sequencing data formatting

newrows2 <- vapply(strsplit(colnames(pab1),"_"), `[`, 1, FUN.VALUE=character(1))
newrows3 <- vapply(strsplit(colnames(pcv1),"_"), `[`, 1, FUN.VALUE=character(1))


colnames(pab1) <- newrows2
colnames(pcv1) <- newrows3



colnames(pab1)[16] <- "1UNCULTURED"
colnames(pab1)[17] <- "65UNCULTURED"
colnames(pab1)[18] <- "110UNCULTURED"
colnames(pab1)[19] <- "111UNCULTURED"

colnames(pcv1)[16] <- "1UNCULTURED"
colnames(pcv1)[17] <- "65UNCULTURED"
colnames(pcv1)[18] <- "110UNCULTURED"
colnames(pcv1)[19] <- "111UNCULTURED" 

pab1.2 <- pab1[!grepl("\\|", row.names(pab1)),]
pcv1.2 <- pcv1[!grepl("\\|", row.names(pcv1)),]


setwd("/Volumes/GoogleDrive/My\ Drive/Primates_Coculture/Primates/MyAnalysis")
primate_meta16s <- read.table("primate_metagenomic_meta.txt", sep = "\t", header = TRUE, row.names ="SampleID")

#deal with metadata formatting first
primate_meta16s$sampname <- rownames(primate_meta16s)
#Get rid of UMGC blank
primate_meta16s <- primate_meta16s[!(primate_meta16s$Microbiome == "BLANK"),]

#deal with sequencing data formatting
###First run the model for pathway abudnances:
primate_genera <- pab1.2
newrows <- vapply(strsplit(rownames(primate_genera),"_"), `[`, 1, FUN.VALUE=character(1))
rownames(primate_genera) <- newrows


#Subset the species so it is between 50-100 using column sums.
gen1 <-primate_genera

#Subset abundance data so that I have only the uncultured samples 
#Get rid of low abundance species
gen1 <-as.data.frame(gen1)
gen1 <- t(gen1)
gen1.1<- gen1[grepl("UNCULTURED",rownames(gen1)),]
#I want to keep species that are represented in at least half the samples (8)
cutoff_nsamples <- colSums(gen1.1>0) > 10
gen1.11<- gen1.1[,cutoff_nsamples]

gen1.2 <-gen1.11[,colSums(gen1.11)>8000]
dim(gen1.2)


gen1.3 <- clr(gen1.2)
gen1.3 <-as.data.frame(gen1.3)
gen1.3$MBID <- rownames(gen1.3)

##Also need to subset the metadata so that it I only have the "uncultured" samples
metagen <- as.data.frame(primate_meta16s)
metagen1 <- metagen[grepl("UNCULTURED", rownames(metagen)),]
#Get rid of  microbiome samples that were not in the colonocyte experiment
metagen1.2 <- metagen1[rownames(metagen1) %in% rownames(gen1.3),]
metagen1.2$MBID <- rownames(metagen1.2)

gen2 <-merge(gen1.3, metagen1.2, by = "MBID", all.y = TRUE)



#These indicies are hardcoded and need to be changed if code is updated!!!!!
colnames(gen2)
gen2[, 2:96][is.na(gen2[, 2:96])] <- 0

gen3 <-gen2[,!grepl("unclassified", colnames(gen2))]


#MAke a new directory for the results
systempath <- "/Volumes/GoogleDrive/My\ Drive/Primates_Coculture/Primates/MyAnalysis_v02"
setwd(systempath)
dir.create(file.path(systempath, format(Sys.time(), "PathAbundanceResults_%F")))
# 
rootpath <- paste0(file.path(systempath, format(Sys.time(), "PathAbundanceResults_%F")))
setwd(paste0(file.path(systempath, format(Sys.time(), "PathAbundanceResults_%F"))))
# 
# 
dir.create(file.path(rootpath, "Plots"))
dir.create(file.path(rootpath, "BaseModel2"))
dir.create(file.path(rootpath, "BaseModel3"))


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


meta <- meta[!(meta$Microbiome %in% c("6")),]
data <- subset(data, select = -CP4.HT2 )



library("gtools")

gen4 <- merge(gen3, meta, by = c("Microbiome", "Species", "ExperimentPlate", "Individual"), all.y = TRUE) 

rownames(gen4) <- gen4$rnaseqlabel
colnames(gen4)
gen4[, 6:100][is.na(gen4[, 6:100])] <- 0
gen4 <- gen4[order(row.names(gen4)),]

gen4 <- gen4[, !(colnames(gen4) %in% c("UNMAPPED", "UNINTEGRATED"))]

##Get all OTUs will test
colnames(gen4)
test <- colnames(gen4)[6:98]


 dataDir <- paste0(rootpath, "/BaseModel3/")
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
   system.time(ddsFull <- DESeq(ddsFull,parallel=TRUE))
   res <- results(ddsFull)
   temp <- data.frame(res@rownames, res$'log2FoldChange', res$'lfcSE', res$'pvalue', res$'padj', stringsAsFactors=FALSE)
   temp2 <- temp[!is.na(temp$res.padj),]
   colnames(temp2) <- c("rownames","log2FoldChange","lfcSE","pvalue","padj")
   save(list=c("res", "ddsFull", "temp2"), file=paste(dataDir, '/DESeq2_',OTU,'.RData', sep=''))
 }
# 
 sapply(test,runmany,stuff=data,cov=gen4)



##Run the analysis

##load data object

rootpath1 <- paste0(rootpath, "/BaseModel/")
rootpath2 <- paste0(rootpath, "/BaseModel2/")
rootpath3 <- "/Volumes/GoogleDrive/My\ Drive/Primates_Coculture/Primates/MyAnalysis_v02/PathAbundanceResults_2018-10-04/BaseModel2/"
rootpath4 <- paste0(rootpath, "/Plots/")
rootpath5 <- paste0(rootpath, "/BaseModel3/")
rootpath6 <- paste0(rootpath, "/Plots3/")





rootpathsaved <- "/Volumes/GoogleDrive/My\ Drive/Primates_Coculture/Primates/MyAnalysis_v02/PathAbundanceResults_2018-10-23/DESeq2_model/"
fileList=dir(rootpathsaved,"*.RData$")
setwd("/Volumes/GoogleDrive/My\ Drive/Primates_Coculture/Primates/MyAnalysis_v02/PathAbundanceResults_2018-10-23/DESeq2_model/")
openRData <- function(file)  {
  x <- gsub("DESeq2_","",file)
  example <- gsub(".RData","",x)
  
  load(paste0(rootpathsaved, file))
  tempe <- data.frame(res@rownames, res$'log2FoldChange', res$'lfcSE', res$'pvalue', res$'padj', stringsAsFactors=FALSE)
  colnames(tempe) <- c("rownames","log2FoldChange","lfcSE","pvalue","padj")
  tempe$OTU <- rep(example,nrow(tempe))
  tempe
}

tempe <- lapply(fileList,openRData)
tempe3 <- do.call(rbind,tempe)


setwd(rootpath5)


colnames(tempe3) <- c("ENSG","log2FoldChange","lfcSE","pvalue","padj","OTU")

##remove those lines with NA in the pval because they are useless in the next analysis
tempe5 <- tempe3[!is.na(tempe3$pvalue),]

##Add ensg


##Recalculate adjusted p
tempe5$newpadj <- p.adjust(tempe5$pvalue,method="BH")


##I think I'll need to do more test corrections but let's just take a look
sig <- tempe5[tempe5$newpadj < 0.1,]
signew <- tempe5[tempe5$newpadj < 0.0001,]
signew <- tempe5[tempe5$newpadj < 0.05,]


##Get important OTUs
otuCount <- plyr::count(signew, vars = "OTU") 
rownames(otuCount) <- otuCount$OTU
sig2 <- merge(signew,otuCount,by="OTU",all=TRUE) ##12820 rows
length(unique(sig2$OTU)) ## 109 pathways
length(unique(sig2$ENSG)) ##389 genes
table(sig2$freq) ##Take each value on the bottom and divide by the top!!!



##Let's get a look at a few examples
##First we'll need the separated gene expression data
library ('tidyr')


##Get sample names so can extract from RData
x <- strsplit(fileList,": ")
z <- sapply(x,function(y) {paste0(y[2])})
samples <- gsub("[.]RData","",z)



DEG <- aux6[aux6$count > 0,]

##Load microbiome data
otu <- t(clr(primate_genera))
#otu <- otu[,1:19]
names <- primate_meta16s
names <- names[names$Treatment=="Uncultured",]

##Focus on OTU stuff that is also DEG with normal individual model
sig3 <- sig2[sig2$ENSG %in% DEG$ENSG,] ##582 ->165
length(unique(sig3$OTU)) ##[1] 89
length(unique(sig3$ENSG)) ##[1] 310

#Getting ready to plot things, need to load baseline data (use COL samples)
###Work on this
otu <-t(otu)

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

setwd(paste0(rootpath5))
pdgspear<- vector()
pdgpear <- vector()
gtrans <- vector()
mic1 <- vector()

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

  setup3$expr <- as.numeric(setup3$expr)
  setup4 <- setup3[setup3$mic != "AA",]
  int <- coef(lm(formula=expr ~ base, data = setup4))[[1]][1]
  sl <- coef(lm(formula=expr ~ base, data = setup4))[[2]][1]
  stuff <- cor.test(setup4$base,setup4$expr,method="spearman")
  stuff2 <- cor.test(setup4$base,setup4$expr,method="pearson")
  
  
  ##Save adjusted pvalues for FDR correction
  pdgspear <<- c(pdgspear, stuff$p.value)
  pdgpear <<- c(pdgpear, stuff$p.value)
  gtrans <<- c(gtrans, symbol)
  
  ##collect stats from correlation test
  spearp <<- c(spearp, stuff$p.value)
  spearrho <<- c(spearrho, stuff$estimate)
  pearp <<- c(pearp, stuff2$p.value)
  pearr <<- c(pearr, stuff2$estimate)
  
  
  ##decide where to annotate
  yval <- ifelse(sl > 0, min(setup2$expr),max(setup2$expr))
  xval <- ifelse(sl > 0, max(setup2$expr),min(setup2$expr))
  numname <- num2[num2$OTU==eval(as.symbol("otusamp")),"num"]
  
  
  zz<- otusamp
  otusamp <- zz

  label <- zz
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
  file1 <- paste0(rootpath5,numname,"_",symbol,"_12132017.pdf")
  
  pdf(file1)
  p<-  ggplot(setup3, aes(x=base, y=expr, color=spec)) +
    geom_point(size=8) +
    labs(title = paste0(transcript,"-",symbol,"\n","r = ",round(as.numeric(stuff2$estimate),digits = 6),"; p-value = ", round(as.numeric(stuff2$p.value), digits = 6),"; \nrho = ",round(as.numeric(stuff$estimate),digits = 6),"; p-value = ",round(as.numeric(stuff$p.value),digits = 6)),x = paste0("CLR Transformed Abundance\n",label),y="Normalized Gene Expression") +
    facet_grid(. ~ facet,scales = "free_x",space="free_x") +
    geom_abline(data=dummy2,aes(intercept = intercept, slope = slope)) +
    scale_x_continuous(breaks = seq(0, max(setup4$base),ifelse(max(setup4$base) > 500, ifelse(max(setup4$base) < 1000, 500, 1000), ifelse(max(setup4$base) < 50, 10, 50)))) +
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
  mic1 <<- c(mic1, numname)
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

sig3.2 <- sig3[sig3$newpadj<0.0001,]
length(unique(sig3.2$OTU)) #[1] 48
length(unique(sig3.2$ENSG)) #[1] 44
sig4 <- merge(sig3.2,def,by="OTU")  
setwd(rootpath)

library(reshape2)
library(plyr)


###Make heatmaps for the correlation data


sensg <- unique(sig4$ENSG)
sotu <-unique(sig4$OTU)
# gg1 <- sensg
# aa1 <- sotu
gg1 <- as.character(sensg)
aa1 <- as.character(sotu)

#Use VST transformed values for gene expression data
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


taxgenlist <-list()
taxplist <- list()
for(i in 1:length(aa1)){
  taxgen <- vector()
  taxpp <-vector()
  for(j in 1:length(bb2)){
    pp <- cor.test(gen6[,aa1[i]], gen6[,bb2[j]], method = "spearman")
    taxgen <- c(taxgen, pp$estimate)
    taxpp <- c(taxpp, pp$p.value)
    
  }
  taxgenlist[[i]] <- taxgen
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
colnames(taxacor) <- bb2$symbol
rownames(taxacor) <- aa1

colnames(taxapval) <-bb2$symbol
rownames(taxapval) <- aa1

write.table(taxacor, file = "PathwaysCorrelations.txt", sep = "\t")

my_palette <- colorRampPalette(c("blue","white","red"))(n = 25)

# Plot heatmap with heatmap.2
#setwd(file.path(rootpath, "Plots"))

Heatmap(taxacor, 
        col=col1,
        top_annotation_height = unit(1, "cm"), 
        km = 1, 
        row_names_gp = gpar(fontsize = 3), 
        column_names_gp = gpar(fontsize=3), width = unit(6, "cm"))


###Now plot the complex heatmap with annotations for which genes are species specific
gs1 <- as.data.frame(colnames(taxacor))

setwd("/Volumes/GoogleDrive/My Drive/Primates_Coculture/Primates/MyAnalysis_v02/Results_2019-02-04/DE_Model2")
Modcombine1 <- read.table("Modcombine_2.csv", sep = ",")

gs2 <- Modcombine1

gs2$symbol <- mapIds(org.Hs.eg.db,  
                     keys=as.character(rownames(gs2)),
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")

toDrop <- c("ENSG00000147996")


gs3 <- gs2[gs2$symbol%in% colnames(taxacor),]
gs3 <- na.omit(gs3)
gs3 <- gs3[ !(rownames(gs3) %in% toDrop), ] 
gs3$Human <- factor(gs3$Human)
gs3$Orangutan <- factor(gs3$Orangutan)
gs3$Gorilla <- factor(gs3$Gorilla)
gs3$Chimp <- factor(gs3$Chimp)
rownames(gs3) <- gs3$symbol

gs4 <- gs3[,1:4]
gs4 <-gs4[order(match(rownames(gs4),colnames(taxacor))),]
col_ha <- HeatmapAnnotation(df = gs4, which="col",
                            col = list(Chimp = c("1"="#F8766D","0"="gray90"),
                                       Gorilla = c("1" = "#7CAE00","0"="gray90"),
                                       Human = c("1" = "#00BFC4" , "0"="gray90"),
                                       Orangutan =c("1" = "#C77Cff", "0"="gray90")),
                            show_legend = TRUE)

draw(col_ha, 1:36) 
ha_rot_cn = HeatmapAnnotation(text = anno_text(colnames(taxacor), rot = 45, just = "right", offset = unit(0, "mm")))



#To draw the barplots on the side for the genera
#I'm going to use gen1.11 because I have removed all of the non UNCULTURED samples and I haven't transformed any of the abundance values
#get averages of each microbial species for each primate species



allchimps <- c("100UNCULTURED", "400UNCULTURED", "300UNCULTURED", "500UNCULTURED")
allgorillas <- c("2UNCULTURED", "4UNCULTURED","14UNCULTURED","26UNCULTURED","40UNCULTURED","19UNCULTURED","24UNCULTURED")
allhumans <- c("1UNCULTURED","65UNCULTURED","110UNCULTURED","111UNCULTURED")
allorangutans <- c("8UNCULTURED","18UNCULTURED","36UNCULTURED","28UNCULTURED")

microc <-gen1.11[rownames(gen1.11) %in% allchimps,]
microg <- gen1.11[rownames(gen1.11) %in% allgorillas,]
microh <- gen1.11[rownames(gen1.11) %in% allhumans,]
microo <- gen1.11[rownames(gen1.11) %in% allorangutans,]

microc1 <- microc[,colnames(microc) %in% rownames(taxacor)]
microg1 <- microg[,colnames(microg) %in% rownames(taxacor)]
microh1 <- microh[,colnames(microh) %in% rownames(taxacor)]
microo1 <- microo[,colnames(microo) %in% rownames(taxacor)]



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

rownames(taxacor) <- unlist(lapply(strsplit(rownames(taxacor), ':', fixed = TRUE), '[', 2))
rownames(taxacor) <- unlist(lapply(strsplit(rownames(taxacor), '(', fixed = TRUE), '[', 1))

rownames(taxapval) <- unlist(lapply(strsplit(rownames(taxapval), ':', fixed = TRUE), '[', 2))
rownames(taxapval) <- unlist(lapply(strsplit(rownames(taxapval), '(', fixed = TRUE), '[', 1))


colnames(microc1) <- unlist(lapply(strsplit(colnames(microc1), ':', fixed = TRUE), '[', 2))
colnames(microc1) <- unlist(lapply(strsplit(colnames(microc1), '(', fixed = TRUE), '[', 1))

colnames(microg1) <- unlist(lapply(strsplit(colnames(microg1), ':', fixed = TRUE), '[', 2))
colnames(microg1) <- unlist(lapply(strsplit(colnames(microg1), '(', fixed = TRUE), '[', 1))

colnames(microh1) <- unlist(lapply(strsplit(colnames(microh1), ':', fixed = TRUE), '[', 2))
colnames(microh1) <- unlist(lapply(strsplit(colnames(microh1), '(', fixed = TRUE), '[', 1))

colnames(microo1) <- unlist(lapply(strsplit(colnames(microo1), ':', fixed = TRUE), '[', 2))
colnames(microo1) <- unlist(lapply(strsplit(colnames(microo1), '(', fixed = TRUE), '[', 1))



rownames(micall1) <- rownames(taxacor)

microc1 <- microc1[,order(match(colnames(microc1),rownames(taxacor)))]
microg1 <- microg1[,order(match(colnames(microg1),rownames(taxacor)))]
microh1 <- microh1[,order(match(colnames(microh1),rownames(taxacor)))]
microo1 <- microo1[,order(match(colnames(microo1),rownames(taxacor)))]

colnames(microc1) <- rownames(taxacor)
colnames(microg1) <- rownames(taxacor)
colnames(microh1) <- rownames(taxacor)
colnames(microo1) <- rownames(taxacor)

row_ha <- rowAnnotation(micall1 = row_anno_barplot(micall1, 
                                                   axis = TRUE, 
                                                   axis_side = "bottom",
                                                   gp = gpar(fill = c(Chimp = "#F8766D", Gorilla ="#7CAE00",Human = "#00BFC4", Orangutan  = "#C77Cff"))), 
                        width = unit(1.5, "cm"))

draw(row_ha, 1:10)
colnames(taxacor) <-bb2$symbol
##Set colors for heatmap correlations
col1 = colorRampPalette(c( "#999999", "#999999","#ffffff","#ef8a62","#ef8a62"))(n = 299)

pdf("CompHeatmap_PathwaysGenes.pdf")
p <- Heatmap(taxacor, 
             top_annotation = col_ha, 
             col=col1,
             top_annotation_height = unit(1, "cm"), 
             km = 1, 
             row_names_gp = gpar(fontsize = 3), 
             column_names_gp = gpar(fontsize=3), width = unit(11, "cm")) + row_ha
p <- draw(p, show_heatmap_legend=FALSE)
print(p)
dev.off()

#Check to see what the heatmap looks like
q <- Heatmap(taxacor, top_annotation = col_ha, top_annotation_height = unit(1, "cm"), km = 1, row_names_gp = gpar(fontsize = 5.5), column_names_gp = gpar(fontsize=5.5)) + row_ha
draw(q, show_heatmap_legend=FALSE)




##Make barplot for chimp abundances


ha = HeatmapAnnotation(points = row_anno)
draw(ha, 1:10)

row_hc <- HeatmapAnnotation(boxplot = row_anno_boxplot(t(microc1), outline = FALSE,
                                                       axis = TRUE, 
                                                       axis_side = "bottom", ylim = c(0,14),
                                                       gp = gpar(fill = c(Chimp = "#F8766D"))), 
                            width = unit(2, "cm"),which = "row" )

draw(row_hc, 1:5)

#Make barplot for gorilla abundances
row_hg <- HeatmapAnnotation(boxplot = row_anno_boxplot(t(microg1), outline = FALSE,
                                                       axis = TRUE, 
                                                       axis_side = "bottom", ylim=c(0,14),
                                                       gp = gpar(fill = c(Gorilla ="#7CAE00"))), 
                            width = unit(2, "cm"),which = "row")

draw(row_hg, 1:59)

#Make barplot for human abundances
row_hh <- HeatmapAnnotation(boxplot = row_anno_boxplot(t(microh1), outline=FALSE,
                                                       axis = TRUE, 
                                                       axis_side = "bottom", ylim=c(0,14),
                                                       gp = gpar(fill = c(Human = "#00BFC4"))), 
                            width = unit(2, "cm"), which="row")

draw(row_hh, 1:10)

#Make barplot for orangutan abundances
row_ho <- HeatmapAnnotation(boxplot = row_anno_boxplot(t(microo1), outline = FALSE,
                                                       axis = TRUE, 
                                                       axis_side = "bottom", ylim=c(0,14),
                                                       gp = gpar(fill = c(Orangutan  = "#C77Cff"))), 
                            width = unit(2, "cm"), which="row")

draw(row_ho, 1:10)

systempath <- "/Volumes/GoogleDrive/My\ Drive/Primates_Coculture/Primates/MyAnalysis_v02"
setwd(systempath)
dir.create(file.path(systempath, format(Sys.time(), "Fig4B_%F")))
rootpath <- paste0(file.path(systempath, format(Sys.time(), "Fig4B_%F")))

setwd(rootpath)
pdf("CompHeatmap_sig.pdf")
p <- Heatmap(taxacor, 
             col=col1,
             top_annotation = col_ha, 
             top_annotation_height = unit(1, "cm"), 
             km = 1, 
             row_names_gp = gpar(fontsize = 10, fontface="italic", base_family="Arial"), 
             column_names_gp = gpar(fontsize=3), width = unit(12.1, "cm")) + row_hc + row_hg+ row_hh + row_ho
p <- draw(p, show_heatmap_legend=FALSE)
print(p)
dev.off()

##Draw heatmap without barplots
setwd(rootpath)
pdf("CompHeatmap_nobars.pdf")
p <- Heatmap(taxacor, 
             col=col1,
             top_annotation = col_ha, 
             top_annotation_height = unit(1, "cm"), 
             km = 1, 
             row_names_gp = gpar(fontsize = 10), 
             column_names_gp = gpar(fontsize=3), width = unit(10.1, "cm"))
p <- draw(p, show_heatmap_legend=FALSE)
print(p)
dev.off()

#Heatmap with only chimp bars
pdf("CompHeatmap_chimpbars.pdf")
p <- Heatmap(taxacor, 
             col=col1,
             top_annotation = col_ha, 
             top_annotation_height = unit(1, "cm"), 
             km = 1, 
             row_names_gp = gpar(fontsize = 4, fontface="italic", base_family="Arial"), 
             column_names_gp = gpar(fontsize=3), width = unit(8, "cm")) + row_hc 
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
             column_names_gp = gpar(fontsize=3), width = unit(8, "cm")) + row_hg
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
             column_names_gp = gpar(fontsize=3), width = unit(8, "cm")) + row_hh
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
             column_names_gp = gpar(fontsize=3), width = unit(8, "cm")) + row_ho
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

pvals <- melt(taxapval)
rhos <- melt(taxacor)


tabs1 <- pvals %>% right_join(rhos, by=c("Var1","Var2"))

colnames(tabs1) <-c("Pathway", "Gene", "P-value", "Correlation")
tabs1$Qvalue <- p.adjust(tabs1$`P-value`, method = "BH")

write.table(tabs1,file="SupplementalTable6_CorrPathways.txt", sep = "\t")


sigprint <- sig4

sigprint$symbol <- mapIds(org.Hs.eg.db,  
                          keys=as.character(sigprint$ENSG),
                          column="SYMBOL",
                          keytype="ENSEMBL",
                          multiVals="first")

#

write.table(sigprint, file = "MicrobesPathsGenesModel.txt", sep = "\t")



