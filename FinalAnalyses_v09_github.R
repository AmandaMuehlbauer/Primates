#This is a curated version of my analyses from primates_RNA_16s_v0X.Rmd
# R version 3.3.3 (2017-03-06)
# Platform: x86_64-apple-darwin13.4.0 (64-bit)
# Running under: macOS  10.13.6
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# attached base packages:
#   [1] grid      parallel  stats4    stats     graphics  grDevices utils     datasets 
# [9] methods   base     
# 
# other attached packages:
#   [1] scales_1.0.0               gridExtra_2.3              gplots_3.0.1              
# [4] Hmisc_4.1-1                Formula_1.2-3              survival_2.40-1           
# [7] tidyr_0.8.1                UpSetR_1.3.3               biomaRt_2.30.0            
# [10] hash_2.2.6                 lattice_0.20-34            org.Hs.eg.db_3.4.0        
# [13] AnnotationDbi_1.36.2       vsn_3.42.3                 dplyr_0.7.6               
# [16] PoiClaClu_1.0.2            RColorBrewer_1.1-2         ggplot2_3.0.0             
# [19] DESeq2_1.14.1              parathyroidSE_1.12.0       SummarizedExperiment_1.4.0
# [22] Biobase_2.34.0             GenomicRanges_1.26.4       GenomeInfoDb_1.10.3       
# [25] IRanges_2.8.2              S4Vectors_0.12.2           BiocGenerics_0.20.0       
# [28] hexbin_1.27.1             
# 
# loaded via a namespace (and not attached):
#   [1] bitops_1.0-6          bit64_0.9-7           tools_3.3.3           backports_1.1.2      
# [5] R6_2.2.2              affyio_1.44.0         rpart_4.1-10          KernSmooth_2.23-15   
# [9] DBI_1.0.0             lazyeval_0.2.1        colorspace_1.3-2      nnet_7.3-12          
# [13] withr_2.1.2           tidyselect_0.2.4      bit_1.1-14            preprocessCore_1.36.0
# [17] htmlTable_1.12        labeling_0.3          caTools_1.17.1.1      checkmate_1.8.5      
# [21] genefilter_1.56.0     affy_1.52.0           stringr_1.3.1         digest_0.6.17        
# [25] foreign_0.8-67        XVector_0.14.1        base64enc_0.1-3       pkgconfig_2.0.2      
# [29] htmltools_0.3.6       limma_3.30.13         htmlwidgets_1.2       rlang_0.2.2          
# [33] rstudioapi_0.7        RSQLite_2.1.1         BiocInstaller_1.24.0  bindr_0.1.1          
# [37] BiocParallel_1.8.2    gtools_3.8.1          acepack_1.4.1         RCurl_1.95-4.11      
# [41] magrittr_1.5          Matrix_1.2-8          Rcpp_0.12.18          munsell_0.5.0        
# [45] stringi_1.2.4         zlibbioc_1.20.0       plyr_1.8.4            blob_1.1.1           
# [49] gdata_2.18.0          crayon_1.3.4          splines_3.3.3         annotate_1.52.1      
# [53] locfit_1.5-9.1        knitr_1.20            pillar_1.3.0          geneplotter_1.52.0   
# [57] XML_3.98-1.16         glue_1.3.0            latticeExtra_0.6-28   data.table_1.11.6    
# [61] gtable_0.2.0          purrr_0.2.5           assertthat_0.2.0      xtable_1.8-3         
# [65] tibble_1.4.2          memoise_1.1.0         bindrcpp_0.2.2        cluster_2.0.5        

###Make sure that the correct version of DESeq2 is being used

library("hexbin")

library("DESeq2")
library ('lmtest')
library('ggplot2')

library("pheatmap")
library("RColorBrewer")
library("PoiClaClu")
library("dplyr")
#library("vsn")
library("AnnotationDbi")
library("org.Hs.eg.db")
#source("https://bioconductor.org/biocLite.R")
#biocLite("goseq")
library("lattice")
library("hash")
library("biomaRt")
library("Cormotif")
library("UpSetR")
library("tidyr")
##install.packages("Hmisc", dependencies=T)
library("Hmisc")
library("gplots")
library("gridExtra")
library("grid")



#set directory and load data
setwd("/Volumes/GoogleDrive/My\ Drive/Primates_Coculture/RNAseq_Gopher/R_workspace")
load("subread_DESeq2-data.RDATA")


# 



############################This next bit was used for filtering for protein coding genes###############
# # Filter for protein coding genes
# ########### Only keep protein-coding genes in the read.counts matrix ############## 
# #First retrieve gene ids
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
##Save the data filtered for protein coding sequences in case the ensembl database changes
write.table(data, file = "Filtered_Data.txt", sep = "\t")
########################end section used for filtering protein coding genes#############################
data_filt <- data[rowSums(data)>1,]
write.table(data_filt, file = "Filtered_Genes_Rep.txt", sep="\t")

#MAke a new directory for the results
systempath <- "/Volumes/GoogleDrive/My\ Drive/Primates_Coculture/Primates/MyAnalysis_v02"
setwd(systempath)
dir.create(file.path(systempath, format(Sys.time(), "Results_%F")))

rootpath <- paste0(file.path(systempath, format(Sys.time(), "Results_%F")))
setwd(paste0(file.path(systempath, format(Sys.time(), "Results_%F"))))
# dim(data)
# [1] 19714    32

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


#remove Gorilla 6 (Schroeder replicate). This is the one that looks like an orangutan. This corresponds to HT2 for the RNA seq data
meta <- meta[!(meta$Microbiome %in% c("6")),]
data <- subset(data, select = -CP4.HT2 )
# #Model 1
# ###DE genes over all samples
# First, work on DE genes over all of the samples compared to the controls without microbiome 
# This will give me a list of genes that are expressed in microbiome co-culture regardless of species at higher levels than the no-microbiome control 
# 
# From the DESeq2 manual: In order to benefit from the default settings of the package, you should put the
# variable of interest at the end of the formula and make sure the control level is the first level.

#meta$Species_All <- ifelse(meta$Species_All=="Primate",1,0)

dds_all <- DESeqDataSetFromMatrix(countData = data, colData = meta, design = ~ ExperimentPlate + Species_All)
dds_all <- dds_all[ rowSums(counts(dds_all)) > 1, ]
dds_all <- DESeq(dds_all)

res_all <- results(dds_all, contrast = c("Species_All", "Primate", "AA"))

res_all$symbol <- mapIds(org.Hs.eg.db,
                         keys=row.names(res_all),
                         column="SYMBOL",
                         keytype="ENSEMBL",
                         multiVals="first")
res_all$entrez <- mapIds(org.Hs.eg.db,
                         keys=row.names(res_all),
                         column="ENTREZID",
                         keytype="ENSEMBL",
                         multiVals="first")





summary(res_all)
# out of 17860 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 1940, 11%
# LFC < 0 (down)     : 1690, 9.5%
# outliers [1]       : 8, 0.045%
# low counts [2]     : 3807, 21%
# (mean count < 4)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

sum(res_all$padj < 0.1, na.rm=TRUE)
# [1] 3676

#First number is number of p-values less than 0.1, second number is number of p-values reported. 
#Third number is all genes with an adjusted p-value less than 0.25

sum(res_all$pvalue < 0.1, na.rm=TRUE)
# [1] 5679
sum(!is.na(res_all$pvalue))
# [1] 17724
sum(res_all$padj < 0.1, na.rm=TRUE)
# [1] 3676

resSig_allprimates <- subset(res_all, padj < 0.1)
#Number of genes up regulated
res_all_up <- resSig_allprimates[resSig_allprimates$log2FoldChange >0.25,]
nrow(res_all_up)
# [1] 1958


#Number of genes down regulated
res_all_down <- resSig_allprimates[resSig_allprimates$log2FoldChange < -0.25,]
nrow(res_all_down)
# [1] 1718

###Make table with different cutoffs:

up1 <- c(nrow(resSig_allprimates[resSig_allprimates$log2FoldChange >0.25,]),
         nrow(resSig_allprimates[resSig_allprimates$log2FoldChange >0.5,]),
         nrow(resSig_allprimates[resSig_allprimates$log2FoldChange >0.75,]),
         nrow(resSig_allprimates[resSig_allprimates$log2FoldChange >1,]))
down1 <- c(nrow(resSig_allprimates[resSig_allprimates$log2FoldChange < (-0.25),]),
           nrow(resSig_allprimates[resSig_allprimates$log2FoldChange < (-0.5),]),
           nrow(resSig_allprimates[resSig_allprimates$log2FoldChange < (-0.75),]),
           nrow(resSig_allprimates[resSig_allprimates$log2FoldChange < (-1),]))
mod1 <-rbind(up1, down1)

tot1 <- colSums(mod1)
mod1 <-rbind(mod1, tot1)

colnames(mod1) <- c("LFC_0.25", "LFC_0.5", "LFC_0.75", "LFC_1.0")
rownames(mod1) <- c("Up_DE", "Down_DE", "Total_DE")

#setwd("~/Google Drive/Primates_Coculture/Primates/MyAnalysis/Results_v09/DE_cutoffTables")
write.table(mod1, file = "DE_CutoffTable_Model1.txt", sep = "\t")
write.table(res_all, file = "DE_allgenesModel1.txt", sep= "\t")


dir.create(file.path(rootpath, "DE_Model1" ))
setwd(paste0(file.path(rootpath, "DE_Model1")))

DE_allprimates <- resSig_allprimates[ order(resSig_allprimates$log2FoldChange, decreasing = TRUE), ]

write.table(resSig_allprimates, file = "Model1_v09.txt", sep = "\t")
write.table(resSig_allprimates, file = "SupplementalTable1_TreatmentControl.txt", sep = "\t")
write.table(res_all, file = "Allgenes.txt", sep = "\t")

resSig_allprimates_sub <- cbind(c(resSig_allprimates$symbol), c(resSig_allprimates$log2FoldChange))
write.table(resSig_allprimates_sub, file = "TreatmentControl_GT.txt", sep="\t")

test3 <- varianceStabilizingTransformation(dds_all,blind=TRUE)
test4 <- as.data.frame(assay(test3))
test4 <- t(test4)



ge1 <- row.names(DE_allprimates)



###Heatmap for model 1
#shows speciifc samples

nm1 <- test4[,ge1]
nm1<- t(nm1)
for (cl in 1:length(unique(as.numeric(meta$ExperimentPlate)))) {
  avg <- apply(nm1[,as.numeric(meta$ExperimentPlate)==cl], 1, mean) #figure out means for each plate
  nm1[,as.numeric(meta$ExperimentPlate)==cl] <-
    nm1[,as.numeric(meta$ExperimentPlate)==cl] - avg
}

#nm1 <-test4[,ge1]
#nm1<- t(nm1)
nm2 <- nm1[,1:19]


colnames(nm2) <- c("Gorilla 1", "Gorilla 2", "Gorilla 3", "Gorilla 4", "Gorilla 5",   "Orangutan 1", "Orangutan 2", "Orangutan 3","Orangutan 4","Gorilla 6", "Chimp 1", "Chimp 2", "Chimp 3", "Gorilla 7", "Chimp 4", "Human 1", "Human 2", "Human 3", "Human 4" )

nm3 <- nm2[ ,order(colnames(nm2))]
#nm3 <-head(nm3)
nm3 <- data.matrix(nm3)


#colors = c(seq(-3,-1,length=100),seq(-1,0,length=100),seq(0,0,length=100),seq(0,1,length=100),seq(1,3,length=100))
my_palette <- colorRampPalette(c("blue","blue", "white","red", "red"))(n = 299)



nm4 <- nm3[,]


###################################
#Model 2

###DE genes each species compared to the control

################################Human compared to control with no microbiome

dds_hc <- DESeqDataSetFromMatrix(countData = data, colData = meta, design = ~ ExperimentPlate + Species )
dds_hc <- dds_hc[ rowSums(counts(dds_hc)) > 1, ]
dds_hc <- DESeq(dds_hc)


res_hc <- results(dds_hc, contrast = c("Species", "Human", "AA"))

res_hc$symbol <- mapIds(org.Hs.eg.db,
                        keys=row.names(res_hc),
                        column="SYMBOL",
                        keytype="ENSEMBL",
                        multiVals="first")
res_hc$entrez <- mapIds(org.Hs.eg.db,
                        keys=row.names(res_hc),
                        column="ENTREZID",
                        keytype="ENSEMBL",
                        multiVals="first")

summary(res_hc)


# out of 17860 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 964, 5.4%
# LFC < 0 (down)     : 482, 2.7%
# outliers [1]       : 6, 0.034%
# low counts [2]     : 4500, 25%
# (mean count < 8)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results


sum(res_hc$pvalue < 0.1, na.rm=TRUE)
# [1] 3847
sum(!is.na(res_hc$pvalue))
# [1] 17854
sum(res_hc$padj < 0.1, na.rm=TRUE)
# [1] 1446

resSig_humancontrol <- subset(res_hc, padj < 0.1)
#head(resSig_h[ order(resSig_h$log2FoldChange), ])

#Number of genes up regulated
hum_up <- resSig_humancontrol[resSig_humancontrol$log2FoldChange >0.25,]
nrow(hum_up)
# upregulated, no threshold : [1] 964
#upregulated with threshold of 0.25: [1] 756

#Number of genes down regulated
hum_down <- resSig_humancontrol[resSig_humancontrol$log2FoldChange <=-0.25,]
nrow(hum_down)
# downregulated, no threshold : [1] 482
# downregulated with threshold of -0.25: [1] 344

DE_humancontrol <- resSig_humancontrol[ order(resSig_humancontrol$log2FoldChange, decreasing = TRUE), ]
dir.create(file.path(rootpath, "DE_Model2" ))
setwd(file.path(rootpath, "DE_Model2"))

write.table(resSig_humancontrol, file = "DE_Human_comparedtoNoMBControl_v09.txt", sep = "\t")
#write.table(resSig_humancontrol, file = "SupplementalTable21_HumanLFC.txt", sep = "\t")

###################Chimpanzee compared to no microbiome control
dds_cc <- DESeqDataSetFromMatrix(countData = data, colData = meta, design = ~ ExperimentPlate + Species )
dds_cc <- dds_cc[ rowSums(counts(dds_cc)) > 1, ]
dds_cc <- DESeq(dds_cc)


res_cc <- results(dds_cc, contrast = c("Species", "Chimp", "AA"))

res_cc$symbol <- mapIds(org.Hs.eg.db,
                        keys=row.names(res_cc),
                        column="SYMBOL",
                        keytype="ENSEMBL",
                        multiVals="first")
res_cc$entrez <- mapIds(org.Hs.eg.db,
                        keys=row.names(res_cc),
                        column="ENTREZID",
                        keytype="ENSEMBL",
                        multiVals="first")

summary(res_cc)

sum(res_cc$padj < 0.1, na.rm=TRUE)
# [1] 1341
sum(res_cc$pvalue < 0.1, na.rm=TRUE)
# [1] 3996
sum(!is.na(res_cc$pvalue))
# [1] 17854

resSig_chimpcontrol <- subset(res_cc, padj < 0.1)

#Number of genes up regulated
chimp_up <- resSig_chimpcontrol[resSig_chimpcontrol$log2FoldChange >0.25,]
nrow(chimp_up)
# [1] 788

#Number of genes down regulated
chimp_down<- resSig_chimpcontrol[resSig_chimpcontrol$log2FoldChange <= -0.25,]
nrow(chimp_down)
# [1] 553




DE_chimpcontrol <- resSig_chimpcontrol[ order(resSig_chimpcontrol$log2FoldChange, decreasing = TRUE), ]


write.table(resSig_chimpcontrol, file = "DE_Chimp_comparedtoNoMBControl_v09.txt", sep = "\t")

####Gorilla compared to no microbiome control
dds_gc <- DESeqDataSetFromMatrix(countData = data, colData = meta, design = ~ ExperimentPlate + Species )
dds_gc <- dds_gc[ rowSums(counts(dds_gc)) > 1, ]
dds_gc <- DESeq(dds_gc)


res_gc <- results(dds_gc, contrast = c("Species", "Gorilla", "AA"))

res_gc$symbol <- mapIds(org.Hs.eg.db,
                        keys=row.names(res_gc),
                        column="SYMBOL",
                        keytype="ENSEMBL",
                        multiVals="first")
res_gc$entrez <- mapIds(org.Hs.eg.db,
                        keys=row.names(res_gc),
                        column="ENTREZID",
                        keytype="ENSEMBL",
                        multiVals="first")

summary(res_gc)
# out of 17860 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 1389, 7.8%
# LFC < 0 (down)     : 1024, 5.7%
# outliers [1]       : 6, 0.034%
# low counts [2]     : 3808, 21%
# (mean count < 4)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

sum(res_gc$padj < 0.1, na.rm=TRUE)
#[1] 2413
sum(res_gc$pvalue < 0.1, na.rm=TRUE)
#[1] 4735
sum(!is.na(res_gc$pvalue))
#[1] 17854



resSig_gorillacontrol <- subset(res_gc, padj < 0.1)
#head(resSig_h[ order(resSig_h$log2FoldChange), ])

#Number of genes up regulated
gor_up <- resSig_gorillacontrol[resSig_gorillacontrol$log2FoldChange > 0.25,]
nrow(gor_up)
# [1] 889
#Number of genes down regulated
gor_down<- resSig_gorillacontrol[resSig_gorillacontrol$log2FoldChange <= -0.25,]
nrow(gor_down)
# [1] 509


DE_gorillacontrol <- resSig_gorillacontrol[ order(resSig_gorillacontrol$log2FoldChange, decreasing = TRUE), ]

setwd(file.path(rootpath, "DE_Model2"))
write.table(resSig_gorillacontrol, file = "DE_Gorilla_comparedtoNoMBControl_v09.txt", sep = "\t")
#write.table(resSig_gorillacontrol, file = "SupplementalTable23_GorillaLFC.txt", sep = "\t")

####Orangutan compared to no microbiome control
dds_oc <- DESeqDataSetFromMatrix(countData = data, colData = meta, design = ~ ExperimentPlate + Species)
dds_oc <- dds_oc[ rowSums(counts(dds_oc)) > 1, ]
dds_oc <- DESeq(dds_oc)

res_oc <- results(dds_oc, contrast = c("Species", "Orangutan", "AA"))

res_oc$symbol <- mapIds(org.Hs.eg.db,
                        keys=row.names(res_oc),
                        column="SYMBOL",
                        keytype="ENSEMBL",
                        multiVals="first")
res_oc$entrez <- mapIds(org.Hs.eg.db,
                        keys=row.names(res_oc),
                        column="ENTREZID",
                        keytype="ENSEMBL",
                        multiVals="first")

summary(res_oc)
# out of 17860 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 849, 4.8%
# LFC < 0 (down)     : 408, 2.3%
# outliers [1]       : 6, 0.034%
# low counts [2]     : 4500, 25%
# (mean count < 8)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results


sum(res_oc$padj < 0.1, na.rm=TRUE)
#[1] 1257
sum(res_oc$pvalue < 0.1, na.rm=TRUE)
#[1] 3826
sum(!is.na(res_oc$pvalue))
#[1] 17854


resSig_orangutancontrol <- subset(res_oc, padj < 0.1)
#head(resSig_h[ order(resSig_h$log2FoldChange), ])

#Number of genes up regulated
ora_up <- resSig_orangutancontrol[resSig_orangutancontrol$log2FoldChange >0.25,]
nrow(ora_up)
# [1] 849
#Number of genes down regulated
ora_down<- resSig_orangutancontrol[resSig_orangutancontrol$log2FoldChange <= -0.25,]
nrow(ora_down)
# [1] 408

DE_orangutancontrol <- resSig_orangutancontrol[ order(resSig_orangutancontrol$log2FoldChange, decreasing = TRUE), ]

setwd(file.path(rootpath, "DE_Model2"))
write.table(resSig_orangutancontrol, file = "DE_Orangutan_comparedtoNoMBControl_v09.txt", sep = "\t")
#write.table(resSig_orangutancontrol, file = "SupplementalTable24_OrangutanLFC.txt", sep = "\t")

####Write results from model 2 to files
setwd(file.path(rootpath, "DE_Model2"))
write.table(res_oc, file = "SupplementalTable24_OrangutanLFC.txt", sep = "\t")
write.table(res_hc, file = "SupplementalTable21_HumanLFC.txt", sep = "\t")
write.table(res_cc, file = "SupplementalTable22_ChimpLFC.txt", sep = "\t")
write.table(res_gc, file = "SupplementalTable23_GorillaLFC.txt", sep = "\t")

##Want to make plots similar to species-specific mean-SE for conserved genes from Model 1


dir.create("Plots")
dir.create("Plots/Subset")
file1 <-file.path(rootpath, "DE_Model2/Plots/Subset/")
setwd(file1)



# res_alltop<-subset(res_alltop,!(is.na(res_alltop["symbol"])))
# rownames(res_alltop) <- res_alltop$symbol
# list_genes <- res_alltop$symbol

#########################################################################################
###This plots the top upregulated and bottom downregulated genes for Model 1 (treatment vs control)

list_genes <- c(head(rownames(DE_allprimates), 200), tail(rownames(DE_allprimates), 200))
ac1 <- data.frame()
title <- character()
p<-list()
#plots <- list() 
for(a in list_genes)
{
  
  ac1 <- data.frame(c("Human", "Chimp", "Gorilla", "Orangutan"), 
                    c(res_hc[a, "log2FoldChange"], res_cc[a, "log2FoldChange"], res_gc[a, "log2FoldChange"], res_oc[a, "log2FoldChange"]),
                    c(res_hc[a, "lfcSE"], res_cc[a, "lfcSE"], res_gc[a, "lfcSE"], res_oc[a, "lfcSE"]))
  
  colnames(ac1) <- c("Species", "LFC", "SE")

  
  p[[a]]<- ggplot(ac1, aes(x=Species, y=LFC, group = Species, color=Species)) +
     geom_pointrange(aes(ymin=LFC-SE, ymax=LFC+SE))+labs(title= a) + 
     geom_errorbar(aes(ymin=LFC-SE, ymax=LFC+SE), width=.2) +
     geom_point(shape = 1, size = 3) + 
     theme_bw() + 
     theme(legend.position="none", axis.title.x = element_text( size=10), axis.text.x  = element_text( size=10), axis.title.y = element_text(size=10), axis.text.y = element_text(size=10), plot.title = element_text(size=10)) 
  }


m1 <-marrangeGrob(grobs = p, nrow = 2, ncol = 2)
ggsave(paste0(file.path(rootpath, "/DE_Model2/Plots/Subset/ConservedGenes.pdf")), m1)

setwd(rootpath)
dir.create("SavedModels")
setwd(file.path(rootpath, "/SavedModels"))
saveRDS(res_hc, "res_hc.RDS")
saveRDS(res_gc, "res_gc.RDS")
saveRDS(res_cc, "res_cc.RDS")
saveRDS(res_oc, "res_oc.RDS")
####End of making plots 


####################################################################################################################
###This is to make single plots for the paper. 

setwd(rootpath)
dir.create("PlotsPaper")
setwd(file.path(rootpath, "/PlotsPaper"))
a <- "ENSG00000197641" #set which gene to plot
ac1 <- data.frame(c("Human", "Chimp", "Gorilla", "Orangutan"), 
                  c(res_hc[a, "log2FoldChange"], res_cc[a, "log2FoldChange"], res_gc[a, "log2FoldChange"], res_oc[a, "log2FoldChange"]),
                  c(res_hc[a, "lfcSE"], res_cc[a, "lfcSE"], res_gc[a, "lfcSE"], res_oc[a, "lfcSE"]))

colnames(ac1) <- c("Species", "LFC", "SE")

pdf("Upreg.pdf")
p <-ggplot(ac1, aes(x=Species, y=LFC, group = Species, color=Species)) +
  geom_pointrange(aes(ymin=LFC-SE, ymax=LFC+SE), size=2)+labs(title=res_hc[a,"symbol"]) +
  geom_errorbar(aes(ymin=LFC-SE, ymax=LFC+SE), width=.2, size=2) +
  geom_point(shape = 1, size = 3) +
  theme_bw() +
  theme(legend.position="none", axis.title.x = element_text( size=10), axis.text.x  = element_text( size=10), axis.title.y = element_text(size=10), axis.text.y = element_text(size=10), plot.title = element_text(size=10)) +
  ylim(0,1.2) + geom_hline(yintercept=0, size=.2)
print(p)
dev.off()


a <- "ENSG00000165264" #set which gene to plot
ac1 <- data.frame(c("Human", "Chimp", "Gorilla", "Orangutan"), 
                  c(res_hc[a, "log2FoldChange"], res_cc[a, "log2FoldChange"], res_gc[a, "log2FoldChange"], res_oc[a, "log2FoldChange"]),
                  c(res_hc[a, "lfcSE"], res_cc[a, "lfcSE"], res_gc[a, "lfcSE"], res_oc[a, "lfcSE"]))

colnames(ac1) <- c("Species", "LFC", "SE")
pdf("Downreg.pdf")
p <-ggplot(ac1, aes(x=Species, y=LFC, group = Species, color=Species)) +
  geom_pointrange(aes(ymin=LFC-SE, ymax=LFC+SE), size = 2)+labs(title=res_hc[a,"symbol"]) +
  geom_errorbar(aes(ymin=LFC-SE, ymax=LFC+SE), width=.2,size=2) +
  geom_point(shape = 1, size = 3) +
  theme_bw() +
  theme(legend.position="none", axis.title.x = element_text( size=10), axis.text.x  = element_text( size=10), axis.title.y = element_text(size=10), axis.text.y = element_text(size=10), plot.title = element_text(size=10)) +
  ylim(-0.75,0) + geom_hline(yintercept=0, size=.2)

print(p)
dev.off()

###No response 

no_resp <- subset(res_all, padj > 0.1)
a <- "ENSG00000119509" #set which gene to plot
ac1 <- data.frame(c("Human", "Chimp", "Gorilla", "Orangutan"), 
                  c(res_hc[a, "log2FoldChange"], res_cc[a, "log2FoldChange"], res_gc[a, "log2FoldChange"], res_oc[a, "log2FoldChange"]),
                  c(res_hc[a, "lfcSE"], res_cc[a, "lfcSE"], res_gc[a, "lfcSE"], res_oc[a, "lfcSE"]))

colnames(ac1) <- c("Species", "LFC", "SE")
pdf("No_resp.pdf")
p <-ggplot(ac1, aes(x=Species, y=LFC, group = Species, color=Species)) +
  geom_pointrange(aes(ymin=LFC-SE, ymax=LFC+SE), size = 2)+labs(title=res_hc[a,"symbol"]) +
  geom_errorbar(aes(ymin=LFC-SE, ymax=LFC+SE), width=.2,size=2) +
  geom_point(shape = 1, size = 3) +
  theme_bw() +
  theme(legend.position="none", axis.title.x = element_text( size=10), axis.text.x  = element_text( size=10), axis.title.y = element_text(size=10), axis.text.y = element_text(size=10), plot.title = element_text(size=10)) +
  geom_hline(yintercept=0, size=.2) + ylim(-0.25, 0.25)

print(p)
dev.off()

p
####################################################################################################################

###Heatmaps for Model 2
res_c_d <- as.data.frame(res_cc)
res_h_d <- as.data.frame(res_hc)
res_g_d <- as.data.frame(res_gc)
res_o_d <- as.data.frame(res_oc)

res_c_d$sample <- "chimp"
res_h_d$sample <- "human"
res_g_d$sample <- "gorilla"
res_o_d$sample <- "orangutan"

res_c_d$ENSG <- rownames(res_c_d)
res_h_d$ENSG <- rownames(res_h_d)
res_g_d$ENSG <- rownames(res_g_d)
res_o_d$ENSG <- rownames(res_o_d)

res_total <-rbind(res_c_d, res_h_d, res_g_d, res_o_d)

rownames(res_total) <- seq(1:nrow(res_total))

#sum(match(rownames(res_c_d), rownames(res_h_d)))


# sum(rownames(res_c_d) ==rownames(res_h_d))
# sum(rownames(res_c_d) ==rownames(res_g_d))
# sum(rownames(res_c_d) ==rownames(res_o_d))

# res_total <- cbind(res_c_d, res_h_d, res_g_d, res_o_d)
# res_tot <- res_total[,c(7,8,2,11,20,29)]
# colnames(res_tot) <- c("symbol", "entrez", "chimp", "human","gorilla","orangutan")




setwd("/Volumes/GoogleDrive/My\ Drive/Primates_Coculture/Primates/MyAnalysis_v02/Results_2019-08-20/DE_Model7")
g_c<-read.table(file = "LRT_chimp.txt", sep= "\t")
g_g <- read.table(file = "LRT_gorilla.txt", sep= "\t")
g_o <- read.table( file = "LRT_orangutan.txt", sep= "\t")
g_h <- read.table( file = "LRT_human.txt", sep= "\t")

g_cg <- read.table( file = "LRT_chimp_gorilla.txt", sep = "\t")
g_ch<- read.table( file = "LRT_chimp_human.txt", sep = "\t")
g_co <- read.table( file = "LRT_chimp_orangutan.txt", sep = "\t")
g_gh <- read.table(file = "LRT_gorilla_human.txt", sep = "\t")
g_go <- read.table( file = "LRT_gorilla_orangutan.txt", sep = "\t")
g_ho <- read.table( file = "LRT_human_orangutan.txt", sep = "\t")


g_cgh <- read.table( file = "LRT_Triple_CGH.txt", sep= "\t")
g_cgo <- read.table(file = "LRT_Triple_CGO.txt", sep= "\t")
g_cho <- read.table(file = "LRT_Triple_CHO.txt", sep= "\t")
g_gho <- read.table( file = "LRT_Triple_GHO.txt", sep= "\t")

g_cgho <- read.table( file = "LRT_chimp_gorilla_human_orangutan.txt", sep = "\t")

g_c$chimp<- 1
g_c$gorilla <-0
g_c$human<-0
g_c$orangutan<-0

g_g$chimp<- 0
g_g$gorilla <-1
g_g$human<-0
g_g$orangutan<-0

g_h$chimp<- 0
g_h$gorilla <-0
g_h$human<-1
g_h$orangutan<-0

g_o$chimp<- 0
g_o$gorilla <-0
g_o$human<-0
g_o$orangutan<-1


##Doubles
g_cg$chimp<- 1
g_cg$gorilla <-1
g_cg$human<-0
g_cg$orangutan<-0


g_ch$chimp<- 1
g_ch$gorilla <-0
g_ch$human<-1
g_ch$orangutan<-0


g_co$chimp<- 1
g_co$gorilla <-0
g_co$human<-0
g_co$orangutan<-1


g_gh$chimp<- 0
g_gh$gorilla <-1
g_gh$human<-1
g_gh$orangutan<-0

g_go$chimp<- 0
g_go$gorilla <-1
g_go$human<-0
g_go$orangutan<-1

g_ho$chimp<- 0
g_ho$gorilla <-0
g_ho$human<-1
g_ho$orangutan<-1

#Triples
g_cgh$chimp<- 1
g_cgh$gorilla <-1
g_cgh$human<-1
g_cgh$orangutan<-0

g_cho$chimp<- 1
g_cho$gorilla <-0
g_cho$human<-1
g_cho$orangutan<-1

g_gho$chimp<- 0
g_gho$gorilla <-1
g_gho$human<-1
g_gho$orangutan<-1

g_cgo$chimp<- 1
g_cgo$gorilla <-1
g_cgo$human<-0
g_cgo$orangutan<-1

##All
g_cgho$chimp<- 1
g_cgho$gorilla <-1
g_cgho$human<-1
g_cgho$orangutan<-1

g_combine <- rbind(g_c, g_g, g_h, g_o, g_cg, g_ch, g_co, g_gh, g_go, g_ho,g_cgh, g_cho,g_gho,g_cgo,g_cgho)


###Arrange data to make the heatmap
lfc <- spread(res_total[,c(7,9,10,8,2)],sample,log2FoldChange)
colnames(lfc)[4:7] <- paste0("logFC_",colnames(lfc)[4:7])
se <- spread(res_total[,c(7,9,10,8,3)],sample,lfcSE)
colnames(se)[4:7] <- paste0("se_",colnames(se)[4:7])
padj <- spread(res_total[,c(7,9,10,8,6)],sample,padj)
colnames(padj)[4:7] <- paste0("padj_",colnames(padj)[4:7])
sep2 <- cbind(lfc,se[,c(4:7)],padj[,c(4:7)])

sep2 <- sep2[]


samples <-c("chimp", "gorilla", "human", "orangutan")

# whichGenes <- function(samp,dat) {
#   padj <- paste0("padj_",samp)
#   logFC <- paste0("logFC_",samp)
#   temp <- ifelse(!is.na(dat[padj]),ifelse(dat[padj]<0.1,ifelse(abs(dat[logFC])>0.3,1,0),0),NA)
#   #temp <- ifelse(!is.na(dat[padj]))
#   colnames(temp) <- samp
#   temp
# }

#aux4 <- lapply(samples,whichGenes,dat=sep2)
#aux5 <- do.call(cbind,aux4)
aux6 <- merge(sep2,g_combine, by.x = "ENSG", by.y = "ENSG")

aux6$count <- rowSums(aux6[,c(16:19)],na.rm=TRUE)
dir.create(file.path(rootpath, "Plots"))

setwd(file.path(rootpath, "/Plots"))

write.table(aux6, file = "aux6.txt", sep="\t")


aux7 <- aux6[aux6$count >0,c(1,2,3,4,5,6,7)]
aux8 <-aux7[,c(4:7)]
rownames(aux8) <- aux7$ENSG


####Heatmap
# 
# n_matrix <- data.matrix(aux8)
# colors = c(seq(-3,-1,length=100),seq(-1,0,length=100),seq(0,0,length=100),seq(0,1,length=100),seq(1,3,length=100))
# my_palette <- colorRampPalette(c("blue","white","red"))(n = 299)
# #plot(rep(1,10),col=my_palette, pch=19,cex=2)

dir.create(file.path(rootpath, "Plots"))
dir.create(file.path(rootpath, "Plots/Model_1"))
dir.create(file.path(rootpath, "Plots/Model_1/Heatmap"))
setwd(file.path(rootpath, "/Plots/Model_1/Heatmap"))
# heatmap.2(n_matrix, scale="none", symm=F,symkey=F,symbreaks=T, trace="none",density.info="none", col=my_palette, dendrogram="both", cexRow=0.5, cexCol = 1, margins = c(5,1), labCol = c("Human", "Gorilla", "Chimp", "Orangutan"))
# 
# pdf ("DEheatmap_SpeciesVcontrol_v10.pdf")
# heatmap.2(n_matrix, scale="none", symm=F,symkey=F,symbreaks=T, trace="none",density.info="none", col=my_palette, dendrogram="both", cexRow=0.5, cexCol = 1, margins = c(5,1), labCol = c("Human", "Gorilla", "Chimp", "Orangutan"))
# dev.off()

n_matrix <- data.matrix(aux8)
#colors = c(seq(-3,-1,length=100),seq(-1,0,length=100),seq(0,0,length=100),seq(0,1,length=100),seq(1,3,length=100))
my_palette <- colorRampPalette(c("blue", "white", "red"))(n=299)
#plot(rep(1,10),col=my_palette, pch=19,cex=2)
n_matrix2 <- as.data.frame(n_matrix)
n_matrix2$col_n <-"black"

##Before, using we used a different model for the heatmap. 
##This was to select for genes such that the lfc was >=|0.3| for at least one hominid species
# n_matrix2$col_n[abs(n_matrix2$logFC_chimp)>0.3 & abs(n_matrix2$logFC_gorilla) <0.3 & abs(n_matrix2$logFC_human )<0.3 & abs(n_matrix2$logFC_orangutan)<0.3] <- "red"
# 
# n_matrix2$col_n[abs(n_matrix2$logFC_chimp)<0.3 & abs(n_matrix2$logFC_gorilla) >0.3 & abs(n_matrix2$logFC_human )<0.3 & abs(n_matrix2$logFC_orangutan)<0.3] <- "green"
# 
# n_matrix2$col_n[abs(n_matrix2$logFC_chimp)<0.3 & abs(n_matrix2$logFC_gorilla) <0.3 & abs(n_matrix2$logFC_human )>0.3 & abs(n_matrix2$logFC_orangutan)<0.3] <- "blue"
# 
# n_matrix2$col_n[abs(n_matrix2$logFC_chimp)<0.3 & abs(n_matrix2$logFC_gorilla) <0.3 & abs(n_matrix2$logFC_human )<0.3 & abs(n_matrix2$logFC_orangutan)>0.3] <- "orange"

# 
# n_matrix2$col_n[abs(n_matrix2$logFC_chimp)>0.3 & abs(n_matrix2$logFC_gorilla) >0.3 & abs(n_matrix2$logFC_human )>0.3 & abs(n_matrix2$logFC_orangutan)>0.3] <- "grey"

# as.factor(n_matrix2$col_n)
n_matrix2 <-as.data.frame(n_matrix2)
n_matrix2 <- n_matrix2[order(n_matrix2[,"col_n"]),]

# goror<- rownames(n_matrix2[order(n_matrix2$logFC_gorilla),])
# oror<- rownames(n_matrix2[order(n_matrix2$logFC_orangutan),])
# chor<- rownames(n_matrix2[order(n_matrix2$logFC_chimp),])
# humor<- rownames(n_matrix2[order(n_matrix2$logFC_human),])

# 
#n_matrix3 <- data.matrix(n_matrix2[1:4])
n_matrix3 <- n_matrix[rownames(n_matrix2),]



rownames(n_matrix3) <- factor(rownames(n_matrix3),levels = rownames(n_matrix3))

setwd(file.path(rootpath, "/Plots/Model_1/Heatmap"))


setwd(file.path(rootpath, "DE_Model2"))
Modcombine1 <- g_combine

rownames(Modcombine1) <- Modcombine1$ENSG


aux <- Modcombine1
#aux <-aux[,2:5]
colnames(aux) <- c("ENSG","Chimp", "Gorilla", "Human", "Orangutan")
# aux <- cbind(aux$ENSG, aux$Human, aux$Orangutan, aux$Gorilla, aux$Chimp)

# colnames(aux) <- c("ENSG","Human", "Orangutan", "Gorilla", "Chimp")
# aux <-as.data.frame(aux)
rownames(aux) <- aux$ENSG
aux$X1[aux$Human == 1 & aux$Chimp == 0 & aux$Gorilla == 0 & aux$Orangutan == 0] <- 1
aux$X2[aux$Human == 0 & aux$Chimp == 1 & aux$Gorilla == 0 & aux$Orangutan == 0] <- 1
aux$X3[aux$Human == 0 & aux$Chimp == 0 & aux$Gorilla == 1 & aux$Orangutan == 0] <- 1
aux$X4[aux$Human == 0 & aux$Chimp == 0 & aux$Gorilla == 0 & aux$Orangutan == 1] <- 1

aux$X5[aux$Human == 1 & aux$Chimp == 1 & aux$Gorilla == 0 & aux$Orangutan == 0] <- 1
aux$X6[aux$Human == 1 & aux$Chimp == 0 & aux$Gorilla == 1 & aux$Orangutan == 0] <- 1
aux$X7[aux$Human == 1 & aux$Chimp == 0 & aux$Gorilla == 0 & aux$Orangutan == 1] <- 1
aux$X8[aux$Human == 0 & aux$Chimp == 1 & aux$Gorilla == 1 & aux$Orangutan == 0] <- 1

aux$X9[aux$Human == 0 & aux$Chimp == 1 & aux$Gorilla == 0 & aux$Orangutan == 1] <- 1
aux$X10[aux$Human == 0 & aux$Chimp == 0 & aux$Gorilla == 1 & aux$Orangutan == 1] <- 1

aux$X11[aux$Human == 1 & aux$Chimp == 1 & aux$Gorilla == 1 & aux$Orangutan == 0] <- 1
aux$X12[aux$Human == 1 & aux$Chimp == 1 & aux$Gorilla == 0 & aux$Orangutan == 1] <- 1
aux$X13[aux$Human == 1 & aux$Chimp == 0 & aux$Gorilla == 1 & aux$Orangutan == 1] <- 1
aux$X14[aux$Human == 0 & aux$Chimp == 1 & aux$Gorilla == 1 & aux$Orangutan == 1] <- 1
aux$X15[aux$Human == 1 & aux$Chimp == 1 & aux$Gorilla == 1 & aux$Orangutan == 1] <- 1

aux[is.na(aux)] <- 0

rownames(aux6) <- aux6$ENSG
aux$Human <- aux6[rownames(aux),"logFC_human"] 
aux$Chimp <- aux6[rownames(aux),"logFC_chimp"] 
aux$Gorilla <- aux6[rownames(aux),"logFC_gorilla"] 
aux$Orangutan <- aux6[rownames(aux),"logFC_orangutan"] 




oct <- c("X1", "X2", "X3", "X4", "X5", "X6", "X7", "X8", "X9", "X10", "X11", "X12", "X13", "X14", "X15")
oot <- c("Human", "Chimp", "Gorilla", "Orangutan", "Human_Chimp", "Human_Gorilla", "Human_Orangutan", "Chimp_Gorilla", "Chimp_Orangutan", "Gorilla_Orangutan", "Human_Chimp_Gorilla", "Human_Chimp_Orangutan", "Human_Gorilla_Orangutan", "Chimp_Gorilla_Orangutan", "Human_Chimp_Gorilla_Orangutan")

ootsimp <-c("Human", "Chimp", "Gorilla", "Orangutan", "Human", "Gorilla", "Orangutan", "Chimp", "Chimp", "Gorilla", "Human", "Orangutan", "Human", "Orangutan", "Human")

ordertable1 <- data.frame(oct, oot)
colnames(ordertable1) <-c("column", "order")

ordertable2 <- data.frame(oct, ootsimp)
colnames(ordertable2) <-c("column", "order")

# 
# column    order
# X1    Human
# X2    Gorilla



orderrow <- function(col,dat,table) {
    #print("1")
  dat2 <- dat[dat[col]==1,]
     #print("2")
  species <- table[table$column == col,"order"]
    #print(species)
  dat2$cat <- paste0(col,"_",species)
  colnames(dat2)[colnames(dat2) == eval(as.symbol("species"))] <- 'example'
  dat3 <- dat2[order(dat2$example),]
  colnames(dat3) <- gsub("example", species, colnames(dat3))
  # print("4")
  #dat3
  print(dat3)
  #dim(dat3)
}
auxd <- aux[,c(2:20)]
auxd <- auxd[c("Human", "Orangutan", "Gorilla", "Chimp", "X1","X2","X3", "X4", "X5", "X6", "X7","X8","X9","X10","X11","X12","X13","X14","X15")]
auxd<-as.data.frame(auxd)
auxv <- colnames(auxd[,c(5:19)])



sep3.1 <- lapply(auxv,orderrow,dat=auxd,table=ordertable2)


sep4 <- do.call(rbind,sep3.1)


sep4$col_n <- "gray" #other/else
sep4$col_n[sep4$X1 == 1] <- "#00BFC4" #human
sep4$col_n[sep4$X2 == 1] <- "#F8766D" #chimp
sep4$col_n[sep4$X3 == 1] <- "#7CAE00" #gorilla 
sep4$col_n[sep4$X4 == 1] <- "#C77Cff" # orangutan
sep4$col_n[sep4$X15 == 1] <- "black" #all

####Finally, make the reordered heatmap.

n_matrix4<-sep4[rownames(sep4) %in% rownames(n_matrix2),c(1:4,21)]
n_matrix5 <-as.matrix(n_matrix4[,c(1:4)])


library(scales)
# n_matrix4<-n_matrix4[order(n_matrix4[,sep4$order],decreasing=T),]
dir.create(file.path(rootpath, "Plots/Model_2"))
setwd("/Volumes/GoogleDrive/My Drive/Primates_Coculture/Primates/MyAnalysis_v02/Results_2018-11-06/")
write.table(n_matrix5, file = "n_matrix5_LFC4Modelsf.txt", sep = "\t")
write.table(n_matrix4, file = "n_matrix4_LFC4Modelsf.txt", sep = "\t")
saveRDS(n_matrix5, "n_matrix5.RDS")
saveRDS(n_matrix4, "n_matrix4.RDS")



setwd(file.path(rootpath, "Plots/Model_2"))
#myBreaks <- seq(-1.5, 1.5, length.out=300)
pdf ("DEheatmap_Reorder_1.v09_.pdf")
heatmap.2(n_matrix5,
          RowSideColors=as.character(n_matrix4[as.character(unlist(rownames(n_matrix4))),]$col_n), 
          scale="none",
          symm=F,symkey=F,symbreaks=T, 
          trace="none",density.info="none", 
          col=my_palette, cexRow=0.5, cexCol = 1, 
          margins = c(5,1), labCol = c("Human", "Gorilla", "Chimp", "Orangutan"), 
          dendrogram = "none", Rowv = FALSE)
legend("bottomleft",
     legend = c("Chimp", "Gorilla", "Human", "Orangutan", "All", "Other"),
     col = c("#F8766D", "#7CAE00", "#00BFC4","#C77Cff","black","gray"),
       lty= 1,
       lwd = 5,
       cex=.7
)
dev.off()


#default colors in ggplot: 
##F8766D = red = chimp
##7CAE00=green = gorilla
##00BFC4=teal = human
##C77Cff=purple = orangutan

 

####Upset plot with the reordered heatmap.
#Make a corresponding upset plot with the subsetted data. 

aux7 <- aux6

aux7 <- aux7[!is.na(aux7$padj_chimp),]
aux7 <- aux7[!is.na(aux7$padj_human),]
aux7 <- aux7[!is.na(aux7$padj_gorilla),]
aux7 <- aux7[!is.na(aux7$padj_orangutan),]

#setwd(file.path(File_folder, "Plots/Model_2Combine"))
#pdf("UpSet_Model2Combine_LFC_0.3_v07.pdf")
#bmp("UpSet_Model1_Model2_20171110_v01.bmp")
upset(aux6,sets = colnames(aux6[,c(16:19)]), sets.bar.color = "#56B4E9",
      order.by = "freq", empty.intersections = "on", point.size=3,text.scale=c(1.5,1.5,1.5,1.5,1.5,1))







##############LRT Model proposed by Roger ##################

meta$Species_Hu[meta$Species == "Human"] <- 1
meta$Species_Hu[meta$Microbiome == "AA"] <- 0
meta$Species_Hu[is.na(meta$Species_Hu)] <- 0
#meta$Species_H[meta$Species_H == "AA"] <- " "

meta$Species_Ch[meta$Species == "Chimp"] <- 1
meta$Species_Ch[meta$Microbiome == "AA"] <- 0
meta$Species_Ch[is.na(meta$Species_Ch)] <- 0
#meta$Species_C[meta$Species_C == "AA"] <- " "

meta$Species_Go[meta$Species == "Gorilla"] <- 1
meta$Species_Go[meta$Microbiome == "AA"] <- 0
meta$Species_Go[is.na(meta$Species_Go)] <- 0
#meta$Species_G[meta$Species_G == "AA"] <- " "

meta$Species_Or[meta$Species == "Orangutan"] <- 1
meta$Species_Or[meta$Microbiome == "AA"] <- 0
meta$Species_Or[is.na(meta$Species_Or)] <- 0

meta$Species_HuCh[meta$Species == "Human" | meta$Species == "Chimp"]<- 1
meta$Species_HuCh[meta$Microbiome == "AA"] <- 0
meta$Species_HuCh[is.na(meta$Species_HuCh)] <- 0

meta$Species_HuGo[meta$Species == "Human" | meta$Species == "Gorilla"]<- 1
meta$Species_HuGo[meta$Microbiome == "AA"] <- 0
meta$Species_HuGo[is.na(meta$Species_HuGo)] <- 0

meta$Species_HuOr[meta$Species == "Human" | meta$Species == "Orangutan"]<- 1
meta$Species_HuOr[meta$Microbiome == "AA"] <- 0
meta$Species_HuOr[is.na(meta$Species_HuOr)] <- 0

meta$Species_GoCh[meta$Species == "Gorilla" | meta$Species == "Chimp"]<- 1
meta$Species_GoCh[meta$Microbiome == "AA"] <- 0
meta$Species_GoCh[is.na(meta$Species_GoCh)] <- 0

meta$Species_GoOr[meta$Species == "Gorilla" | meta$Species == "Orangutan"]<- 1
meta$Species_GoOr[meta$Microbiome == "AA"] <- 0
meta$Species_GoOr[is.na(meta$Species_GoOr)] <- 0

meta$Species_OrCh[meta$Species == "Chimp" | meta$Species == "Orangutan"]<- 1
meta$Species_OrCh[meta$Microbiome == "AA"] <- 0
meta$Species_OrCh[is.na(meta$Species_OrCh)] <- 0

meta$Species_HuChGo[meta$Species == "Human" | meta$Species == "Chimp" | meta$Species =="Gorilla"]<- 1
meta$Species_HuChGo[meta$Microbiome == "AA"] <- 0
meta$Species_HuChGo[is.na(meta$Species_HuChGo)] <- 0

meta$Species_HuChOr[meta$Species == "Human" | meta$Species == "Chimp" | meta$Species =="Orangutan"]<- 1
meta$Species_HuChOr[meta$Microbiome == "AA"] <- 0
meta$Species_HuChOr[is.na(meta$Species_HuChOr)] <- 0

meta$Species_HuGoOr[meta$Species == "Human" | meta$Species == "Orangutan" | meta$Species =="Gorilla"]<- 1
meta$Species_HuGoOr[meta$Microbiome == "AA"] <- 0
meta$Species_HuGoOr[is.na(meta$Species_HuGoOr)] <- 0

meta$Species_ChGoOr[meta$Species == "Chimp" | meta$Species == "Orangutan" | meta$Species =="Gorilla"]<- 1
meta$Species_ChGoOr[meta$Microbiome == "AA"] <- 0
meta$Species_ChGoOr[is.na(meta$Species_ChGoOr)] <- 0

meta$Species_All[meta$Species == "AA"] <- "AA"
meta$Species_All[meta$Species!= "AA"] <- "Primate"
meta$Species_All[is.na(meta$Species_All)] <- "AA"

meta$Species_All1[meta$Species == "AA"] <- 0
meta$Species_All1[meta$Species!= "AA"] <- 1
meta$Species_All1[is.na(meta$Species_All1)] <- 0

meta$Species_All1 <-factor(meta$Species_All1)
meta$Species_All <-factor(meta$Species_All)
meta$Species_Hu <-factor(meta$Species_Hu)
meta$Species_Ch <-factor(meta$Species_Ch)
meta$Species_Go <-factor(meta$Species_Go)
meta$Species_Or <-factor(meta$Species_Or)
meta$Species_HuCh <-factor(meta$Species_HuCh)
meta$Species_HuGo <-factor(meta$Species_HuGo)
meta$Species_HuOr <-factor(meta$Species_HuOr)
meta$Species_HuChGo <-factor(meta$Species_HuChGo)
meta$Species_HuChOr <-factor(meta$Species_HuChOr)
meta$Species_HuGoOr <-factor(meta$Species_HuGoOr)
meta$Species_GoCh <-factor(meta$Species_GoCh)
meta$Species_GoOr <-factor(meta$Species_GoOr)
meta$Species_OrCh <-factor(meta$Species_OrCh)
meta$Species_ChGoOr <-factor(meta$Species_ChGoOr)

meta$ExperimentPlate <-factor(meta$ExperimentPlate)


##All
dds_cgho <- DESeqDataSetFromMatrix(countData = data,  colData = meta,  design = ~ ExperimentPlate + Species_All)
dds_cgho <- DESeq(dds_cgho, test = "LRT", reduced =~ ExperimentPlate)

##Single species
dds_Hu <- DESeqDataSetFromMatrix(countData =data,  colData = meta,  design = ~ ExperimentPlate + Species_Hu)
dds_Hu <- DESeq(dds_Hu, test = "LRT", reduced =~ ExperimentPlate)

dds_Ch <- DESeqDataSetFromMatrix(countData =data, colData = meta,   design = ~ ExperimentPlate + Species_Ch)
dds_Ch <- DESeq(dds_Ch, test = "LRT", reduced =~ ExperimentPlate)

dds_Go <- DESeqDataSetFromMatrix(countData =data,  colData = meta,  design = ~ ExperimentPlate + Species_Go)
dds_Go <- DESeq(dds_Go, test = "LRT", reduced =~ ExperimentPlate)

dds_Or <- DESeqDataSetFromMatrix(countData =data,  colData = meta,  design = ~ ExperimentPlate + Species_Or)
dds_Or <- DESeq(dds_Or, test = "LRT", reduced =~ ExperimentPlate)

#Two species
dds_HuCh <- DESeqDataSetFromMatrix(countData =data, colData = meta,   design = ~ ExperimentPlate + Species_HuCh)
dds_HuCh <- DESeq(dds_HuCh, test = "LRT", reduced =~ ExperimentPlate)

dds_HuGo <- DESeqDataSetFromMatrix(countData =data, colData = meta,   design = ~ ExperimentPlate + Species_HuGo)
dds_HuGo <- DESeq(dds_HuGo, test = "LRT", reduced =~ ExperimentPlate)

dds_HuOr <- DESeqDataSetFromMatrix(countData =data,  colData = meta,  design = ~ ExperimentPlate + Species_HuOr)
dds_HuOr <- DESeq(dds_HuOr, test = "LRT", reduced =~ ExperimentPlate)

dds_GoOr <- DESeqDataSetFromMatrix(countData =data,  colData = meta,  design = ~ ExperimentPlate + Species_GoOr)
dds_GoOr <- DESeq(dds_GoOr, test = "LRT", reduced =~ ExperimentPlate)

dds_GoCh <- DESeqDataSetFromMatrix(countData =data,  colData = meta,  design = ~ ExperimentPlate + Species_GoCh)
dds_GoCh <- DESeq(dds_GoCh, test = "LRT", reduced =~ ExperimentPlate)

dds_OrCh <- DESeqDataSetFromMatrix(countData =data,  colData = meta,  design = ~ ExperimentPlate + Species_OrCh)
dds_OrCh <- DESeq(dds_OrCh, test = "LRT", reduced =~ ExperimentPlate)

##Three species
dds_HuChGo <- DESeqDataSetFromMatrix(countData =data,  colData = meta,  design = ~ ExperimentPlate + Species_HuChGo)
dds_HuChGo <- DESeq(dds_HuChGo, test = "LRT", reduced =~ ExperimentPlate)

dds_HuChOr <- DESeqDataSetFromMatrix(countData =data,  colData = meta,  design = ~ ExperimentPlate + Species_HuChOr)
dds_HuChOr <- DESeq(dds_HuChOr, test = "LRT", reduced =~ ExperimentPlate)

dds_HuGoOr<- DESeqDataSetFromMatrix(countData =data, colData = meta,   design = ~ ExperimentPlate + Species_HuGoOr)
dds_HuGoOr <- DESeq(dds_HuGoOr, test = "LRT", reduced =~ ExperimentPlate)

dds_ChGoOr<- DESeqDataSetFromMatrix(countData =data,  colData = meta,  design = ~ ExperimentPlate + Species_ChGoOr)
dds_ChGoOr <- DESeq(dds_ChGoOr, test = "LRT", reduced =~ ExperimentPlate)


#single 
dds_Hu <- dds_Hu[ rowSums(counts(dds_Hu)) > 1, ]
rdds_Hu <- results(dds_Hu)
sum(rdds_Hu$padj < 0.1, na.rm=TRUE)
sum(rdds_Hu$pvalue < 0.1, na.rm=TRUE)
sum(!is.na(rdds_Hu$pvalue))


dds_Ch <- dds_Ch[ rowSums(counts(dds_Ch)) > 1, ]
rdds_Ch <- results(dds_Ch)
sum(rdds_Ch$padj < 0.1, na.rm=TRUE)
sum(rdds_Ch$pvalue < 0.1, na.rm=TRUE)
sum(!is.na(rdds_Ch$pvalue))

dds_Go <- dds_Go[ rowSums(counts(dds_Go)) > 1, ]
rdds_Go <- results(dds_Go)
sum(rdds_Go$padj < 0.1, na.rm=TRUE)
sum(rdds_Go$pvalue < 0.1, na.rm=TRUE)
sum(!is.na(rdds_Go$pvalue))

dds_Or <- dds_Or[ rowSums(counts(dds_Or)) > 1, ]
rdds_Or <- results(dds_Or)
sum(rdds_Or$padj < 0.1, na.rm=TRUE)
sum(rdds_Or$pvalue < 0.1, na.rm=TRUE)
sum(!is.na(rdds_Or$pvalue))

#pairs
dds_HuCh <- dds_HuCh[ rowSums(counts(dds_HuCh)) > 1, ]
rdds_HuCh <- results(dds_HuCh)

dds_HuGo <- dds_HuGo[ rowSums(counts(dds_HuGo)) > 1, ]
rdds_HuGo <- results(dds_HuGo)

dds_HuOr <- dds_HuOr[ rowSums(counts(dds_HuOr)) > 1, ]
rdds_HuOr <- results(dds_HuOr)

dds_GoCh <- dds_GoCh[ rowSums(counts(dds_GoCh)) > 1, ]
rdds_GoCh <- results(dds_GoCh)

dds_GoOr <- dds_GoOr[ rowSums(counts(dds_GoOr)) > 1, ]
rdds_GoOr <- results(dds_GoOr)

dds_OrCh <- dds_OrCh[ rowSums(counts(dds_OrCh)) > 1, ]
rdds_OrCh <- results(dds_OrCh)

#triples
dds_HuChGo <- dds_HuChGo[ rowSums(counts(dds_HuChGo)) > 1, ]
rdds_HuChGo <-results(dds_HuChGo)

dds_HuChOr <- dds_HuChOr[ rowSums(counts(dds_HuChOr)) > 1, ]
rdds_HuChOr <- results(dds_HuChOr)

dds_HuGoOr <- dds_HuGoOr[ rowSums(counts(dds_HuGoOr)) > 1, ]
rdds_HuGoOr <- results(dds_HuGoOr)

dds_ChGoOr <- dds_ChGoOr[ rowSums(counts(dds_ChGoOr)) > 1, ]
rdds_ChGoOr <- results(dds_ChGoOr)

##All
dds_cgho <- dds_cgho[rowSums(counts(dds_cgho)) >1,]
rdds_cgho <- results(dds_cgho)


#Rename columns
#All
names(rdds_cgho) <-  c("baseMean_chimp_gorilla_human_orangutan", "log2FoldChange_chimp_gorilla_human_orangutan", "lfcSE_chimp_gorilla_human_orangutan", "stat_chimp_gorilla_human_orangutan", "pvalue_chimp_gorilla_human_orangutan", "padj_chimp_gorilla_human_orangutan")

#singles
names(rdds_Hu) <-  c("baseMean_human", "log2FoldChange_human", "lfcSE_human", "stat_human", "pvalue_human", "padj_human")

names(rdds_Ch) <-  c("baseMean_chimp", "log2FoldChange_chimp", "lfcSE_chimp", "stat_chimp", "pvalue_chimp", "padj_chimp")

names(rdds_Go) <-  c("baseMean_gorilla", "log2FoldChange_gorilla", "lfcSE_gorilla", "stat_gorilla", "pvalue_gorilla", "padj_gorilla")

names(rdds_Or) <-  c("baseMean_orangutan", "log2FoldChange_orangutan", "lfcSE_orangutan", "stat_orangutan", "pvalue_orangutan", "padj_orangutan")

#doubles
names(rdds_HuCh) <-  c("baseMean_human_chimp", "log2FoldChange_human_chimp", "lfcSE_human_chimp", "stat_human_chimp", "pvalue_human_chimp", "padj_human_chimp")

names(rdds_HuGo) <-  c("baseMean_human_gorilla", "log2FoldChange_human_gorilla", "lfcSE_human_gorilla", "stat_human_gorilla", "pvalue_human_gorilla", "padj_human_gorilla")

names(rdds_HuOr) <-  c("baseMean_human_orangutan", "log2FoldChange_human_orangutan", "lfcSE_human_orangutan", "stat_human_orangutan", "pvalue_human_orangutan", "padj_human_orangutan")

names(rdds_GoOr) <-  c("baseMean_gorilla_orangutan", "log2FoldChange_gorilla_orangutan", "lfcSE_gorilla_orangutan", "stat_gorilla_orangutan", "pvalue_gorilla_orangutan", "padj_gorilla_orangutan")

names(rdds_GoCh) <-  c("baseMean_gorilla_chimp", "log2FoldChange_gorilla_chimp", "lfcSE_gorilla_chimp", "stat_gorilla_chimp", "pvalue_gorilla_chimp", "padj_gorilla_chimp")

names(rdds_OrCh) <-  c("baseMean_orangutan_chimp", "log2FoldChange_orangutan_chimp", "lfcSE_orangutan_chimp", "stat_orangutan_chimp", "pvalue_orangutan_chimp", "padj_orangutan_chimp")


#triples
names(rdds_HuChGo) <-  c("baseMean_human_chimp_gorilla", "log2FoldChange_human_chimp_gorilla", "lfcSE_human_chimp_gorilla", "stat_human_chimp_gorilla", "pvalue_human_chimp_gorilla", "padj_human_chimp_gorilla")

names(rdds_HuChOr) <-  c("baseMean_human_chimp_orangutan", "log2FoldChange_human_chimp_orangutan", "lfcSE_human_chimp_orangutan", "stat_human_chimp_orangutan", "pvalue_human_chimp_orangutan", "padj_human_chimp_orangutan")

names(rdds_HuGoOr) <-  c("baseMean_human_gorilla_orangutan", "log2FoldChange_human_gorilla_orangutan", "lfcSE_human_gorilla_orangutan", "stat_human_gorilla_orangutan", "pvalue_human_gorilla_orangutan", "padj_human_gorilla_orangutan")

names(rdds_ChGoOr) <-  c("baseMean_chimp_gorilla_orangutan", "log2FoldChange_chimp_gorilla_orangutan", "lfcSE_chimp_gorilla_orangutan", "stat_chimp_gorilla_orangutan", "pvalue_chimp_gorilla_orangutan", "padj_chimp_gorilla_orangutan")



lrtall<- cbind( rdds_Ch, rdds_Go, rdds_Hu, rdds_Or, rdds_HuCh, rdds_HuGo, rdds_HuOr, rdds_GoOr, rdds_GoCh, rdds_OrCh, rdds_HuChGo, rdds_HuChOr, rdds_HuGoOr, rdds_ChGoOr, rdds_cgho)

lrtall$symbol <- mapIds(org.Hs.eg.db,
                        keys=row.names(lrtall),
                        column="SYMBOL",
                        keytype="ENSEMBL",
                        multiVals="first")
lrtall$entrez <- mapIds(org.Hs.eg.db,
                        keys=row.names(lrtall),
                        column="ENTREZID",
                        keytype="ENSEMBL",
                        multiVals="first")

dir.create(file.path(rootpath, "DE_Model7"))
setwd(file.path(rootpath, "DE_Model7"))
write.table(lrtall, "Results_LRT_15_v02.txt", sep = "\t")
write.table(lrtall, "SupplementalTable4_LRT.txt", sep = "\t")

#setwd("~/Google Drive/Primates_Coculture/Primates/MyAnalysis/Results_v08/LRT/")
#singles

###This next part came from Roger

library(tidyverse)

library("readr")
library("tidyr")
library("dplyr")
library(magrittr)
library(data.table)
#This next line is for reading in the data if the data hasn't been already loaded
#The tidyverse package is weird and requires me to restart the session for it to work, 
#so I'm just going to use the model data generated above

d <- read_tsv(file.path( "/Volumes/GoogleDrive/My\ Drive/Primates_Coculture/Primates/MyAnalysis_v02/Results_2018-11-06/DE_Model7/Results_LRT_15_v02.txt"))

#d <-lrtall
## to convert stat to pvalue 1-pchisq(d$stat_chimp*2,df=1))
## this assumes from Wiki LRT page:
# A convenient result by Samuel S. Wilks, says that as the sample size n
# approaches  \infty , the test statistic -2 \log(\Lambda) for a nested model 
# asymptotically will be chi-squared distributed \chi ^{2} with degrees of
# freedom equal to the difference in dimensionality of  \Theta  and  \Theta_0
# (here 1) when  H_{0} holds true.[6] This means that for a great variety of
# hypotheses, a practitioner can compute the likelihood ratio \Lambda  for the
# data and compare -2\log(\Lambda) to the \chi ^{2} value corresponding to a
# desired statistical significance as an approximate statistical test.


## Pick all the colums starting with stat
ds <- d %>% dplyr::select(starts_with("stat_")) %>% as.matrix
ds <- ds * 2
rownames(ds) <- d$X1
saveRDS(ds, "ds.RDS")
saveRDS(d, "d.RDS")
## All the values should be >0, but they are not, this could be just numerical precisions for genes that are probably not differentially express across any condition. 
## I think we could remove all the rows that do not generate a pvalue for at least one of the tests. 

## Identify those genes that are at least significant for one of the models
padj <- d %>% dplyr::select(starts_with("padj_")) %>% as.matrix
sig <- rowSums(padj<0.1,na.rm=T)>0

## subset the singificant ones

aux <- exp(ds[sig,])

## now we have the matrix of likelihood ratios P(D/Hi)/P(D/H0) for i being each of the 15 configurations 

## Bayes theorem P(H_i|D) \propto P(D|H_i) * P(H_i) 
## we consisder all the configurations equally likely a priory. 
## so to calculate the posteriors we need to normalize the data. 
## the posterior for all the configurations needs to add to 1

rs <- rowSums(aux)
post <- aux/rs
colnames(post) <- gsub("stat","post",colnames(post))

## this makes a plot with the expected number of signals in each configuration. 
psum <- colSums(post)
mydf <- tibble(conf=names(psum),psum=psum)
mydf %>%
  ggplot(aes(x=reorder(conf,psum),y=psum)) +
  geom_col() + coord_flip() + geom_text(aes(label=round(psum,1)))



pt <- as.tibble(post) %>% mutate(ensg=rownames(post))

best_conf <- gather(pt,key="conf",value="post",-ensg) %>% group_by(ensg) %>% top_n(1,post) 

write.table(post,file=gzfile("post.txt.gz"),sep="\t",quote=F)

write_tsv(best_conf,"best_conf.txt")

mem1<-table(best_conf$conf)



#singles
g_c <- as.data.frame(best_conf %>% filter(conf=="post_chimp") %>% dplyr::select(ensg) %>% unlist)
g_g <-as.data.frame(best_conf %>% filter(conf=="post_gorilla") %>% dplyr::select(ensg) %>% unlist)
g_o <- as.data.frame(best_conf %>% filter(conf=="post_orangutan") %>% dplyr::select(ensg) %>% unlist)
g_h <- as.data.frame(best_conf %>% filter(conf=="post_human") %>% dplyr::select(ensg) %>% unlist)

#doubles
g_cg <-  as.data.frame(best_conf %>% filter(conf=="post_gorilla_chimp") %>% dplyr::select(ensg) %>% unlist)
g_ch <-  as.data.frame(best_conf %>% filter(conf=="post_human_chimp") %>% dplyr::select(ensg) %>% unlist)
g_co <-  as.data.frame(best_conf %>% filter(conf=="post_orangutan_chimp") %>% dplyr::select(ensg) %>% unlist)
g_gh <- as.data.frame(best_conf %>% filter(conf=="post_human_gorilla") %>% dplyr::select(ensg) %>% unlist)
g_go <- as.data.frame(best_conf %>% filter(conf=="post_gorilla_orangutan") %>% dplyr::select(ensg) %>% unlist)
g_ho <- as.data.frame(best_conf %>% filter(conf=="post_human_orangutan") %>% dplyr::select(ensg) %>% unlist)


#triples
g_cgh<- as.data.frame(best_conf %>% filter(conf=="post_human_chimp_gorilla") %>% dplyr::select(ensg) %>% unlist)
g_cgo<- as.data.frame(best_conf %>% filter(conf=="post_chimp_gorilla_orangutan") %>% dplyr::select(ensg) %>% unlist)
g_cho<- as.data.frame(best_conf %>% filter(conf=="post_human_chimp_orangutan") %>% dplyr::select(ensg) %>% unlist)
g_gho<- as.data.frame(best_conf %>% filter(conf=="post_human_gorilla_orangutan") %>% dplyr::select(ensg) %>% unlist)

#all
g_cgho <- as.data.frame(best_conf %>% filter(conf=="post_chimp_gorilla_human_orangutan") %>% dplyr::select(ensg) %>% unlist)

colnames(g_c) <- "ENSG"
colnames(g_g) <- "ENSG"
colnames(g_h) <- "ENSG"
colnames(g_o) <- "ENSG"

colnames(g_cg) <- "ENSG"
colnames(g_ch) <- "ENSG"
colnames(g_co) <- "ENSG"
colnames(g_gh) <- "ENSG"
colnames(g_go) <- "ENSG"
colnames(g_ho) <- "ENSG"

colnames(g_cgh) <- "ENSG"
colnames(g_cgo) <- "ENSG"
colnames(g_cho) <- "ENSG"
colnames(g_gho) <- "ENSG"

colnames(g_cgho) <- "ENSG"


setwd(file.path(rootpath, "DE_Model7"))
write.table(g_c, file = "LRT_chimp.txt", sep= "\t")
write.table(g_g, file = "LRT_gorilla.txt", sep= "\t")
write.table(g_o, file = "LRT_orangutan.txt", sep= "\t")
write.table(g_h, file = "LRT_human.txt", sep= "\t")

write.table(g_cg, file = "LRT_chimp_gorilla.txt", sep = "\t")
write.table(g_ch, file = "LRT_chimp_human.txt", sep = "\t")
write.table(g_co, file = "LRT_chimp_orangutan.txt", sep = "\t")
write.table(g_gh, file = "LRT_gorilla_human.txt", sep = "\t")
write.table(g_go, file = "LRT_gorilla_orangutan.txt", sep = "\t")
write.table(g_ho, file = "LRT_human_orangutan.txt", sep = "\t")


write.table(g_cgh, file = "LRT_Triple_CGH.txt", sep= "\t")
write.table(g_cgo, file = "LRT_Triple_CGO.txt", sep= "\t")
write.table(g_cho, file = "LRT_Triple_CHO.txt", sep= "\t")
write.table(g_gho, file = "LRT_Triple_GHO.txt", sep= "\t")

write.table(g_cgho, file = "LRT_chimp_gorilla_human_orangutan.txt", sep = "\t")

dim(g_c)
dim(g_g)
dim(g_h)
dim(g_o)

dim(g_cg)
dim(g_ch)
dim(g_co)
dim(g_gh)
dim(g_go)
dim(g_ho)

dim(g_cgh)
dim(g_cgo)
dim(g_cho)
dim(g_gho)

dim(g_cgho)


library("UpSetR")
mem1 <-as.data.frame(mem1)
mem2 <- c(chimp=24, gorilla = 57, human = 71, orangutan=12, `chimp&gorilla&human&orangutan`=2261, `chimp&gorilla&orangutan`=495, `gorilla&chimp`=113, `gorilla&orangutan`=184, `human&chimp`=57,`human&chimp&gorilla`=309, `human&chimp&orangutan`=144, `human&gorilla`=161, `human&gorilla&orangutan`=365, `human&orangutan`=36, `orangutan&chimp`=40)

pdf("Upset_LRT.pdf")
upset(fromExpression(mem2), order.by = "freq",point.size=3,text.scale=c(1.5,1.5,1.5,1.5,1.5,1))

dev.off()


g_c <- read.table("LRT_chimp.txt", sep="\t")
g_g <- read.table("LRT_gorilla.txt", sep="\t")
g_h <- read.table("LRT_human.txt", sep="\t")
g_o <- read.table("LRT_orangutan.txt", sep="\t")

#doubles
g_cg <- read.table("LRT_chimp_gorilla.txt", sep ="\t")
g_ch <- read.table("LRT_chimp_human.txt", sep ="\t")
g_co <- read.table("LRT_chimp_orangutan.txt", sep = "\t")

g_gh <- read.table("LRT_gorilla_human.txt", sep ="\t")
g_go <- read.table("LRT_gorilla_orangutan.txt", sep = "\t")
g_ho <- read.table("LRT_human_orangutan.txt", sep = "\t")

#triples
g_cgh <- read.table("LRT_Triple_CGH.txt", sep = "\t")
g_cgo <-read.table("LRT_Triple_CGO.txt", sep = "\t")
g_cho <-read.table("LRT_Triple_CHO.txt", sep = "\t")
g_gho <-read.table("LRT_Triple_GHO.txt", sep = "\t")

#all

g_cgho <- read.table("LRT_chimp_gorilla_human_orangutan.txt", sep = "\t")




#singles
rownames(g_c) <- g_c$ENSG
rownames(g_g) <- g_g$ENSG
rownames(g_h) <- g_h$ENSG
rownames(g_o) <- g_o$ENSG

#doubles
rownames(g_cg) <- g_cg$ENSG
rownames(g_ch) <- g_ch$ENSG
rownames(g_co) <- g_co$ENSG
rownames(g_gh) <- g_gh$ENSG
rownames(g_go) <- g_go$ENSG
rownames(g_ho) <- g_ho$ENSG

#triples
rownames(g_cgh) <- g_cgh$ENSG
rownames(g_cgo) <- g_cgo$ENSG
rownames(g_cho) <- g_cho$ENSG
rownames(g_gho) <- g_gho$ENSG

#all
rownames(g_cgho) <-g_cgho$ENSG

##adding on gene names
g_c$symbol <- mapIds(org.Hs.eg.db,
                     keys=row.names(g_c),
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")
g_c$entrez <- mapIds(org.Hs.eg.db,
                     keys=row.names(g_c),
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")

g_g$symbol <- mapIds(org.Hs.eg.db,
                     keys=row.names(g_g),
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")
g_g$entrez <- mapIds(org.Hs.eg.db,
                     keys=row.names(g_g),
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")

g_h$symbol <- mapIds(org.Hs.eg.db,
                     keys=row.names(g_h),
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")
g_h$entrez <- mapIds(org.Hs.eg.db,
                     keys=row.names(g_h),
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")

g_o$symbol <- mapIds(org.Hs.eg.db,
                     keys=row.names(g_o),
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")
g_o$entrez <- mapIds(org.Hs.eg.db,
                     keys=row.names(g_o),
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")

#doubles
g_cg$symbol <- mapIds(org.Hs.eg.db,
                      keys=row.names(g_cg),
                      column="SYMBOL",
                      keytype="ENSEMBL",
                      multiVals="first")
g_cg$entrez <- mapIds(org.Hs.eg.db,
                      keys=row.names(g_cg),
                      column="ENTREZID",
                      keytype="ENSEMBL",
                      multiVals="first")

g_ch$symbol <- mapIds(org.Hs.eg.db,
                      keys=row.names(g_ch),
                      column="SYMBOL",
                      keytype="ENSEMBL",
                      multiVals="first")
g_ch$entrez <- mapIds(org.Hs.eg.db,
                      keys=row.names(g_ch),
                      column="ENTREZID",
                      keytype="ENSEMBL",
                      multiVals="first")

g_co$symbol <- mapIds(org.Hs.eg.db,
                      keys=row.names(g_co),
                      column="SYMBOL",
                      keytype="ENSEMBL",
                      multiVals="first")
g_co$entrez <- mapIds(org.Hs.eg.db,
                      keys=row.names(g_co),
                      column="ENTREZID",
                      keytype="ENSEMBL",
                      multiVals="first")

g_gh$symbol <- mapIds(org.Hs.eg.db,
                      keys=row.names(g_gh),
                      column="SYMBOL",
                      keytype="ENSEMBL",
                      multiVals="first")
g_gh$entrez <- mapIds(org.Hs.eg.db,
                      keys=row.names(g_gh),
                      column="ENTREZID",
                      keytype="ENSEMBL",
                      multiVals="first")

g_go$symbol <- mapIds(org.Hs.eg.db,
                      keys=row.names(g_go),
                      column="SYMBOL",
                      keytype="ENSEMBL",
                      multiVals="first")
g_go$entrez <- mapIds(org.Hs.eg.db,
                      keys=row.names(g_go),
                      column="ENTREZID",
                      keytype="ENSEMBL",
                      multiVals="first")

g_ho$symbol <- mapIds(org.Hs.eg.db,
                      keys=row.names(g_ho),
                      column="SYMBOL",
                      keytype="ENSEMBL",
                      multiVals="first")
g_ho$entrez <- mapIds(org.Hs.eg.db,
                      keys=row.names(g_ho),
                      column="ENTREZID",
                      keytype="ENSEMBL",
                      multiVals="first")




#triples
g_cgh$symbol <- mapIds(org.Hs.eg.db,
                       keys=row.names(g_cgh),
                       column="SYMBOL",
                       keytype="ENSEMBL",
                       multiVals="first")
g_cgh$entrez <- mapIds(org.Hs.eg.db,
                       keys=row.names(g_cgh),
                       column="ENTREZID",
                       keytype="ENSEMBL",
                       multiVals="first")

g_cgo$symbol <- mapIds(org.Hs.eg.db,
                       keys=row.names(g_cgo),
                       column="SYMBOL",
                       keytype="ENSEMBL",
                       multiVals="first")
g_cgo$entrez <- mapIds(org.Hs.eg.db,
                       keys=row.names(g_cgo),
                       column="ENTREZID",
                       keytype="ENSEMBL",
                       multiVals="first")

g_cho$symbol <- mapIds(org.Hs.eg.db,
                       keys=row.names(g_cho),
                       column="SYMBOL",
                       keytype="ENSEMBL",
                       multiVals="first")
g_cho$entrez <- mapIds(org.Hs.eg.db,
                       keys=row.names(g_cho),
                       column="ENTREZID",
                       keytype="ENSEMBL",
                       multiVals="first")

g_gho$symbol <- mapIds(org.Hs.eg.db,
                       keys=row.names(g_gho),
                       column="SYMBOL",
                       keytype="ENSEMBL",
                       multiVals="first")
g_gho$entrez <- mapIds(org.Hs.eg.db,
                       keys=row.names(g_gho),
                       column="ENTREZID",
                       keytype="ENSEMBL",
                       multiVals="first")

#All

g_cgho$symbol <- mapIds(org.Hs.eg.db,
                        keys=row.names(g_cgho),
                        column="SYMBOL",
                        keytype="ENSEMBL",
                        multiVals="first")
g_cgho$entrez <- mapIds(org.Hs.eg.db,
                        keys=row.names(g_cgho),
                        column="ENTREZID",
                        keytype="ENSEMBL",
                        multiVals="first")

###Get lists of genes for enrichment analysis/IPA:

res_cc1 <- res_cc
res_cc1$ensg <- rownames(res_cc1)
res_gc1 <- res_gc
res_gc1$ensg <- rownames(res_gc1)
res_hc1 <- res_hc
res_hc1$ensg <- rownames(res_hc1)
res_oc1 <- res_oc
res_oc1$ensg <- rownames(res_oc1)

listc <- g_c$ENSG
listg <- g_g$ENSG
listh <- g_h$ENSG
listo <- g_o$ENSG

listcg <- g_cg$ENSG
listch <- g_ch$ENSG
listco <- g_co$ENSG
listgh <- g_gh$ENSG
listgo <- g_go$ENSG
listho <- g_ho$ENSG

listcgh <- g_cgh$ENSG
listcgo <- g_cgo$ENSG
listcho <- g_cho$ENSG
listgho <- g_gho$ENSG

listca <- c(as.character(listc), as.character(listgho))
listga <- c(as.character(listg), as.character(listcho))
listha <- c(as.character(listh), as.character(listcgo))
listoa <- c(as.character(listo), as.character(listcgh))

chimpinc <- as.data.frame(res_cc1[res_cc1$ensg %in% listca,])
gorillainc <- as.data.frame(res_gc1[res_gc1$ensg %in% listga,])
humaninc <- as.data.frame(res_hc1[res_hc1$ensg %in% listha,])
orangutaninc <- as.data.frame(res_hc1[res_hc1$ensg %in% listoa,])

setwd(file.path(rootpath, "DE_Model7"))

write.table(chimpinc, file = "Chimp_allgenesIPA.txt", sep = "\t")
write.table(gorillainc, file = "Gorilla_allgenesIPA.txt", sep = "\t")
write.table(humaninc, file = "Human_allgenesIPA.txt", sep = "\t")
write.table(orangutaninc, file = "Orangutan_allgenesIPA.txt", sep = "\t")

#################END LRT Model##############

library("gridExtra")
library("grid")
library("lattice")

dir.create(file.path(rootpath, "Plots/Model_7"))
dir.create(file.path(rootpath, "Plots/Model_7/Human"))
dir.create(file.path(rootpath, "Plots/Model_7/Chimp"))
dir.create(file.path(rootpath, "Plots/Model_7/Gorilla"))
dir.create(file.path(rootpath, "Plots/Model_7/Orangutan"))

##Make plots for human-chimp combo
dir.create(file.path(rootpath, "Plots/Model_7/HumanChimp"))

#list_genes <- g_g$best_conf.....filter.conf.....post_gorilla.......select.ensg......unlist
#list_genes <- g_cgo$best_conf.....filter.conf.....post_chimp_gorilla_orangutan.......select.ensg......unlist
#list_genes <- g_cgh$best_conf.....filter.conf.....post_human_chimp_gorilla.......select.ensg......unlist
#list_genes<- g_gho$best_conf.....filter.conf.....post_human_gorilla_orangutan.......select.ensg......unlist
#list_genes <- gorillalist$ENSG
#list_genes <- humanlist$ENSG
#list_genes <- orangutanlist$ENSG
list_genes <- g_h$ENSG
list_genes <- g_g$ENSG
list_genes <- g_c$ENSG
list_genes <- g_o$ENSG
list_genes <- g_ch$ENSG


ac1 <- data.frame()
title <- character()
p<-list()
#plots <- list() 

file1 <-file.path(rootpath, "Plots/Model_7/Orangutan/")
setwd(file1)
for(a in list_genes)
{
  
  ac1 <- data.frame(c("Human", "Chimp", "Gorilla", "Orangutan"), 
                    c(res_hc[a, "log2FoldChange"], res_cc[a, "log2FoldChange"], res_gc[a, "log2FoldChange"], res_oc[a, "log2FoldChange"]),
                    c(res_hc[a, "lfcSE"], res_cc[a, "lfcSE"], res_gc[a, "lfcSE"], res_oc[a, "lfcSE"]))
  
  colnames(ac1) <- c("Species", "LFC", "SE")
  
  
  p[[a]]<- ggplot(ac1, aes(x=Species, y=LFC, group = Species, color=Species)) +
    geom_pointrange(aes(ymin=LFC-SE, ymax=LFC+SE))+labs(title=res_hc[a,"symbol"]) + 
    geom_errorbar(aes(ymin=LFC-SE, ymax=LFC+SE), width=.2) +
    geom_point(shape = 1, size = 3) + 
    theme_bw() + 
    theme(legend.position="none", axis.title.x = element_text( size=10), axis.text.x  = element_text( size=10), axis.title.y = element_text(size=10), axis.text.y = element_text(size=10), plot.title = element_text(size=10)) 
  
  # pdf(res_hc[a, "symbol"], ".pdf")
  # print(p)
  # dev.off()
}

m1 <-marrangeGrob(grobs = p, nrow = 2, ncol = 2)


ggsave(paste0(file.path(rootpath, "Plots/Model_7/Orangutan.pdf")), m1)

setwd(file.path(rootpath, "Plots/Model_7/"))
a <- "ENSG00000112309"
ac1 <- data.frame(c("Human", "Chimp", "Gorilla", "Orangutan"), 
                  c(res_hc[a, "log2FoldChange"], res_cc[a, "log2FoldChange"], res_gc[a, "log2FoldChange"], res_oc[a, "log2FoldChange"]),
                  c(res_hc[a, "lfcSE"], res_cc[a, "lfcSE"], res_gc[a, "lfcSE"], res_oc[a, "lfcSE"]))

colnames(ac1) <- c("Species", "LFC", "SE")
pdf("B3GAT2.pdf")
p <-ggplot(ac1, aes(x=Species, y=LFC, group = Species, color=Species)) +
  geom_pointrange(aes(ymin=LFC-SE, ymax=LFC+SE), size = 2)+labs(title=res_hc[a,"symbol"]) +
  geom_errorbar(aes(ymin=LFC-SE, ymax=LFC+SE), width=.2,size=2) +
  geom_point(shape = 1, size = 3) +
  theme_bw() +
  theme(legend.position="none", axis.title.x = element_text( size=10), axis.text.x  = element_text( size=10), axis.title.y = element_text(size=10), axis.text.y = element_text(size=10), plot.title = element_text(size=10)) 

print(p)
dev.off()

a<- "ENSG00000117593"
ac1 <- data.frame(c("Human", "Chimp", "Gorilla", "Orangutan"), 
                  c(res_hc[a, "log2FoldChange"], res_cc[a, "log2FoldChange"], res_gc[a, "log2FoldChange"], res_oc[a, "log2FoldChange"]),
                  c(res_hc[a, "lfcSE"], res_cc[a, "lfcSE"], res_gc[a, "lfcSE"], res_oc[a, "lfcSE"]))

colnames(ac1) <- c("Species", "LFC", "SE")
pdf("DARS2.pdf")
p <-ggplot(ac1, aes(x=Species, y=LFC, group = Species, color=Species)) +
  geom_pointrange(aes(ymin=LFC-SE, ymax=LFC+SE), size = 2)+labs(title=res_hc[a,"symbol"]) +
  geom_errorbar(aes(ymin=LFC-SE, ymax=LFC+SE), width=.2,size=2) +
  geom_point(shape = 1, size = 3) +
  theme_bw() +
  theme(legend.position="none", axis.title.x = element_text( size=10), axis.text.x  = element_text( size=10), axis.title.y = element_text(size=10), axis.text.y = element_text(size=10), plot.title = element_text(size=10)) 

print(p)
dev.off()


a<- "ENSG00000144048"
ac1 <- data.frame(c("Human", "Chimp", "Gorilla", "Orangutan"), 
                  c(res_hc[a, "log2FoldChange"], res_cc[a, "log2FoldChange"], res_gc[a, "log2FoldChange"], res_oc[a, "log2FoldChange"]),
                  c(res_hc[a, "lfcSE"], res_cc[a, "lfcSE"], res_gc[a, "lfcSE"], res_oc[a, "lfcSE"]))

colnames(ac1) <- c("Species", "LFC", "SE")
pdf("DUSP11.pdf")
p <-ggplot(ac1, aes(x=Species, y=LFC, group = Species, color=Species)) +
  geom_pointrange(aes(ymin=LFC-SE, ymax=LFC+SE), size = 2)+labs(title=res_hc[a,"symbol"]) +
  geom_errorbar(aes(ymin=LFC-SE, ymax=LFC+SE), width=.2,size=2) +
  geom_point(shape = 1, size = 3) +
  theme_bw() +
  theme(legend.position="none", axis.title.x = element_text( size=10), axis.text.x  = element_text( size=10), axis.title.y = element_text(size=10), axis.text.y = element_text(size=10), plot.title = element_text(size=10)) 

print(p)
dev.off()

a<- "ENSG00000138771" 
ac1 <- data.frame(c("Human", "Chimp", "Gorilla", "Orangutan"), 
                  c(res_hc[a, "log2FoldChange"], res_cc[a, "log2FoldChange"], res_gc[a, "log2FoldChange"], res_oc[a, "log2FoldChange"]),
                  c(res_hc[a, "lfcSE"], res_cc[a, "lfcSE"], res_gc[a, "lfcSE"], res_oc[a, "lfcSE"]))

colnames(ac1) <- c("Species", "LFC", "SE")
pdf("SHROOM3.pdf")
p <-ggplot(ac1, aes(x=Species, y=LFC, group = Species, color=Species)) +
  geom_pointrange(aes(ymin=LFC-SE, ymax=LFC+SE), size = 2)+labs(title=res_hc[a,"symbol"]) +
  geom_errorbar(aes(ymin=LFC-SE, ymax=LFC+SE), width=.2,size=2) +
  geom_point(shape = 1, size = 3) +
  theme_bw() +
  theme(legend.position="none", axis.title.x = element_text( size=10), axis.text.x  = element_text( size=10), axis.title.y = element_text(size=10), axis.text.y = element_text(size=10), plot.title = element_text(size=10)) 


print(p)
dev.off()


list_genes <- c(head(rownames(g_ch), 20), tail(rownames(g_ch), 20))
ac1 <- data.frame()
title <- character()
p<-list()
#plots <- list() 
for(a in list_genes)
{
  
  ac1 <- data.frame(c("Human", "Chimp", "Gorilla", "Orangutan"), 
                    c(res_hc[a, "log2FoldChange"], res_cc[a, "log2FoldChange"], res_gc[a, "log2FoldChange"], res_oc[a, "log2FoldChange"]),
                    c(res_hc[a, "lfcSE"], res_cc[a, "lfcSE"], res_gc[a, "lfcSE"], res_oc[a, "lfcSE"]))
  
  colnames(ac1) <- c("Species", "LFC", "SE")
  
  
  p[[a]]<- ggplot(ac1, aes(x=Species, y=LFC, group = Species, color=Species)) +
    geom_pointrange(aes(ymin=LFC-SE, ymax=LFC+SE))+labs(title= a) + 
    geom_errorbar(aes(ymin=LFC-SE, ymax=LFC+SE), width=.2) +
    geom_point(shape = 1, size = 3) + 
    theme_bw() + 
    theme(legend.position="none", axis.title.x = element_text( size=10), axis.text.x  = element_text( size=10), axis.title.y = element_text(size=10), axis.text.y = element_text(size=10), plot.title = element_text(size=10)) 
}


m1 <-marrangeGrob(grobs = p, nrow = 2, ncol = 2)
ggsave(paste0(file.path(rootpath, "/Plots/Model_7/HumanChimp.pdf")), m1)

#/Volumes/GoogleDrive/My\ Drive/Primates_Coculture/Primates/MyAnalysis_v02/Results_2019-01-30/Plots/Model_7/HumanChimp
#+





##########GWAS enrichment############
##Test for gwas enrichment for Model6 results. Model 2 is 4 models, species vs control

setwd("~/Google Drive/GWAS_Catalog_2018/")
library("readr")
library("AnnotationDbi")
gwas <- read_tsv("gwas_catalog_v1.0.1-associations_e91_r2018-01-16.tsv")

dir.create(file.path(rootpath, "GWAS"))
setwd(file.path(rootpath, "GWAS"))


###This is to determine GWAS enrichment for genes identified in the LRT model
#Including non-expressed genes for each species
library("data.table")
hall <- as.data.frame(c(g_h$symbol, g_cgo$symbol))
colnames(hall) <- "symbol"
call <- as.data.frame(c(g_c$symbol, g_gho$symbol))
colnames(call) <- "symbol"
gall <- as.data.frame(c(g_g$symbol, g_cho$symbol))
colnames(gall) <- "symbol"
oall <- as.data.frame(c(g_o$symbol, g_cgh$symbol))
colnames(oall) <- "symbol"

gwas <- subset(gwas, !is.na(MAPPED_GENE))

cgwas2 <- intersect(call$symbol, gwas$MAPPED_GENE)
cgwas21 <- gwas[gwas$MAPPED_GENE %in% cgwas2,]

ggwas2 <- intersect(gall$symbol, gwas$MAPPED_GENE)
ggwas21 <- gwas[gwas$MAPPED_GENE %in% ggwas2,]

hgwas2 <- intersect(hall$symbol, gwas$MAPPED_GENE)
hgwas21 <- gwas[gwas$MAPPED_GENE %in% hgwas2,]

ogwas2 <- intersect(oall$symbol, gwas$MAPPED_GENE)
ogwas21 <- gwas[gwas$MAPPED_GENE %in% ogwas2,]

hl1 <- c("ileocolitis, Crohn's disease", "inflammatory bowel disease", "Crohn's disease", 
         "gastric adenocarcinoma", "gastric carcinoma", "type I diabetes mellitus",
         "obesity", "Parkinson's disease", "gout", "colorectal cancer")

hinterest <- hgwas21[hgwas21$MAPPED_TRAIT %in% hl1,]

p<-ggplot(data=hinterest, aes(x=MAPPED_TRAIT, y=PVALUE_MLOG)) +
  geom_bar(stat="identity")
p


###Add gene symbols 
setwd(file.path(rootpath))
bgall <- read.table("Filtered_Data.txt", sep="\t")
bgall$symbol <- mapIds(org.Hs.eg.db,
                       keys=row.names(bgall),
                       column="SYMBOL",
                       keytype="ENSEMBL",
                       multiVals="first")
bgall$entrez <- mapIds(org.Hs.eg.db,
                       keys=row.names(bgall),
                       column="ENTREZID",
                       keytype="ENSEMBL",
                       multiVals="first")

bggwas2 <- intersect(bgall$symbol, gwas$MAPPED_GENE)
bggwas21 <- gwas[gwas$MAPPED_GENE %in% bggwas2,]


h1 <- data.frame(c(gwas = (length(hgwas2)), notgwas =(nrow(hall) - length(hgwas2))), 
                 c(length(bggwas2), nrow(bgall)-length(bggwas2)))
colnames(h1) <- c("Human", "Background")


c1 <- data.frame(c(gwas = (length(cgwas2)), notgwas =(nrow(call) - length(cgwas2))), 
                 c(length(bggwas2), nrow(bgall)-length(bggwas2)))
colnames(c1) <- c("Chimp", "Background")               

g1 <- data.frame(c(gwas = (length(ggwas2)), notgwas =(nrow(gall) - length(ggwas2))),
                 c(length(bggwas2), nrow(bgall)-length(bggwas2)))
colnames(g1) <- c("Gorilla", "Background") 


o1 <- data.frame(c(gwas = (length(ogwas2)), notgwas =(nrow(oall) - length(ogwas2))),
                 c(length(bggwas2), nrow(bgall)-length(bggwas2)))
colnames(o1) <- c("Orangutan", "Background") 

###Run fisher tests to determine enrichment
fisher.test(h1, alternative = "greater")
# Fisher's Exact Test for Count Data
# 
# data:  h1
# p-value = 8.962e-09
# alternative hypothesis: true odds ratio is greater than 1
# 95 percent confidence interval:
#  1.250746      Inf
# sample estimates:
# odds ratio 
#   1.372531 

fisher.test(c1, alternative = "greater")
# Fisher's Exact Test for Count Data
# 
# data:  c1
# p-value = 0.0002849
# alternative hypothesis: true odds ratio is greater than 1
# 95 percent confidence interval:
#  1.11152     Inf
# sample estimates:
# odds ratio 
#   1.225317 

fisher.test(g1, alternative = "greater")
# Fisher's Exact Test for Count Data
# 
# data:  g1
# p-value = 1.28e-12
# alternative hypothesis: true odds ratio is greater than 1
# 95 percent confidence interval:
#  1.235679      Inf
# sample estimates:
# odds ratio 
#   1.319059 

fisher.test(o1, alternative = "greater")
# Fisher's Exact Test for Count Data
# 
# data:  o1
# p-value = 0.0009011
# alternative hypothesis: true odds ratio is greater than 1
# 95 percent confidence interval:
#  1.090195      Inf
# sample estimates:
# odds ratio 
#   1.201355 


#########################

setwd("~/Google Drive/Primates_Coculture/Primates/MyAnalysis/Results_v08/LRT/")
read.table()


# get genes that are exclusively associated with a bunch of diseases from gwas catalog
##Crohn's disease
#subset the GWAS database
crohntrait <- c("Crohn's disease", "perianal Crohn's disease" )
crohns <- gwas[gwas$MAPPED_TRAIT %in% crohntrait,]

ibdtrait <- c("inflammatory bowel disease" )
ibd <- gwas[gwas$MAPPED_TRAIT %in% ibdtrait,]

colorectalcancertrait <- c("colorectal cancer")
colorectalcancer <- gwas[gwas$MAPPED_TRAIT %in% colorectalcancertrait,]

uctrait<- c ("ulcerative colitis")
uc <- gwas[gwas$MAPPED_TRAIT %in%  uctrait,]

typeitrait <- c("type I diabetes mellitus")
typei <- gwas[gwas$MAPPED_TRAIT %in% typeitrait,]

typeiitrait <- c("type II diabetes mellitus")
typeii <- gwas[gwas$MAPPED_TRAIT %in% typeiitrait,]

crptrait <- c("C-reactive protein measurement")
crp <- gwas[gwas$MAPPED_TRAIT %in% crptrait,]

cdd <- gwas[gwas$MAPPED_TRAIT %in% "inflammatory response",]

bmitrait <- c("body mass index")
bmi <- gwas[gwas$MAPPED_TRAIT %in% bmitrait,]

choltrait <- c("total cholesterol measurement")
chol <- gwas[gwas$MAPPED_TRAIT %in% choltrait,]

alstrait <- c("sporadic amyotrophic lateral sclerosis")
als <- gwas[gwas$MAPPED_TRAIT %in% alstrait,]

igg <-gwas[gwas$MAPPED_TRAIT %in%"serum IgG glycosylation measurement" ,] 

schizophrenia <- gwas[gwas$MAPPED_TRAIT %in%"schizophrenia" ,] 

bc <-  gwas[gwas$MAPPED_TRAIT %in%"breast carcinoma" ,] 

#make a list of dataframes for specific gwas traits
gtrait <- list(crohns, ibd, colorectalcancer, uc, typei, typeii, crp, bmi, chol, als, igg, schizophrenia, bc)


trait <- c("type I diabetes mellitus","C-reactive protein measurement" ,  "Crohn's disease", "inflammatory bowel disease" ,
           "ulcerative colitis" , "inflammatory response", "type II diabetes mellitus, autoantibody measurement","type II diabetes mellitus, metabolic syndrome" ,
           "psoriasis, Crohn's disease" , "ulcerative colitis, Crohn's disease" , "gout", "type II diabetes mellitus, obesity", 
           "stricture, Crohn's disease, complicated disease course","erythema nodosum, Crohn's disease", "ileocolitis, Crohn's disease",                                                                                                                                                                                                                                                                                                                  
           "Crohn's disease, mild disease course", "perianal Crohn's disease, complicated disease course",                                                                                                                                                                                                                                                                                             
           "perianal Crohn's disease","colorectal cancer, microsatellite instability measurement, overall survival",                                                                                                                                                                                                                                                                 
"colorectal cancer, microsatellite instability measurement", "Crohn's disease, disease prognosis measurement" 
)


###Cycle through this list of diseases

nm <- c("Crohn's disease", "Inflammatory bowel disease", "Colorectal cancer", "Ulcerative colitis", 
        "Type I diabetes", "Type II diabetes", "C-reactive protein", "Body mass index", "Total Cholesterol Measurement", "ALS", "IGG", "Schi", "Breast carcinoma")

hp <- vector()
cp <- vector()
gp <- vector()
op <- vector()

for(i in 1:length(gtrait))
{
  #First get the background genes
  bggwas2 <- intersect(bgall$symbol, gtrait[[i]]$MAPPED_GENE)
  bggwas21 <- gtrait[[i]][gtrait[[i]]$MAPPED_GENE %in% bggwas2,]
  
  #human enrichment
  hgwas2 <- intersect(hall$symbol, gtrait[[i]]$MAPPED_GENE)
  hgwas21 <- gtrait[[i]][gtrait[[i]]$MAPPED_GENE %in% hgwas2,]
  

  print(length(unique(hgwas21$MAPPED_GENE)))
  print(length(unique(bggwas21$MAPPED_GENE)))
  
  h1 <- data.frame(c(gwas = (length(hgwas2)), notgwas =(nrow(hall) - length(hgwas2))), 
                   c(length(bggwas2), nrow(bgall)-length(bggwas2)))
  colnames(h1) <- c("Human", "Background")
  htest <- fisher.test(h1, alternative = "greater")
  hp <- c(hp, htest$p.value)
  print(h1)
  
  
  #chimp enrichment
  cgwas2 <- intersect(call$symbol, gtrait[[i]]$MAPPED_GENE)
  cgwas21 <- gtrait[[i]][gtrait[[i]]$MAPPED_GENE %in% cgwas2,]
  
  c1 <- data.frame(c(gwas = (length(cgwas2)), notgwas =(nrow(call) - length(cgwas2))), 
                   c(length(bggwas2), nrow(bgall)-length(bggwas2)))
  colnames(c1) <- c("Chimp", "Background")          
  ctest <- fisher.test(c1)
  cp <- c(cp, ctest$p.value)
  
  print(c1)
  #gorilla enrichment
  ggwas2 <- intersect(gall$symbol, gtrait[[i]]$MAPPED_GENE)
  ggwas21 <- gtrait[[i]][gtrait[[i]]$MAPPED_GENE %in% ggwas2,]
  
  g1 <- data.frame(c(gwas = (length(ggwas2)), notgwas =(nrow(gall) - length(ggwas2))),
                   c(length(bggwas2), nrow(bgall)-length(bggwas2)))
  colnames(g1) <- c("Gorilla", "Background") 
  gtest <- fisher.test(g1)
  gp <- c(gp, gtest$p.value)
  
  print(g1)
  #orangutan enrichment
  ogwas2 <- intersect(oall$symbol, gtrait[[i]]$MAPPED_GENE)
  ogwas21 <- gtrait[[i]][gtrait[[i]]$MAPPED_GENE %in% ogwas2,]
  
  o1 <- data.frame(c(gwas = (length(ogwas2)), notgwas =(nrow(oall) - length(ogwas2))),
                   c(length(bggwas2), nrow(bgall)-length(bggwas2)))
  colnames(o1) <- c("Orangutan", "Background") 
  otest <- fisher.test(o1)
  op <- c(op, otest$p.value)
  
  print(o1)
}

h2 <- as.data.frame(cbind(nm, hp))
h2$hp <- as.numeric(as.character(h2$hp))
colnames(h2) <- c("Trait", "P-value")

c2 <- as.data.frame(cbind(nm, cp))
c2$cp <- as.numeric(as.character(c2$cp))
colnames(c2) <- c("Trait", "P-value")

g2 <- as.data.frame(cbind(nm, gp))
g2$gp <- as.numeric(as.character(g2$gp))
colnames(g2) <- c("Trait", "P-value")

o2 <- as.data.frame(cbind(nm, op))
o2$op <- as.numeric(as.character(o2$op))
colnames(o2) <- c("Trait", "P-value")

h2$Species <- "Human"
c2$Species <- "Chimp"
g2$Species <- "Gorilla"
o2$Species <- "Orangutan"

h2$Color <- "#00BFC4"
c2$Color <- "#F8766D"
g2$Color <- "#7CAE00"
o2$Color <- "#C77Cff"  


h2$log <- -log10(h2$`P-value`)
c2$log <- -log10(c2$`P-value`)
g2$log <- -log10(g2$`P-value`)
o2$log <- -log10(o2$`P-value`)

ggplot(data = h2, aes(x = h2$Trait, y = h2$log)) + coord_flip() + geom_bar(stat="identity", width=0.5) + theme_bw() + ylab("-log(p-value)") + xlab ("Human") + geom_hline(yintercept = -log10(0.05))
ggplot(data = c2, aes(x = c2$Trait, y = c2$log)) + coord_flip() + geom_bar(stat="identity", width=0.5) + theme_bw() + ylab("-log(p-value)") + xlab("Chimp")+ geom_hline(yintercept = -log10(0.05))
ggplot(data = g2, aes(x = g2$Trait, y = g2$log)) + coord_flip() + geom_bar(stat="identity", width=0.5) + theme_bw() + ylab("-log(p-value)") + xlab("Gorilla")+ geom_hline(yintercept = -log10(0.05))
ggplot(data = o2, aes(x = o2$Trait, y = o2$log)) + coord_flip() + geom_bar(stat="identity", width=0.5) + theme_bw() + ylab("-log(p-value)") + xlab("Orangutan")+ geom_hline(yintercept = -log10(0.05))

siglist <- rbind(h2[h2$`P-value`<0.05,], c2[c2$`P-value` <0.05,], g2[g2$`P-value`<0.05,], o2[o2$`P-value`<0.05,])


  
siglist <- transform(siglist, Trait = reorder(Trait, log))


setwd(file.path(rootpath, "GWAS"))
pdf("Gwas_enrichment.pdf")

p <- ggplot(data = siglist, aes(x = siglist$Trait, y = siglist$log)) + 
  coord_flip()+ 
  geom_bar(stat="identity", width=0.5, aes(fill = as.factor(Species))) +
  theme_bw() +
  ylab("-log(p-value)") + 
  xlab ("Trait") + 
  geom_hline(yintercept = -log10(0.05)) +
  scale_fill_manual(values=c("#F8766D",  "#7CAE00", "#00BFC4"), name = "Species")
print(p)

dev.off()
#default colors in ggplot: 
##F8766D = red = chimp
##7CAE00=green = gorilla
##00BFC4=teal = human
##C77Cff=purple = orangutan

table(gwas$MAPPED_TRAIT)



