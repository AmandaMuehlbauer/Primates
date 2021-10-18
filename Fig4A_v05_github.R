library("gplots")
library("ggplot2")
library("ComplexHeatmap")
library("circlize")
library("colorspace")
library("GetoptLong")
library("org.Hs.eg.db")
library("compositions")
#Get the data from the Pathway-Genes correlations model.
setwd("/Volumes/GoogleDrive/My Drive/Primates_Coculture/Primates/MyAnalysis_v02/Fig4Aii_2019-02-26")

sig5 <- read.table("significant_spearmanSpeciesGenes.txt",sep="\t", header = TRUE)

gen4 <- read.table("MetadataPathways.txt", sep="\t", header = TRUE, check.names = FALSE)
data <- read.table("Genes.txt", sep = "\t", header = TRUE)

###Heatmaps
##First need to make 
#Using pearson correlations

length(unique(sig4$OTU))
#30
length(unique(sig4$ENSG))
#127
sig4.1 <- sig4[sig4$pearsonpvalue<0.1,]
length(unique(sig4.1$OTU))
#25
length(unique(sig4.1$ENSG))
#96

dim(sig4.1)
# [1] 183  15
sig4.2 <- as.data.frame(unique(sig4.1$ENSG))
write.table(sig4.2, file ="taxagenes_unique.txt", sep = "\t")

#Get unique taxa and genes
sensg <- unique(sig4.1$ENSG)
sotu <- unique(sig4.1$OTU)

gg1 <- as.character(sensg)
aa1 <- as.character(sotu)

dat1 <- t(data)
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
colnames(taxapval) <- bb2$bb2
rownames(taxapval) <- aa1

library("reshape2")
pvals <- melt(taxapval)
rvals <- melt(taxacor)

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
#rownames(gs4) <- gs4$symbol
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

#metagen1.2 <- metagen1[rownames(metagen1) %in% rownames(gen1.3),]

allchimps <- c("100UNCULTURED", "400UNCULTURED", "300UNCULTURED", "500UNCULTURED")
allgorillas <- c("2UNCULTURED", "4UNCULTURED","14UNCULTURED","26UNCULTURED","40UNCULTURED","19UNCULTURED","24UNCULTURED")
allhumans <- c("1UNCULTURED","65UNCULTURED","110UNCULTURED","111UNCULTURED")
allorangutans <- c("8UNCULTURED","18UNCULTURED","36UNCULTURED","28UNCULTURED")

setwd("/Volumes/GoogleDrive/My\ Drive/Primates_Coculture/Primates/MyAnalysis")
primate_meta16s <- read.table("primate_metagenomic_meta.txt", sep = "\t", header = TRUE, row.names ="SampleID")
setwd("/Volumes/GoogleDrive/My\ Drive/Primates_Coculture/Sequencing20180530/ShotgunSequencing/Tables_Combined")
primate_genera <- t(read.table("Merged_Species.txt", sep = "\t", header = TRUE, row.names = 1, check.names = FALSE))
#primate_genera
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

rownames(taxacor) <-str_replace_all(rownames(taxacor), "[^[:alnum:]]", " ")

row_ha <- rowAnnotation(micall1 = row_anno_boxplot(micall1, 
                                                   axis = TRUE, 
                                                   axis_side = "bottom",
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
                                                       axis_side = "bottom", ylim = c(-9,5),
                                                       gp = gpar(fill = c(Chimp = "#F8766D"))), 
                            width = unit(2, "cm"),which = "row" )

draw(row_hc, 1:10)

#Make barplot for gorilla abundances
row_hg <- HeatmapAnnotation(boxplot = row_anno_boxplot(microg1, outline = FALSE,
                                                       axis = TRUE, 
                                                       axis_side = "bottom", ylim=c(-9,5),
                                                       gp = gpar(fill = c(Gorilla ="#7CAE00"))), 
                            width = unit(2, "cm"),which = "row")


draw(row_hg, 1:10)

#Make barplot for human abundances
row_hh <- HeatmapAnnotation(boxplot = row_anno_boxplot(microh1, outline=FALSE,
                                                       axis = TRUE, 
                                                       axis_side = "bottom", ylim=c(-9,5),
                                                       gp = gpar(fill = c(Human = "#00BFC4"))), 
                            width = unit(2, "cm"), which="row")

draw(row_hh, 1:10)

#Make barplot for orangutan abundances
row_ho <- HeatmapAnnotation(boxplot = row_anno_boxplot(microo1, outline = FALSE,
                                                       axis = TRUE, 
                                                       axis_side = "bottom", ylim=c(-9,5),
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
             top_annotation_height = unit(1, "cm"), 
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
             top_annotation_height = unit(1, "cm"), 
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
             top_annotation_height = unit(1, "cm"), 
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

