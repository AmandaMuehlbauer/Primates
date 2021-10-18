# this script runs {GO, KEGG, Reactome} enrichment analyses on DEGs from DESeq2, using all 17000 background genes instead of the filtered 3000
library(tidyverse)
library(knitr)
#library(DESeq2)
# the following 2 packages work only in R/3.5.2 in new environment
#install.packages("devtools")
#devtools::install_github("stephenturner/annotables")
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("clusterProfiler")

##Gene ratio is number of annotated genes in a particular node/ number of genes of interest. 
##(This has nothing to do with background list)
library(extrafont)
library(annotables)
library(clusterProfiler)
#library(qqman)

#BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)
library(ReactomePA)
library(pathview)
library("enrichplot")
library("DOSE")
library("egg")
library(reshape2)
library(gridExtra)
library("ggforce")

font_import()
#set directory and load data
setwd("/Volumes/GoogleDrive/My\ Drive/Primates_Coculture/RNAseq_Gopher/R_workspace")
load("subread_DESeq2-data.RDATA")

###This is the protein coding filtered dataset file
###The list of genes was finalized on 10/09/2018. 
###There should be 19714 genes and 32 samples. 
###One sample will be removed because it seems to be an outlier
# setwd("/Volumes/GoogleDrive/My\ Drive/Primates_Coculture/Primates/MyAnalysis/Results_2018-05-02")
# data <- read.table("Filtered_Data.txt", sep="\t", check.names = FALSE)
# dim(data)
###Make directories that change each time I run this script
# x <- getwd()  ## I'm going to demo in a tempdir
# setwd(tempdir())
# list.dirs()
# # [1] "."
# # [2] "./downloaded_packages"
# # [3] "./rs-graphics-16e13b20-59b3-4ef3-bdcd-02852b1ea576"
# newdir <- paste("Primates", Sys.Date(), sep = "_")
# dir.create(newdir)
# setwd(newdir)
# list.dirs()
# # [1] "."
# # [2] "./downloaded_packages"
# # [3] "./rs-graphics-16e13b20-59b3-4ef3-bdcd-02852b1ea576"
# # [4] "./Test_2014-09-03"
# setwd(x)     ## Reset to original working directory
# 


setwd("/Volumes/GoogleDrive/My\ Drive/Primates_Coculture/Primates/MyAnalysis_v02/Results_2019-07-22/")
background <-read.table("DE_allgenesModel1.txt", sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)

####Read in the Treatment vs control genes
###Using all ~17000 genes that are detected for the background genes
# setwd("/Volumes/GoogleDrive/My\ Drive/Primates_Coculture/Primates/MyAnalysis_v02/Results_2018-11-06/DE_Model1")
# background <-read.table("Model1_v09.txt", sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
background <- as.character(background$symbol)
background.df <- bitr(background, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
background <- as.character(background.df[,2])


setwd("/Volumes/GoogleDrive/My\ Drive/Primates_Coculture/Primates/MyAnalysis_v02/PrimateComparev01_2019-07-01")
high_list <- read.table("HighGenes_All.txt", sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)

high_list <- as.character(high_list$symbol)
DEGs.df <- bitr(high_list, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
DEGs <- as.character(DEGs.df[,2]) 


# GO over-representation test (takes a couple minutes):
GO_h <- enrichGO(DEGs, 
               OrgDb="org.Hs.eg.db", 
               keyType = "ENTREZID", 
               ont = "BP",
               pvalueCutoff = 0.1, 
               pAdjustMethod = "BH", 
               universe=background,
               qvalueCutoff = 0.1, 
               minGSSize = 10, 
               maxGSSize = 500,
               readable = FALSE, 
               pool = FALSE)
# KEGG over-representation test:
KEGG_h <- enrichKEGG(gene=DEGs, 
                   organism='hsa', 
                   keyType = "kegg",
                   pvalueCutoff = 0.1, 
                   pAdjustMethod = "BH", 
                   universe=background,
                   qvalueCutoff = 0.1, 
                   minGSSize = 10) 

# REACTOME analysis:
REACT_h <- enrichPathway(gene=DEGs, 
                       organism="human",
                       pvalueCutoff = 0.1, 
                       pAdjustMethod = "BH", 
                       universe=background,
                       qvalueCutoff=0.1, 
                       minGSSize=10)


systempath <- "/Volumes/GoogleDrive/My\ Drive/Primates_Coculture/Primates/MyAnalysis_v02"
setwd(systempath)
dir.create(file.path(systempath, format(Sys.time(), "ClusterProfiler_17000_%F")))

rootpath <- paste0(file.path(systempath, format(Sys.time(), "ClusterProfiler_17000_%F")))
setwd(paste0(file.path(systempath, format(Sys.time(), "ClusterProfiler_17000_%F"))))

#save the results:
write.table(GO_h@result, file=paste0("GO_high.txt"), sep="\t", quote=FALSE)
write.table(KEGG_h@result, file=paste0("KEGG_high.txt"), sep="\t", quote=FALSE)
write.table(REACT_h@result, file=paste0("REACTOM_high.txt"), sep="\t", quote=FALSE)

# combined + shortened:
newGO_h <- GO_h@result[GO_h@result$Count>1 & GO_h@result$qvalue<0.1,]
newKEGG_h <- KEGG_h@result[KEGG_h@result$Count>1 & KEGG_h@result$qvalue<0.1,]
newREACT_h <- REACT_h@result[REACT_h@result$Count>1 & REACT_h@result$qvalue<0.1,]

merged <- rbind(newGO_h[0:min(10,nrow(newGO_h)),], newKEGG_h[0:min(10,nrow(newKEGG_h)),], newREACT_h[0:min(10,nrow(newREACT_h)),])
# write.table(merged, file=paste0("./results/shortened/", 
#                                 myvar, "_enrichment_summary_all.txt"), 
#             sep="\t",
#             quote=FALSE, 
#             row.names=FALSE)



de <- names(DEGs)[abs(DEGs) > 2]


if (nrow(head(GO_h))>0){
  pdf("./GO_dotplot_highDEGs.pdf")   
  a<-dotplot(GO_h)
  print(a)
  dev.off()
  pdf("./GO_emapplot_highDEGs.pdf")
  b<-emapplot(GO_h)
  print(b)
  dev.off()
  pdf("./GO_cnetplot_highDEGs.pdf")
  c<-cnetplot(GO_h, categorySize="pvalue", foldChange=DEGs)
  print(c)
  dev.off()
  pdf("./GO_goplot_highDEGs.pdf")
  d<-goplot(GO_h)
  print(d)
  dev.off()
}

if (nrow(head(KEGG_h))>0){
  pdf("./KEGG_dotplot_highDEGs.pdf")
  e<-dotplot(KEGG_h)
  print(e)
  dev.off()
  pdf("./KEGG_emapplot_highDEGs.pdf")
  f<-emapplot(KEGG_h)
  print(f)
  dev.off()
  pdf("./KEGG_cnetplot_highDEGs.pdf")
  g<-cnetplot(KEGG_h, categorySize="pvalue", foldChange=DEGs)
  print(g)
  dev.off()
  for (i in 1:nrow(newKEGG_h)){
    pathview(gene.data  =DEGs, pathway.id =newKEGG_h[i,1], species    = "hsa", limit= list(gene=1, cpd=1))
  }
}

if (nrow(head(REACT_h))>0){
  pdf("./REACT_dotplot_highDEGs.pdf")
  h<-dotplot(REACT_h)
  print(h)
  dev.off()
  pdf("./REACT_emapplot_highDEGs.pdf")
  i<-emapplot(REACT_h)
  print(i)
  dev.off()
  pdf("./REACT_cnetplot_highDEGs.pdf")
  j<-cnetplot(REACT_h, categorySize="pvalue", foldChange=DEGs)
  print(j)
  dev.off()
}

pdf("testdotplot.pdf")
p <-clusterProfiler::dotplot(GO_h, showCategory=20)
print(p)
dev.off()


clusterProfiler::dotplot(GO_h)


##Capitalize the first letter of the categories
CapStr <- function(y) {
  c <- strsplit(y, " ")[[1]]
  if(c =="mRNA"){
    
  }
  paste(toupper(substring(c, 1,1)), substring(c, 2),
        sep="", collapse=" ")
}


GO_h@result$Description_cap<-sapply(GO_h@result$Description, CapStr)
KEGG_h@result$Description_cap<-sapply(KEGG_h@result$Description, CapStr)
REACT_h@result$Description_cap<-sapply(REACT_h@result$Description, CapStr)


####Remake the figures so they are nicer
pdf("GO_high_fig.pdf")
gohigh <-ggplot(head(GO_h, 10), aes(x =GeneRatio, y=reorder(Description_cap, Count))) + 
  theme_bw()+ 
  geom_point(aes(size = Count, color = p.adjust))+  
  labs(color = "Q-Value", size = "Count", x="Gene Ratio", y = " ") + 
  theme(axis.text.y = element_text(size = 8)) +
  scale_y_discrete(position="right") + 
  scale_size(guide=guide_legend(reverse=TRUE))  +
  scale_color_gradientn(colours = c("#7aa6c2", "#004c6d")) +
scale_size_continuous(name = "Count",
                      breaks = c(0, 10, 20, 30))
print(gohigh)
dev.off()


##Plot KEGG high
pdf("KEGG_high_fig.pdf")
kegghigh <-ggplot(head(KEGG_h, 10), aes(x =GeneRatio, y=reorder(Description_cap, Count))) + 
  theme_bw()+ 
  geom_point(aes(size = Count, color = p.adjust))+  
  labs(color = "Q-Value", size = "Count", x="Gene Ratio", y = " ") + 
  theme(axis.text.y = element_text(size = 8)) +
  scale_y_discrete(position="right") + 
  scale_size(guide=guide_legend(reverse=TRUE))  +
  scale_color_gradientn(colours = c("#7aa6c2", "#004c6d"))
print(kegghigh)
dev.off()

###Plot REACTOME
pdf("REACT_high_fig.pdf")
reacthigh <-ggplot(head(REACT_h, 10), aes(x =GeneRatio, y=reorder(Description_cap, Count))) + 
  theme_bw()+ 
  geom_point(aes(size = Count, color = p.adjust))+  
  labs(color = "Q-Value", size = "Count", x="Gene Ratio", y = " ") + 
  theme(axis.text.y = element_text(size = 8)) +
  scale_y_discrete(position="right") + 
  scale_size(guide=guide_legend(reverse=TRUE))  +
  scale_color_gradientn(colours = c("#7aa6c2", "#004c6d"))
print(reacthigh)
dev.off()




dotplot(KEGG_h) + scale_color_gradient(low="dodgerblue", high="salmon") 

dotplot(REACT_h) + scale_color_gradient(low="darkblue", high="salmon") 






###Do same analysis for the low genes
setwd("/Volumes/GoogleDrive/My\ Drive/Primates_Coculture/Primates/MyAnalysis_v02/PrimateComparev01_2019-07-01")
low_list <- read.table("LowGenes_All.txt", sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)

low_list <- as.character(low_list$symbol)
DEGs.df <- bitr(low_list, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
DEGs <- as.character(DEGs.df[,2]) 
# GO over-representation test (takes a couple minutes):
GO_l <- enrichGO(DEGs, 
               OrgDb="org.Hs.eg.db", 
               keyType = "ENTREZID", 
               ont = "BP",
               pvalueCutoff = 0.1, 
               pAdjustMethod = "BH", 
               universe=background,qvalueCutoff = 0.1, 
               minGSSize = 10, 
               maxGSSize = 500,
               readable = FALSE, 
               pool = FALSE)
# KEGG over-representation test:
KEGG_l <- enrichKEGG(gene=DEGs, 
                   organism='hsa', 
                   keyType = "kegg",
                   pvalueCutoff = 0.1, 
                   pAdjustMethod = "BH", 
                   universe=background,
                   qvalueCutoff = 0.1, 
                   minGSSize = 10) 

# REACTOME analysis:
REACT_l <- enrichPathway(gene=DEGs, 
                       organism="human",
                       pvalueCutoff = 0.1, 
                       pAdjustMethod = "BH", 
                       universe=background,qvalueCutoff=0.1, minGSSize=10)



setwd(paste0(file.path(systempath, format(Sys.time(), "ClusterProfiler_17000_%F"))))

#save the results:
write.table(GO_l@result, file=paste0("GO_low.txt"), sep="\t", quote=FALSE)
write.table(KEGG_l@result, file=paste0("KEGG_low.txt"), sep="\t", quote=FALSE)
write.table(REACT_l@result, file=paste0("REACTOM_low.txt"), sep="\t", quote=FALSE)

# combined + shortened:
newGO_l <- GO_l@result[GO_l@result$Count>1 & GO_l@result$qvalue<0.1,]
newKEGG_l <- KEGG_l@result[KEGG_l@result$Count>1 & KEGG_l@result$qvalue<0.1,]
newREACT_l <- REACT_l@result[REACT_l@result$Count>1 & REACT_l@result$qvalue<0.1,]

merged <- rbind(newGO_l[0:min(10,nrow(newGO_l)),], newKEGG_l[0:min(10,nrow(newKEGG_l)),], newREACT_l[0:min(10,nrow(newREACT_l)),])
# write.table(merged, file=paste0("./results/shortened/", 
#                                 myvar, "_enrichment_summary_all.txt"), 
#             sep="\t",
#             quote=FALSE, 
#             row.names=FALSE)



de <- names(DEGs)[abs(DEGs) > 2]


if (nrow(head(GO_l))>0){
  pdf("./GO_dotplot_lowDEGs.pdf")   
  a<-dotplot(GO_l) 
  print(a)
  dev.off()
  pdf("./GO_emapplot_lowDEGs.pdf")
  b<-emapplot(GO_l)
  print(b)
  dev.off()
  pdf("./GO_cnetplot_lowDEGs.pdf")
  c<-cnetplot(GO_l, categorySize="pvalue", foldChange=DEGs)
  print(c)
  dev.off()
  pdf("./GO_goplot_lowDEGs.pdf")
  d<-goplot(GO_l)
  print(d)
  dev.off()
}

if (nrow(head(KEGG_l))>0){
  pdf("./KEGG_dotplot_lowDEGs.pdf")
  e<-dotplot(KEGG_l)
  print(e)
  dev.off()
  pdf("./KEGG_emapplot_lowDEGs.pdf")
  f<-emapplot(KEGG_l)
  print(f)
  dev.off()
  pdf("./KEGG_cnetplot_lowDEGs.pdf")
  g<-cnetplot(KEGG_l, categorySize="pvalue", foldChange=DEGs)
  print(g)
  dev.off()
  for (i in 1:nrow(newKEGG_l)){
    pathview(gene.data  =DEGs, pathway.id =newKEGG_l[i,1], species    = "hsa", limit= list(gene=1, cpd=1))
  }
}

if (nrow(head(REACT_l))>0){
  pdf("./REACT_dotplot_lowDEGs.pdf")
  h<-dotplot(REACT_l)
  print(h)
  dev.off()
  pdf("./REACT_emapplot_lowDEGs.pdf")
  i<-emapplot(REACT_l)
  print(i)
  dev.off()
  pdf("./REACT_cnetplot_lowDEGs.pdf")
  j<-cnetplot(REACT_l, categorySize="pvalue", foldChange=DEGs)
  print(j)
  dev.off()
}

pdf("testdotplot.pdf")
p <-clusterProfiler::dotplot(GO_l, showCategory=20)
print(p)
dev.off()


GO_l@result$Description_cap<-sapply(GO_l@result$Description, CapStr)
KEGG_l@result$Description_cap<-sapply(KEGG_l@result$Description, CapStr)
REACT_l@result$Description_cap<-sapply(REACT_l@result$Description, CapStr)



####Remake the figures so they are nicer
pdf("GO_low_fig.pdf")
golow <-ggplot(head(GO_l, 10), aes(x =GeneRatio, y=reorder(Description_cap, Count))) + 
  theme_bw()+ 
  geom_point(aes(size = Count, color = p.adjust))+  
  labs(color = "Q-Value", size = "Count", x="Gene Ratio", y = " ") + 
  theme(axis.text.y = element_text(size = 8)) +
  scale_y_discrete(position="right") + 
  scale_size(guide=guide_legend(reverse=TRUE))  +
  scale_color_gradientn(colours = c("#7aa6c2", "#004c6d"))
print(golow)
dev.off()


##Plot KEGG high
pdf("KEGG_low_fig.pdf")
kegglow <-ggplot(head(KEGG_l, 10), aes(x =GeneRatio, y=reorder(Description_cap, Count))) + 
  theme_bw()+ 
  geom_point(aes(size = Count, color = p.adjust))+  
  labs(color = "Q-Value", size = "Count", x="Gene Ratio", y = " ") + 
  theme(axis.text.y = element_text(size = 8)) +
  scale_y_discrete(position="right") + 
  scale_size(guide=guide_legend(reverse=TRUE))  +
  scale_color_gradientn(colours = c("#7aa6c2", "#004c6d"))
print(kegglow)
dev.off()

###Plot REACTOME
pdf("REACT_low_fig.pdf")
reactlow <-ggplot(head(REACT_l, 10), aes(x =GeneRatio, y=reorder(Description_cap, Count))) + 
  theme_bw()+ 
  geom_point(aes(size = Count, color = p.adjust))+  
  labs(color = "Q-Value", size = "Count", x="Gene Ratio", y = " ") + 
  theme(axis.text.y = element_text(size = 8)) +
  scale_y_discrete(position="right") + 
  scale_size(guide=guide_legend(reverse=TRUE))  +
  scale_color_gradientn(colours = c("#7aa6c2", "#004c6d"))
print(reactlow)
dev.off()

grid.arrange(gohigh, golow, ncol=1)

g2 <- ggplotGrob(gohigh)
g3 <- ggplotGrob(golow)
g <- rbind(g2, g3, size = "first")
g$widths <- unit.pmax(g2$widths, g3$widths)
grid.newpage()
grid.draw(g)

g2 <- ggplotGrob(kegghigh)
g3 <- ggplotGrob(kegglow)
g <- rbind(g2, g3, size = "first")
g$widths <- unit.pmax(g2$widths, g3$widths)
grid.newpage()
grid.draw(g)

GO_h@result$Div <- "High"
GO_l@result$Div <- "Low"

KEGG_h@result$Div <- "High"
KEGG_l@result$Div <- "Low"

REACT_h@result$Div <- "High"
REACT_l@result$Div <- "Low"


##Make plots so that they are somewhat standard

GO_20 <- rbind(head(GO_h, 10), head(GO_l, 10) )
KEGG_20 <- rbind(head(KEGG_h, 10), head(KEGG_l, 10) )
REACT_20 <- rbind(head(REACT_h, 10), head(REACT_l, 10) )


GO_10 <- rbind(head(GO_h, 5), head(GO_l, 5) )
KEGG_10 <- rbind(head(KEGG_h, 5), head(KEGG_l, 5) )
REACT_10 <- rbind(head(REACT_h, 5), head(REACT_l, 5) )

GO_20$GeneRatioDec <- sapply(GO_20$GeneRatio, function(x) eval(parse(text=x)))
KEGG_20$GeneRatioDec <- sapply(KEGG_20$GeneRatio, function(x) eval(parse(text=x)))
REACT_20$GeneRatioDec <- sapply(REACT_20$GeneRatio, function(x) eval(parse(text=x)))

GO_10$GeneRatioDec <- sapply(GO_10$GeneRatio, function(x) eval(parse(text=x)))
KEGG_10$GeneRatioDec <- sapply(KEGG_10$GeneRatio, function(x) eval(parse(text=x)))
REACT_10$GeneRatioDec <- sapply(REACT_10$GeneRatio, function(x) eval(parse(text=x)))

GO_20$Database <- "GO"
KEGG_20$Database <- "KEGG"
REACT_20$Database <- "REACTOME"

GO_10$Database <- "GO"
KEGG_10$Database <- "KEGG"
REACT_10$Database <- "REACTOME"

all60 <- rbind(GO_20, KEGG_20, REACT_20)

all30 <- rbind(GO_10, KEGG_10, REACT_10)

pdf("GO_enrich.pdf")
go <-ggplot(GO_20, aes(x =GeneRatioDec, y=reorder(Description_cap, Count))) + 
  theme_bw()+ 
  geom_point(aes(size = Count, color = p.adjust))+  
  labs(color = "Q-Value", size = "Count", x="Gene Ratio", y = " ") + 
  theme(axis.text.y = element_text(size = 8)) +
  scale_y_discrete(position="right") + 
  scale_size(guide=guide_legend(reverse=TRUE))  +
  scale_color_gradientn(colours = c("#7aa6c2", "#004c6d"))+ facet_grid(Div ~ .,scales="free")
print(go)
dev.off()

pdf("GO_enrich_nolabel.pdf")
go <-ggplot(GO_20, aes(x =GeneRatioDec, y=reorder(Description_cap, Count))) + 
  theme_bw()+ 
  geom_point(aes(size = Count, color = p.adjust))+  
  labs(color = "Q-Value", size = "Count", x="Gene Ratio", y = " ") + 
  theme(axis.text.y = element_text(size = 8)) +
  scale_y_discrete(position="right") + 
  scale_size(guide=guide_legend(reverse=TRUE))  +
  scale_color_gradientn(colours = c("#008080", "#ccffff"))+ facet_grid(Div ~ .,scales="free") +
  theme(axis.title.x=element_blank(),
        axis.text.y = element_blank())
print(go)
dev.off()

pdf("KEGG_enrich.pdf")
kegg <- ggplot(KEGG_20, aes(x =GeneRatioDec, y=reorder(Description_cap, Count))) + 
  theme_bw()+ 
  geom_point(aes(size = Count, color = p.adjust))+  
  labs(color = "Q-Value", size = "Count", x="Gene Ratio", y = " ") + 
  theme(axis.text.y = element_text(size = 8)) +
  scale_y_discrete(position="right") + 
  scale_size(guide=guide_legend(reverse=TRUE))  +
  scale_color_gradientn(colours = c("#7aa6c2", "#004c6d"))+ facet_grid(Div ~ .,scales="free")
print(kegg)
dev.off()


pdf("KEGG_enrich_nolabel.pdf")
kegg <- ggplot(KEGG_20, aes(x =GeneRatioDec, y=reorder(Description_cap, Count))) + 
  theme_bw()+ 
  geom_point(aes(size = Count, color = p.adjust))+  
  labs(color = "Q-Value", size = "Count", x="Gene Ratio", y = " ") + 
  theme(axis.text.y = element_text(size = 8)) +
  scale_y_discrete(position="right") + 
  scale_size(guide=guide_legend(reverse=TRUE))  +
  scale_color_gradientn(colours = c("#b34700", "#fff0e6"))+ facet_grid(Div ~ .,scales="free") +
  theme(axis.title.x=element_blank(),
        axis.text.y = element_blank())
print(kegg)
dev.off()

pdf("REACT_enrich.pdf")
react <- ggplot(REACT_20, aes(x =GeneRatioDec, y=reorder(Description_cap, Count))) + 
  theme_bw()+ 
  geom_point(aes(size = Count, color = p.adjust))+  
  labs(color = "Q-Value", size = "Count", x="Gene Ratio", y = " ") + 
  theme(axis.text.y = element_text(size = 8)) +
  scale_y_discrete(position="right") + 
  scale_size(guide=guide_legend(reverse=TRUE))  +
  scale_color_gradientn(colours = c("#7aa6c2", "#004c6d"))+ facet_grid(Div ~ .,scales="free")
print(react)
dev.off()



pdf("REACT_enrich_nolabel.pdf")
react <- ggplot(REACT_20, aes(x =GeneRatioDec, y=reorder(Description_cap, Count))) + 
  theme_bw()+ 
  geom_point(aes(size = Count, color = p.adjust))+  
  labs(color = "Q-Value", size = "Count", x="Gene Ratio", y = " ") + 
  theme(axis.text.y = element_text(size = 8)) +
  scale_y_discrete(position="right") + 
  scale_size(guide=guide_legend(reverse=TRUE))  +
  scale_color_gradientn(colours = c("#191966", "#eaeafa"))+ facet_grid(Div ~ .,scales="free") +
  theme(axis.title.x=element_blank(),
        axis.text.y = element_blank())
print(react)
dev.off()

pdf("all_enrich.pdf")
all <- ggplot(all60, aes(x =GeneRatioDec, y=reorder(Description_cap, Count))) + 
  theme_bw()+ 
  geom_point(aes(size = Count, color = p.adjust))+  
  labs(color = "Q-Value", size = "Count", x="Gene Ratio", y = " ") + 
  theme(axis.text.y = element_text(size = 8)) +
  scale_y_discrete(position="right") + 
  scale_size(guide=guide_legend(reverse=TRUE))  +
  scale_color_gradientn(colours = c("#154360", "#A9CCE3", "#D5D8DC"))+ facet_grid(Database + Div ~ .,scales="free")
print(all)
dev.off()

pdf("all_nolab.pdf")
all_nl <-ggplot(all60, aes(x =GeneRatioDec, y=reorder(Description_cap, Count))) + 
  theme_bw()+ 
  geom_point(aes(size = Count, color = p.adjust))+  
  labs(color = "Q-Value", size = "Count", x="Gene Ratio", y = " ") + 
  scale_y_discrete(position="right") + 
  scale_size(guide=guide_legend(reverse=TRUE))  +
  scale_color_gradientn(colours = c("#154360", "#A9CCE3", "#D5D8DC"))+ facet_grid(Database + Div ~ .,scales="free")+ 
  theme(axis.title.x=element_blank(),
        axis.text.y=element_blank())
print(all_nl)
dev.off()


pdf("all30_nolab.pdf")
all_nl <-ggplot(all30, aes(x =GeneRatioDec, y=reorder(Description_cap, Count))) + 
  theme_bw()+ 
  geom_point(aes(size = Count, color = p.adjust))+  
  labs(color = "Q-Value", size = "Count", x="Gene Ratio", y = " ") + 
  scale_y_discrete(position="right") + 
  scale_size(guide=guide_legend(reverse=TRUE))  +
  scale_color_gradientn(colours = c("#154360", "#A9CCE3", "#D5D8DC"))+  
  facet_wrap_paginate(Database + Div ~ .,scales="free", ncol = 1,page=3, strip.position = "right") +
  theme(axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank())
print(all_nl)
dev.off()

pdf("all30_lab.pdf")
all_nl <-ggplot(all30, aes(x =GeneRatioDec, y=reorder(Description_cap, Count))) + 
  theme_bw()+ 
  geom_point(aes(size = Count, color = p.adjust))+  
  labs(color = "Q-Value", size = "Count", x="Gene Ratio", y = " ") + 
  scale_y_discrete(position="right") + 
  scale_size(guide=guide_legend(reverse=TRUE))  +
  scale_color_gradientn(colours = c("#154360", "#A9CCE3", "#D5D8DC"))+  
  facet_wrap_paginate(Database + Div ~ .,scales="free", ncol = 1,page=3, strip.position = "right") +
  theme(axis.title.x=element_blank())
print(all_nl)
dev.off()

kr1 <- rbind(KEGG_20, REACT_20)

# pdf("keggreact_lab.pdf")
# all_nl <-ggplot(kr1, aes(x =GeneRatioDec, y=reorder(Description_cap, Count))) + 
#   theme_bw()+ 
#   geom_point(aes(size = Count, color = p.adjust))+  
#   labs(color = "Q-Value", size = "Count", x="Gene Ratio", y = " ") + 
#   scale_y_discrete(position="right") + 
#   scale_size(guide=guide_legend(reverse=TRUE))  +
#   scale_color_gradientn(colours = c("#154360", "#A9CCE3", "#D5D8DC"))+  
#   facet_wrap_paginate(Database + Div ~ .,scales="free", ncol = 1,page=3, strip.position = "right") +
#   theme(axis.title.x=element_blank(),axis.ticks.x=element_blank())
# print(all_nl)
# dev.off()




react <- ggplot(REACT_20, aes(x =GeneRatioDec, y=reorder(Description_cap, Count))) + 
  theme_bw()+ 
  geom_point(aes(size = Count, color = p.adjust))+  
  labs(color = "Q-Value", size = "Count", x="Gene Ratio", y = " ") + 
  theme(axis.text.y = element_text(size = 8), panel.border = element_rect(colour = "black", size=2)) +
  scale_y_discrete(position="right") + 
  scale_size(guide=guide_legend(reverse=TRUE))  +
  scale_color_gradientn(colours = c("#008080", "#ccffff"))+ facet_grid(Div ~ .,scales="free") +
  theme(axis.title.x=element_blank(), 
        axis.text.y=element_blank())


kegg <- ggplot(KEGG_20, aes(x =GeneRatioDec, y=reorder(Description_cap, Count))) + 
  theme_bw()+ 
  geom_point(aes(size = Count, color = p.adjust))+  
  labs(color = "Q-Value", size = "Count", x="Gene Ratio", y = " ") + 
  theme(axis.text.y = element_text(size = 8),panel.border = element_rect(colour = "black", size=2)) +
  scale_y_discrete(position="right") + 
  scale_size(guide=guide_legend(reverse=TRUE))  +
  scale_color_gradientn(colours = c("#991f00","#ffd6cc"))+ facet_grid(Div ~ .,scales="free", space = "free_y")+
  theme(axis.title.x=element_blank(), 
        axis.text.y=element_blank())

pdf("kegg_react_v2.pdf")
keggreact <- grid.arrange(grobs = lapply(list(kegg, react), set_panel_size,
                       width = unit(2.5, "in"),height = unit(1.5, "in")))
print(keggreact)
dev.off()

pdf("kegg_react_test.pdf")
kr3 <-ggarrange(kegg, react, widths = c(0.5))
print(kr3)
dev.off()


kegg <- ggplot_gtable(ggplot_build(kegg))
react <- ggplot_gtable(ggplot_build(react))


maxWidth = unit.pmax(kegg$widths[2:3], react$widths[2:3])

kegg$widths[2:3] <- maxWidth
react$widths[2:3] <- maxWidth

grid.arrange(kegg, react, heights = c(3, 2))


ggplot(KEGG_20, aes(x =-log10(qvalue), y=reorder(Description_cap, Count), colour = factor(Div))) + 
  theme_bw()+ 
  geom_point(aes(size = Count))+  
  labs(color = "Q-Value", size = "Count", x="Gene Ratio", y = " ") + 
  theme(axis.text.y = element_text(size = 8),panel.border = element_rect(colour = "black", size=2)) +
  scale_y_discrete(position="right") + 
  scale_size(guide=guide_legend(reverse=TRUE))  +
  facet_grid(Div ~ .,scales="free", space = "free_y")+
  theme(axis.title.x=element_blank(), 
        axis.text.y=element_blank())


react <- ggplot(REACT_20, aes(x =-log10(qvalue), y=reorder(Description_cap, Count), colour = factor(Div))) + 
  theme_bw()+ 
  geom_point(aes(size = Count ))+  
  labs(color = "Q-Value", size = "Count", x="Gene Ratio", y = " ") + 
  theme(axis.text.y = element_text(size = 8), 
        panel.border = element_rect(colour = "black", size=2)) +
  scale_y_discrete(position="right") + 
  scale_size(guide=guide_legend(reverse=TRUE))  +
  facet_grid(Div ~ .,scales="free", space = "free") +
  scale_color_manual(values = c("#991f00", "#008080"))+
  theme(axis.title.x=element_blank(),
        axis.text.y=element_blank()) + 
  guides(colour = "none")


kegg <- ggplot(KEGG_20, aes(x =-log10(qvalue), y=reorder(Description_cap, Count), colour = factor(Div))) + 
  theme_bw()+ 
  geom_point(aes(size = Count ))+  
  labs(color = "Q-Value", size = "Count", x="Gene Ratio", y = " ") + 
  theme(axis.text.y = element_text(size = 8), panel.border = element_rect(colour = "black", size=2)) +
  scale_y_discrete(position="right") + 
  scale_size(guide=guide_legend(reverse=TRUE))  +
  facet_grid(Div ~ .,scales="free", space = "free") +
  scale_color_manual(values = c("#991f00", "#008080"))+
  theme(axis.title.x=element_blank(), 
        axis.text.y=element_blank())

kegg <- ggplot_gtable(ggplot_build(kegg))
react <- ggplot_gtable(ggplot_build(react))



maxWidth = unit.pmax(kegg$widths[2:3], react$widths[2:3])

kegg$widths[2:3] <- maxWidth
react$widths[2:3] <- maxWidth

grid.arrange(kegg, react, heights = c(1))

grid.arrange(kegg, react, ncol = 1)


pdf("test1.pdf")
p <-  ggplot(kr1, aes(x =-log10(qvalue),
                        y=reorder(Description_cap, Count), 
                        colour = factor(Div))) + 
  theme_bw()+ 
  geom_point(aes(size = Count ))+  
  labs(color = "Q-Value", size = "Count", x="Gene Ratio", y = " ") + 
  theme(axis.text.y = element_text(size = 8), 
        panel.border = element_rect(colour = "black", size=2)) +
  scale_y_discrete(position="right") + 
  scale_size(guide=guide_legend(reverse=TRUE))  +
  facet_grid(Database~.,scales="free", 
             space = "free") +
  scale_color_manual(values = c("#991f00", "#008080"))+
  theme(axis.title.x=element_blank(),
        text = element_text(family = "Arial"), 
        axis.text.y =element_blank())
print(p)
dev.off()

divs <- c("test", "Low", "High", "Low")
names(divs)<- c("High", "Low", "High", "Low")

dats <- c(" "," "," ", " " )
names(dats) <- c("KEGG", "KEGG", "REACTOME", "REACTOME")

ggplot(kr1, aes(x =-log10(qvalue),
                y=reorder(Description_cap, Count), 
                colour = factor(Div))) + 
  theme_bw()+ 
  geom_point(aes(size = Count ))+  
  labs(color = "Q-Value", size = "Count", x="Gene Ratio", y = " ") + 
  theme(axis.text.y = element_text(size = 8), 
        panel.border = element_rect(colour = "black", size=2)) +
  scale_y_discrete(position="right") + 
  scale_size(guide=guide_legend(reverse=TRUE))  +
  facet_grid(Database + Div~.,scales="free", 
             space = "free",  
             label_value(divs, multi_line = FALSE)) +
  scale_color_manual(values = c("#991f00", "#008080"))+
  theme(axis.title.x=element_blank(),
        axis.text.y =element_blank()) 





  pdf("test2.pdf")
p <- ggplot(kr1, aes(x =-log10(qvalue),
                y=reorder(Description_cap, Count), 
                colour = factor(Div))) + 
  theme_bw()+ 
  geom_point(aes(size = Count ))+  
  labs(color = "Q-Value", size = "Count", x="Gene Ratio", y = " ") + 
  theme(axis.text.y = element_text(size = 8), 
        panel.border = element_rect(colour = "black", size=2)) +
  scale_y_discrete(position="right") + 
  scale_size(guide=guide_legend(reverse=TRUE))  +
  facet_grid(Database + Div~.,scales="free", 
             space = "free",  labeller = divs) +
  scale_color_manual(values = c("#991f00", "#008080"))+
  theme(axis.title.x=element_blank(),
        axis.text.y =element_blank())

print(p)
dev.off()


pdf("test_3.pdf", width = 12, height = 13)
krp2 <-ggplot(kr1, aes(x =-log10(qvalue),
                y=reorder(Description_cap, Count), 
                colour = factor(Div))) + 
  theme_bw()+ 
  geom_point(aes(size = Count ))+  
  labs(color = "Divergence", size = "Count", x="-log(Q-value)", y = " ") + 
  theme(axis.text.y = element_text(size = 20), 
        panel.border = element_rect(colour = "black", size=2)) +
  scale_y_discrete(position="right") + 
  scale_size(guide=guide_legend(reverse=TRUE))  +
  facet_grid(Database + Div~.,scales="free", 
             space = "free") +
  scale_color_manual(values = c("#991f00", "#008080"))+
  theme(axis.title.x=element_text(size=20), 
        strip.text.y = element_blank(), 
        axis.text.x = element_text(size=20)) 
print(krp2)
dev.off()

ggplot(kr1, aes(x =-log10(qvalue),
                y=reorder(Description, Count), 
                colour = factor(Div))) + 
  theme_bw()+ 
  geom_point(aes(size = Count ))+  
  labs(color = "Divergence", size = "Count", x="-log(Q-value)", y = " ") + 
  theme(axis.text.y = element_text(size = 20), 
        panel.border = element_rect(colour = "black", size=2),
        panel.spacing = unit(0.7, "lines")) +
  scale_y_discrete(position="right") + 
  scale_size(guide=guide_legend(reverse=TRUE))  +
  facet_grid(Database + Div~.,scales="free", 
             space = "free") +
  scale_color_manual(values = c("#991f00", "#008080"), guide = "none")+
  theme(axis.title.x=element_text(size=20), 
        strip.text.y = element_blank(), 
        axis.text.x = element_text(size=20)) 

pdf("test4.pdf", width = 12, height = 13)

ppp1 <- ggplot(kr1, aes(x =-log10(qvalue),
                y=reorder(Description, Count), 
                colour = factor(Div))) + 
  theme_bw()+ 
  geom_point(aes(size = Count ))+  
  labs(color = "Divergence", size = "Count", x="-log(Q-value)", y = " ") + 
  theme(axis.text.y = element_text(size = 20), 
        panel.border = element_rect(colour = "black", size=2),
        panel.spacing = unit(0.7, "lines")) +
  scale_y_discrete(position="right") + 
  scale_size(guide=guide_legend(reverse=TRUE))  +
  facet_grid(Database + Div~.,scales="free", 
             space = "free") +
  scale_color_manual(values = c("#991f00", "#008080"), guide = "none")+
  theme(axis.title.x=element_text(size=20), 
        strip.text.y = element_blank(), 
        axis.text.x = element_text(size=20)) 

print(ppp1)
dev.off()

pdf("test5.pdf",  width = 12, height = 13)
pp2 <-ggplot(kr1, aes(x =-log10(qvalue),
                y=reorder(Description, -log10(qvalue)), 
                colour = factor(Div))) + 
  theme_bw()+ 
  geom_point(aes(size = Count ))+  
  labs(color = "Divergence", size = "Count", x="-log(Q-value)", y = " ") + 
  theme(axis.text.y = element_text(size = 20), 
        panel.border = element_rect(colour = "black", size=2),
        panel.spacing = unit(0.7, "lines")) +
  scale_y_discrete(position="right") + 
  scale_size(guide=guide_legend(reverse=TRUE))  +
  facet_grid(Database + Div~.,scales="free", 
             space = "free") +
  scale_color_manual(values = c("#991f00", "#008080"), guide = "none")+
  theme(axis.title.x=element_text(size=20), 
        strip.text.y = element_blank(), 
        axis.text.x = element_text(size=20)) 
print(pp2)
dev.off()


Keggtest <- KEGG_20
Keggtest$logqval <- -log10(Keggtest$qvalue)
Keggtest$logpadjval <- -log10(Keggtest$p.adjust)




###Deal with the duplicate value so that the plot orders the categories correctly

kr1["hsa04657","Description"] <- "IL-17 signaling pathway "


pdf("Fig3C.pdf", width = 12, height = 13)
pp1<-ggplot(kr1, aes(x =-log10(p.adjust),
                y=reorder(Description, -log10(p.adjust), FUN=sum), 
                colour = factor(Div))) + 
  theme_bw()+ 
  geom_point(aes(size = Count ))+  
  labs(color = "Divergence", size = "Number of Genes", x="-log(Q-value)", y = " ") + 
  theme(axis.text.y = element_text(size = 20), 
        panel.border = element_rect(colour = "black", size=2),
        panel.spacing = unit(0.7, "lines")) +
  scale_y_discrete(position="right") + 
  scale_size(guide=guide_legend(reverse=TRUE))  +
  facet_grid(Database + Div~., scales="free", 
             space = "free") +
  scale_color_manual(values = c("#991f00", "#008080"), guide = "none")+
  theme(axis.title.x=element_text(size=20), 
        strip.text.y = element_blank(), 
        axis.text.x = element_text(size=20)) + scale_size_continuous(range = c(1, 10))
print(pp1)

dev.off()

