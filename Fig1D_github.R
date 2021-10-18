###To make Fig. 1D, plots for differences in abundances 

library('ggplot2')
library('vegan')
library("dplyr")
library("gplots")
library("DESeq2")
library("compositions")


##Create a directory
systempath <- "/Volumes/GoogleDrive/My\ Drive/Primates_Coculture/Primates/MyAnalysis_v02"
setwd(systempath)
dir.create(file.path(systempath, format(Sys.time(), "Fig1D_v01_%F")))
# 
rootpath <- paste0(file.path(systempath, format(Sys.time(), "Fig1D_v01_%F")))
#setwd(paste0(file.path(systempath, format(Sys.time(), "ShotgunMicrobeAbundance_%F"))))
setwd(rootpath)

#Read in metagenomic shotgun data (species level)

setwd("/Volumes/GoogleDrive/My\ Drive/Primates_Coculture/Sequencing20180530/ShotgunSequencing/Tables_Combined")
spec <- t(read.table("Merged_Species.txt", sep = "\t", header = TRUE, row.names = 1, check.names = FALSE))

#deal with sequencing data formatting
newrows <- vapply(strsplit(rownames(spec),"_"), `[`, 1, FUN.VALUE=character(1))
rownames(spec) <- newrows
#Rename the luca samples
rownames(spec)[16] <- "1UNCULTURED"
rownames(spec)[17] <- "65UNCULTURED"
rownames(spec)[18] <- "110UNCULTURED"
rownames(spec)[19] <- "111UNCULTURED"

#Filter the species table
cut <- colSums(spec>0) > 9
spec1<- spec[,cut]
dim(spec1)
# [1] 19 37


spec2 <- spec1
spec2 <- as.data.frame(spec2)
spec2$X <- rownames(spec2)


#load metadata
setwd("/Volumes/GoogleDrive/My\ Drive/Primates_Coculture/Metagenomics_Gopher/Run019")
met <- read.table("meta_hum.txt", sep = "\t", header = TRUE)
rownames(met) <- met$X
met1 <- met[!(rownames(met)) %in%c("UMGC-H20-BLANK", "UMGC-H2O-Blank", "UMGC.H2O.Blank", 
                                   "16UNCULTURED", "16UNPREPARED",
                                   "20UNCULTURED", "20UNPREPARED",
                                   "32UNCULTURED", "32UNPREPARED",
                                   "38UNCULTURED", "38UNPREPARED",
                                   "42UNCULTURED", "42UNPREPARED",
                                   "44UNCULTURED", "44UNPREPARED",
                                   "46UNCULTURED", "46UNPREPARED",
                                   "600UNCULTURED", "600UNPREPARED",
                                   "6UNCULTURED", "6UNPREPARED",
                                   "6COL", "6ONLY"),]
#Filter metadata 
met2 <- met1[grepl("UNCULTURED", rownames(met1)),]

#combine metadata and microbe abundance table
spec3 <- merge( met2, spec2, by="X")

colnames(spec3)
test <- colnames(spec3[,44:79])
#131

#####Make taxa plots
library(reshape2)
library(plyr)

tt1 <- as.data.frame(spec1)
tt1 <- tt1/100
cut2 <- colSums(tt1)>0.01
tt1 <- tt1[,cut2]
dim(tt1)
# [1] 19 36
tt1$SampleID <- rownames(tt1)
tt1 <- melt(tt1, id.vars = "SampleID",
            variable.name = "Taxa",
            value.name = "RelativeAbundance")
met3 <- met2
met3$SampleID <-rownames(met2)
tt2 <- merge(tt1, met3, by="SampleID")
library(ggplot2)

library(extrafont)
font_import()
loadfonts()


dataDir <- paste0(file.path(systempath, format(Sys.time(), "Fig1D_v01_%F")))



######Plot non transformed abundances
spec4 <- as.data.frame(spec1)
spec4$X <- rownames(spec4)
spec5 <- merge( met2, spec4, by="X")

# ##Create a full panel for 4 select microbial species:
# ###Bacteroides ovatus
# ###Prevotella copri
# ###Faecalbacterium prausnitzii
# ###Phascolarctobacterium succinatutens
# 

spec5$Species <- factor(spec5$Species, levels = c("Human", "Chimp", "Gorilla", "Orangutan"))

###b ovatus
bovatus <- ggplot(spec5, aes(x=Species, y=spec5$s__Bacteroides_ovatus, fill = Species)) + 
  geom_boxplot(outlier.shape=NA, alpha=0.25,lwd=1)+ 
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1) + 
  scale_fill_manual(values=c( "#00BFC4","#F8766D", "#7CAE00", "#C77CFF"))+
  theme_bw() +
  theme(axis.title.x = element_blank(), 
        axis.title.y= element_text(size=18),
        panel.border = element_rect(colour = "black", size=2),
        text=element_text(family="Arial", size=18),
        axis.text.y=element_text(size=22, color= "black"),
        axis.text.x = element_blank(),
        axis.title=element_text(size=26,face="bold")) + 
  labs(y = paste0("Relative Abundance (%)"), title = expression(italic("Bacteroides ovatus"))) + 
  guides(fill=FALSE) 

print(bovatus)

###Prevotella
pcopri <- ggplot(spec5, aes(x=Species, y=spec5$s__Prevotella_copri, fill = Species)) + 
  geom_boxplot(outlier.shape=NA, alpha=0.25,lwd=1)+ 
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1) + 
  scale_fill_manual(values=c( "#00BFC4","#F8766D", "#7CAE00", "#C77CFF"))+
  theme_bw() +
  theme(axis.title.x = element_blank(), 
        axis.title.y= element_blank(),
        panel.border = element_rect(colour = "black", size=2),
        text=element_text(family="Arial", size=18),
        axis.text.y=element_text(size=22, color= "black"),
        axis.text.x = element_blank(),
        axis.title=element_text(size=26,face="bold")) + 
  labs( title = expression(italic("Prevotella copri"))) + 
  guides(fill=FALSE) 

print(pcopri)

###Faecal blahblah prausnitzi
fprau <- ggplot(spec5, aes(x=Species, y=spec5$s__Faecalibacterium_prausnitzii, fill = Species)) + 
  geom_boxplot(outlier.shape=NA, alpha=0.25,lwd=1)+ 
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1) + 
  scale_fill_manual(values=c( "#00BFC4","#F8766D", "#7CAE00", "#C77CFF"))+
  theme_bw() +
  theme(axis.title.x = element_blank(), 
        axis.title.y= element_text(size=18),
        panel.border = element_rect(colour = "black", size=2),
        text=element_text(family="Arial", size = 18),
        axis.text.y=element_text(size=22, color= "black"),
        axis.text.x = element_blank(),
        axis.title=element_text(size=26,face="bold"))  + 
  labs(y = paste0("Relative Abundance (%)"), title = expression(italic("Faecalibacterium prausnitzii"))) + 
  guides(fill=FALSE) 


print(fprau)

###P succinatuens 

psuc <-  ggplot(spec5, aes(x=Species, y=spec5$s__Phascolarctobacterium_succinatutens, fill = Species)) + 
  geom_boxplot(outlier.shape=NA, alpha=0.25,lwd=1)+ 
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1) + 
  scale_fill_manual(values=c( "#00BFC4","#F8766D", "#7CAE00", "#C77CFF"))+
  theme_bw() +
  theme(axis.title.x = element_blank(), 
        axis.title.y= element_blank(),
        panel.border = element_rect(colour = "black", size=2),
        text=element_text(family="Arial", size = 18),
        axis.text.y=element_text(size=22, color= "black"),
        axis.text.x = element_blank(),
        axis.title=element_text(size=26,face="bold")) + 
  labs( title = expression(italic(atop("Phascolarctobacterium","succinatutens")))) + 
  guides(fill=FALSE) 
print(psuc)

library("gridExtra")


library(egg)
library(grid)



library("cowplot")


setwd(rootpath)
pdf("Fig1D_labels.pdf")
plot_grid(bovatus, psuc, fprau, pcopri, ncol = 2, align= "hv",rel_widths = c(1, 1), rel_heights = c(1,1))
dev.off()

################################################
###Make plots without labels for the final figure (fill in the labels in Keynote)


bovatus1 <- ggplot(spec5, aes(x=Species, y=spec5$s__Bacteroides_ovatus, fill = Species)) + 
  geom_boxplot(outlier.shape=NA, alpha=0.25,lwd=1)+ 
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1) + 
  scale_fill_manual(values=c( "#00BFC4","#F8766D", "#7CAE00", "#C77CFF"))+
  theme_bw() +
  theme(axis.title.x = element_blank(), 
        axis.title.y= element_blank(),
        panel.border = element_rect(colour = "black", size=2),
        text=element_text(family="Arial", size=18),
        axis.text.y=element_text(size=22, color= "black"),
        axis.text.x = element_blank(),
        axis.title=element_text(size=26,face="bold")) + 
  guides(fill=FALSE) 

print(bovatus1)

###Prevotella
pcopri1 <- ggplot(spec5, aes(x=Species, y=spec5$s__Prevotella_copri, fill = Species)) + 
  geom_boxplot(outlier.shape=NA, alpha=0.25,lwd=1)+ 
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1) + 
  scale_fill_manual(values=c( "#00BFC4","#F8766D", "#7CAE00", "#C77CFF"))+
  theme_bw() +
  theme(axis.title.x = element_blank(), 
        axis.title.y= element_blank(),
        panel.border = element_rect(colour = "black", size=2),
        text=element_text(family="Arial", size=18),
        axis.text.y=element_text(size=22, color= "black"),
        axis.text.x = element_blank(),
        axis.title=element_text(size=26,face="bold")) + 
  guides(fill=FALSE) 

print(pcopri1)

###Faecal blahblah prausnitzi
fprau1 <- ggplot(spec5, aes(x=Species, y=spec5$s__Faecalibacterium_prausnitzii, fill = Species)) + 
  geom_boxplot(outlier.shape=NA, alpha=0.25,lwd=1)+ 
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1) + 
  scale_fill_manual(values=c( "#00BFC4","#F8766D", "#7CAE00", "#C77CFF"))+
  theme_bw() +
  theme(axis.title.x = element_blank(), 
        axis.title.y= element_blank(),
        panel.border = element_rect(colour = "black", size=2),
        text=element_text(family="Arial", size=18),
        axis.text.y=element_text(size=22, color= "black"),
        axis.text.x = element_blank(),
        axis.title=element_text(size=26,face="bold"))  + 
  guides(fill=FALSE) 


print(fprau1)

###P succinatuens 

psuc1 <-  ggplot(spec5, aes(x=Species, y=spec5$s__Phascolarctobacterium_succinatutens, fill = Species)) + 
  geom_boxplot(outlier.shape=NA, alpha=0.25,lwd=1)+ 
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1) + 
  scale_fill_manual(values=c( "#00BFC4","#F8766D", "#7CAE00", "#C77CFF"))+
  theme_bw() +
  theme(axis.title.x = element_blank(), 
        axis.title.y= element_blank(),
        panel.border = element_rect(colour = "black", size=2),
        text=element_text(family="Arial", size=18),
        axis.text.y=element_text(size=22, color= "black"),
        axis.text.x = element_blank(),
        axis.title=element_text(size=26,face="bold")) + 
  guides(fill=FALSE) 
print(psuc1)


pdf("Fig1D_nolabels.pdf")
plot_grid(bovatus1, psuc1, fprau1, pcopri1, ncol = 2, align= "hv",rel_widths = c(1.3, 1.3), rel_heights = c(1.3,1.3), scale = 0.85)
dev.off()

pdf("bovatus.pdf")
print(bovatus)
dev.off()

