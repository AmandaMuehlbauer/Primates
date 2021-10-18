##This script makes the individual gene plots for Fig. 2C. 

library(extrafont)
font_import()
loadfonts()


library("gridExtra")

library("ggplot2")
library(egg)
library(grid)



library("cowplot")

setwd("/Volumes/GoogleDrive/My Drive/Primates_Coculture/Primates/MyAnalysis_v02/Results_2019-01-30/SavedModels")
res_hc <- readRDS("res_hc.RDS")
res_cc <- readRDS("res_cc.RDS")
res_gc <- readRDS("res_gc.RDS")
res_oc <- readRDS("res_oc.RDS")

res_hc <- as.data.frame(res_hc)
res_cc <- as.data.frame(res_cc)
res_gc <- as.data.frame(res_gc)
res_oc <- as.data.frame(res_oc)

systempath <- "/Volumes/GoogleDrive/My\ Drive/Primates_Coculture/Primates/MyAnalysis_v02"
setwd(systempath)
dir.create(file.path(systempath, format(Sys.time(), "Fig2C_%F")))

rootpath <- paste0(file.path(systempath, format(Sys.time(), "Fig2C_%F")))
setwd(rootpath)
a <- "ENSG00000197641" #set which gene to plot
ac1 <- data.frame(c("Human", "Chimp", "Gorilla", "Orangutan"), 
                  c(res_hc[a, "log2FoldChange"], res_cc[a, "log2FoldChange"], res_gc[a, "log2FoldChange"], res_oc[a, "log2FoldChange"]),
                  c(res_hc[a, "lfcSE"], res_cc[a, "lfcSE"], res_gc[a, "lfcSE"], res_oc[a, "lfcSE"]))

colnames(ac1) <- c("Species", "LFC", "SE")
ac1$Species <- factor(ac1$Species, levels = c("Human", "Chimp", "Gorilla", "Orangutan"))

pdf("Upreg.pdf")
p1 <-ggplot(ac1, aes(x=Species, y=LFC, group = Species, color=Species)) +
  geom_pointrange(aes(ymin=LFC-SE, ymax=LFC+SE), size=2)+labs(title=res_hc[a,"symbol"]) +
  geom_errorbar(aes(ymin=LFC-SE, ymax=LFC+SE), width=.2, size=2) +
  geom_point(shape = 1, size = 3) +
  theme_bw() +
   theme(text=element_text(family="Arial"), 
        title = element_text(size=27, face= "italic"),
        axis.text=element_text(size=36, color= "black"),
        axis.title=element_blank(), 
        axis.text.x=element_blank(),
        axis.text.y = element_text(size=30, color = "black"),
        panel.border = element_rect(size=2, color="black")) +
  theme(legend.position="none") + ylim(0,0.75) + scale_color_manual(values=c( "#00BFC4","#F8766D", "#7CAE00", "#C77CFF")) 
print(p1)
dev.off()
p

#theme(legend.position="none", axis.title.x = element_text( size=10), axis.text.x  = element_text( size=10), axis.title.y = element_text(size=10), axis.text.y = element_text(size=10), plot.title = element_text(size=10)) +
 # ylim(0,1.2) +

a <- "ENSG00000165264" #set which gene to plot
ac1 <- data.frame(c("Human", "Chimp", "Gorilla", "Orangutan"), 
                  c(res_hc[a, "log2FoldChange"], res_cc[a, "log2FoldChange"], res_gc[a, "log2FoldChange"], res_oc[a, "log2FoldChange"]),
                  c(res_hc[a, "lfcSE"], res_cc[a, "lfcSE"], res_gc[a, "lfcSE"], res_oc[a, "lfcSE"]))

colnames(ac1) <- c("Species", "LFC", "SE")
ac1$Species <- factor(ac1$Species, levels = c("Human", "Chimp", "Gorilla", "Orangutan"))

pdf("Downreg.pdf")
p2 <-ggplot(ac1, aes(x=Species, y=LFC, group = Species, color=Species)) +
  geom_pointrange(aes(ymin=LFC-SE, ymax=LFC+SE), size=2)+labs(title=res_hc[a,"symbol"]) +
  geom_errorbar(aes(ymin=LFC-SE, ymax=LFC+SE), width=.2, size=2) +
  geom_point(shape = 1, size = 3) +
  theme_bw() +
  theme(text=element_text(family="Arial"), title = element_text(size=27, face = "italic"),
        axis.text=element_text(size=6, color= "black"),
        axis.title=element_blank(), axis.text.x=element_blank(),axis.text.y = element_text(size=30, color = "black"),
        panel.border = element_rect(size=2, color="black")) +
  theme(legend.position="none") + ylim(-0.6,0) + scale_color_manual(values=c( "#00BFC4","#F8766D", "#7CAE00", "#C77CFF")) 
  
print(p2)
dev.off()

###No response 

no_resp <- subset(res_all, padj > 0.1)
a <- "ENSG00000119509" #set which gene to plot
ac1 <- data.frame(c("Human", "Chimp", "Gorilla", "Orangutan"), 
                  c(res_hc[a, "log2FoldChange"], res_cc[a, "log2FoldChange"], res_gc[a, "log2FoldChange"], res_oc[a, "log2FoldChange"]),
                  c(res_hc[a, "lfcSE"], res_cc[a, "lfcSE"], res_gc[a, "lfcSE"], res_oc[a, "lfcSE"]))

colnames(ac1) <- c("Species", "LFC", "SE")
ac1$Species <- factor(ac1$Species, levels = c("Human", "Chimp", "Gorilla", "Orangutan"))

pdf("No_resp.pdf")
p3 <-ggplot(ac1, aes(x=Species, y=LFC, group = Species, color=Species)) +
  geom_pointrange(aes(ymin=LFC-SE, ymax=LFC+SE), size=2)+labs(title=res_hc[a,"symbol"]) +
  geom_errorbar(aes(ymin=LFC-SE, ymax=LFC+SE), width=.2, size=2) +
  geom_point(shape = 1, size = 3) +
  theme_bw() +
  theme(text=element_text(family="Arial"), title = element_text(size=27, face = "italic"),
        axis.text=element_text(size=30, color= "black"),
        axis.title=element_blank(), axis.text.x=element_blank(),axis.text.y = element_text(size=30, color = "black"),
        panel.border = element_rect(size=2, color="black")) +
  theme(legend.position="none") + ylim(-0.3,0.3)  + scale_color_manual(values=c( "#00BFC4","#F8766D", "#7CAE00", "#C77CFF")) 

print(p3)
dev.off()


a <- "ENSG00000112309"
ac1 <- data.frame(c("Human", "Chimp", "Gorilla", "Orangutan"), 
                  c(res_hc[a, "log2FoldChange"], res_cc[a, "log2FoldChange"], res_gc[a, "log2FoldChange"], res_oc[a, "log2FoldChange"]),
                  c(res_hc[a, "lfcSE"], res_cc[a, "lfcSE"], res_gc[a, "lfcSE"], res_oc[a, "lfcSE"]))

colnames(ac1) <- c("Species", "LFC", "SE")
ac1$Species <- factor(ac1$Species, levels = c("Human", "Chimp", "Gorilla", "Orangutan"))

pdf("B3GAT2.pdf")
p4 <-ggplot(ac1, aes(x=Species, y=LFC, group = Species, color=Species)) +
  geom_pointrange(aes(ymin=LFC-SE, ymax=LFC+SE), size=2)+labs(title=res_hc[a,"symbol"]) +
  geom_errorbar(aes(ymin=LFC-SE, ymax=LFC+SE), width=.2, size=2) +
  geom_point(shape = 1, size = 3) +
  theme_bw() +
  theme(text=element_text(family="Arial"), title = element_text(size=27, face = "italic"),
        axis.text=element_text(size=30, color= "black"),
        axis.title=element_blank(), axis.text.x=element_blank(),axis.text.y = element_text(size=30, color = "black"),
        panel.border = element_rect(size=2, color="black")) +
  theme(legend.position="none") + scale_color_manual(values=c( "#00BFC4","#F8766D", "#7CAE00", "#C77CFF")) 
print(p4)
dev.off()

a<- "ENSG00000117593"
ac1 <- data.frame(c("Human", "Chimp", "Gorilla", "Orangutan"), 
                  c(res_hc[a, "log2FoldChange"], res_cc[a, "log2FoldChange"], res_gc[a, "log2FoldChange"], res_oc[a, "log2FoldChange"]),
                  c(res_hc[a, "lfcSE"], res_cc[a, "lfcSE"], res_gc[a, "lfcSE"], res_oc[a, "lfcSE"]))

colnames(ac1) <- c("Species", "LFC", "SE")
ac1$Species <- factor(ac1$Species, levels = c("Human", "Chimp", "Gorilla", "Orangutan"))

pdf("DARS2.pdf")
p5 <-ggplot(ac1, aes(x=Species, y=LFC, group = Species, color=Species)) +
  geom_pointrange(aes(ymin=LFC-SE, ymax=LFC+SE), size=2)+labs(title=res_hc[a,"symbol"]) +
  geom_errorbar(aes(ymin=LFC-SE, ymax=LFC+SE), width=.2, size=2) +
  geom_point(shape = 1, size = 3) +
  theme_bw() +
  theme(text=element_text(family="Arial"), title = element_text(size=27, face = "italic"),
        axis.text=element_text(size=30, color= "black"),
        axis.title=element_blank(), axis.text.x=element_blank(),axis.text.y = element_text(size=30, color = "black"),
        panel.border = element_rect(size=2, color="black")) +
  theme(legend.position="none") + scale_color_manual(values=c( "#00BFC4","#F8766D", "#7CAE00", "#C77CFF")) 
print(p5)
dev.off()


a<- "ENSG00000144048"
ac1 <- data.frame(c("Human", "Chimp", "Gorilla", "Orangutan"), 
                  c(res_hc[a, "log2FoldChange"], res_cc[a, "log2FoldChange"], res_gc[a, "log2FoldChange"], res_oc[a, "log2FoldChange"]),
                  c(res_hc[a, "lfcSE"], res_cc[a, "lfcSE"], res_gc[a, "lfcSE"], res_oc[a, "lfcSE"]))

colnames(ac1) <- c("Species", "LFC", "SE")
ac1$Species <- factor(ac1$Species, levels = c("Human", "Chimp", "Gorilla", "Orangutan"))

pdf("DUSP11.pdf")
p6 <-ggplot(ac1, aes(x=Species, y=LFC, group = Species, color=Species)) +
  geom_pointrange(aes(ymin=LFC-SE, ymax=LFC+SE), size=2)+labs(title=res_hc[a,"symbol"]) +
  geom_errorbar(aes(ymin=LFC-SE, ymax=LFC+SE), width=.2, size=2) +
  geom_point(shape = 1, size = 3) +
  theme_bw() +
  theme(text=element_text(family="Arial"), title = element_text(size=27, face = "italic"),
        axis.text=element_text(size=30, color= "black"),
        axis.title=element_blank(), axis.text.x=element_blank(),axis.text.y = element_text(size=30, color = "black"),
        panel.border = element_rect(size=2, color="black")) +
  theme(legend.position="none") + scale_color_manual(values=c( "#00BFC4","#F8766D", "#7CAE00", "#C77CFF")) 
print(p6)
dev.off()

a<- "ENSG00000138771" 
ac1 <- data.frame(c("Human", "Chimp", "Gorilla", "Orangutan"), 
                  c(res_hc[a, "log2FoldChange"], res_cc[a, "log2FoldChange"], res_gc[a, "log2FoldChange"], res_oc[a, "log2FoldChange"]),
                  c(res_hc[a, "lfcSE"], res_cc[a, "lfcSE"], res_gc[a, "lfcSE"], res_oc[a, "lfcSE"]))

colnames(ac1) <- c("Species", "LFC", "SE")
ac1$Species <- factor(ac1$Species, levels = c("Human", "Chimp", "Gorilla", "Orangutan"))

pdf("SHROOM3.pdf")
p7 <-ggplot(ac1, aes(x=Species, y=LFC, group = Species, color=Species)) +
  geom_pointrange(aes(ymin=LFC-SE, ymax=LFC+SE), size=2)+labs(title=res_hc[a,"symbol"]) +
  geom_errorbar(aes(ymin=LFC-SE, ymax=LFC+SE), width=.2, size=2) +
  geom_point(shape = 1, size = 3) +
  theme_bw() +
  theme(text=element_text(family="Arial"), title = element_text(size=27, face = "italic"),
        axis.text=element_text(size=30, color= "black"),
        axis.title=element_blank(), axis.text.x=element_blank(),axis.text.y = element_text(size=30, color = "black"),
        panel.border = element_rect(size=2, color="black")) +
  theme(legend.position="none") + scale_color_manual(values=c( "#00BFC4","#F8766D", "#7CAE00", "#C77CFF")) 

print(p7)
dev.off()



###Responds to two of the species (Human and Chimp)
a<- "ENSG00000159228" 
ac1 <- data.frame(c("Human", "Chimp", "Gorilla", "Orangutan"), 
                  c(res_hc[a, "log2FoldChange"], res_cc[a, "log2FoldChange"], res_gc[a, "log2FoldChange"], res_oc[a, "log2FoldChange"]),
                  c(res_hc[a, "lfcSE"], res_cc[a, "lfcSE"], res_gc[a, "lfcSE"], res_oc[a, "lfcSE"]))

colnames(ac1) <- c("Species", "LFC", "SE")
ac1$Species <- factor(ac1$Species, levels = c("Human", "Chimp", "Gorilla", "Orangutan"))


##differenitally Expressed in human and chimp
pdf("CBR1.pdf")
p8 <-ggplot(ac1, aes(x=Species, y=LFC, group = Species, color=Species)) +
  geom_pointrange(aes(ymin=LFC-SE, ymax=LFC+SE), size=2)+labs(title=res_hc[a,"symbol"]) +
  geom_errorbar(aes(ymin=LFC-SE, ymax=LFC+SE), width=.2, size=2) +
  geom_point(shape = 1, size = 3) +
  theme_bw() +
  theme(text=element_text(family="Arial"), title = element_text(size=27, face = "italic"),
        axis.text=element_text(size=30, color= "black"),
        axis.title=element_blank(), axis.text.x=element_blank(),axis.text.y = element_text(size=30, color = "black"),
        panel.border = element_rect(size=2, color="black")) +
  theme(legend.position="none") + scale_color_manual(values=c( "#00BFC4","#F8766D", "#7CAE00", "#C77CFF")) 

print(p8)
dev.off()

pdf("Fig2C_1.pdf")
plot_grid(p1, p2,p3,p4, ncol = 2, align= "hv",rel_widths = c(1.3, 1.3), rel_heights = c(1.3,1.3), scale = 1)
dev.off()

pdf("Fig2C_2.pdf")
plot_grid(p5,p6,p7,p8, ncol = 2, align= "hv",rel_widths = c(1.3, 1.3), rel_heights = c(1.3,1.3), scale = 1)
dev.off()

pdf("Fig2C_all.pdf")
plot_grid(p1, p2,p3,p4,p5,p6,p7,p8, ncol = 4, align= "hv",rel_widths = c(1.3, 1.3), rel_heights = c(1.3,1.3), scale = 0.5)
dev.off()

plotall <- list(p1,p2, p3, p4, p5, p6, p7, p8)
pdf("Fig2C_two.pdf")
marrangeGrob(plotall, nrow = 2, ncol=2)
dev.off()


library(grid)
grid.newpage()
pdf("Fig2C_four.pdf")
grid.layout(nrow=2, ncol=2, widths=c(1,2), heights=c(2,1), respect=TRUE)
grid.draw(
  cbind(
    rbind(ggplotGrob(p3), ggplotGrob(p4), size = "first"),
    rbind(ggplotGrob(p1), ggplotGrob(p2), size = "first"), 
    size = "first"))
dev.off()

pdf("Fig2C_eight.pdf")
grid.layout(nrow=2, ncol=2, widths=c(1,2), heights=c(2,1), respect=TRUE)
grid.draw(
  cbind(
    rbind(ggplotGrob(p7), ggplotGrob(p8), size = "first"),
    rbind(ggplotGrob(p5), ggplotGrob(p6), size = "first"), 
    size = "first"))
dev.off()


#The order I want in the final plot is: 
#Top: 
##Shroom3, B3Gat2, DUSP11, DARS2
#Bottom: 
#INVS(no response), NDUFB6 (down), SERPBIN13 (up), CBR1

###That ends up being: Page1: P7, P4, P3, P2
##Page2: P6, P5, P1, P8

pdf("Fig2C_right.pdf")
grid.layout(nrow=2, ncol=2, widths=c(1,2), heights=c(2,1), respect=TRUE)
grid.draw(
  rbind(
    cbind(ggplotGrob(p7), ggplotGrob(p4), size = "first"),
    cbind(ggplotGrob(p3), ggplotGrob(p2), size = "first"), 
    size = "first"))
dev.off()

pdf("Fig2C_left.pdf")
grid.layout(nrow=2, ncol=2, widths=c(1,2), heights=c(2,1), respect=TRUE)
grid.draw(
  rbind(
    cbind(ggplotGrob(p6), ggplotGrob(p5), size = "first"),
    cbind(ggplotGrob(p1), ggplotGrob(p8), size = "first"), 
    size = "first"))
dev.off()

pdf("Fig2C_all_label.pdf",paper='A4r')
plot_grid(p7, p4, p3, p2, p6, p5, p1, p8, ncol = 4, align= "hv",rel_widths = c(1, 1), rel_heights = c(1,1), scale=1)
dev.off()


######################################################

setwd(rootpath)
a <- "ENSG00000197641" #set which gene to plot
ac1 <- data.frame(c("Human", "Chimp", "Gorilla", "Orangutan"), 
                  c(res_hc[a, "log2FoldChange"], res_cc[a, "log2FoldChange"], res_gc[a, "log2FoldChange"], res_oc[a, "log2FoldChange"]),
                  c(res_hc[a, "lfcSE"], res_cc[a, "lfcSE"], res_gc[a, "lfcSE"], res_oc[a, "lfcSE"]))

colnames(ac1) <- c("Species", "LFC", "SE")
ac1$Species <- factor(ac1$Species, levels = c("Human", "Chimp", "Gorilla", "Orangutan"))

pdf("Upreg_nolabel.pdf")
p1 <-ggplot(ac1, aes(x=Species, y=LFC, group = Species, color=Species)) +
  geom_pointrange(aes(ymin=LFC-SE, ymax=LFC+SE), size=2)+labs(title=res_hc[a,"symbol"]) +
  geom_errorbar(aes(ymin=LFC-SE, ymax=LFC+SE), width=.2, size=2) +
  geom_point(shape = 1, size = 3) +
  theme_bw() +
  theme(text=element_text(family="Arial"), 
        title = element_blank(),
        axis.text=element_blank(),
        axis.title=element_blank(), 
        axis.text.x=element_blank(),
        axis.text.y = element_blank(),
        panel.border = element_rect(size=2, color="black")) +
  theme(legend.position="none") + ylim(0,0.75) + scale_color_manual(values=c( "#00BFC4","#F8766D", "#7CAE00", "#C77CFF")) 
print(p1)
dev.off()
p

#theme(legend.position="none", axis.title.x = element_text( size=10), axis.text.x  = element_text( size=10), axis.title.y = element_text(size=10), axis.text.y = element_text(size=10), plot.title = element_text(size=10)) +
# ylim(0,1.2) +

a <- "ENSG00000165264" #set which gene to plot
ac1 <- data.frame(c("Human", "Chimp", "Gorilla", "Orangutan"), 
                  c(res_hc[a, "log2FoldChange"], res_cc[a, "log2FoldChange"], res_gc[a, "log2FoldChange"], res_oc[a, "log2FoldChange"]),
                  c(res_hc[a, "lfcSE"], res_cc[a, "lfcSE"], res_gc[a, "lfcSE"], res_oc[a, "lfcSE"]))

colnames(ac1) <- c("Species", "LFC", "SE")
ac1$Species <- factor(ac1$Species, levels = c("Human", "Chimp", "Gorilla", "Orangutan"))

pdf("Downreg_nolabel.pdf")
p2 <-ggplot(ac1, aes(x=Species, y=LFC, group = Species, color=Species)) +
  geom_pointrange(aes(ymin=LFC-SE, ymax=LFC+SE), size=2)+labs(title=res_hc[a,"symbol"]) +
  geom_errorbar(aes(ymin=LFC-SE, ymax=LFC+SE), width=.2, size=2) +
  geom_point(shape = 1, size = 3) +
  theme_bw() +
  theme(text=element_text(family="Arial"), 
        title = element_blank(),
        axis.text=element_blank(),
        axis.title=element_blank(), 
        axis.text.x=element_blank(),
        axis.text.y = element_blank(),
        panel.border = element_rect(size=2, color="black")) +
  theme(legend.position="none") + ylim(-0.6,0) + scale_color_manual(values=c( "#00BFC4","#F8766D", "#7CAE00", "#C77CFF")) 

print(p2)
dev.off()

###No response 

no_resp <- subset(res_all, padj > 0.1)
a <- "ENSG00000119509" #set which gene to plot
ac1 <- data.frame(c("Human", "Chimp", "Gorilla", "Orangutan"), 
                  c(res_hc[a, "log2FoldChange"], res_cc[a, "log2FoldChange"], res_gc[a, "log2FoldChange"], res_oc[a, "log2FoldChange"]),
                  c(res_hc[a, "lfcSE"], res_cc[a, "lfcSE"], res_gc[a, "lfcSE"], res_oc[a, "lfcSE"]))

colnames(ac1) <- c("Species", "LFC", "SE")
ac1$Species <- factor(ac1$Species, levels = c("Human", "Chimp", "Gorilla", "Orangutan"))

pdf("No_resp_nolabel.pdf")
p3 <-ggplot(ac1, aes(x=Species, y=LFC, group = Species, color=Species)) +
  geom_pointrange(aes(ymin=LFC-SE, ymax=LFC+SE), size=2)+labs(title=res_hc[a,"symbol"]) +
  geom_errorbar(aes(ymin=LFC-SE, ymax=LFC+SE), width=.2, size=2) +
  geom_point(shape = 1, size = 3) +
  theme_bw() +
  theme(text=element_text(family="Arial"), 
        title = element_blank(),
        axis.text=element_blank(),
        axis.title=element_blank(), 
        axis.text.x=element_blank(),
        axis.text.y = element_blank(),
        panel.border = element_rect(size=2, color="black")) +
  theme(legend.position="none") + ylim(-0.3,0.3)  + scale_color_manual(values=c( "#00BFC4","#F8766D", "#7CAE00", "#C77CFF")) 

print(p3)
dev.off()


a <- "ENSG00000112309"
ac1 <- data.frame(c("Human", "Chimp", "Gorilla", "Orangutan"), 
                  c(res_hc[a, "log2FoldChange"], res_cc[a, "log2FoldChange"], res_gc[a, "log2FoldChange"], res_oc[a, "log2FoldChange"]),
                  c(res_hc[a, "lfcSE"], res_cc[a, "lfcSE"], res_gc[a, "lfcSE"], res_oc[a, "lfcSE"]))

colnames(ac1) <- c("Species", "LFC", "SE")
ac1$Species <- factor(ac1$Species, levels = c("Human", "Chimp", "Gorilla", "Orangutan"))

pdf("B3GAT2_nolabel.pdf")
p4 <-ggplot(ac1, aes(x=Species, y=LFC, group = Species, color=Species)) +
  geom_pointrange(aes(ymin=LFC-SE, ymax=LFC+SE), size=2)+labs(title=res_hc[a,"symbol"]) +
  geom_errorbar(aes(ymin=LFC-SE, ymax=LFC+SE), width=.2, size=2) +
  geom_point(shape = 1, size = 3) +
  theme_bw() +
  theme(text=element_text(family="Arial"), 
        title = element_blank(),
        axis.text=element_blank(),
        axis.title=element_blank(), 
        axis.text.x=element_blank(),
        axis.text.y = element_blank(),
        panel.border = element_rect(size=2, color="black")) +
  theme(legend.position="none") + scale_color_manual(values=c( "#00BFC4","#F8766D", "#7CAE00", "#C77CFF")) 
print(p4)
dev.off()

a<- "ENSG00000117593"
ac1 <- data.frame(c("Human", "Chimp", "Gorilla", "Orangutan"), 
                  c(res_hc[a, "log2FoldChange"], res_cc[a, "log2FoldChange"], res_gc[a, "log2FoldChange"], res_oc[a, "log2FoldChange"]),
                  c(res_hc[a, "lfcSE"], res_cc[a, "lfcSE"], res_gc[a, "lfcSE"], res_oc[a, "lfcSE"]))

colnames(ac1) <- c("Species", "LFC", "SE")
ac1$Species <- factor(ac1$Species, levels = c("Human", "Chimp", "Gorilla", "Orangutan"))

pdf("DARS2_nolabel.pdf")
p5 <-ggplot(ac1, aes(x=Species, y=LFC, group = Species, color=Species)) +
  geom_pointrange(aes(ymin=LFC-SE, ymax=LFC+SE), size=2)+labs(title=res_hc[a,"symbol"]) +
  geom_errorbar(aes(ymin=LFC-SE, ymax=LFC+SE), width=.2, size=2) +
  geom_point(shape = 1, size = 3) +
  theme_bw() +
  theme(text=element_text(family="Arial"), 
        title = element_blank(),
        axis.text=element_blank(),
        axis.title=element_blank(), 
        axis.text.x=element_blank(),
        axis.text.y = element_blank(),
        panel.border = element_rect(size=2, color="black")) +
  theme(legend.position="none") + scale_color_manual(values=c( "#00BFC4","#F8766D", "#7CAE00", "#C77CFF")) 
print(p5)
dev.off()


a<- "ENSG00000144048"
ac1 <- data.frame(c("Human", "Chimp", "Gorilla", "Orangutan"), 
                  c(res_hc[a, "log2FoldChange"], res_cc[a, "log2FoldChange"], res_gc[a, "log2FoldChange"], res_oc[a, "log2FoldChange"]),
                  c(res_hc[a, "lfcSE"], res_cc[a, "lfcSE"], res_gc[a, "lfcSE"], res_oc[a, "lfcSE"]))

colnames(ac1) <- c("Species", "LFC", "SE")
ac1$Species <- factor(ac1$Species, levels = c("Human", "Chimp", "Gorilla", "Orangutan"))

pdf("DUSP11_nolabel.pdf")
p6 <-ggplot(ac1, aes(x=Species, y=LFC, group = Species, color=Species)) +
  geom_pointrange(aes(ymin=LFC-SE, ymax=LFC+SE), size=2)+labs(title=res_hc[a,"symbol"]) +
  geom_errorbar(aes(ymin=LFC-SE, ymax=LFC+SE), width=.2, size=2) +
  geom_point(shape = 1, size = 3) +
  theme_bw() +
  theme(text=element_text(family="Arial"), 
        title = element_blank(),
        axis.text=element_blank(),
        axis.title=element_blank(), 
        axis.text.x=element_blank(),
        axis.text.y = element_blank(),
        panel.border = element_rect(size=2, color="black")) +
  theme(legend.position="none") + scale_color_manual(values=c( "#00BFC4","#F8766D", "#7CAE00", "#C77CFF")) 
print(p6)
dev.off()

a<- "ENSG00000138771" 
ac1 <- data.frame(c("Human", "Chimp", "Gorilla", "Orangutan"), 
                  c(res_hc[a, "log2FoldChange"], res_cc[a, "log2FoldChange"], res_gc[a, "log2FoldChange"], res_oc[a, "log2FoldChange"]),
                  c(res_hc[a, "lfcSE"], res_cc[a, "lfcSE"], res_gc[a, "lfcSE"], res_oc[a, "lfcSE"]))

colnames(ac1) <- c("Species", "LFC", "SE")
ac1$Species <- factor(ac1$Species, levels = c("Human", "Chimp", "Gorilla", "Orangutan"))

pdf("SHROOM3_nolabel.pdf")
p7 <-ggplot(ac1, aes(x=Species, y=LFC, group = Species, color=Species)) +
  geom_pointrange(aes(ymin=LFC-SE, ymax=LFC+SE), size=2)+labs(title=res_hc[a,"symbol"]) +
  geom_errorbar(aes(ymin=LFC-SE, ymax=LFC+SE), width=.2, size=2) +
  geom_point(shape = 1, size = 3) +
  theme_bw() +
  theme(text=element_text(family="Arial"), 
        title = element_blank(),
        axis.text=element_blank(),
        axis.title=element_blank(), 
        axis.text.x=element_blank(),
        axis.text.y = element_blank(),
        panel.border = element_rect(size=2, color="black")) +
  theme(legend.position="none") + scale_color_manual(values=c( "#00BFC4","#F8766D", "#7CAE00", "#C77CFF")) 

print(p7)
dev.off()



###Responds to two of the species (Human and Chimp)
a<- "ENSG00000159228" 
ac1 <- data.frame(c("Human", "Chimp", "Gorilla", "Orangutan"), 
                  c(res_hc[a, "log2FoldChange"], res_cc[a, "log2FoldChange"], res_gc[a, "log2FoldChange"], res_oc[a, "log2FoldChange"]),
                  c(res_hc[a, "lfcSE"], res_cc[a, "lfcSE"], res_gc[a, "lfcSE"], res_oc[a, "lfcSE"]))

colnames(ac1) <- c("Species", "LFC", "SE")
ac1$Species <- factor(ac1$Species, levels = c("Human", "Chimp", "Gorilla", "Orangutan"))


##differenitally Expressed in human and chimp
pdf("CBR1_nolabel.pdf")
p8 <-ggplot(ac1, aes(x=Species, y=LFC, group = Species, color=Species)) +
  geom_pointrange(aes(ymin=LFC-SE, ymax=LFC+SE), size=2)+labs(title=res_hc[a,"symbol"]) +
  geom_errorbar(aes(ymin=LFC-SE, ymax=LFC+SE), width=.2, size=2) +
  geom_point(shape = 1, size = 3) +
  theme_bw() +
  theme(text=element_text(family="Arial"), 
        title = element_blank(),
        axis.text=element_blank(),
        axis.title=element_blank(), 
        axis.text.x=element_blank(),
        axis.text.y = element_blank(),
        panel.border = element_rect(size=2, color="black")) +
  theme(legend.position="none") + scale_color_manual(values=c( "#00BFC4","#F8766D", "#7CAE00", "#C77CFF")) 

print(p8)
dev.off()




pdf("Fig2C_left_nolabel.pdf")
plot_grid(p7,p4,p3,p2, ncol = 2, align= "hv",rel_widths = c(1, 1), rel_heights = c(1,1), scale= 0.8)
dev.off()


pdf("Fig2C_right_nolabel.pdf")
plot_grid(p6, p5, p1, p8, ncol = 2, align= "hv",rel_widths = c(1, 1), rel_heights = c(1,1), scale= 0.8)
dev.off()

pdf("Fig2C_all_nolabels.pdf",paper='A4r')
plot_grid(p7,p4,p3,p2,p6, p5, p1, p8, ncol = 4, align= "hv",rel_widths = c(1, 1), rel_heights = c(1,1), scale=1)
dev.off()

