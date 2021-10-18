###This makes the heatmap for Fig2A
library("ggplot2")
library("RColorBrewer")
library("gplots")

##Create directory
systempath <- "/Volumes/GoogleDrive/My\ Drive/Primates_Coculture/Primates/MyAnalysis_v02"
setwd(systempath)
dir.create(file.path(systempath, format(Sys.time(), "Fig2A_%F")))

rootpath <- paste0(file.path(systempath, format(Sys.time(), "Fig2A_%F")))

#Read in data from DESeq2 analysis
setwd("/Volumes/GoogleDrive/My Drive/Primates_Coculture/Primates/MyAnalysis_v02/Results_2018-11-06/")
n_matrix5.1 <- readRDS("n_matrix5.RDS")
n_matrix4.1 <-readRDS("n_matrix4.RDS")
n_matrix5.2 <-as.data.frame(n_matrix5.1)
n_matrix4.2 <- as.data.frame(n_matrix4.1[c("Chimp", "Gorilla", "Human", "Orangutan", "col_n")])
n_matrix4.3 <- n_matrix4.2
n_matrix4.3$col_n[n_matrix4.3$col_n == "gray" ] <- "wheat3"
n_matrix5.2 <- as.data.frame(n_matrix5.2[c("Chimp", "Gorilla", "Human", "Orangutan")])
n_matrix5.2 <- as.matrix(n_matrix5.2)

my_palette <- colorRampPalette(c("#67a9cf","#67a9cf", "#f7f7f7","#ef8a62", "#ef8a62"))(n = 299)


specorder <- c("Human", "Chimp", "Gorilla", "Orangutan")
specorder1 <- factor(specorder, levels = c("Human", "Chimp", "Gorilla", "Orangutan"))



setwd(paste0(file.path(systempath, format(Sys.time(), "Fig2A_%F"))))


my_palette <- colorRampPalette(c("#999999", "#ffffff","#ef8a62"  ))(n = 9)

colors = c(seq(min(n_matrix5.2),-0.1,length=100), seq(-0.10000001,0,length=100), seq(0.00000001,max(n_matrix5.2),length=100))
pdf ("test.pdf")
heatmap.2(head(n_matrix5.2,20),
          scale="none",
          symm=F,symkey=F,symbreaks=T,
          trace="none",density.info="none",
          col=my_palette, cexRow=0.5, cexCol = 1,
          margins = c(5,1), labCol = c("Chimp", "Gorilla", "Human", "Orangutan"),
          dendrogram = "none", Rowv = F,Colv = c(specorder1),
          lhei=c(2,10), lwid=c(2,3.5), breaks = colors)
            
          


dev.off()


colors = unique(c(seq(min(n_matrix5.2),-0.1,length=10), seq(-0.1,0.1,length=10), seq(0,max(n_matrix5.2),length=10)))


my_palette <- colorRampPalette(c("#999999", "#ffffff","#ef8a62"  ))(n = 28)

pdf ("Heatmap5.pdf")
heatmap.2(n_matrix5.2,
          RowSideColors=as.character(n_matrix4.3[as.character(unlist(rownames(n_matrix4.3))),]$col_n), 
          scale="none",
          symm=F,symkey=F,symbreaks=T, 
          trace="none",density.info="none", 
          col=my_palette, cexRow=0.5, cexCol = 1, 
          margins = c(5,1), labCol = c("Chimp", "Gorilla", "Human", "Orangutan"), 
          dendrogram = "none", Rowv = F,Colv = c(specorder1),
          lhei=c(2,10), lwid=c(2,3.5), breaks = colors)


dev.off()


pdf ("test3.pdf")
heatmap.2(head(n_matrix5.2, 20),
          scale="none",
          symm=F,symkey=F,symbreaks=T, 
          trace="none",density.info="none", 
          col=my_palette, cexRow=0.5, cexCol = 1, 
          margins = c(5,1), labCol = c("Chimp", "Gorilla", "Human", "Orangutan"), 
          dendrogram = "none", Rowv = F,Colv = c(specorder1),
          lhei=c(2,10), lwid=c(2,3.5), breaks = colors)


dev.off()




# centering with 'scale()'
center_scale <- function(x) {
  scale(x, scale = FALSE)
}

# apply it
nmat1 <- center_scale(n_matrix5.2)
myBreaks <- seq(min(nmat1), max(nmat1), length.out=10)
my_palette <- colorRampPalette(c("#999999", "#ffffff","#ef8a62"  ))(n = 9)
pdf ("test2.pdf")
heatmap.2(nmat1,
          RowSideColors=as.character(n_matrix4.3[as.character(unlist(rownames(n_matrix4.3))),]$col_n), 
          scale="none",
          symm=F,symkey=F,symbreaks=T, 
          trace="none",density.info="none", 
          col=my_palette, cexRow=0.5, cexCol = 1, 
          margins = c(5,1), labCol = c("Chimp", "Gorilla", "Human", "Orangutan"), 
          dendrogram = "none", Rowv = F,Colv = c(specorder1),
          lhei=c(2,10), lwid=c(2,3.5), breaks = myBreaks)


dev.off()

my_palette <- colorRampPalette(c("dodgerblue2", "#f7f7f7","darkorange2"))(n = 299)


pdf ("Heatmap6.pdf")
heatmap.2(n_matrix5.2,
          RowSideColors=as.character(n_matrix4.3[as.character(unlist(rownames(n_matrix4.3))),]$col_n), 
          scale="none",
          symm=F,symkey=F,symbreaks=T, 
          trace="none",density.info="none", 
          col=my_palette, cexRow=0.5, cexCol = 1, 
          margins = c(5,1), labCol = c("Chimp", "Gorilla", "Human", "Orangutan"), 
          dendrogram = "none", Rowv = F,Colv = c(specorder1),
          lhei=c(2,10), lwid=c(2,3.5))

dev.off()


my_palette <- colorRampPalette(c("#ef8a62", "#ffffff","#999999"))(n = 299)

pdf ("Heatmap7.pdf")
heatmap.2(n_matrix5.2,
          RowSideColors=as.character(n_matrix4.3[as.character(unlist(rownames(n_matrix4.3))),]$col_n), 
          scale="none",
          symm=F,symkey=F,symbreaks=T, 
          trace="none",density.info="none", 
          col=my_palette, cexRow=0.5, cexCol = 1, 
          margins = c(5,1), labCol = c("Chimp", "Gorilla", "Human", "Orangutan"), 
          dendrogram = "none", Rowv = FALSE,  lhei=c(2,10), lwid=c(2,3.5))

dev.off()

