##This script makes the upset plot, Fig. 2B. This requires the analysis from the LRT tests in the Final
#Analysis script. 

#Create new directory
systempath <- "/Volumes/GoogleDrive/My\ Drive/Primates_Coculture/Primates/MyAnalysis_v02"
setwd(systempath)
dir.create(file.path(systempath, format(Sys.time(), "Fig2B_%F")))

rootpath <- paste0(file.path(systempath, format(Sys.time(), "Fig2B_%F")))
setwd(rootpath)

#Make the upset plot.
#This uses hard coded numbers to make the plot. 

library("UpSetR")
mem1 <-as.data.frame(mem1)

mem2 <- c(Chimp=24, Gorilla = 57, Human = 71, Orangutan=12, `Chimp&Gorilla&Human&Orangutan`=2261, `Chimp&Gorilla&Orangutan`=495, `Gorilla&Chimp`=113, `Gorilla&Orangutan`=184, `Human&Chimp`=57,`Human&Chimp&Gorilla`=309, `Human&Chimp&Orangutan`=144, `Human&Gorilla`=161, `Human&Gorilla&Orangutan`=365, `Human&Orangutan`=36, `Orangutan&Chimp`=40)

pdf("Upset_LRT.pdf")
upset(fromExpression(mem2), order.by = "freq",point.size=3,text.scale=c(1.5,1.5,1.5,1.5,1.5,1))

dev.off()
