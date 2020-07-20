#!/usr/bin/env Rscript

#load libraries
library(ggplot2)
library(gridExtra)


#load libraries
library(ggplot2)
library(gridExtra)

#check if data exist and load them
#downstream
if (file.exists("../ParasiTE_output/neighbor_TEs-length/RLC-neighbor-TE-downstream.bed"))
RLCi <- read.delim("../ParasiTE_output/neighbor_TEs-length/RLC-neighbor-TE-downstream.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/neighbor_TEs-length/RLG-neighbor-TE-downstream.bed"))
RLGi <- read.delim("../ParasiTE_output/neighbor_TEs-length/RLG-neighbor-TE-downstream.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/neighbor_TEs-length/Caulimoviridae-neighbor-TE-downstream.bed"))
Caulimoviridaei <- read.delim("../ParasiTE_output/neighbor_TEs-length/Caulimoviridae-neighbor-TE-downstream.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/neighbor_TEs-length/RIX-neighbor-TE-downstream.bed"))
RIXi <- read.delim("../ParasiTE_output/neighbor_TEs-length/RIX-neighbor-TE-downstream.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/neighbor_TEs-length/RLX-neighbor-TE-downstream.bed"))
RLXi <- read.delim("../ParasiTE_output/neighbor_TEs-length/RLX-neighbor-TE-downstream.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/neighbor_TEs-length/RYX-neighbor-TE-downstream.bed"))
RYXi <- read.delim("../ParasiTE_output/neighbor_TEs-length/RYX-neighbor-TE-downstream.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/neighbor_TEs-length/RXX-neighbor-TE-downstream.bed"))
RXXi <- read.delim("../ParasiTE_output/neighbor_TEs-length/RXX-neighbor-TE-downstream.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/neighbor_TEs-length/RSX-neighbor-TE-downstream.bed"))
RSXi <- read.delim("../ParasiTE_output/neighbor_TEs-length/RSX-neighbor-TE-downstream.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/neighbor_TEs-length/DTT-neighbor-TE-downstream.bed"))
DTTi <- read.delim("../ParasiTE_output/neighbor_TEs-length/DTT-neighbor-TE-downstream.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/neighbor_TEs-length/DTM-neighbor-TE-downstream.bed"))
DTMi <- read.delim("../ParasiTE_output/neighbor_TEs-length/DTM-neighbor-TE-downstream.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/neighbor_TEs-length/DTA-neighbor-TE-downstream.bed"))
DTAi <- read.delim("../ParasiTE_output/neighbor_TEs-length/DTA-neighbor-TE-downstream.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/neighbor_TEs-length/DTC-neighbor-TE-downstream.bed"))
DTCi <- read.delim("../ParasiTE_output/neighbor_TEs-length/DTC-neighbor-TE-downstream.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/neighbor_TEs-length/DTH-neighbor-TE-downstream.bed"))
DTHi <- read.delim("../ParasiTE_output/neighbor_TEs-length/DTH-neighbor-TE-downstream.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/neighbor_TEs-length/DTX-neighbor-TE-downstream.bed"))
DTXi <- read.delim("../ParasiTE_output/neighbor_TEs-length/DTX-neighbor-TE-downstream.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/neighbor_TEs-length/DHX-neighbor-TE-downstream.bed"))
DHXi <- read.delim("../ParasiTE_output/neighbor_TEs-length/DHX-neighbor-TE-downstream.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/neighbor_TEs-length/TXX-neighbor-TE-downstream.bed"))
TXXi <- read.delim("../ParasiTE_output/neighbor_TEs-length/TXX-neighbor-TE-downstream.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/neighbor_TEs-length/TOTAL-neighbor-TE-downstream.bed"))
TOTALi <- read.delim("../ParasiTE_output/neighbor_TEs-length/TOTAL-neighbor-TE-downstream.bed", header = FALSE, sep = "\t", dec = ".")

#upstream
if (file.exists("../ParasiTE_output/neighbor_TEs-length/RLC-neighbor-TE-upstream.bed"))
RLCe <- read.delim("../ParasiTE_output/neighbor_TEs-length/RLC-neighbor-TE-upstream.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/neighbor_TEs-length/RLG-neighbor-TE-upstream.bed"))
RLGe <- read.delim("../ParasiTE_output/neighbor_TEs-length/RLG-neighbor-TE-upstream.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/neighbor_TEs-length/Caulimoviridae-neighbor-TE-upstream.bed"))
Caulimoviridaee <- read.delim("../ParasiTE_output/neighbor_TEs-length/Caulimoviridae-neighbor-TE-upstream.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/neighbor_TEs-length/RYX-neighbor-TE-upstream.bed"))
RYXe <- read.delim("../ParasiTE_output/neighbor_TEs-length/RYX-neighbor-TE-upstream.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/neighbor_TEs-length/RIX-neighbor-TE-upstream.bed"))
RIXe <- read.delim("../ParasiTE_output/neighbor_TEs-length/RIX-neighbor-TE-upstream.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/neighbor_TEs-length/RLX-neighbor-TE-upstream.bed"))
RLXe <- read.delim("../ParasiTE_output/neighbor_TEs-length/RLX-neighbor-TE-upstream.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/neighbor_TEs-length/RXX-neighbor-TE-upstream.bed"))
RXXe <- read.delim("../ParasiTE_output/neighbor_TEs-length/RXX-neighbor-TE-upstream.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/neighbor_TEs-length/RSX-neighbor-TE-upstream.bed"))
RSXe <- read.delim("../ParasiTE_output/neighbor_TEs-length/RSX-neighbor-TE-upstream.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/neighbor_TEs-length/DTT-neighbor-TE-upstream.bed"))
DTTe <- read.delim("../ParasiTE_output/neighbor_TEs-length/DTT-neighbor-TE-upstream.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/neighbor_TEs-length/DTM-neighbor-TE-upstream.bed"))
DTMe <- read.delim("../ParasiTE_output/neighbor_TEs-length/DTM-neighbor-TE-upstream.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/neighbor_TEs-length/DTA-neighbor-TE-upstream.bed"))
DTAe <- read.delim("../ParasiTE_output/neighbor_TEs-length/DTA-neighbor-TE-upstream.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/neighbor_TEs-length/DTC-neighbor-TE-upstream.bed"))
DTCe <- read.delim("../ParasiTE_output/neighbor_TEs-length/DTC-neighbor-TE-upstream.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/neighbor_TEs-length/DTH-neighbor-TE-upstream.bed"))
DTHe <- read.delim("../ParasiTE_output/neighbor_TEs-length/DTH-neighbor-TE-upstream.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/neighbor_TEs-length/DTX-neighbor-TE-upstream.bed"))
DTXe <- read.delim("../ParasiTE_output/neighbor_TEs-length/DTX-neighbor-TE-upstream.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/neighbor_TEs-length/DHX-neighbor-TE-upstream.bed"))
DHXe <- read.delim("../ParasiTE_output/neighbor_TEs-length/DHX-neighbor-TE-upstream.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/neighbor_TEs-length/TXX-neighbor-TE-upstream.bed"))
TXXe <- read.delim("../ParasiTE_output/neighbor_TEs-length/TXX-neighbor-TE-upstream.bed", header = FALSE, sep = "\t", dec = ".")


if (file.exists("../ParasiTE_output/neighbor_TEs-length/TOTAL-neighbor-TE-upstream.bed"))
TOTALe <- read.delim("../ParasiTE_output/neighbor_TEs-length/TOTAL-neighbor-TE-upstream.bed", header = FALSE, sep = "\t", dec = ".")


#create identifiers for downstream
empty <- data.frame(matrix(ncol = 1, nrow = 1))


#RLC#
#RLC downstream
if (exists("RLCi"))
  { 
RLCinb <- nrow(RLCi)
RLCirep <- rep("RLC", length(RLCinb)) 
RLCidata <- data.frame(RLCi, RLCirep)
  } else {RLCidata <- data.frame(empty, rep("RLC", 1))
  }
colnames(RLCidata) <- c("length","TEs")

#RLC upstream
if (exists("RLCe"))
  {
RLCenb <- nrow(RLCe) 
RLCerep <- rep("RLC", length(RLCenb)) 
RLCedata <- data.frame(RLCe, RLCerep)
  } else {RLCedata <- data.frame(empty, rep("RLC", 1))
  }
colnames(RLCedata) <- c("length","TEs")


#RLG#
#RLG downstream
if (exists("RLGi"))
  {
RLGinb <- nrow(RLGi)  
RLGirep <- rep("RLG", length(RLGinb)) 
RLGidata <- data.frame(RLGi, RLGirep)
  } else {RLGidata <- data.frame(empty, rep("RLG", 1))
  }
colnames(RLGidata) <- c("length","TEs")

#RLG upstream
if (exists("RLGe"))
  {
RLGenb <- nrow(RLGe)    
RLGerep <- rep("RLG", length(RLGenb)) 
RLGedata <- data.frame(RLGe, RLGerep)
  } else {RLGedata <- data.frame(empty, rep("RLG", 1))
  }
colnames(RLGedata) <- c("length","TEs")

#Caulimoviridae#
#Caulimoviridae downstream
if (exists("Caulimoviridaei"))
  { 
Caulimoviridaeinb <- nrow(Caulimoviridaei)      
Caulimoviridaeirep <- rep("Cau", length(Caulimoviridaeinb)) 
Caulimoviridaeidata <- data.frame(Caulimoviridaei, Caulimoviridaeirep)
  } else {Caulimoviridaeidata <- data.frame(empty, rep("Cau", 1))
  }
colnames(Caulimoviridaeidata) <- c("length","TEs")

#Caulimoviridae upstream
if (exists("Caulimoviridaee"))
  { 
Caulimoviridaeenb <- nrow(Caulimoviridaee)    
Caulimoviridaeerep <- rep("Cau", length(Caulimoviridaeenb)) 
Caulimoviridaeedata <- data.frame(Caulimoviridaee, Caulimoviridaeerep)
  } else {Caulimoviridaeedata <- data.frame(empty, rep("Cau", 1))
  }
colnames(Caulimoviridaeedata) <- c("length","TEs")

#RYX#
#RYX downstream
if (exists("RYXi"))
  { 
RYXinb <- nrow(RYXi)
RYXirep <- rep("RYX", length(RYXinb)) 
RYXidata <- data.frame(RYXi, RYXirep)
  } else {RYXidata <- data.frame(empty, rep("RYX", 1))
  }
colnames(RYXidata) <- c("length","TEs")

#RYX upstream
if (exists("RYXe"))
  {
RYXenb <- nrow(RYXe)  
RYXerep <- rep("RYX", length(RYXenb)) 
RYXedata <- data.frame(RYXe, RYXerep)
  } else {RYXedata <- data.frame(empty, rep("RYX", 1))
  }
colnames(RYXedata) <- c("length","TEs")


#RIX#
#RIX downstream
if (exists("RIXi"))
  { 
RIXinb <- nrow(RIXi)
RIXirep <- rep("RIX", length(RIXinb)) 
RIXidata <- data.frame(RIXi, RIXirep)
  } else {RIXidata <- data.frame(empty, rep("RIX", 1))
  }
colnames(RIXidata) <- c("length","TEs")

#RIX upstream
if (exists("RIXe"))
  {
RIXenb <- nrow(RIXe)  
RIXerep <- rep("RIX", length(RIXenb)) 
RIXedata <- data.frame(RIXe, RIXerep)
  } else {RIXedata <- data.frame(empty, rep("RIX", 1))
  }
colnames(RIXedata) <- c("length","TEs")

#RLX#
#RLX downstream
if (exists("RLXi"))
  { 
RLXinb <- nrow(RLXi)   
RLXirep <- rep("RLX", length(RLXinb)) 
RLXidata <- data.frame(RLXi, RLXirep)
  } else {RLXidata <- data.frame(empty, rep("RLX", 1))
  }
colnames(RLXidata) <- c("length","TEs")

#RLX upstream
if (exists("RLXe"))
  { 
RLXenb <- nrow(RLXe)
RLXerep <- rep("RLX", length(RLXenb)) 
RLXedata <- data.frame(RLXe, RLXerep)
  } else {RLXedata <- data.frame(empty, rep("RLX", 1))
  }
colnames(RLXedata) <- c("length","TEs")

#RXX#
#RXX downstream
if (exists("RXXi"))
  { 
RXXinb <- nrow(RXXi)
RXXirep <- rep("RXX", length(RXXinb)) 
RXXidata <- data.frame(RXXi, RXXirep)
  } else {RXXidata <- data.frame(empty, rep("RXX", 1))
  }
colnames(RXXidata) <- c("length","TEs")

#RXX upstream
if (exists("RXXe"))
  { 
RXXenb <- nrow(RXXe)
RXXerep <- rep("RXX", length(RXXenb)) 
RXXedata <- data.frame(RXXe, RXXerep)
  } else {RXXedata <- data.frame(empty, rep("RXX", 1))
  }
colnames(RXXedata) <- c("length","TEs")

#RSX#
#RSX downstream
if (exists("RSXi"))
  { 
RSXinb <- nrow(RSXi)
RSXirep <- rep("RSX", length(RSXinb)) 
RSXidata <- data.frame(RSXi, RSXirep)
  } else {RSXidata <- data.frame(empty, rep("RSX", 1))
  }
colnames(RSXidata) <- c("length","TEs")

#RSX upstream
if (exists("RSXe"))
  { 
RSXenb <- nrow(RSXe)
RSXerep <- rep("RSX", length(RSXenb)) 
RSXedata <- data.frame(RSXe, RSXerep)
  } else {RSXedata <- data.frame(empty, rep("RSX", 1))
  }
colnames(RSXedata) <- c("length","TEs")

#DTM#
#DTM downstream
if (exists("DTMi"))
  { 
DTMinb <- nrow(DTMi)  
DTMirep <- rep("DTM", length(DTMinb)) 
DTMidata <- data.frame(DTMi, DTMirep)
  } else {DTMidata <- data.frame(empty, rep("DTM", 1))
  }
colnames(DTMidata) <- c("length","TEs")

#DTM upstream
if (exists("DTMe"))
  { 
DTMenb <- nrow(DTMe) 
DTMerep <- rep("DTM", length(DTMenb)) 
DTMedata <- data.frame(DTMe, DTMerep)
  } else {DTMedata <- data.frame(empty, rep("DTM", 1))
  }
colnames(DTMedata) <- c("length","TEs")

#DTA#
#DTA downstream
if (exists("DTAi"))
  { 
DTAinb <- nrow(DTAi) 
DTAirep <- rep("DTA", length(DTAinb)) 
DTAidata <- data.frame(DTAi, DTAirep)
  } else {DTAidata <- data.frame(empty, rep("DTA", 1))
  }
colnames(DTAidata) <- c("length","TEs")

#DTA upstream
if (exists("DTAe"))
  { 
DTAenb <- nrow(DTAe) 
DTAerep <- rep("DTA", length(DTAenb)) 
DTAedata <- data.frame(DTAe, DTAerep)
  } else {DTAedata <- data.frame(empty, rep("DTA", 1))
  }
colnames(DTAedata) <- c("length","TEs")

#DTC#
#DTC downstream
if (exists("DTCi"))
  { 
DTCinb <- nrow(DTCi)
DTCirep <- rep("DTC", length(DTCinb)) 
DTCidata <- data.frame(DTCi, DTCirep)
  } else {DTCidata <- data.frame(empty, rep("DTC", 1))
  }
colnames(DTCidata) <- c("length","TEs")

#DTC upstream
if (exists("DTCe"))
  { 
DTCenb <- nrow(DTCe)
DTCerep <- rep("DTC", length(DTCenb)) 
DTCedata <- data.frame(DTCe, DTCerep)
  } else {DTCedata <- data.frame(empty, rep("DTC", 1))
  }
colnames(DTCedata) <- c("length","TEs")

#DTH#
#DTH downstream
if (exists("DTHi"))
  { 
DTHinb <- nrow(DTHi)
DTHirep <- rep("DTH", length(DTHinb)) 
DTHidata <- data.frame(DTHi, DTHirep)
  } else {DTHidata <- data.frame(empty, rep("DTH", 1))
  }
colnames(DTHidata) <- c("length","TEs")

#DTH upstream
if (exists("DTHe"))
  { 
DTHenb <- nrow(DTHe)  
DTHerep <- rep("DTH", length(DTHenb)) 
DTHedata <- data.frame(DTHe, DTHerep)
  } else {DTHedata <- data.frame(empty, rep("DTH", 1))
  }
colnames(DTHedata) <- c("length","TEs")

#DTT#
#DTT downstream
if (exists("DTTi"))
  { 
DTTinb <- nrow(DTTi)  
DTTirep <- rep("DTT", length(DTTinb)) 
DTTidata <- data.frame(DTTi, DTTirep)
  } else {DTTidata <- data.frame(empty, rep("DTT", 1))
  }
colnames(DTTidata) <- c("length","TEs")

#DTT upstream
if (exists("DTTe"))
  {
DTTenb <- nrow(DTTe)    
DTTerep <- rep("DTT", length(DTTenb)) 
DTTedata <- data.frame(DTTe, DTTerep)
  } else {DTTedata <- data.frame(empty, rep("DTT", 1))
  }
colnames(DTTedata) <- c("length","TEs")

#DTX#
#DTX downstream
if (exists("DTXi"))
  { 
DTXinb <- nrow(DTXi)  
DTXirep <- rep("DTX", length(DTXinb)) 
DTXidata <- data.frame(DTXi, DTXirep)
  } else {DTXidata <- data.frame(empty, rep("DTX", 1))
  }
colnames(DTXidata) <- c("length","TEs")

#DTX upstream
if (exists("DTXe"))
  { 
DTXenb <- nrow(DTXe)  
DTXerep <- rep("DTX", length(DTXenb)) 
DTXedata <- data.frame(DTXe, DTXerep)
  } else {DTXedata <- data.frame(empty, rep("DTX", 1))
  }
colnames(DTXedata) <- c("length","TEs")

#DHX#
#DHX downstream
if (exists("DHXi"))
  {
DHXinb <- nrow(DHXi)    
DHXirep <- rep("DHX", length(DHXinb)) 
DHXidata <- data.frame(DHXi, DHXirep)
  } else {DHXidata <- data.frame(empty, rep("DHX", 1))
  }
colnames(DHXidata) <- c("length","TEs")

#DHX upstream
if (exists("DHXe"))
  {
DHXenb <- nrow(DHXe)    
DHXerep <- rep("DHX", length(DHXenb)) 
DHXedata <- data.frame(DHXe, DHXerep)
  } else {DHXedata <- data.frame(empty, rep("DHX", 1))
  }
colnames(DHXedata) <- c("length","TEs")

#TXX#
#TXX downstream
if (exists("TXXi"))
  {
TXXinb <- nrow(TXXi)    
TXXirep <- rep("TXX", length(TXXinb)) 
TXXidata <- data.frame(TXXi, TXXirep)
  } else {TXXidata <- data.frame(empty, rep("TXX", 1))
  }
colnames(TXXidata) <- c("length","TEs")

#TXX upstream
if (exists("TXXe"))
  {
TXXenb <- nrow(TXXe)    
TXXerep <- rep("TXX", length(TXXenb)) 
TXXedata <- data.frame(TXXe, TXXerep)
  } else {TXXedata <- data.frame(empty, rep("TXX", 1))
  }
colnames(TXXedata) <- c("length","TEs")

#Plot downstream
iClassITEs <- rbind(RLCidata, RLGidata, Caulimoviridaeidata, RIXidata, RYXidata, RLXidata, RSXidata, RXXidata)
iClassIITEs <- rbind(DTTidata, DTMidata, DTAidata, DTCidata, DTHidata, DTXidata, DHXidata, TXXidata)

iClassITEsplot  <- ggplot(iClassITEs, aes(x=length, color=TEs)) + geom_freqpoly() + theme_grey() + ylab("Number of TEs") + xlab("Length") + theme(legend.title = element_blank()) + ggtitle("Downstream")
iClassIITEsplot  <- ggplot(iClassIITEs, aes(x=length, color=TEs)) + geom_freqpoly() + theme_grey() + ylab("Number of TEs") + xlab("Length") + theme(legend.title = element_blank())+ ggtitle(" ")


#Plot upstream
eClassITEs <- rbind(RLCedata, RLGedata, Caulimoviridaeedata, RIXedata, RYXedata, RLXedata, RSXedata, RXXedata) 
eClassIITEs <- rbind(DTTedata, DTMedata, DTAedata, DTCedata, DTHedata, DTXedata, DHXedata, TXXedata)

eClassITEsplot  <- ggplot(eClassITEs, aes(x=length, color=TEs)) + geom_freqpoly() + theme_grey() + ylab("Number of TEs") + xlab("Length") + theme(legend.title = element_blank()) + ggtitle("Upstream")
eClassIITEsplot  <- ggplot(eClassIITEs, aes(x=length, color=TEs)) + geom_freqpoly() + theme_grey() + ylab("Number of TEs") + xlab("Length") + theme(legend.title = element_blank())+ ggtitle(" ")


pdf(file = "../ParasiTE_output/Results/Plots/neighbor_TEs_2.pdf")
grid.arrange(iClassITEsplot,iClassIITEsplot, eClassITEsplot, eClassIITEsplot, top="Frequency of the length of neighbor TE events")
dev.off()


