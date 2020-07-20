#!/usr/bin/env Rscript

#load libraries
library(ggplot2)
library(cowplot)


#Global data intragenic TEs #########################################################################################################################################################################################################################################################################################################################################################################################################################################################################

if (file.exists("../ParasiTE_output/Results/Annotations_TEs/Intragenic_intronic_TEs.bed"))
intronic <- read.delim("../ParasiTE_output/Results/Annotations_TEs/Intragenic_intronic_TEs.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/Results/Annotations_TEs/Intragenic_exonic_TEs.bed"))
exonic <- read.delim("../ParasiTE_output/Results/Annotations_TEs/Intragenic_exonic_TEs.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/Results/Annotations_TEs/Intragenic_ambigous_TEs.bed"))
ambigous <- read.delim("../ParasiTE_output/Results/Annotations_TEs/Intragenic_ambigous_TEs.bed", header = FALSE, sep = "\t", dec = ".")


#create identifiers intragenic

emptydata <- data.frame(matrix(ncol = 1, nrow = 1))

intragenicrep <- rep("Intragenic", length(1))
intronicrep <- rep("Intronic", length(1))
exonicrep <- rep("Exonic", length(1))
ambigousrep <- rep("Ambigous", length(1))

#intronic
if (exists("intronic"))
  { 
intronicnb <- nrow(intronic)
intronic_ident <- data.frame(intragenicrep, intronicrep, intronicnb)  
  } else {intronic_ident <- data.frame(intragenicrep, intronicrep, emptydata)
  }
colnames(intronic_ident) <- c("class", "group","number" )

#exonic
if (exists("exonic"))
  { 
exonicnb <- nrow(exonic)
exonic_ident <- data.frame(intragenicrep, exonicrep, exonicnb)  
  } else {exonic_ident <- data.frame(intragenicrep, exonicrep, emptydata)
  }
colnames(exonic_ident) <- c("class", "group","number" )

#ambigous
if (exists("ambigous"))
  { 
ambigousnb <- nrow(ambigous)
ambigous_ident <- data.frame(intragenicrep, ambigousrep, ambigousnb)  
  } else {ambigous_ident <- data.frame(intragenicrep, ambigousrep, emptydata)
  }
colnames(ambigous_ident) <- c("class", "group","number" )


tot_intragenic <- rbind(ambigous_ident, exonic_ident, intronic_ident)



#Proportion of Intronic TEs ###################################################################################################################################


#check if data exist and load them
#intragenic
if (file.exists("../ParasiTE_output/intragenic_TEs/RLC-intronic_TEs.bed"))
RLCi <- read.delim("../ParasiTE_output/intragenic_TEs/RLC-intronic_TEs.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/intragenic_TEs/RLG-intronic_TEs.bed"))
RLGi <- read.delim("../ParasiTE_output/intragenic_TEs/RLG-intronic_TEs.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/intragenic_TEs/Caulimoviridae-intronic_TEs.bed"))
Caulimoviridaei <- read.delim("../ParasiTE_output/intragenic_TEs/Caulimoviridae-intronic_TEs.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/intragenic_TEs/RYX-intronic_TEs.bed"))
RYXi <- read.delim("../ParasiTE_output/intragenic_TEs/RYX-intronic_TEs.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/intragenic_TEs/RIX-intronic_TEs.bed"))
RIXi <- read.delim("../ParasiTE_output/intragenic_TEs/RIX-intronic_TEs.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/intragenic_TEs/RLX-intronic_TEs.bed"))
RLXi <- read.delim("../ParasiTE_output/intragenic_TEs/RLX-intronic_TEs.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/intragenic_TEs/RXX-intronic_TEs.bed"))
RXXi <- read.delim("../ParasiTE_output/intragenic_TEs/RXX-intronic_TEs.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/intragenic_TEs/RSX-intronic_TEs.bed"))
RSXi <- read.delim("../ParasiTE_output/intragenic_TEs/RSX-intronic_TEs.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/intragenic_TEs/DTT-intronic_TEs.bed"))
DTTi <- read.delim("../ParasiTE_output/intragenic_TEs/DTT-intronic_TEs.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/intragenic_TEs/DTM-intronic_TEs.bed"))
DTMi <- read.delim("../ParasiTE_output/intragenic_TEs/DTM-intronic_TEs.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/intragenic_TEs/DTA-intronic_TEs.bed"))
DTAi <- read.delim("../ParasiTE_output/intragenic_TEs/DTA-intronic_TEs.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/intragenic_TEs/DTC-intronic_TEs.bed"))
DTCi <- read.delim("../ParasiTE_output/intragenic_TEs/DTC-intronic_TEs.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/intragenic_TEs/DTH-intronic_TEs.bed"))
DTHi <- read.delim("../ParasiTE_output/intragenic_TEs/DTH-intronic_TEs.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/intragenic_TEs/DTX-intronic_TEs.bed"))
DTXi <- read.delim("../ParasiTE_output/intragenic_TEs/DTX-intronic_TEs.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/intragenic_TEs/DHX-intronic_TEs.bed"))
DHXi <- read.delim("../ParasiTE_output/intragenic_TEs/DHX-intronic_TEs.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/intragenic_TEs/TXX-intronic_TEs.bed"))
TXXi <- read.delim("../ParasiTE_output/intragenic_TEs/TXX-intronic_TEs.bed", header = FALSE, sep = "\t", dec = ".")

#if (file.exists("../ParasiTE_output/intragenic_TEs/TOTAL-intronic_TEs.bed"))
#TOTALi <- read.delim("../ParasiTE_output/intragenic_TEs/TOTAL-intronic_TEs.bed", header = FALSE, sep = "\t", dec = ".")

#exonic
if (file.exists("../ParasiTE_output/intragenic_TEs/RLC-exonic_TEs.bed"))
RLCe <- read.delim("../ParasiTE_output/intragenic_TEs/RLC-exonic_TEs.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/intragenic_TEs/RLG-exonic_TEs.bed"))
RLGe <- read.delim("../ParasiTE_output/intragenic_TEs/RLG-exonic_TEs.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/intragenic_TEs/Caulimoviridae-exonic_TEs.bed"))
Caulimoviridaee <- read.delim("../ParasiTE_output/intragenic_TEs/Caulimoviridae-exonic_TEs.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/intragenic_TEs/RYX-exonic_TEs.bed"))
RYXe <- read.delim("../ParasiTE_output/intragenic_TEs/RYX-exonic_TEs.bed", header = FALSE, sep = "\t", dec = ".")


if (file.exists("../ParasiTE_output/intragenic_TEs/RIX-exonic_TEs.bed"))
RIXe <- read.delim("../ParasiTE_output/intragenic_TEs/RIX-exonic_TEs.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/intragenic_TEs/RLX-exonic_TEs.bed"))
RLXe <- read.delim("../ParasiTE_output/intragenic_TEs/RLX-exonic_TEs.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/intragenic_TEs/RXX-exonic_TEs.bed"))
RXXe <- read.delim("../ParasiTE_output/intragenic_TEs/RXX-exonic_TEs.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/intragenic_TEs/RSX-exonic_TEs.bed"))
RSXe <- read.delim("../ParasiTE_output/intragenic_TEs/RSX-exonic_TEs.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/intragenic_TEs/DTT-exonic_TEs.bed"))
DTTe <- read.delim("../ParasiTE_output/intragenic_TEs/DTT-exonic_TEs.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/intragenic_TEs/DTM-exonic_TEs.bed"))
DTMe <- read.delim("../ParasiTE_output/intragenic_TEs/DTM-exonic_TEs.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/intragenic_TEs/DTA-exonic_TEs.bed"))
DTAe <- read.delim("../ParasiTE_output/intragenic_TEs/DTA-exonic_TEs.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/intragenic_TEs/DTC-exonic_TEs.bed"))
DTCe <- read.delim("../ParasiTE_output/intragenic_TEs/DTC-exonic_TEs.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/intragenic_TEs/DTH-exonic_TEs.bed"))
DTHe <- read.delim("../ParasiTE_output/intragenic_TEs/DTH-exonic_TEs.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/intragenic_TEs/DTX-exonic_TEs.bed"))
DTXe <- read.delim("../ParasiTE_output/intragenic_TEs/DTX-exonic_TEs.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/intragenic_TEs/DHX-exonic_TEs.bed"))
DHXe <- read.delim("../ParasiTE_output/intragenic_TEs/DHX-exonic_TEs.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/intragenic_TEs/TXX-exonic_TEs.bed"))
TXXe <- read.delim("../ParasiTE_output/intragenic_TEs/TXX-exonic_TEs.bed", header = FALSE, sep = "\t", dec = ".")

#if (file.exists("../ParasiTE_output/intragenic_TEs/TOTAL-exonic_TEs.bed"))
#TOTALe <- read.delim("../ParasiTE_output/intragenic_TEs/TOTAL-exonic_TEs.bed", header = FALSE, sep = "\t", dec = ".")

#ambigous
if (file.exists("../ParasiTE_output/intragenic_TEs/RLC-ambigous_TEs.bed"))
RLCa <- read.delim("../ParasiTE_output/intragenic_TEs/RLC-ambigous_TEs.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/intragenic_TEs/RLG-ambigous_TEs.bed"))
RLGa <- read.delim("../ParasiTE_output/intragenic_TEs/RLG-ambigous_TEs.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/intragenic_TEs/Caulimoviridae-ambigous_TEs.bed"))
Caulimoviridaea <- read.delim("../ParasiTE_output/intragenic_TEs/Caulimoviridae-ambigous_TEs.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/intragenic_TEs/RYX-ambigous_TEs.bed"))
RYXa <- read.delim("../ParasiTE_output/intragenic_TEs/RYX-ambigous_TEs.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/intragenic_TEs/RIX-ambigous_TEs.bed"))
RIXa <- read.delim("../ParasiTE_output/intragenic_TEs/RIX-ambigous_TEs.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/intragenic_TEs/RLX-ambigous_TEs.bed"))
RLXa <- read.delim("../ParasiTE_output/intragenic_TEs/RLX-ambigous_TEs.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/intragenic_TEs/RXX-ambigous_TEs.bed"))
RXXa <- read.delim("../ParasiTE_output/intragenic_TEs/RXX-ambigous_TEs.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/intragenic_TEs/RSX-ambigous_TEs.bed"))
RSXa <- read.delim("../ParasiTE_output/intragenic_TEs/RSX-ambigous_TEs.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/intragenic_TEs/DTT-ambigous_TEs.bed"))
DTTa <- read.delim("../ParasiTE_output/intragenic_TEs/DTT-ambigous_TEs.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/intragenic_TEs/DTM-ambigous_TEs.bed"))
DTMa <- read.delim("../ParasiTE_output/intragenic_TEs/DTM-ambigous_TEs.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/intragenic_TEs/DTA-ambigous_TEs.bed"))
DTAa <- read.delim("../ParasiTE_output/intragenic_TEs/DTA-ambigous_TEs.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/intragenic_TEs/DTC-ambigous_TEs.bed"))
DTCa <- read.delim("../ParasiTE_output/intragenic_TEs/DTC-ambigous_TEs.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/intragenic_TEs/DTH-ambigous_TEs.bed"))
DTHa <- read.delim("../ParasiTE_output/intragenic_TEs/DTH-ambigous_TEs.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/intragenic_TEs/DTX-ambigous_TEs.bed"))
DTXa <- read.delim("../ParasiTE_output/intragenic_TEs/DTX-ambigous_TEs.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/intragenic_TEs/DHX-ambigous_TEs.bed"))
DHXa <- read.delim("../ParasiTE_output/intragenic_TEs/DHX-ambigous_TEs.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/intragenic_TEs/TXX-ambigous_TEs.bed"))
TXXa <- read.delim("../ParasiTE_output/intragenic_TEs/TXX-ambigous_TEs.bed", header = FALSE, sep = "\t", dec = ".")

#if (file.exists("../ParasiTE_output/intragenic_TEs/TOTAL-ambigous_TEs.bed"))
#TOTALa <- read.delim("../ParasiTE_output/intragenic_TEs/TOTAL-ambigous_TEs.bed", header = FALSE, sep = "\t", dec = ".")

# 1st step calculate the proportion of the intragenic TE familly

intronicrep <- rep("Intronic", 1)
exonicrep <- rep("Exonic", 1)
ambigousrep <- rep("Ambigous", 1)

RLCrep <- rep("RLC", 1)
RLGrep <- rep("RLG", 1)
Caulimoviridaerep <- rep("Cau", 1)
RYXrep <- rep("RYX", 1)
RIXrep <- rep("RIX", 1)
RLXrep <- rep("RLX", 1)
RXXrep <- rep("RXX", 1)
RSXrep <- rep("RSX", 1)

DTTrep <- rep("DTT", 1)
DTMrep <- rep("DTM", 1)
DTArep <- rep("DTA", 1)
DTCrep <- rep("DTC", 1)
DTHrep <- rep("DTH", 1)
DTXrep <- rep("DTX", 1)
DHXrep <- rep("DHX", 1)
TXXrep <- rep("TXX", 1)
TOTALrep <- rep("TOTAL", 1)
empty <- data.frame(matrix(ncol = 1, nrow = 1))

#RLC#
#RLC intronic
if (exists("RLCi"))
  { 
RLCinb <- nrow(RLCi)
RLCidata <- data.frame(RLCinb, RLCrep, intronicrep)  
  } else {RLCidata <- data.frame(empty, RLCrep, intronicrep)
  }
colnames(RLCidata) <- c("number","TEs", "group")

#RLC exonic
if (exists("RLCe"))
  { 
RLCenb <- nrow(RLCe)
RLCedata <- data.frame(RLCenb, RLCrep, exonicrep)
  } else {RLCedata <- data.frame(empty, RLCrep, exonicrep)
  }
colnames(RLCedata) <- c("number","TEs", "group")


#RLC ambigous
if (exists("RLCa"))
  { 
RLCanb <- nrow(RLCa)
RLCadata <- data.frame(RLCanb, RLCrep, ambigousrep)
  } else {RLCadata <- data.frame(empty, RLCrep, ambigousrep)
  }
colnames(RLCadata) <- c("number","TEs", "group")

#RLG
#RLG intronic
if (exists("RLGi"))
  { 
RLGinb <- nrow(RLGi)
RLGidata <- data.frame(RLGinb, RLGrep, intronicrep)
  } else {RLGidata <- data.frame(empty, RLGrep, intronicrep)
  }
colnames(RLGidata) <- c("number","TEs", "group")

#RLG exonic
if (exists("RLGe"))
  { 
RLGenb <- nrow(RLGe)
RLGedata <- data.frame(RLGenb, RLGrep, exonicrep)
  } else {RLGedata <- data.frame(empty, RLGrep, exonicrep)
  }
colnames(RLGedata) <- c("number","TEs", "group")

#RLG ambigous
if (exists("RLGa"))
  { 
RLGanb <- nrow(RLGa)
RLGadata <- data.frame(RLGanb, RLGrep, ambigousrep)
  } else {RLGadata <- data.frame(empty, RLGrep, ambigousrep)
  }
colnames(RLGadata) <- c("number","TEs", "group")

#Caulimoviridae#
#Caulimoviridae intronic

if (exists("Caulimoviridaei"))
  { 
Caulimoviridaeinb <- nrow(Caulimoviridaei)
Caulimoviridaeidata <- data.frame(Caulimoviridaeinb, Caulimoviridaerep, intronicrep)
  } else {Caulimoviridaeidata <- data.frame(empty, Caulimoviridaerep, intronicrep)
  }
colnames(Caulimoviridaeidata) <- c("number","TEs", "group")

#Caulimoviridae exonic

if (exists("Caulimoviridaee"))
  { 
Caulimoviridaeenb <- nrow(Caulimoviridaee)
Caulimoviridaeedata <- data.frame(Caulimoviridaeenb, Caulimoviridaerep, exonicrep)
  } else {Caulimoviridaeedata <- data.frame(empty, Caulimoviridaerep, exonicrep)
  }
colnames(Caulimoviridaeedata) <- c("number","TEs", "group")

#Caulimoviridae ambigous

if (exists("Caulimoviridaea"))
  { 
Caulimoviridaeanb <- nrow(Caulimoviridaea)
Caulimoviridaeadata <- data.frame(Caulimoviridaeanb, Caulimoviridaerep, ambigousrep)
  } else {Caulimoviridaeadata <- data.frame(empty, Caulimoviridaerep, ambigousrep)
  }
colnames(Caulimoviridaeadata) <- c("number","TEs", "group")

#RYX#
#RYX intronic
if (exists("RYXi"))
  { 
RYXinb <- nrow(RYXi)
RYXidata <- data.frame(RYXinb, RYXrep, intronicrep)  
  } else {RYXidata <- data.frame(empty, RYXrep, intronicrep)
  }
colnames(RYXidata) <- c("number","TEs", "group")

#RYX exonic
if (exists("RYXe"))
  { 
RYXenb <- nrow(RYXe)
RYXedata <- data.frame(RYXenb, RYXrep, exonicrep)
  } else {RYXedata <- data.frame(empty, RYXrep, exonicrep)
  }
colnames(RYXedata) <- c("number","TEs", "group")


#RYX ambigous
if (exists("RYXa"))
  { 
RYXanb <- nrow(RYXa)
RYXadata <- data.frame(RYXanb, RYXrep, ambigousrep)
  } else {RYXadata <- data.frame(empty, RYXrep, ambigousrep)
  }
colnames(RYXadata) <- c("number","TEs", "group")

#RIX#
#RIX intronic
if (exists("RIXi"))
  { 
RIXinb <- nrow(RIXi)
RIXidata <- data.frame(RIXinb, RIXrep, intronicrep)  
  } else {RIXidata <- data.frame(empty, RIXrep, intronicrep)
  }
colnames(RIXidata) <- c("number","TEs", "group")

#RIX exonic
if (exists("RIXe"))
  { 
RIXenb <- nrow(RIXe)
RIXedata <- data.frame(RIXenb, RIXrep, exonicrep)
  } else {RIXedata <- data.frame(empty, RIXrep, exonicrep)
  }
colnames(RIXedata) <- c("number","TEs", "group")


#RIX ambigous
if (exists("RIXa"))
  { 
RIXanb <- nrow(RIXa)
RIXadata <- data.frame(RIXanb, RIXrep, ambigousrep)
  } else {RIXadata <- data.frame(empty, RIXrep, ambigousrep)
  }
colnames(RIXadata) <- c("number","TEs", "group")

#RLX#
#RLX intronic
if (exists("RLXi"))
  { 
RLXinb <- nrow(RLXi)
RLXidata <- data.frame(RLXinb, RLXrep, intronicrep)  
  } else {RLXidata <- data.frame(empty, RLXrep, intronicrep)
  }
colnames(RLXidata) <- c("number","TEs", "group")

#RLX exonic
if (exists("RLXe"))
  { 
RLXenb <- nrow(RLXe)
RLXedata <- data.frame(RLXenb, RLXrep, exonicrep)
  } else {RLXedata <- data.frame(empty, RLXrep, exonicrep)
  }
colnames(RLXedata) <- c("number","TEs", "group")


#RLX ambigous
if (exists("RLXa"))
  { 
RLXanb <- nrow(RLXa)
RLXadata <- data.frame(RLXanb, RLXrep, ambigousrep)
  } else {RLXadata <- data.frame(empty, RLXrep, ambigousrep)
  }
colnames(RLXadata) <- c("number","TEs", "group")


#RXX#
#RXX intronic
if (exists("RXXi"))
  { 
RXXinb <- nrow(RXXi)
RXXidata <- data.frame(RXXinb, RXXrep, intronicrep)  
  } else {RXXidata <- data.frame(empty, RXXrep, intronicrep)
  }
colnames(RXXidata) <- c("number","TEs", "group")

#RXX exonic
if (exists("RXXe"))
  { 
RXXenb <- nrow(RXXe)
RXXedata <- data.frame(RXXenb, RXXrep, exonicrep)
  } else {RXXedata <- data.frame(empty, RXXrep, exonicrep)
  }
colnames(RXXedata) <- c("number","TEs", "group")


#RXX ambigous
if (exists("RXXa"))
  { 
RXXanb <- nrow(RXXa)
RXXadata <- data.frame(RXXanb, RXXrep, ambigousrep)
  } else {RXXadata <- data.frame(empty, RXXrep, ambigousrep)
  }
colnames(RXXadata) <- c("number","TEs", "group")
  
#RSX#
#RSX intronic
if (exists("RSXi"))
  { 
RSXinb <- nrow(RSXi)
RSXidata <- data.frame(RSXinb, RSXrep, intronicrep)  
  } else {RSXidata <- data.frame(empty, RSXrep, intronicrep)
  }
colnames(RSXidata) <- c("number","TEs", "group")

#RSX exonic
if (exists("RSXe"))
  { 
RSXenb <- nrow(RSXe)
RSXedata <- data.frame(RSXenb, RSXrep, exonicrep)
  } else {RSXedata <- data.frame(empty, RSXrep, exonicrep)
  }
colnames(RSXedata) <- c("number","TEs", "group")


#RSX ambigous
if (exists("RSXa"))
  { 
RSXanb <- nrow(RSXa)
RSXadata <- data.frame(RSXanb, RSXrep, ambigousrep)
  } else {RSXadata <- data.frame(empty, RSXrep, ambigousrep)
  }
colnames(RSXadata) <- c("number","TEs", "group")


#DTT#
#DTT intronic
if (exists("DTTi"))
  { 
DTTinb <- nrow(DTTi)
DTTidata <- data.frame(DTTinb, DTTrep, intronicrep)  
  } else {DTTidata <- data.frame(empty, DTTrep, intronicrep)
  }
colnames(DTTidata) <- c("number","TEs", "group")

#DTT exonic
if (exists("DTTe"))
  { 
DTTenb <- nrow(DTTe)
DTTedata <- data.frame(DTTenb, DTTrep, exonicrep)
  } else {DTTedata <- data.frame(empty, DTTrep, exonicrep)
  }
colnames(DTTedata) <- c("number","TEs", "group")


#DTT ambigous
if (exists("DTTa"))
  { 
DTTanb <- nrow(DTTa)
DTTadata <- data.frame(DTTanb, DTTrep, ambigousrep)
  } else {DTTadata <- data.frame(empty, DTTrep, ambigousrep)
  }
colnames(DTTadata) <- c("number","TEs", "group")


#DTM#
#DTM intronic
if (exists("DTMi"))
  { 
DTMinb <- nrow(DTMi)
DTMidata <- data.frame(DTMinb, DTMrep, intronicrep)  
  } else {DTMidata <- data.frame(empty, DTMrep, intronicrep)
  }
colnames(DTMidata) <- c("number","TEs", "group")

#DTM exonic
if (exists("DTMe"))
  { 
DTMenb <- nrow(DTMe)
DTMedata <- data.frame(DTMenb, DTMrep, exonicrep)
  } else {DTMedata <- data.frame(empty, DTMrep, exonicrep)
  }
colnames(DTMedata) <- c("number","TEs", "group")


#DTM ambigous
if (exists("DTMa"))
  { 
DTManb <- nrow(DTMa)
DTMadata <- data.frame(DTManb, DTMrep, ambigousrep)
  } else {DTMadata <- data.frame(empty, DTMrep, ambigousrep)
  }
colnames(DTMadata) <- c("number","TEs", "group")

#DTA#
#DTA intronic
if (exists("DTAi"))
  { 
DTAinb <- nrow(DTAi)
DTAidata <- data.frame(DTAinb, DTArep, intronicrep)  
  } else {DTAidata <- data.frame(empty, DTArep, intronicrep)
  }
colnames(DTAidata) <- c("number","TEs", "group")

#DTA exonic
if (exists("DTAe"))
  { 
DTAenb <- nrow(DTAe)
DTAedata <- data.frame(DTAenb, DTArep, exonicrep)
  } else {DTAedata <- data.frame(empty, DTArep, exonicrep)
  }
colnames(DTAedata) <- c("number","TEs", "group")


#DTA ambigous
if (exists("DTAa"))
  { 
DTAanb <- nrow(DTAa)
DTAadata <- data.frame(DTAanb, DTArep, ambigousrep)
  } else {DTAadata <- data.frame(empty, DTArep, ambigousrep)
  }
colnames(DTAadata) <- c("number","TEs", "group")

#DTC
#DTC#
#DTC intronic
if (exists("DTCi"))
  { 
DTCinb <- nrow(DTCi)
DTCidata <- data.frame(DTCinb, DTCrep, intronicrep)  
  } else {DTCidata <- data.frame(empty, DTCrep, intronicrep)
  }
colnames(DTCidata) <- c("number","TEs", "group")

#DTC exonic
if (exists("DTCe"))
  { 
DTCenb <- nrow(DTCe)
DTCedata <- data.frame(DTCenb, DTCrep, exonicrep)
  } else {DTCedata <- data.frame(empty, DTCrep, exonicrep)
  }
colnames(DTCedata) <- c("number","TEs", "group")


#DTC ambigous
if (exists("DTCa"))
  { 
DTCanb <- nrow(DTCa)
DTCadata <- data.frame(DTCanb, DTCrep, ambigousrep)
  } else {DTCadata <- data.frame(empty, DTCrep, ambigousrep)
  }
colnames(DTCadata) <- c("number","TEs", "group")


#DTH
#DTH intronic
if (exists("DTHi"))
  { 
DTHinb <- nrow(DTHi)
DTHidata <- data.frame(DTHinb, DTHrep, intronicrep)  
  } else {DTHidata <- data.frame(empty, DTHrep, intronicrep)
  }
colnames(DTHidata) <- c("number","TEs", "group")

#DTH exonic
if (exists("DTHe"))
  { 
DTHenb <- nrow(DTHe)
DTHedata <- data.frame(DTHenb, DTHrep, exonicrep)
  } else {DTHedata <- data.frame(empty, DTHrep, exonicrep)
  }
colnames(DTHedata) <- c("number","TEs", "group")


#DTH ambigous
if (exists("DTHa"))
  { 
DTHanb <- nrow(DTHa)
DTHadata <- data.frame(DTHanb, DTHrep, ambigousrep)
  } else {DTHadata <- data.frame(empty, DTHrep, ambigousrep)
  }
colnames(DTHadata) <- c("number","TEs", "group")


#DTX#
#DTX intronic
if (exists("DTXi"))
  { 
DTXinb <- nrow(DTXi)
DTXidata <- data.frame(DTXinb, DTXrep, intronicrep)  
  } else {DTXidata <- data.frame(empty, DTXrep, intronicrep)
  }
colnames(DTXidata) <- c("number","TEs", "group")

#DTX exonic
if (exists("DTXe"))
  { 
DTXenb <- nrow(DTXe)
DTXedata <- data.frame(DTXenb, DTXrep, exonicrep)
  } else {DTXedata <- data.frame(empty, DTXrep, exonicrep)
  }
colnames(DTXedata) <- c("number","TEs", "group")


#DTX ambigous
if (exists("DTXa"))
  { 
DTXanb <- nrow(DTXa)
DTXadata <- data.frame(DTXanb, DTXrep, ambigousrep)
  } else {DTXadata <- data.frame(empty, DTXrep, ambigousrep)
  }
colnames(DTXadata) <- c("number","TEs", "group")

#DHX#
#DHX intronic
if (exists("DHXi"))
  { 
DHXinb <- nrow(DHXi)
DHXidata <- data.frame(DHXinb, DHXrep, intronicrep)  
  } else {DHXidata <- data.frame(empty, DHXrep, intronicrep)
  }
colnames(DHXidata) <- c("number","TEs", "group")

#DHX exonic
if (exists("DHXe"))
  { 
DHXenb <- nrow(DHXe)
DHXedata <- data.frame(DHXenb, DHXrep, exonicrep)
  } else {DHXedata <- data.frame(empty, DHXrep, exonicrep)
  }
colnames(DHXedata) <- c("number","TEs", "group")


#DHX ambigous
if (exists("DHXa"))
  { 
DHXanb <- nrow(DHXa)
DHXadata <- data.frame(DHXanb, DHXrep, ambigousrep)
  } else {DHXadata <- data.frame(empty, DHXrep, ambigousrep)
  }
colnames(DHXadata) <- c("number","TEs", "group")


#TXX#
#TXX intronic
if (exists("TXXi"))
  { 
TXXinb <- nrow(TXXi)
TXXidata <- data.frame(TXXinb, TXXrep, intronicrep)  
  } else {TXXidata <- data.frame(empty, TXXrep, intronicrep)
  }
colnames(TXXidata) <- c("number","TEs", "group")

#TXX exonic
if (exists("TXXe"))
  { 
TXXenb <- nrow(TXXe)
TXXedata <- data.frame(TXXenb, TXXrep, exonicrep)
  } else {TXXedata <- data.frame(empty, TXXrep, exonicrep)
  }
colnames(TXXedata) <- c("number","TEs", "group")


#TXX ambigous
if (exists("TXXa"))
  { 
TXXanb <- nrow(TXXa)
TXXadata <- data.frame(TXXanb, TXXrep, ambigousrep)
  } else {TXXadata <- data.frame(empty, TXXrep, ambigousrep)
  }
colnames(TXXadata) <- c("number","TEs", "group")


#TOTAL
#TOTAL intronic
#if (exists("TOTALi"))
#  { 
#TOTALinb <- nrow(TOTALi)
#TOTALidata <- data.frame(TOTALinb, TOTALrep, intronicrep)  
#  } else {TOTALidata <- data.frame(empty, TOTALrep, intronicrep)
#  }
#colnames(TOTALidata) <- c("number","TEs", "group")

#TOTAL exonic
#if (exists("TOTALe"))
#  { 
#TOTALenb <- nrow(TOTALe)
#TOTALedata <- data.frame(TOTALenb, TOTALrep, exonicrep)
#  } else {TOTALedata <- data.frame(empty, TOTALrep, exonicrep)
#  }
#colnames(TOTALedata) <- c("number","TEs", "group")


#TOTAL ambigous
#if (exists("TOTALa"))
#  { 
#TOTALanb <- nrow(TOTALa)
#TOTALadata <- data.frame(TOTALanb, TOTALrep, ambigousrep)
#  } else {TOTALadata <- data.frame(empty, TOTALrep, ambigousrep)
#  }
#colnames(TOTALadata) <- c("number","TEs", "group")


#Merge
MergeTEs <- rbind(TXXadata, DHXadata, DTHadata, DTCadata, DTAadata, DTMadata, DTTadata, RSXadata, RXXadata, RIXadata, RYXadata, Caulimoviridaeadata, RLXadata, RLGadata, RLCadata, TXXedata, DHXedata, DTHedata, DTCedata, DTAedata, DTMedata, DTTedata, RSXedata, RXXedata,  RIXedata,  RYXedata, Caulimoviridaeedata, RLXedata, RLGedata, RLCedata, TXXidata, DHXidata, DTHidata, DTCidata, DTAidata, DTMidata,  DTTidata, RSXidata, RXXidata, RIXidata, RYXidata, Caulimoviridaeidata, RLXidata, RLGidata, RLCidata )

#Export txt file
Export <- MergeTEs[c("group", "TEs", "number")]
write.table(Export, "../ParasiTE_output/Results/Plots/Intragenic_TEs_3.txt", sep="\t", row.names = F)

blank_theme <- theme_minimal()+
  theme(
  axis.title.x = element_text(size=15),
  axis.title.y = element_blank(),
  panel.border = element_blank(),
  panel.grid=element_blank(),
  axis.ticks = element_blank(),
  plot.title=element_text(size=15),
  legend.text= element_text(size=15), 
  legend.title = element_blank(), 
  axis.text.x=element_blank(),
  axis.text.y=element_blank()
  )

TEsplot  <- ggplot(data=MergeTEs, aes(x=TEs, y=number, fill=group)) + geom_bar(stat="identity")  + coord_flip() + theme_grey(base_size=20) + scale_color_manual(values=c( "blue","green", "red")) + ggtitle("Proportion of intragenic TE families")+ ylab("Number of TEs") + theme(legend.position="none", axis.title.x =element_text(size=15), axis.title.y = element_blank(), plot.title=element_text(size=15))
 
  
plot_intragenic <- ggplot(tot_intragenic, aes(x= class, y = number, fill= group)) + geom_col() + geom_text(aes(label = number), position = position_stack(vjust = 0.5), size=6 ) 
pie <- plot_intragenic+ coord_polar("y", start=0) + blank_theme + ggtitle("Proportion of intragenic TEs") + ylab("Number of TEs")
Global <- plot_grid(pie, TEsplot, align= "hv", rel_widths= c(2,3))
save_plot("../ParasiTE_output/Results/Plots/Intragenic_TEs_1.pdf", Global, ncol=2)
dev.off()


