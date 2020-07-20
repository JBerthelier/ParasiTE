#!/usr/bin/env Rscript

#load libraries
library(ggplot2)
library(gridExtra)

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

if (file.exists("../ParasiTE_output/intragenic_TEs/TOTAL-intronic_TEs.bed"))
Totali <- read.delim("../ParasiTE_output/intragenic_TEs/TOTAL-intronic_TEs.bed", header = FALSE, sep = "\t", dec = ".")

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

if (file.exists("../ParasiTE_output/intragenic_TEs/TOTAL-exonic_TEs.bed"))
Totale <- read.delim("../ParasiTE_output/intragenic_TEs/TOTAL-exonic_TEs.bed", header = FALSE, sep = "\t", dec = ".")

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

if (file.exists("../ParasiTE_output/intragenic_TEs/TOTAL-ambigous_TEs.bed"))
Totala <- read.delim("../ParasiTE_output/intragenic_TEs/TOTAL-ambigous_TEs.bed", header = FALSE, sep = "\t", dec = ".")

#create identifiers for intronic
empty <- data.frame(matrix(ncol = 1, nrow = 1))


#RLC#
#RLC intronic
if (exists("RLCi"))
  { 
RLCinbl <- nrow(RLCi)
RLCirepl <- rep("RLC", length(RLCinbl)) 
RLCidatal <- data.frame(RLCi, RLCirepl)
  } else {RLCidatal <- data.frame(empty, rep("RLC", 1))
  }
colnames(RLCidatal) <- c("length","TEs")

#RLC exonic
if (exists("RLCe"))
  {
RLCenbl <- nrow(RLCe) 
RLCerepl <- rep("RLC", length(RLCenbl)) 
RLCedatal <- data.frame(RLCe, RLCerepl)
  } else {RLCedatal <- data.frame(empty, rep("RLC", 1))
  }
colnames(RLCedatal) <- c("length","TEs")

#RLC ambigous
if (exists("RLCa"))
  {
RLCanbl <- nrow(RLCa)  
RLCarepl <- rep("RLC", length(RLCanbl)) 
RLCadatal <- data.frame(RLCa, RLCarepl)
  } else {RLCadatal <- data.frame(empty, rep("RLC", 1))
  }
colnames(RLCadatal) <- c("length","TEs")

#RLG#
#RLG intronic
if (exists("RLGi"))
  {
RLGinbl <- nrow(RLGi)  
RLGirepl <- rep("RLG", length(RLGinbl)) 
RLGidatal <- data.frame(RLGi, RLGirepl)
  } else {RLGidatal <- data.frame(empty, rep("RLG", 1))
  }
colnames(RLGidatal) <- c("length","TEs")

#RLG exonic
if (exists("RLGe"))
  {
RLGenbl <- nrow(RLGe)    
RLGerepl <- rep("RLG", length(RLGenbl)) 
RLGedatal <- data.frame(RLGe, RLGerepl)
  } else {RLGedatal <- data.frame(empty, rep("RLG", 1))
  }
colnames(RLGedatal) <- c("length","TEs")

#RLG ambigous
if (exists("RLGa"))
  { 
RLGanbl <- nrow(RLGa)    
RLGarepl <- rep("RLG", length(RLGanbl)) 
RLGadatal <- data.frame(RLGa, RLGarepl)
  } else {RLGadatal <- data.frame(empty, rep("RLG", 1))
  }
colnames(RLGadatal) <- c("length","TEs")

#Caulimoviridae#
#Caulimoviridae intronic
if (exists("Caulimoviridaei"))
  { 
Caulimoviridaeinbl <- nrow(Caulimoviridaei)      
Caulimoviridaeirepl <- rep("Cau", length(Caulimoviridaeinbl)) 
Caulimoviridaeidatal <- data.frame(Caulimoviridaei, Caulimoviridaeirepl)
  } else {Caulimoviridaeidatal <- data.frame(empty, rep("Cau", 1))
  }
colnames(Caulimoviridaeidatal) <- c("length","TEs")

#Caulimoviridae exonic
if (exists("Caulimoviridaee"))
  { 
Caulimoviridaeenbl <- nrow(Caulimoviridaee)    
Caulimoviridaeerepl <- rep("Cau", length(Caulimoviridaeenbl)) 
Caulimoviridaeedatal <- data.frame(Caulimoviridaee, Caulimoviridaeerepl)
  } else {Caulimoviridaeedatal <- data.frame(empty, rep("Cau", 1))
  }
colnames(Caulimoviridaeedatal) <- c("length","TEs")

#Caulimoviridae ambigous
if (exists("Caulimoviridaea"))
  { 
Caulimoviridaeanbl <- nrow(Caulimoviridaea)  
Caulimoviridaearepl <- rep("Cau", length(Caulimoviridaeanbl)) 
Caulimoviridaeadatal <- data.frame(Caulimoviridaea, Caulimoviridaearepl)
  } else {Caulimoviridaeadatal <- data.frame(empty, rep("Cau", 1))
  }
colnames(Caulimoviridaeadatal) <- c("length","TEs")


#RYX#
#RYX intronic
if (exists("RYXi"))
  { 
RYXinbl <- nrow(RYXi)
RYXirepl <- rep("RYX", length(RYXinbl)) 
RYXidatal <- data.frame(RYXi, RYXirepl)
  } else {RYXidatal <- data.frame(empty, rep("RYX", 1))
  }
colnames(RYXidatal) <- c("length","TEs")

#RYX exonic
if (exists("RYXe"))
  {
RYXenbl <- nrow(RYXe)  
RYXerepl <- rep("RYX", length(RYXenbl)) 
RYXedatal <- data.frame(RYXe, RYXerepl)
  } else {RYXedatal <- data.frame(empty, rep("RYX", 1))
  }
colnames(RYXedatal) <- c("length","TEs")

#RYX ambigous
if (exists("RYXa"))
  { 
RYXanbl <- nrow(RYXa) 
RYXarepl <- rep("RYX", length(RYXanbl)) 
RYXadatal <- data.frame(RYXa, RYXarepl)
  } else {RYXadatal <- data.frame(empty, rep("RYX", 1))
  }
colnames(RYXadatal) <- c("length","TEs")



#RIX#
#RIX intronic
if (exists("RIXi"))
  { 
RIXinbl <- nrow(RIXi)
RIXirepl <- rep("RIX", length(RIXinbl)) 
RIXidatal <- data.frame(RIXi, RIXirepl)
  } else {RIXidatal <- data.frame(empty, rep("RIX", 1))
  }
colnames(RIXidatal) <- c("length","TEs")

#RIX exonic
if (exists("RIXe"))
  {
RIXenbl <- nrow(RIXe)  
RIXerepl <- rep("RIX", length(RIXenbl)) 
RIXedatal <- data.frame(RIXe, RIXerepl)
  } else {RIXedatal <- data.frame(empty, rep("RIX", 1))
  }
colnames(RIXedatal) <- c("length","TEs")

#RIX ambigous
if (exists("RIXa"))
  { 
RIXanbl <- nrow(RIXa) 
RIXarepl <- rep("RIX", length(RIXanbl)) 
RIXadatal <- data.frame(RIXa, RIXarepl)
  } else {RIXadatal <- data.frame(empty, rep("RIX", 1))
  }
colnames(RIXadatal) <- c("length","TEs")


#RLX#
#RLX intronic
if (exists("RLXi"))
  { 
RLXinbl <- nrow(RLXi)   
RLXirepl <- rep("RLX", length(RLXinbl)) 
RLXidatal <- data.frame(RLXi, RLXirepl)
  } else {RLXidatal <- data.frame(empty, rep("RLX", 1))
  }
colnames(RLXidatal) <- c("length","TEs")

#RLX exonic
if (exists("RLXe"))
  { 
RLXenbl <- nrow(RLXe)
RLXerepl <- rep("RLX", length(RLXenbl)) 
RLXedatal <- data.frame(RLXe, RLXerepl)
  } else {RLXedatal <- data.frame(empty, rep("RLX", 1))
  }
colnames(RLXedatal) <- c("length","TEs")

#RLX ambigous
if (exists("RLXa"))
  { 
RLXanbl <- nrow(RLXa)
RLXarepl <- rep("RLX", length(RLXanbl)) 
RLXadatal <- data.frame(RLXa, RLXarepl)
  } else {RLXadatal <- data.frame(empty, rep("RLX", 1))
  }
colnames(RLXadatal) <- c("length","TEs")

#RXX#
#RXX intronic
if (exists("RXXi"))
  { 
RXXinbl <- nrow(RXXi)
RXXirepl <- rep("RXX", length(RXXinbl)) 
RXXidatal <- data.frame(RXXi, RXXirepl)
  } else {RXXidatal <- data.frame(empty, rep("RXX", 1))
  }
colnames(RXXidatal) <- c("length","TEs")

#RXX exonic
if (exists("RXXe"))
  { 
RXXenbl <- nrow(RXXe)
RXXerepl <- rep("RXX", length(RXXenbl)) 
RXXedatal <- data.frame(RXXe, RXXerepl)
  } else {RXXedatal <- data.frame(empty, rep("RXX", 1))
  }
colnames(RXXedatal) <- c("length","TEs")

#RXX ambigous
if (exists("RXXa"))
  { 
RXXanbl <- nrow(RXXa)
RXXarepl <- rep("RXX", length(RXXanbl)) 
RXXadatal <- data.frame(RXXa, RXXarepl)
  } else {RXXadatal <- data.frame(empty, rep("RXX", 1))
  }
colnames(RXXadatal) <- c("length","TEs")

#RSX#
#RSX intronic
if (exists("RSXi"))
  { 
RSXinbl <- nrow(RSXi)
RSXirepl <- rep("RSX", length(RSXinbl)) 
RSXidatal <- data.frame(RSXi, RSXirepl)
  } else {RSXidatal <- data.frame(empty, rep("RSX", 1))
  }
colnames(RSXidatal) <- c("length","TEs")

#RSX exonic
if (exists("RSXe"))
  { 
RSXenbl <- nrow(RSXe)
RSXerepl <- rep("RSX", length(RSXenbl)) 
RSXedatal <- data.frame(RSXe, RSXerepl)
  } else {RSXedatal <- data.frame(empty, rep("RSX", 1))
  }
colnames(RSXedatal) <- c("length","TEs")

#RSX ambigous
if (exists("RSXa"))
  { 
RSXanbl <- nrow(RSXa)
RSXarepl <- rep("RSX", length(RSXanbl)) 
RSXadatal <- data.frame(RSXa, RSXarepl)
  } else {RSXadatal <- data.frame(empty, rep("RSX", 1))
  }
colnames(RSXadatal) <- c("length","TEs")

#DTM#
#DTM intronic
if (exists("DTMi"))
  { 
DTMinbl <- nrow(DTMi)  
DTMirepl <- rep("DTM", length(DTMinbl)) 
DTMidatal <- data.frame(DTMi, DTMirepl)
  } else {DTMidatal <- data.frame(empty, rep("DTM", 1))
  }
colnames(DTMidatal) <- c("length","TEs")

#DTM exonic
if (exists("DTMe"))
  { 
DTMenbl <- nrow(DTMe) 
DTMerepl <- rep("DTM", length(DTMenbl)) 
DTMedatal <- data.frame(DTMe, DTMerepl)
  } else {DTMedatal <- data.frame(empty, rep("DTM", 1))
  }
colnames(DTMedatal) <- c("length","TEs")

#DTM ambigous
if (exists("DTMa"))
  { 
DTManbl <- nrow(DTMa) 
DTMarepl <- rep("DTM", length(DTManbl)) 
DTMadatal <- data.frame(DTMa, DTMarepl)
  } else {DTMadatal <- data.frame(empty, rep("DTM", 1))
  }
colnames(DTMadatal) <- c("length","TEs")

#DTA#
#DTA intronic
if (exists("DTAi"))
  { 
DTAinbl <- nrow(DTAi) 
DTAirepl <- rep("DTA", length(DTAinbl)) 
DTAidatal <- data.frame(DTAi, DTAirepl)
  } else {DTAidatal <- data.frame(empty, rep("DTA", 1))
  }
colnames(DTAidatal) <- c("length","TEs")

#DTA exonic
if (exists("DTAe"))
  { 
DTAenbl <- nrow(DTAe) 
DTAerepl <- rep("DTA", length(DTAenbl)) 
DTAedatal <- data.frame(DTAe, DTAerepl)
  } else {DTAedatal <- data.frame(empty, rep("DTA", 1))
  }
colnames(DTAedatal) <- c("length","TEs")

#DTA ambigous
if (exists("DTAa"))
  { 
DTAanbl <- nrow(DTAa)
DTAarepl <- rep("DTA", length(DTAanbl)) 
DTAadatal <- data.frame(DTAa, DTAarepl)
  } else {DTAadatal <- data.frame(empty, rep("DTA", 1))
  }
colnames(DTAadatal) <- c("length","TEs")


#DTC#
#DTC intronic
if (exists("DTCi"))
  { 
DTCinbl <- nrow(DTCi)
DTCirepl <- rep("DTC", length(DTCinbl)) 
DTCidatal <- data.frame(DTCi, DTCirepl)
  } else {DTCidatal <- data.frame(empty, rep("DTC", 1))
  }
colnames(DTCidatal) <- c("length","TEs")

#DTC exonic
if (exists("DTCe"))
  { 
DTCenbl <- nrow(DTCe)
DTCerepl <- rep("DTC", length(DTCenbl)) 
DTCedatal <- data.frame(DTCe, DTCerepl)
  } else {DTCedatal <- data.frame(empty, rep("DTC", 1))
  }
colnames(DTCedatal) <- c("length","TEs")

#DTC ambigous
if (exists("DTCa"))
  {
DTCanbl <- nrow(DTCa)  
DTCarepl <- rep("DTC", length(DTCanbl)) 
DTCadatal <- data.frame(DTCa, DTCarepl)
  } else {DTCadatal <- data.frame(empty, rep("DTC", 1))
  }
colnames(DTCadatal) <- c("length","TEs")

#DTH#
#DTH intronic
if (exists("DTHi"))
  { 
DTHinbl <- nrow(DTHi)
DTHirepl <- rep("DTH", length(DTHinbl)) 
DTHidatal <- data.frame(DTHi, DTHirepl)
  } else {DTHidatal <- data.frame(empty, rep("DTH", 1))
  }
colnames(DTHidatal) <- c("length","TEs")

#DTH exonic
if (exists("DTHe"))
  { 
DTHenbl <- nrow(DTHe)  
DTHerepl <- rep("DTH", length(DTHenbl)) 
DTHedatal <- data.frame(DTHe, DTHerepl)
  } else {DTHedatal <- data.frame(empty, rep("DTH", 1))
  }
colnames(DTHedatal) <- c("length","TEs")

#DTH ambigous
if (exists("DTHa"))
  { 
DTHanbl <- nrow(DTHa)  
DTHarepl <- rep("DTH", length(DTHanbl)) 
DTHadatal <- data.frame(DTHa, DTHarepl)
  } else {DTHadatal <- data.frame(empty, rep("DTH", 1))
  }
colnames(DTHadatal) <- c("length","TEs")

#DTT#
#DTT intronic
if (exists("DTTi"))
  { 
DTTinbl <- nrow(DTTi)  
DTTirepl <- rep("DTT", length(DTTinbl)) 
DTTidatal <- data.frame(DTTi, DTTirepl)
  } else {DTTidatal <- data.frame(empty, rep("DTT", 1))
  }
colnames(DTTidatal) <- c("length","TEs")

#DTT exonic
if (exists("DTTe"))
  {
DTTenbl <- nrow(DTTe)    
DTTerepl <- rep("DTT", length(DTTenbl)) 
DTTedatal <- data.frame(DTTe, DTTerepl)
  } else {DTTedatal <- data.frame(empty, rep("DTT", 1))
  }
colnames(DTTedatal) <- c("length","TEs")

#DTT ambigous
if (exists("DTTa"))
  {
DTTanbl <- nrow(DTTa)  
DTTarepl <- rep("DTT", length(DTTanbl)) 
DTTadatal <- data.frame(DTTa, DTTarepl)
  } else {DTTadatal <- data.frame(empty, rep("DTT", 1))
  }
colnames(DTTadatal) <- c("length","TEs")

#DTX#
#DTX intronic
if (exists("DTXi"))
  { 
DTXinbl <- nrow(DTXi)  
DTXirepl <- rep("DTX", length(DTXinbl)) 
DTXidatal <- data.frame(DTXi, DTXirepl)
  } else {DTXidatal <- data.frame(empty, rep("DTX", 1))
  }
colnames(DTXidatal) <- c("length","TEs")

#DTX exonic
if (exists("DTXe"))
  { 
DTXenbl <- nrow(DTXe)  
DTXerepl <- rep("DTX", length(DTXenbl)) 
DTXedatal <- data.frame(DTXe, DTXerepl)
  } else {DTXedatal <- data.frame(empty, rep("DTX", 1))
  }
colnames(DTXedatal) <- c("length","TEs")

#DTX ambigous
if (exists("DTXa"))
  { 
DTXanbl <- nrow(DTXa)  
DTXarepl <- rep("DTX", length(DTXanbl)) 
DTXadatal <- data.frame(DTXa, DTXarepl)
  } else {DTXadatal <- data.frame(empty, rep("DTX", 1))
  }
colnames(DTXadatal) <- c("length","TEs")

#DHX#
#DHX intronic
if (exists("DHXi"))
  {
DHXinbl <- nrow(DHXi)    
DHXirepl <- rep("DHX", length(DHXinbl)) 
DHXidatal <- data.frame(DHXi, DHXirepl)
  } else {DHXidatal <- data.frame(empty, rep("DHX", 1))
  }
colnames(DHXidatal) <- c("length","TEs")

#DHX exonic
if (exists("DHXe"))
  {
DHXenbl <- nrow(DHXe)    
DHXerepl <- rep("DHX", length(DHXenbl)) 
DHXedatal <- data.frame(DHXe, DHXerepl)
  } else {DHXedatal <- data.frame(empty, rep("DHX", 1))
  }
colnames(DHXedatal) <- c("length","TEs")

#DHX ambigous
if (exists("DHXa"))
  { 
DHXanbl <- nrow(DHXa)  
DHXarepl <- rep("DHX", length(DHXanbl)) 
DHXadatal <- data.frame(DHXa, DHXarepl)
  } else {DHXadatal <- data.frame(empty, rep("DHX", 1))
  }
colnames(DHXadatal) <- c("length","TEs")

#TXX#
#TXX intronic
if (exists("TXXi"))
  {
TXXinbl <- nrow(TXXi)    
TXXirepl <- rep("TXX", length(TXXinbl)) 
TXXidatal <- data.frame(TXXi, TXXirepl)
  } else {TXXidatal <- data.frame(empty, rep("TXX", 1))
  }
colnames(TXXidatal) <- c("length","TEs")

#TXX exonic
if (exists("TXXe"))
  {
TXXenbl <- nrow(TXXe)    
TXXerepl <- rep("TXX", length(TXXenbl)) 
TXXedatal <- data.frame(TXXe, TXXerepl)
  } else {TXXedatal <- data.frame(empty, rep("TXX", 1))
  }
colnames(TXXedatal) <- c("length","TEs")

#TXX ambigous
if (exists("TXXa"))
  { 
TXXanbl <- nrow(TXXa)  
TXXarepl <- rep("TXX", length(TXXanbl)) 
TXXadatal <- data.frame(TXXa, TXXarepl)
  } else {TXXadatal <- data.frame(empty, rep("TXX", 1))
  }
colnames(TXXadatal) <- c("length","TEs")

#Total#
#Total intronic
#if (exists("Totali"))
#  {
#Totalinbl <- nrow(Totali)    
#Totalirepl <- rep("Totali", length(Totalinbl)) 
#Totalidatal <- data.frame(Totali, Totalirepl)
#  } else {Totalidatal <- data.frame(empty, rep("Totali", 1))
#  }
#colnames(Totalidatal) <- c("length","TEs")
#
#Total exonic
#if (exists("Totale"))
#  {
#Totalenbl <- nrow(Totale)    
#Totalerepl <- rep("Totale", length(Totalenbl)) 
#Totaledatal <- data.frame(Totale, Totalerepl)
#  } else {Totaledatal <- data.frame(empty, rep("Totale", 1))
#  }
#colnames(Totaledatal) <- c("length","TEs")
#
#Total ambigous
#if (exists("Totala"))
#  { 
#Totalanbl <- nrow(Totala)  
#Totalarepl <- rep("Totala", length(Totalanbl)) 
#Totaladatal <- data.frame(Totala, Totalarepl)
#  } else {Totaladatal <- data.frame(empty, rep("Totala", 1))
#  }
#colnames(Totaladatal) <- c("length","TEs")


#Plot intragenic
#AllRLC <- rbind(RLCidatal, RLCadatal, RLCedatal)
#AllRLG <- rbind(RLGidatal, RLGadatal, RLGedatal)
#AllRLX <- rbind(RLXidatal, RLXadatal, RLXedatal)
#AllRYX <- rbind(RYXidatal, RYXadatal, RYXedatal)
#AllRIX <- rbind(RIXidatal, RIXadatal, RIXedatal)
#AllCaulimoviridae <- rbind(Caulimoviridaeidatal, Caulimoviridaeadatal, Caulimoviridaeedatal)
#AllRSX <- rbind(RSXidatal, RSXadatal, RSXedatal)
#AllRXX <- rbind(RXXidatal, RXXadatal, RXXedatal)
#AllDTT <- rbind(DTTidatal, DTTadatal, DTTedatal)
#AllDTM <- rbind(DTMidatal, DTMadatal, DTMedatal)
#AllDTA <- rbind(DTAidatal, DTAadatal, DTAedatal)
#AllDTC <- rbind(DTCidatal, DTCadatal, DTCedatal)
#AllDTH <- rbind(DTHidatal, DTHadatal, DTHedatal)
#AllDTX <- rbind(DTXidatal, DTXadatal, DTXedatal)
#AllDHX <- rbind(DHXidatal, DHXadatal, DHXedatal)
#AllTXX <- rbind(TXXidatal, TXXadatal, TXXedatal)

#AllClassITEs <- rbind(AllRLC, AllRLG, AllCaulimoviridae, AllRYX, AllRIX, AllRLX, AllRSX, AllRXX)
#AllClassIITEs <- rbind(AllDTT, AllDTM, AllDTA, AllDTC, AllDTH, AllDTX, AllDHX, AllTXX)

#AllClassITEsplot  <- ggplot(AllClassITEs, aes(x=length, color=TEs)) + geom_freqpoly() + theme_grey() + ylab("Number of TEs") + xlab("Length") + theme(legend.title = element_blank()) + ggtitle("Intragenic")
#AllClassIITEsplot  <- ggplot(AllClassIITEs, aes(x=length, color=TEs)) + geom_freqpoly() + theme_grey() + ylab("Number of TEs") + xlab("Length") + theme(legend.title = element_blank()) + ggtitle(" ")

#Plot intronic
iClassITEs <- rbind(RLCidatal, RLGidatal, Caulimoviridaeidatal, RYXidatal, RIXidatal, RLXidatal, RSXidatal, RXXidatal)
iClassIITEs <- rbind(DTTidatal, DTMidatal, DTAidatal, DTCidatal, DTHidatal, DTXidatal, DHXidatal, TXXidatal)

iClassITEsplot  <- ggplot(iClassITEs, aes(x=length, color=TEs)) + geom_freqpoly() + theme_grey() + ylab("Number of TEs") + xlab("Length") + theme(legend.title = element_blank()) + ggtitle("Intronic")
iClassIITEsplot  <- ggplot(iClassIITEs, aes(x=length, color=TEs)) + geom_freqpoly() + theme_grey() + ylab("Number of TEs") + xlab("Length") + theme(legend.title = element_blank()) + ggtitle(" ")

#Plot exonic
eClassITEs <- rbind(RLCedatal, RLGedatal, Caulimoviridaeedatal, RYXedatal, RIXedatal, RLXedatal, RSXedatal, RXXedatal)
eClassIITEs <- rbind(DTTedatal, DTMedatal, DTAedatal, DTCedatal, DTHedatal, DTXedatal, DHXedatal, TXXedatal)

eClassITEsplot  <- ggplot(eClassITEs, aes(x=length, color=TEs)) + geom_freqpoly() + theme_grey() + ylab("Number of TEs") + xlab("Length") + theme(legend.title = element_blank()) + ggtitle("Exonic")
eClassIITEsplot  <- ggplot(eClassIITEs, aes(x=length, color=TEs)) + geom_freqpoly() + theme_grey() + ylab("Number of TEs") + xlab("Length") + theme(legend.title = element_blank()) + ggtitle(" ")

#Plot ambigous
aClassITEs <- rbind(RLCadatal, RLGadatal, Caulimoviridaeadatal, RYXadatal, RIXadatal, RLXadatal, RSXadatal, RXXadatal)
aClassIITEs <- rbind(DTTadatal, DTMadatal, DTAadatal, DTCadatal, DTHadatal, DTXadatal, DHXadatal, TXXadatal)

aClassITEsplot  <- ggplot(aClassITEs, aes(x=length, color=TEs)) + geom_freqpoly() + theme_grey() + ylab("Number of TEs") + xlab("Length") + theme(legend.title = element_blank()) + ggtitle("Ambigous")
aClassIITEsplot  <- ggplot(aClassIITEs, aes(x=length, color=TEs)) + geom_freqpoly() + theme_grey() + ylab("Number of TEs") + xlab("Length") + theme(legend.title = element_blank()) + ggtitle(" ")


pdf(file = "../ParasiTE_output/Results/Plots/Intragenic_TEs_2.pdf")
#grid.arrange(AllClassITEsplot, AllClassIITEsplot, iClassITEsplot,iClassIITEsplot, eClassITEsplot, eClassIITEsplot, aClassITEsplot, aClassIITEsplot, ncol=2, top="Frequency of the length of intragenic TEs")
grid.arrange(iClassITEsplot,iClassIITEsplot, eClassITEsplot, eClassIITEsplot, aClassITEsplot, aClassIITEsplot, ncol=2, top="Frequency of the length of intragenic TEs")

dev.off()


