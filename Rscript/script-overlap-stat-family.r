#!/usr/bin/env Rscript

#load libraries
library(ggplot2)
library(cowplot)
library(scales)

#load data

#downstream
#intergenic

#!/usr/bin/env Rscript

#load libraries
library(ggplot2)
library(gridExtra)


#Overlap#######################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################


Overlap_upsteam <- read.delim("../ParasiTE_output/Results/Annotations_TEs/Overlap_TEs_upstream.bed", header = FALSE, sep = "\t", dec = ".")
Overlap_downstream <- read.delim("../ParasiTE_output/Results/Annotations_TEs/Overlap_TEs_downstream.bed", header = FALSE, sep = "\t", dec = ".")
Full_Overlap <- read.delim("../ParasiTE_output/Results/Annotations_TEs/Overlap_TEs_full.bed", header = FALSE, sep = "\t", dec = ".")


#create identifiers overlap

Overlaprep <- rep("Overlap", length(1))


Overlap_upsteamnb <- nrow(Overlap_upsteam)
Overlap_upsteamrep <- rep("Upstream", length(1))
Overlap_upsteam_ident <- data.frame(Overlaprep, Overlap_upsteamrep, Overlap_upsteamnb )
colnames(Overlap_upsteam_ident) <- c("class", "group","number" )

Overlap_downsteamnb <- nrow(Overlap_downstream)
Overlap_downsteamrep <- rep("Downstream", length(1))
Overlap_downsteam_ident <- data.frame(Overlaprep, Overlap_downsteamrep, Overlap_downsteamnb )
colnames(Overlap_downsteam_ident) <- c("class", "group","number" )

Full_Overlapnb <- nrow(Full_Overlap)
Full_Overlaprep <- rep("Ambigous", length(1))
Full_Overlap_ident <- data.frame(Overlaprep, Full_Overlaprep, Full_Overlapnb )
colnames(Full_Overlap_ident) <- c("class", "group","number" )


tot_Overlap <- rbind(Full_Overlap_ident, Overlap_upsteam_ident, Overlap_downsteam_ident)

#check if data exist and load them
#downstream
if (file.exists("../ParasiTE_output/overlap_TEs-length/RLC-overlap-TE-downstream.bed"))
RLCi <- read.delim("../ParasiTE_output/overlap_TEs-length/RLC-overlap-TE-downstream.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/overlap_TEs-length/RLG-overlap-TE-downstream.bed"))
RLGi <- read.delim("../ParasiTE_output/overlap_TEs-length/RLG-overlap-TE-downstream.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/overlap_TEs-length/Caulimoviridae-overlap-TE-downstream.bed"))
Caulimoviridaei <- read.delim("../ParasiTE_output/overlap_TEs-length/Caulimoviridae-overlap-TE-downstream.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/overlap_TEs-length/RYX-overlap-TE-downstream.bed"))
RYXi <- read.delim("../ParasiTE_output/overlap_TEs-length/RYX-overlap-TE-downstream.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/overlap_TEs-length/RIX-overlap-TE-downstream.bed"))
RIXi <- read.delim("../ParasiTE_output/overlap_TEs-length/RIX-overlap-TE-downstream.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/overlap_TEs-length/RLX-overlap-TE-downstream.bed"))
RLXi <- read.delim("../ParasiTE_output/overlap_TEs-length/RLX-overlap-TE-downstream.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/overlap_TEs-length/RXX-overlap-TE-downstream.bed"))
RXXi <- read.delim("../ParasiTE_output/overlap_TEs-length/RXX-overlap-TE-downstream.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/overlap_TEs-length/RSX-overlap-TE-downstream.bed"))
RSXi <- read.delim("../ParasiTE_output/overlap_TEs-length/RSX-overlap-TE-downstream.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/overlap_TEs-length/DTT-overlap-TE-downstream.bed"))
DTTi <- read.delim("../ParasiTE_output/overlap_TEs-length/DTT-overlap-TE-downstream.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/overlap_TEs-length/DTM-overlap-TE-downstream.bed"))
DTMi <- read.delim("../ParasiTE_output/overlap_TEs-length/DTM-overlap-TE-downstream.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/overlap_TEs-length/DTA-overlap-TE-downstream.bed"))
DTAi <- read.delim("../ParasiTE_output/overlap_TEs-length/DTA-overlap-TE-downstream.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/overlap_TEs-length/DTC-overlap-TE-downstream.bed"))
DTCi <- read.delim("../ParasiTE_output/overlap_TEs-length/DTC-overlap-TE-downstream.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/overlap_TEs-length/DTH-overlap-TE-downstream.bed"))
DTHi <- read.delim("../ParasiTE_output/overlap_TEs-length/DTH-overlap-TE-downstream.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/overlap_TEs-length/DTX-overlap-TE-downstream.bed"))
DTXi <- read.delim("../ParasiTE_output/overlap_TEs-length/DTX-overlap-TE-downstream.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/overlap_TEs-length/DHX-overlap-TE-downstream.bed"))
DHXi <- read.delim("../ParasiTE_output/overlap_TEs-length/DHX-overlap-TE-downstream.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/overlap_TEs-length/TXX-overlap-TE-downstream.bed"))
TXXi <- read.delim("../ParasiTE_output/overlap_TEs-length/TXX-overlap-TE-downstream.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/overlap_TEs-length/TOTAL-overlap-TE-downstream.bed"))
TOTALi <- read.delim("../ParasiTE_output/overlap_TEs-length/TOTAL-overlap-TE-downstream.bed", header = FALSE, sep = "\t", dec = ".")

#upstream
if (file.exists("../ParasiTE_output/overlap_TEs-length/RLC-overlap-TE-upstream.bed"))
RLCe <- read.delim("../ParasiTE_output/overlap_TEs-length/RLC-overlap-TE-upstream.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/overlap_TEs-length/RLG-overlap-TE-upstream.bed"))
RLGe <- read.delim("../ParasiTE_output/overlap_TEs-length/RLG-overlap-TE-upstream.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/overlap_TEs-length/Caulimoviridae-overlap-TE-upstream.bed"))
Caulimoviridaee <- read.delim("../ParasiTE_output/overlap_TEs-length/Caulimoviridae-overlap-TE-upstream.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/overlap_TEs-length/RYX-overlap-TE-upstream.bed"))
RYXe <- read.delim("../ParasiTE_output/overlap_TEs-length/RYX-overlap-TE-upstream.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/overlap_TEs-length/RIX-overlap-TE-upstream.bed"))
RIXe <- read.delim("../ParasiTE_output/overlap_TEs-length/RIX-overlap-TE-upstream.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/overlap_TEs-length/RLX-overlap-TE-upstream.bed"))
RLXe <- read.delim("../ParasiTE_output/overlap_TEs-length/RLX-overlap-TE-upstream.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/overlap_TEs-length/RXX-overlap-TE-upstream.bed"))
RXXe <- read.delim("../ParasiTE_output/overlap_TEs-length/RXX-overlap-TE-upstream.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/overlap_TEs-length/RSX-overlap-TE-upstream.bed"))
RSXe <- read.delim("../ParasiTE_output/overlap_TEs-length/RSX-overlap-TE-upstream.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/overlap_TEs-length/DTT-overlap-TE-upstream.bed"))
DTTe <- read.delim("../ParasiTE_output/overlap_TEs-length/DTT-overlap-TE-upstream.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/overlap_TEs-length/DTM-overlap-TE-upstream.bed"))
DTMe <- read.delim("../ParasiTE_output/overlap_TEs-length/DTM-overlap-TE-upstream.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/overlap_TEs-length/DTA-overlap-TE-upstream.bed"))
DTAe <- read.delim("../ParasiTE_output/overlap_TEs-length/DTA-overlap-TE-upstream.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/overlap_TEs-length/DTC-overlap-TE-upstream.bed"))
DTCe <- read.delim("../ParasiTE_output/overlap_TEs-length/DTC-overlap-TE-upstream.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/overlap_TEs-length/DTH-overlap-TE-upstream.bed"))
DTHe <- read.delim("../ParasiTE_output/overlap_TEs-length/DTH-overlap-TE-upstream.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/overlap_TEs-length/DTX-overlap-TE-upstream.bed"))
DTXe <- read.delim("../ParasiTE_output/overlap_TEs-length/DTX-overlap-TE-upstream.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/overlap_TEs-length/DHX-overlap-TE-upstream.bed"))
DHXe <- read.delim("../ParasiTE_output/overlap_TEs-length/DHX-overlap-TE-upstream.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/overlap_TEs-length/TXX-overlap-TE-upstream.bed"))
TXXe <- read.delim("../ParasiTE_output/overlap_TEs-length/TXX-overlap-TE-upstream.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/overlap_TEs-length/TOTAL-overlap-TE-upstream.bed"))
TOTALe <- read.delim("../ParasiTE_output/overlap_TEs-length/TOTAL-overlap-TE-upstream.bed", header = FALSE, sep = "\t", dec = ".")

#Full
if (file.exists("../ParasiTE_output/overlap_TEs-length/RLC-overlap-TE-full.bed"))
RLCa <- read.delim("../ParasiTE_output/overlap_TEs-length/RLC-overlap-TE-full.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/overlap_TEs-length/RLG-overlap-TE-full.bed"))
RLGa <- read.delim("../ParasiTE_output/overlap_TEs-length/RLG-overlap-TE-full.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/overlap_TEs-length/Caulimoviridae-overlap-TE-full.bed"))
Caulimoviridaea <- read.delim("../ParasiTE_output/overlap_TEs-length/Caulimoviridae-overlap-TE-full.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/overlap_TEs-length/RYX-overlap-TE-full.bed"))
RYXa <- read.delim("../ParasiTE_output/overlap_TEs-length/RYX-overlap-TE-full.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/overlap_TEs-length/RIX-overlap-TE-full.bed"))
RIXa <- read.delim("../ParasiTE_output/overlap_TEs-length/RIX-overlap-TE-full.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/overlap_TEs-length/RLX-overlap-TE-full.bed"))
RLXa <- read.delim("../ParasiTE_output/overlap_TEs-length/RLX-overlap-TE-full.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/overlap_TEs-length/RXX-overlap-TE-full.bed"))
RXXa <- read.delim("../ParasiTE_output/overlap_TEs-length/RXX-overlap-TE-full.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/overlap_TEs-length/RSX-overlap-TE-full.bed"))
RSXa <- read.delim("../ParasiTE_output/overlap_TEs-length/RSX-overlap-TE-full.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/overlap_TEs-length/DTT-overlap-TE-full.bed"))
DTTa <- read.delim("../ParasiTE_output/overlap_TEs-length/DTT-overlap-TE-full.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/overlap_TEs-length/DTM-overlap-TE-full.bed"))
DTMa <- read.delim("../ParasiTE_output/overlap_TEs-length/DTM-overlap-TE-full.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/overlap_TEs-length/DTA-overlap-TE-full.bed"))
DTAa <- read.delim("../ParasiTE_output/overlap_TEs-length/DTA-overlap-TE-full.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/overlap_TEs-length/DTC-overlap-TE-full.bed"))
DTCa <- read.delim("../ParasiTE_output/overlap_TEs-length/DTC-overlap-TE-full.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/overlap_TEs-length/DTH-overlap-TE-full.bed"))
DTHa <- read.delim("../ParasiTE_output/overlap_TEs-length/DTH-overlap-TE-full.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/overlap_TEs-length/DTX-overlap-TE-full.bed"))
DTXa <- read.delim("../ParasiTE_output/overlap_TEs-length/DTX-overlap-TE-full.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/overlap_TEs-length/DHX-overlap-TE-full.bed"))
DHXa <- read.delim("../ParasiTE_output/overlap_TEs-length/DHX-overlap-TE-full.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/overlap_TEs-length/TXX-overlap-TE-full.bed"))
TXXa <- read.delim("../ParasiTE_output/overlap_TEs-length/TXX-overlap-TE-full.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/overlap_TEs-length/TOTAL-overlap-TE-full.bed"))
TOTALa <- read.delim("../ParasiTE_output/overlap_TEs-length/TOTAL-overlap-TE-full.bed", header = FALSE, sep = "\t", dec = ".")

#create identifiers
downrep <- rep("Downstream", 1)
uprep <- rep("Upstream", 1)
fullrep <- rep("Full", 1)

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


#RLC

#RLCinb <- nrow(RLCi)
#RLCidata <- data.frame(RLCinb, RLCrep, downrep)
#colnames(RLCidata) <- c("number","TEs", "group")

#RLC#
#RLC downstream
if (exists("RLCi"))
  { 
RLCinb <- nrow(RLCi)
RLCidata <- data.frame(RLCinb, RLCrep, downrep)  
  } else {RLCidata <- data.frame(empty, RLCrep, downrep)
  }
colnames(RLCidata) <- c("number","TEs", "group")

#RLC upstream
if (exists("RLCe"))
  { 
RLCenb <- nrow(RLCe)
RLCedata <- data.frame(RLCenb, RLCrep, uprep)
  } else {RLCedata <- data.frame(empty, RLCrep, uprep)
  }
colnames(RLCedata) <- c("number","TEs", "group")


#RLC full
if (exists("RLCa"))
  { 
RLCanb <- nrow(RLCa)
RLCadata <- data.frame(RLCanb, RLCrep, fullrep)
  } else {RLCadata <- data.frame(empty, RLCrep, fullrep)
  }
colnames(RLCadata) <- c("number","TEs", "group")

#RLG
#RLG downstream
if (exists("RLGi"))
  { 
RLGinb <- nrow(RLGi)
RLGidata <- data.frame(RLGinb, RLGrep, downrep)
  } else {RLGidata <- data.frame(empty, RLGrep, downrep)
  }
colnames(RLGidata) <- c("number","TEs", "group")

#RLG upstream
if (exists("RLGe"))
  { 
RLGenb <- nrow(RLGe)
RLGedata <- data.frame(RLGenb, RLGrep, uprep)
  } else {RLGedata <- data.frame(empty, RLGrep, uprep)
  }
colnames(RLGedata) <- c("number","TEs", "group")


#RLG full
if (exists("RLGa"))
  { 
RLGanb <- nrow(RLGa)
RLGadata <- data.frame(RLGanb, RLGrep, fullrep)
  } else {RLGadata <- data.frame(empty, RLGrep, fullrep)
  }
colnames(RLGadata) <- c("number","TEs", "group")

#Caulimoviridae#
#Caulimoviridae downstream

if (exists("Caulimoviridaei"))
  { 
Caulimoviridaeinb <- nrow(Caulimoviridaei)
Caulimoviridaeidata <- data.frame(Caulimoviridaeinb, Caulimoviridaerep, downrep)
  } else {Caulimoviridaeidata <- data.frame(empty, Caulimoviridaerep, downrep)
  }
colnames(Caulimoviridaeidata) <- c("number","TEs", "group")

#Caulimoviridae upstream

if (exists("Caulimoviridaee"))
  { 
Caulimoviridaeenb <- nrow(Caulimoviridaee)
Caulimoviridaeedata <- data.frame(Caulimoviridaeenb, Caulimoviridaerep, uprep)
  } else {Caulimoviridaeedata <- data.frame(empty, Caulimoviridaerep, uprep)
  }
colnames(Caulimoviridaeedata) <- c("number","TEs", "group")

#Caulimoviridae full

if (exists("Caulimoviridaea"))
  { 
Caulimoviridaeanb <- nrow(Caulimoviridaea)
Caulimoviridaeadata <- data.frame(Caulimoviridaeanb, Caulimoviridaerep, fullrep)
  } else {Caulimoviridaeadata <- data.frame(empty, Caulimoviridaerep, fullrep)
  }
colnames(Caulimoviridaeadata) <- c("number","TEs", "group")


#RYX#
#RYX downstream
if (exists("RYXi"))
  { 
RYXinb <- nrow(RYXi)
RYXidata <- data.frame(RYXinb, RYXrep, downrep)  
  } else {RYXidata <- data.frame(empty, RYXrep, downrep)
  }
colnames(RYXidata) <- c("number","TEs", "group")

#RYX upstream
if (exists("RYXe"))
  { 
RYXenb <- nrow(RYXe)
RYXedata <- data.frame(RYXenb, RYXrep, uprep)
  } else {RYXedata <- data.frame(empty, RYXrep, uprep)
  }
colnames(RYXedata) <- c("number","TEs", "group")


#RYX full
if (exists("RYXa"))
  { 
RYXanb <- nrow(RYXa)
RYXadata <- data.frame(RYXanb, RYXrep, fullrep)
  } else {RYXadata <- data.frame(empty, RYXrep, fullrep)
  }
colnames(RYXadata) <- c("number","TEs", "group")


#RIX#
#RIX downstream
if (exists("RIXi"))
  { 
RIXinb <- nrow(RIXi)
RIXidata <- data.frame(RIXinb, RIXrep, downrep)  
  } else {RIXidata <- data.frame(empty, RIXrep, downrep)
  }
colnames(RIXidata) <- c("number","TEs", "group")

#RIX upstream
if (exists("RIXe"))
  { 
RIXenb <- nrow(RIXe)
RIXedata <- data.frame(RIXenb, RIXrep, uprep)
  } else {RIXedata <- data.frame(empty, RIXrep, uprep)
  }
colnames(RIXedata) <- c("number","TEs", "group")


#RIX full
if (exists("RIXa"))
  { 
RIXanb <- nrow(RIXa)
RIXadata <- data.frame(RIXanb, RIXrep, fullrep)
  } else {RIXadata <- data.frame(empty, RIXrep, fullrep)
  }
colnames(RIXadata) <- c("number","TEs", "group")

#RLX#
#RLX downstream
if (exists("RLXi"))
  { 
RLXinb <- nrow(RLXi)
RLXidata <- data.frame(RLXinb, RLXrep, downrep)  
  } else {RLXidata <- data.frame(empty, RLXrep, downrep)
  }
colnames(RLXidata) <- c("number","TEs", "group")

#RLX upstream
if (exists("RLXe"))
  { 
RLXenb <- nrow(RLXe)
RLXedata <- data.frame(RLXenb, RLXrep, uprep)
  } else {RLXedata <- data.frame(empty, RLXrep, uprep)
  }
colnames(RLXedata) <- c("number","TEs", "group")


#RLX full
if (exists("RLXa"))
  { 
RLXanb <- nrow(RLXa)
RLXadata <- data.frame(RLXanb, RLXrep, fullrep)
  } else {RLXadata <- data.frame(empty, RLXrep, fullrep)
  }
colnames(RLXadata) <- c("number","TEs", "group")


#RXX#
#RXX downstream
if (exists("RXXi"))
  { 
RXXinb <- nrow(RXXi)
RXXidata <- data.frame(RXXinb, RXXrep, downrep)  
  } else {RXXidata <- data.frame(empty, RXXrep, downrep)
  }
colnames(RXXidata) <- c("number","TEs", "group")

#RXX upstream
if (exists("RXXe"))
  { 
RXXenb <- nrow(RXXe)
RXXedata <- data.frame(RXXenb, RXXrep, uprep)
  } else {RXXedata <- data.frame(empty, RXXrep, uprep)
  }
colnames(RXXedata) <- c("number","TEs", "group")


#RXX full
if (exists("RXXa"))
  { 
RXXanb <- nrow(RXXa)
RXXadata <- data.frame(RXXanb, RXXrep, fullrep)
  } else {RXXadata <- data.frame(empty, RXXrep, fullrep)
  }
colnames(RXXadata) <- c("number","TEs", "group")
  
#RSX#
#RSX downstream
if (exists("RSXi"))
  { 
RSXinb <- nrow(RSXi)
RSXidata <- data.frame(RSXinb, RSXrep, downrep)  
  } else {RSXidata <- data.frame(empty, RSXrep, downrep)
  }
colnames(RSXidata) <- c("number","TEs", "group")

#RSX upstream
if (exists("RSXe"))
  { 
RSXenb <- nrow(RSXe)
RSXedata <- data.frame(RSXenb, RSXrep, uprep)
  } else {RSXedata <- data.frame(empty, RSXrep, uprep)
  }
colnames(RSXedata) <- c("number","TEs", "group")


#RSX full
if (exists("RSXa"))
  { 
RSXanb <- nrow(RSXa)
RSXadata <- data.frame(RSXanb, RSXrep, fullrep)
  } else {RSXadata <- data.frame(empty, RSXrep, fullrep)
  }
colnames(RSXadata) <- c("number","TEs", "group")


#DTT#
#DTT downstream
if (exists("DTTi"))
  { 
DTTinb <- nrow(DTTi)
DTTidata <- data.frame(DTTinb, DTTrep, downrep)  
  } else {DTTidata <- data.frame(empty, DTTrep, downrep)
  }
colnames(DTTidata) <- c("number","TEs", "group")

#DTT upstream
if (exists("DTTe"))
  { 
DTTenb <- nrow(DTTe)
DTTedata <- data.frame(DTTenb, DTTrep, uprep)
  } else {DTTedata <- data.frame(empty, DTTrep, uprep)
  }
colnames(DTTedata) <- c("number","TEs", "group")


#DTT full
if (exists("DTTa"))
  { 
DTTanb <- nrow(DTTa)
DTTadata <- data.frame(DTTanb, DTTrep, fullrep)
  } else {DTTadata <- data.frame(empty, DTTrep, fullrep)
  }
colnames(DTTadata) <- c("number","TEs", "group")


#DTM#
#DTM downstream
if (exists("DTMi"))
  { 
DTMinb <- nrow(DTMi)
DTMidata <- data.frame(DTMinb, DTMrep, downrep)  
  } else {DTMidata <- data.frame(empty, DTMrep, downrep)
  }
colnames(DTMidata) <- c("number","TEs", "group")

#DTM upstream
if (exists("DTMe"))
  { 
DTMenb <- nrow(DTMe)
DTMedata <- data.frame(DTMenb, DTMrep, uprep)
  } else {DTMedata <- data.frame(empty, DTMrep, uprep)
  }
colnames(DTMedata) <- c("number","TEs", "group")


#DTM full
if (exists("DTMa"))
  { 
DTManb <- nrow(DTMa)
DTMadata <- data.frame(DTManb, DTMrep, fullrep)
  } else {DTMadata <- data.frame(empty, DTMrep, fullrep)
  }
colnames(DTMadata) <- c("number","TEs", "group")

#DTA downstream
if (exists("DTAi"))
  { 
DTAinb <- nrow(DTAi)
DTAidata <- data.frame(DTAinb, DTArep, downrep)  
  } else {DTAidata <- data.frame(empty, DTArep, downrep)
  }
colnames(DTAidata) <- c("number","TEs", "group")

#DTA upstream
if (exists("DTAe"))
  { 
DTAenb <- nrow(DTAe)
DTAedata <- data.frame(DTAenb, DTArep, uprep)
  } else {DTAedata <- data.frame(empty, DTArep, uprep)
  }
colnames(DTAedata) <- c("number","TEs", "group")


#DTA full
if (exists("DTAa"))
  { 
DTAanb <- nrow(DTAa)
DTAadata <- data.frame(DTAanb, DTArep, fullrep)
  } else {DTAadata <- data.frame(empty, DTArep, fullrep)
  }
colnames(DTAadata) <- c("number","TEs", "group")

#DTC
#DTC#
#DTC downstream
if (exists("DTCi"))
  { 
DTCinb <- nrow(DTCi)
DTCidata <- data.frame(DTCinb, DTCrep, downrep)  
  } else {DTCidata <- data.frame(empty, DTCrep, downrep)
  }
colnames(DTCidata) <- c("number","TEs", "group")

#DTC upstream
if (exists("DTCe"))
  { 
DTCenb <- nrow(DTCe)
DTCedata <- data.frame(DTCenb, DTCrep, uprep)
  } else {DTCedata <- data.frame(empty, DTCrep, uprep)
  }
colnames(DTCedata) <- c("number","TEs", "group")


#DTC full
if (exists("DTCa"))
  { 
DTCanb <- nrow(DTCa)
DTCadata <- data.frame(DTCanb, DTCrep, fullrep)
  } else {DTCadata <- data.frame(empty, DTCrep, fullrep)
  }
colnames(DTCadata) <- c("number","TEs", "group")


#DTH
#DTH downstream
if (exists("DTHi"))
  { 
DTHinb <- nrow(DTHi)
DTHidata <- data.frame(DTHinb, DTHrep, downrep)  
  } else {DTHidata <- data.frame(empty, DTHrep, downrep)
  }
colnames(DTHidata) <- c("number","TEs", "group")

#DTH upstream
if (exists("DTHe"))
  { 
DTHenb <- nrow(DTHe)
DTHedata <- data.frame(DTHenb, DTHrep, uprep)
  } else {DTHedata <- data.frame(empty, DTHrep, uprep)
  }
colnames(DTHedata) <- c("number","TEs", "group")


#DTH full
if (exists("DTHa"))
  { 
DTHanb <- nrow(DTHa)
DTHadata <- data.frame(DTHanb, DTHrep, fullrep)
  } else {DTHadata <- data.frame(empty, DTHrep, fullrep)
  }
colnames(DTHadata) <- c("number","TEs", "group")


#DTX#
#DTX downstream
if (exists("DTXi"))
  { 
DTXinb <- nrow(DTXi)
DTXidata <- data.frame(DTXinb, DTXrep, downrep)  
  } else {DTXidata <- data.frame(empty, DTXrep, downrep)
  }
colnames(DTXidata) <- c("number","TEs", "group")

#DTX upstream
if (exists("DTXe"))
  { 
DTXenb <- nrow(DTXe)
DTXedata <- data.frame(DTXenb, DTXrep, uprep)
  } else {DTXedata <- data.frame(empty, DTXrep, uprep)
  }
colnames(DTXedata) <- c("number","TEs", "group")


#DTX full
if (exists("DTXa"))
  { 
DTXanb <- nrow(DTXa)
DTXadata <- data.frame(DTXanb, DTXrep, fullrep)
  } else {DTXadata <- data.frame(empty, DTXrep, fullrep)
  }
colnames(DTXadata) <- c("number","TEs", "group")

#DHX#
#DHX downstream
if (exists("DHXi"))
  { 
DHXinb <- nrow(DHXi)
DHXidata <- data.frame(DHXinb, DHXrep, downrep)  
  } else {DHXidata <- data.frame(empty, DHXrep, downrep)
  }
colnames(DHXidata) <- c("number","TEs", "group")

#DHX upstream
if (exists("DHXe"))
  { 
DHXenb <- nrow(DHXe)
DHXedata <- data.frame(DHXenb, DHXrep, uprep)
  } else {DHXedata <- data.frame(empty, DHXrep, uprep)
  }
colnames(DHXedata) <- c("number","TEs", "group")


#DHX full
if (exists("DHXa"))
  { 
DHXanb <- nrow(DHXa)
DHXadata <- data.frame(DHXanb, DHXrep, fullrep)
  } else {DHXadata <- data.frame(empty, DHXrep, fullrep)
  }
colnames(DHXadata) <- c("number","TEs", "group")

#TXX#
#TXX downstream
if (exists("TXXi"))
  { 
TXXinb <- nrow(TXXi)
TXXidata <- data.frame(TXXinb, TXXrep, downrep)  
  } else {TXXidata <- data.frame(empty, TXXrep, downrep)
  }
colnames(TXXidata) <- c("number","TEs", "group")

#TXX upstream
if (exists("TXXe"))
  { 
TXXenb <- nrow(TXXe)
TXXedata <- data.frame(TXXenb, TXXrep, uprep)
  } else {TXXedata <- data.frame(empty, TXXrep, uprep)
  }
colnames(TXXedata) <- c("number","TEs", "group")


#TXX full
if (exists("TXXa"))
  { 
TXXanb <- nrow(TXXa)
TXXadata <- data.frame(TXXanb, TXXrep, fullrep)
  } else {TXXadata <- data.frame(empty, TXXrep, fullrep)
  }
colnames(TXXadata) <- c("number","TEs", "group")


#TOTAL
#TOTAL downstream
if (exists("TOTALi"))
  { 
TOTALinb <- nrow(TOTALi)
TOTALidata <- data.frame(TOTALinb, TOTALrep, downrep)  
  } else {TOTALidata <- data.frame(empty, TOTALrep, downrep)
  }
colnames(TOTALidata) <- c("number","TEs", "group")

#TOTAL upstream
if (exists("TOTALe"))
  { 
TOTALenb <- nrow(TOTALe)
TOTALedata <- data.frame(TOTALenb, TOTALrep, uprep)
  } else {TOTALedata <- data.frame(empty, TOTALrep, uprep)
  }
colnames(TOTALedata) <- c("number","TEs", "group")


#TOTAL full
if (exists("TOTALa"))
  { 
TOTALanb <- nrow(TOTALa)
TOTALadata <- data.frame(TOTALanb, TOTALrep, fullrep)
  } else {TOTALadata <- data.frame(empty, TOTALrep, fullrep)
  }
colnames(TOTALadata) <- c("number","TEs", "group")

#Merge 
MergeTEs <- rbind(TXXadata, DHXadata, DTHadata, DTCadata, DTAadata,  DTMadata, DTTadata, RSXadata, RXXadata, RIXadata, RYXadata, Caulimoviridaeadata, RLXadata, RLGadata, RLCadata, TXXedata, DHXedata, DTHedata, DTCedata, DTAedata, DTMedata, DTTedata, RSXedata, RXXedata, RIXedata, RYXedata, Caulimoviridaeedata, RLXedata, RLGedata, RLCedata, TXXidata, DHXidata, DTHidata, DTCidata, DTAidata, DTMidata,  DTTidata,  RSXidata, RXXidata, RIXidata, RYXidata, Caulimoviridaeidata, RLXidata, RLGidata, RLCidata )

#Export txt file
Export <- MergeTEs[c("group", "TEs", "number")]
write.table(Export, "../ParasiTE_output/Results/Plots/overlap_TEs_3.txt", sep="\t", row.names = F)


TEsplot  <- ggplot() + geom_bar(data=MergeTEs, aes(x=TEs, y=number, fill=group), stat="identity") + coord_flip() + theme_grey(base_size=20) + ylab("Number of TEs") + ggtitle("Proportion of overlap TE families") +  theme(legend.position="none", axis.title.y =element_blank(), axis.title.x =element_text(size=15),plot.title=element_text(size=15) )

#grid.arrange(RLCplot, RLGplot, RIXplot, RLXplot, RXXplot, ClassITEsplot, ncol=2)

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

plot_overlap <- ggplot(tot_Overlap, aes(x= class, y = number, fill= group)) + geom_col() + geom_text(aes(label = number), position = position_stack(vjust = 0.5), size=6 ) 
pie <- plot_overlap + coord_polar("y", start=0) + blank_theme + ggtitle("Proportion of overlap TEs")+ ylab("Number of TEs") 

Global <- plot_grid(pie, TEsplot, align= "h", rel_widths= c(2,3))
save_plot("../ParasiTE_output/Results/Plots/overlap_TEs_1.pdf", Global, ncol=2)
dev.off()
