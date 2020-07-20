#!/usr/bin/env Rscript

#load libraries
library(ggplot2)
library(gridExtra)

#load data
if (file.exists("../ParasiTE_output/neighbor_TEs/RLC-neighbor_TEs-downstream.bed"))
RLC <- read.delim("../ParasiTE_output/neighbor_TEs/RLC-neighbor_TEs-downstream.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/neighbor_TEs/RLG-neighbor_TEs-downstream.bed"))
RLG <- read.delim("../ParasiTE_output/neighbor_TEs/RLG-neighbor_TEs-downstream.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/neighbor_TEs/Caulimoviridae-neighbor_TEs-downstream.bed"))
Caulimoviridae <- read.delim("../ParasiTE_output/neighbor_TEs/Caulimoviridae-neighbor_TEs-downstream.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/neighbor_TEs/RIX-neighbor_TEs-downstream.bed"))
RIX <- read.delim("../ParasiTE_output/neighbor_TEs/RIX-neighbor_TEs-downstream.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/neighbor_TEs/RLX-neighbor_TEs-downstream.bed"))
RLX <- read.delim("../ParasiTE_output/neighbor_TEs/RLX-neighbor_TEs-downstream.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/neighbor_TEs/RXX-neighbor_TEs-downstream.bed"))
RXX <- read.delim("../ParasiTE_output/neighbor_TEs/RXX-neighbor_TEs-downstream.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/neighbor_TEs/RSX-neighbor_TEs-downstream.bed"))
RSX <- read.delim("../ParasiTE_output/neighbor_TEs/RSX-neighbor_TEs-downstream.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/neighbor_TEs/DTM-neighbor_TEs-downstream.bed"))
DTM <- read.delim("../ParasiTE_output/neighbor_TEs/DTM-neighbor_TEs-downstream.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/neighbor_TEs/DTA-neighbor_TEs-downstream.bed"))
DTA <- read.delim("../ParasiTE_output/neighbor_TEs/DTA-neighbor_TEs-downstream.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/neighbor_TEs/DTC-neighbor_TEs-downstream.bed"))
DTC <- read.delim("../ParasiTE_output/neighbor_TEs/DTC-neighbor_TEs-downstream.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/neighbor_TEs/DTH-neighbor_TEs-downstream.bed"))
DTH <- read.delim("../ParasiTE_output/neighbor_TEs/DTH-neighbor_TEs-downstream.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/neighbor_TEs/DTT-neighbor_TEs-downstream.bed"))
DTT <- read.delim("../ParasiTE_output/neighbor_TEs/DTT-neighbor_TEs-downstream.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/neighbor_TEs/DTX-neighbor_TEs-downstream.bed"))
DTX <- read.delim("../ParasiTE_output/neighbor_TEs/DTX-neighbor_TEs-downstream.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/neighbor_TEs/DHX-neighbor_TEs-downstream.bed"))
DHX <- read.delim("../ParasiTE_output/neighbor_TEs/DHX-neighbor_TEs-downstream.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/neighbor_TEs/TXX-neighbor_TEs-downstream.bed"))
TXX <- read.delim("../ParasiTE_output/neighbor_TEs/TXX-neighbor_TEs-downstream.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/neighbor_TEs/TOTAL-neighbor_TEs-downstream.bed"))
TOTAL <- read.delim("../ParasiTE_output/neighbor_TEs/TOTAL-neighbor_TEs-downstream.bed", header = FALSE, sep = "\t", dec = ".")

#create dataifiers

empty <- data.frame(matrix(ncol = 1, nrow = 1))

if (exists("RLC"))
  { 
RLCnb <- nrow(RLC)
RLCrep <- rep("RLC", length(RLCnb)) 
RLCdata <- data.frame(RLC, RLCrep)
colnames(RLCdata) <- c("length","TEs")
  } else {RLCdata <- data.frame(empty, rep("RLC", 1))
  }
colnames(RLCdata) <- c("length","TEs")


if (exists("RLG"))
  { 
RLGnb <- nrow(RLG)
RLGrep <- rep("RLG", length(RLGnb)) 
RLGdata <- data.frame(RLG, RLGrep)
colnames(RLGdata) <- c("length","TEs")
  } else {RLGdata <- data.frame(empty, rep("RLG", 1))
  }
colnames(RLGdata) <- c("length","TEs")


if (exists("Cau"))
  { 
Caulimoviridaenb <- nrow(Caulimoviridae)
Caulimoviridaerep <- rep("Cau", length(Caulimoviridaenb)) 
Caulimoviridaedata <- data.frame(Caulimoviridae, Caulimoviridaerep)
colnames(Caulimoviridaedata) <- c("length","TEs")
  } else {Caulimoviridaedata <- data.frame(empty, rep("Cau", 1))
  }
colnames(Caulimoviridaedata) <- c("length","TEs")


if (exists("RIX"))
  { 
RIXnb <- nrow(RIX)
RIXrep <- rep("RIX", length(RIXnb)) 
RIXdata <- data.frame(RIX, RIXrep)
colnames(RIXdata) <- c("length","TEs")
  } else {RIXdata <- data.frame(empty, rep("RIX", 1))
  }
colnames(RIXdata) <- c("length","TEs")


if (exists("RLX"))
  { 
RLXnb <- nrow(RLX)
RLXrep <- rep("RLX", length(RLXnb)) 
RLXdata <- data.frame(RLX, RLXrep)
colnames(RLXdata) <- c("length","TEs")
  } else {RLXdata <- data.frame(empty, rep("RLX", 1))
  }
colnames(RLXdata) <- c("length","TEs")


if (exists("RXX"))
  { 
RXXnb <- nrow(RXX)
RXXrep <- rep("RXX", length(RXXnb)) 
RXXdata <- data.frame(RXX, RXXrep)
colnames(RXXdata) <- c("length","TEs")
  } else {RXXdata <- data.frame(empty, rep("RXX", 1))
  }
colnames(RXXdata) <- c("length","TEs")


if (exists("RSX"))
  { 
RSXnb <- nrow(RSX)
RSXrep <- rep("RSX", length(RSXnb)) 
RSXdata <- data.frame(RSX, RSXrep)
colnames(RSXdata) <- c("length","TEs")
  } else {RSXdata <- data.frame(empty, rep("RSX", 1))
  }
colnames(RSXdata) <- c("length","TEs")


if (exists("DTM"))
  { 
DTMnb <- nrow(DTM)
DTMrep <- rep("DTM", length(DTMnb)) 
DTMdata <- data.frame(DTM, DTMrep)
colnames(DTMdata) <- c("length","TEs")
  } else {DTMdata <- data.frame(empty, rep("DTM", 1))
  }
colnames(DTMdata) <- c("length","TEs")


if (exists("DTA"))
  { 
DTAnb <- nrow(DTA)
DTArep <- rep("DTA", length(DTAnb)) 
DTAdata <- data.frame(DTA, DTArep)
colnames(DTAdata) <- c("length","TEs")
  } else {DTAdata <- data.frame(empty, rep("DTA", 1))
  }
colnames(DTAdata) <- c("length","TEs")


if (exists("DTC"))
  { 
DTCnb <- nrow(DTC)
DTCrep <- rep("DTC", length(DTCnb)) 
DTCdata <- data.frame(DTC, DTCrep)
colnames(DTCdata) <- c("length","TEs")
  } else {DTCdata <- data.frame(empty, rep("DTC", 1))
  }
colnames(DTCdata) <- c("length","TEs")


if (exists("DTH"))
  { 
DTHnb <- nrow(DTH)
DTHrep <- rep("DTH", length(DTHnb)) 
DTHdata <- data.frame(DTH, DTHrep)
colnames(DTHdata) <- c("length","TEs")
  } else {DTHdata <- data.frame(empty, rep("DTH", 1))
  }
colnames(DTHdata) <- c("length","TEs")


if (exists("DTT"))
  { 
DTTnb <- nrow(DTT)
DTTrep <- rep("DTT", length(DTTnb)) 
DTTdata <- data.frame(DTT, DTTrep)
colnames(DTTdata) <- c("length","TEs")
  } else {DTTdata <- data.frame(empty, rep("DTT", 1))
  }
colnames(DTTdata) <- c("length","TEs")


if (exists("DTX"))
  { 
DTXnb <- nrow(DTX)
DTXrep <- rep("DTX", length(DTXnb)) 
DTXdata <- data.frame(DTX, DTXrep)
colnames(DTXdata) <- c("length","TEs")
  } else {DTXdata <- data.frame(empty, rep("DTX", 1))
  }
colnames(DTXdata) <- c("length","TEs")


if (exists("DHX"))
  { 
DHXnb <- nrow(DHX)
DHXrep <- rep("DHX", length(DHXnb)) 
DHXdata <- data.frame(DHX, DHXrep)
colnames(DHXdata) <- c("length","TEs")
  } else {DHXdata <- data.frame(empty, rep("DHX", 1))
  }
colnames(DHXdata) <- c("length","TEs")

if (exists("TXX"))
  { 
TXXnb <- nrow(TXX)
TXXrep <- rep("TXX", length(TXXnb)) 
TXXdata <- data.frame(TXX, TXXrep)
colnames(TXXdata) <- c("length","TEs")
  } else {TXXdata <- data.frame(empty, rep("TXX", 1))
  }
colnames(TXXdata) <- c("length","TEs")

if (exists("TOTAL"))
  { 
TOTALnb <- nrow(TOTAL)
TOTALrep <- rep("TOTAL", length(TOTALnb)) 
TOTALdata <- data.frame(TOTAL, TOTALrep)
colnames(TOTALdata) <- c("length","TEs")
  } else {TOTALdata <- data.frame(empty, rep("TOTAL", 1))
  }
colnames(TOTALdata) <- c("length","TEs")

#Merge

ClassITEs <- rbind(RLCdata, RLGdata, RLXdata, Caulimoviridaedata, RIXdata, RXXdata, RSXdata)
ClassIITEs <- rbind(DTMdata, DTAdata, DTCdata, DTHdata, DTXdata, DHXdata, TXXdata)

#fonction from https://medium.com/@gscheithauer/how-to-add-number-of-observations-to-a-ggplot2-boxplot-b22710f7ef80 to add text
stat_box_data <- function(y, upper_limit = max(ClassITEs$length) * 1.15) {
  return( 
    data.frame(
      y = 0.95 * upper_limit,
      label = paste('count =', length(y), '\n',
                    'mean =', round(mean(y), 1), '\n')
    )
  )
}


plot1 <- ggplot(ClassITEs, aes(x = TEs, y = length)) + geom_violin(trim = FALSE) + geom_boxplot( width=0.2) + stat_summary(fun.data = stat_box_data, geom = "text",hjust = 0.5,vjust = 0.2, size=3)

plot2 <- ggplot(ClassIITEs, aes(x = TEs, y = length)) + geom_violin(trim = FALSE) + geom_boxplot(width=0.2) + stat_summary(fun.data = stat_box_data, geom = "text",hjust = 0.5,vjust = 0.2, size=3)

pdf(file = "../ParasiTE_output/Distance_neighbor_TEs_downstream.pdf")
grid.arrange(plot1, plot2, nrow=2)
dev.off()
