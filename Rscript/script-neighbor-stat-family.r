#!/usr/bin/env Rscript

#load libraries
library(ggplot2)
library(cowplot)

#load data

if (file.exists("../ParasiTE_output/Results/Annotations_TEs/Neighbor_TEs_upstream.bed"))
neighbor_TEs_upstream <- read.delim("../ParasiTE_output/Results/Annotations_TEs/Neighbor_TEs_upstream.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/Results/Annotations_TEs/Neighbor_TEs_downstream.bed"))
neighbor_TEs_downstream <- read.delim("../ParasiTE_output/Results/Annotations_TEs/Neighbor_TEs_downstream.bed", header = FALSE, sep = "\t", dec = ".")

#create identifiers neighbor
neighbor_TEs_upstreamnb <- nrow(neighbor_TEs_upstream)
neighbor_TEsrep <- rep("neighbor", length(1))
neighbor_TEs_upstreamrep <- rep("upstream", length(1))
neighbor_TEs_upstream_ident <- data.frame(neighbor_TEsrep, neighbor_TEs_upstreamrep, neighbor_TEs_upstreamnb )
colnames(neighbor_TEs_upstream_ident) <- c("class", "group","number" )


neighbor_TEs_downstreamnnb <- nrow(neighbor_TEs_downstream)
neighbor_TEs_downstreamnrep <- rep("downstream", length(1))
neighbor_TEs_downstreamn_ident <- data.frame(neighbor_TEsrep, neighbor_TEs_downstreamnrep, neighbor_TEs_downstreamnnb )
colnames(neighbor_TEs_downstreamn_ident) <- c("class", "group","number" )


tot_neighbor <- rbind(neighbor_TEs_downstreamn_ident, neighbor_TEs_upstream_ident)

#Familly 
#check if data exist and load them

if (file.exists("../ParasiTE_output/neighbor_TEs-length/RLC-neighbor-TE-downstream.bed"))
RLCi <- read.delim("../ParasiTE_output/neighbor_TEs-length/RLC-neighbor-TE-downstream.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/neighbor_TEs-length/RLG-neighbor-TE-downstream.bed"))
RLGi <- read.delim("../ParasiTE_output/neighbor_TEs-length/RLG-neighbor-TE-downstream.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/neighbor_TEs-length/Caulimoviridae-neighbor-TE-downstream.bed"))
Caulimoviridaei <- read.delim("../ParasiTE_output/neighbor_TEs-length/Caulimoviridae-neighbor-TE-downstream.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/neighbor_TEs-length/RYX-neighbor-TE-downstream.bed"))
RYXi <- read.delim("../ParasiTE_output/neighbor_TEs-length/RYX-neighbor-TE-downstream.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/neighbor_TEs-length/RIX-neighbor-TE-downstream.bed"))
RIXi <- read.delim("../ParasiTE_output/neighbor_TEs-length/RIX-neighbor-TE-downstream.bed", header = FALSE, sep = "\t", dec = ".")

if (file.exists("../ParasiTE_output/neighbor_TEs-length/RLX-neighbor-TE-downstream.bed"))
RLXi <- read.delim("../ParasiTE_output/neighbor_TEs-length/RLX-neighbor-TE-downstream.bed", header = FALSE, sep = "\t", dec = ".")

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

#create identifiers
downrep <- rep("downstream", 1)
uprep <- rep("upstream", 1)

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

#DTA#
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

MergeTEs <- rbind( TXXedata, DHXedata, DTHedata, DTCedata, DTAedata, DTMedata, DTTedata, RSXedata, RXXedata, RIXedata, RYXedata, Caulimoviridaeedata, RLXedata, RLGedata, RLCedata,
TXXidata, DHXidata, DTHidata, DTCidata, DTAidata, DTMidata, DTTidata, RSXidata, RXXidata, RIXidata, RYXidata, Caulimoviridaeidata, RLXidata, RLGidata, RLCidata )

#Export txt file
Export <- MergeTEs[c("group", "TEs", "number")]
write.table(Export, "../ParasiTE_output/Results/Plots/neighbor_TEs_3.txt", sep="\t", row.names = F)



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

  
TEsplot  <- ggplot(data=MergeTEs, aes(x=TEs, y=number, fill=group)) + geom_bar(stat="identity")  + coord_flip() + theme_grey(base_size=20) + scale_color_manual(values=c( "blue","green", "red")) + ggtitle("Proportion of neighbor TE families")+ ylab("Number of TEs") + theme(legend.position="none", axis.title.x =element_text(size=15), axis.title.y = element_blank(), plot.title=element_text(size=15))
 
plot_neighbor <- ggplot(tot_neighbor, aes(x= class, y = number, fill= group)) + geom_col() + geom_text(aes(label = number), position = position_stack(vjust = 0.5), size=6 ) 
pie <- plot_neighbor+ coord_polar("y", start=0) + blank_theme + ggtitle("Proportion of neighbor TE events") + ylab("Number of TEs")
Global <- plot_grid(pie, TEsplot, align= "h", rel_widths= c(2,3))
save_plot("../ParasiTE_output/Results/Plots/neighbor_TEs_1.pdf", Global, ncol=2)
dev.off()
