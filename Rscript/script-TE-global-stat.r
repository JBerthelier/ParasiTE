#!/usr/bin/env Rscript

#load libraries
library(ggplot2)
library(gridExtra)

#load data

#intragenic
#intergenic

if (file.exists("../ParasiTE_output/Results/Annotations_TEs/All_TEs.bed"))
total_TEs <- read.delim("../ParasiTE_output/Results/Annotations_TEs/All_TEs.bed", header = FALSE, sep = "\t", dec = ".")
if (file.exists("../ParasiTE_output/Results/Annotations_TEs/All_Intragenic_TEs.bed"))
intragenic <- read.delim("../ParasiTE_output/Results/Annotations_TEs/All_Intragenic_TEs.bed", header = FALSE, sep = "\t", dec = ".")
if (file.exists("../ParasiTE_output/Results/Annotations_TEs/All_Intergenic_TEs.bed"))
intergenic <- read.delim("../ParasiTE_output/Results/Annotations_TEs/All_Intergenic_TEs.bed", header = FALSE, sep = "\t", dec = ".")

#create identifiers intragenic
total_TEsnb <- nrow(total_TEs)
total_TEsrep <- rep("Total", length(1))
total_TEs_ident <- data.frame(total_TEsrep, total_TEsnb )
colnames(total_TEs_ident) <- c("group","number" )


intragenicnb <- nrow(intragenic)
intragenicrep <- rep("Intragenic", length(1))
intragenic_ident <- data.frame(intragenicrep, intragenicnb )
colnames(intragenic_ident) <- c("group","number" )


intergenicnb <- nrow(intergenic)
intergenicrep <- rep("Intergenic", length(1))
intergenic_ident <- data.frame(intergenicrep, intergenicnb )
colnames(intergenic_ident) <- c("group","number" )
print(intergenic_ident)

tot_TEs <- rbind(total_TEs_ident, intergenic_ident, intragenic_ident)
print(tot_TEs)

#plot
pdf(file="../ParasiTE_output/Results/Plots/global-stat.pdf")
plot <- ggplot(tot_TEs, aes(x=group, y=number)) + geom_bar(stat="identity", fill="black") + geom_text(aes(label=number), vjust=-0.5, color="black", size=5) + theme_grey(base_size=22) +theme(axis.title.x = element_blank()) + ggtitle("Proportion of TEs") + ylab("Number of TEs")

grid.arrange(plot, nrow=1)
dev.off()
