library(stringr)
library(data.table) ## v 1.9.6+ 
library(dplyr)
library(splitstackshape)
library(tidyr)


####### 1/ Loading of input data
## 1.1/ Stringtie transcript annotation is loaded
print("loading input data")
All_transcript <-  read.table("../ParasiTE_output/transcript_annotation_wo_FP.bed", header = FALSE, sep = '\t', check.names=FALSE)
colnames(All_transcript)<-c("chromosome","start", "end", "poin", "other", "strand","tool", "type", "point","id")

## 1.2/ TE-related transcripts candidates are loaded

transposons <- read.table("../ParasiTE_output/Results/Annotations_chimerick_TE_events/Candidate_TEs.bed", header = FALSE, sep = '\t', check.names=FALSE)
colnames(transposons)<-c("techromosome","testart", "teend", "tename", "other", "strand","localisation","method", "point", "idte","exon_chromosome","exon_start", "exon_end", "exonpoint", "exonother", "exonstrand","exontool","exonmatch", "exonpoint2", "idexon", "TE_vs_exonlenght")
#check if the cov is indicated for exon, beacause Stringtie do not give the coverage if we use the option -R in Stringtie. Therefore this cause problem for Parasite.
#Here, the script check if the coverage is indicated, if not it add a false one 0.000000, to make the script work.
#FOR ANNOTATION OF ARAPORT11 
#transposons$idexon <- ifelse(grepl("cov",transposons$idexon),transposons$idexon,paste0(transposons$idexon," cov 0.000000;"))
#FOR ANNOTATION OF STRINGTIE R 
transposons$idexon <- ifelse(grepl("cov",transposons$idexon),transposons$idexon,paste0(transposons$idexon,"cov0.000000;"))


## 1.3/ Stringtie exon annotation is loaded

All_exons <- read.table("../ParasiTE_output/exon_transcript_annotation_wo_FP.bed", header = FALSE, sep = '\t', check.names=FALSE)
colnames(All_exons)<-c("chromosome","start", "end", "poin", "other", "strand","tool", "type", "point","id")
#check if the cov is indicated for exon, beacause Stringtie do not give the coverage if we use the option -R in Stringtie. Therefore this cause problem for Parasite.
#Here, the script check if the coverage is indicated, if not it add a false one 0.000000, to make the script work.
#FOR ANNOTATION OF ARAPORT11 
#All_exons$id <- ifelse(grepl("cov",All_exons$id),All_exons$id,paste0(All_exons$id," cov 0.000000;"))
#FOR ANNOTATION OF STRINGTIE R 
transposons$idexon <- ifelse(grepl("cov",transposons$idexon),transposons$idexon,paste0(transposons$idexon,"cov0.000000;"))

####### 2/ Important impormation are extracted
##2.1/ Extraction of important information from the Stringtie transcript annotation
head(All_transcript)
All_transcript <- setDT(All_transcript)[, paste0("id", 1:2) := tstrsplit(All_transcript$id, "; cov")]
#All_transcript <- setDT(All_transcript)[, paste0("indice", 1:3) := tstrsplit(All_transcript$id2, ";")] error A
All_transcript <- cSplit(All_transcript,"id2",";")
head(All_transcript) 
All_transcript <- setDT(All_transcript)[, paste0("COV", 1:2) := tstrsplit(All_transcript$id2_1, " ")]
head(All_transcript) 
All_transcript <- setDT(All_transcript)[, paste0("FPKM", 1:2) := tstrsplit(All_transcript$id2_2, " ")]
#All_transcript <- setDT(All_transcript)[, paste0("FPKM", 1:3) := tstrsplit(All_transcript$id2_2, " ")] error A
All_transcript <- setDT(All_transcript)[, paste0("TPM", 1:2) := tstrsplit(All_transcript$id2_3, " ")]
#All_transcript <- setDT(All_transcript)[, paste0("TPM", 1:3) := tstrsplit(All_transcript$id2_3, " ")] error A
All_transcript <- cSplit(All_transcript,"id1",";")
All_transcript <- setDT(All_transcript)[, paste0("id_Gene", 1:2) := tstrsplit(All_transcript$id1_1, " ")]
All_transcript$id_iso2 <-  as.character(All_transcript$id_iso2)
All_transcript <- setDT(All_transcript)[, paste0("id_iso", 1:2) := tstrsplit(All_transcript$id1_2, " ")]
All_transcript$COV2 <-  as.character(All_transcript$COV2)
All_transcript$FPKM3 <-  as.character(All_transcript$FPKM3)
All_transcript$TPM3 <-  as.character(All_transcript$TPM3)

##2.2/ Extraction of important information from the Stringtie exon annotation
All_exons <- cSplit(All_exons,"id",";")
All_exons <- setDT(All_exons)[, paste0("id_iso", 1:2) := tstrsplit(All_exons$id_2, " ")]

##2.2/ Extraction of important information from the TE annotation
head(transposons)
transposons <- setDT(transposons)[, paste0("cov", 1:2) := tstrsplit(transposons$idexon, "; cov")]
transposons <- setDT(transposons)[, paste0("ref", 1:2) := tstrsplit(transposons$cov1, "reference_id")]
transposons <- cSplit(transposons,"idexon",";")
transposons <- setDT(transposons)[, paste0("id_gene", 1:2) := tstrsplit(transposons$idexon_1, " ")]
transposons <- setDT(transposons)[, paste0("id_iso", 1:2) := tstrsplit(transposons$idexon_2, " ")]
transposons <- setDT(transposons)[, paste0("cov_exon", 1:2) := tstrsplit(transposons$cov2, ";")]
transposons <- setDT(transposons)[, paste0("idexon_numb", 1:2) := tstrsplit(transposons$idexon_3, " ")]
transposons$cov2 <-  as.character(transposons$cov2)
transposons$ref2 <-  as.character(transposons$ref2)
transposons$localisation <-  as.character(transposons$localisation)
transposons$method <-  as.character(transposons$method)
transposons$tename <-  as.character(transposons$tename)
transposons$exonstrand <-  as.character(transposons$exonstrand)
transposons$idexon1 <-  as.character(transposons$idexon1)
transposons$exon_start <-  as.character(transposons$exon_start)
transposons$exon_end <-  as.character(transposons$exon_end)
transposons$exon_chromosome  <-  as.character(transposons$exon_chromosome)
transposons$TE_vs_exonlenght <-  as.character(transposons$TE_vs_exonlenght)

##2.3/Filters 
###transcript with exon coverage < 1x can be removed here - OFF
#transposons$cov_exon1 <-  as.numeric(transposons$cov_exon1)
#transposons <- transposons[!( transposons$cov_exon1 < 1) ,]

####### 3/ Calculation of the frequence of TE in transcripts of every genes 
##3.1/ Calculate the number of exons for all transcripts
All_exons_list <- list(All_exons$id_iso2)
freq_All_exons <- as.data.frame(table(comb = do.call(paste, All_exons_list )))

##3.2/ Calculate the total number of transcripts for each genes
isoform_cov <- cbind(All_transcript$id_Gene2, All_transcript$id_iso2, All_transcript$COV2, All_transcript$FPKM3, All_transcript$TPM3)
colnames(isoform_cov)<-c("gene", "comb","coverage", "FPKM", "TPM")
isoform_cov <- as.data.frame(isoform_cov)
Genelist <- list(All_transcript$id_Gene2)
freqTranscript <- as.data.frame(table(comb = do.call(paste, Genelist )))
colnames(freqTranscript)<-c("gene", "Total_number_transcript")

##mute warning
options(warn=-1)
isoform_cov <- left_join(isoform_cov, freqTranscript, by = c('gene'))
##restor warning
options(warn=0)


##3.3/ Calculate the number of transcripts with a specific TE

TE_transcript_count <- transposons[,c("id_gene2", "id_iso2", "tename")]
TE_transcript_count$TEtranscript <- paste( TE_transcript_count$tename, TE_transcript_count$id_gene2)
TE_transcript_count <- TE_transcript_count[!duplicated(TE_transcript_count), ]
transcripts_TEs_list <- list(TE_transcript_count$TEtranscript)
freqTranscripts_TEs_list  <- as.data.frame(table(comb = do.call(paste, transcripts_TEs_list )))
colnames(freqTranscripts_TEs_list)<-c("TE_gene", "Number_transcript_with_TE")

#join interesting data
transposons$TEtranscript <- paste(transposons$tename, transposons$id_gene2)
TEexon <- cbind(transposons$TEtranscript, transposons$id_iso2, transposons$id_gene2, transposons$ref2, transposons$cov_exon1, transposons$exon_chromosome, transposons$exon_start, transposons$exon_end, transposons$TE_vs_exonlenght , transposons$tename, transposons$techromosome, transposons$testart, transposons$teend, transposons$localisation, transposons$method, transposons$exonic, transposons$exonstrand, transposons$idexon_numb2 )
colnames(TEexon)<-c("TE_gene","comb", "gene_name", "gene_ref", "exon_coverage", "exon_chromosome", "exon_start", "exon_end", "TE_vs_exonlenght" , "te_name", "te_chromosome", "te_start", "te_end", "te_localisation", "method","exon_strand", "exon_number")
TEexon <- as.data.frame(TEexon)

options(warn=-1)
join <- left_join(TEexon, freqTranscripts_TEs_list, by = c('TE_gene'))
options(warn=0)

join <- join[,c("comb", "TE_gene", "gene_name", "gene_ref", "exon_coverage", "exon_chromosome", "exon_start", "exon_end", "TE_vs_exonlenght" , "te_name", "te_chromosome", "te_start", "te_end", "te_localisation", "method", "exon_strand", "exon_number", "Number_transcript_with_TE")]

options(warn=-1)
join <- left_join(join, freq_All_exons, by = c('comb'))
join <- left_join(join, isoform_cov, by = c('comb'))
options(warn=0)

colnames(join)<-c("Transcript_id", "TE_gene", "Gene_id", "Gene_ref_id", "Coverage_exon", "exon_chromosome", "exon_start", "exon_end", "TE_vs_exonlenght", "TE_id", "TE_chromosome", "TE_start", "TE_end", "TE_localisation", "method", "exon_strand", "exon_number", "Total_transcripts_with_this_TEs", "total_exon_number", "nepasprendre", "Coverage_transcript", "FPKM_transcript", "TPM_transcript", "Total_number_transcripts")

join <- join[,c("Transcript_id", "TE_gene", "Gene_id", "Gene_ref_id", "Total_number_transcripts", "Total_transcripts_with_this_TEs", "FPKM_transcript", "TPM_transcript", "Coverage_transcript", "Coverage_exon", "exon_chromosome", "exon_start", "exon_end", "TE_vs_exonlenght", "TE_id", "TE_chromosome", "TE_start", "TE_end" ,"TE_localisation", "method", "exon_strand", "exon_number", "total_exon_number")]

join$newcolumn <- 0
join$newcolumn <- with(join, ifelse(exon_number == "1" & total_exon_number == "1", "single" , "0"))
join$newcolumn2 <- with(join, ifelse(newcolumn == "0" & exon_strand == "+" & exon_number == total_exon_number, "ending", newcolumn ))
join$newcolumn3 <- with(join, ifelse(newcolumn2 == "0" & exon_strand == "-" & exon_number == total_exon_number, "ending", newcolumn2 ))
join$newcolumn4 <- with(join, ifelse(newcolumn3 == "0" & exon_strand == "+" & exon_number == "1" & total_exon_number > 1, "begining", newcolumn3 ))
join$newcolumn5 <- with(join, ifelse(newcolumn4 == "0" & exon_strand == "-" & exon_number == "1" & total_exon_number > 1, "begining", newcolumn4 ))
join$newcolumn6 <- with(join, ifelse(newcolumn5 == "0" & exon_number != total_exon_number & total_exon_number > 1, "middle", newcolumn5 ))

join <-subset( join, select = -newcolumn)
join <-subset( join, select = -newcolumn2)
join <-subset( join, select = -newcolumn3)
join <-subset( join, select = -newcolumn4)
join <-subset( join, select = -newcolumn5)

colnames(join)<-c("Transcript_id", "TE_gene", "Gene_id", "Gene_ref_id", "Total_number_transcripts", "Total_transcripts_with_this_TEs", "FPKM_transcript", "TPM_transcript", "Coverage_transcript", "Coverage_exon", "exon_chromosome", "exon_start", "exon_end", "TE_vs_exonlenght", "TE_id", "TE_chromosome", "TE_start", "TE_end" ,"TE_localisation", "method", "exon_strand", "exon_number", "total_exon_number", "TE_position_isoform" )

join_gene <- join[,c("Gene_id", "Transcript_id", "TE_gene", "Gene_ref_id", "Total_number_transcripts", "Total_transcripts_with_this_TEs", "FPKM_transcript", "TPM_transcript", "Coverage_transcript", "Coverage_exon", "exon_chromosome", "exon_start", "exon_end", "TE_vs_exonlenght", "TE_id", "TE_chromosome", "TE_start", "TE_end" ,"TE_localisation", "method", "exon_strand", "exon_number", "total_exon_number", "TE_position_isoform" )]

stat <- do.call(data.frame, aggregate(Gene_id ~ TE_position_isoform, data=join, FUN=table))
stat <- as.data.frame(t(stat))
stat <- setDT(stat, keep.rownames = TRUE )[]
stat = stat[-1,]

colnames(stat)<-c("id","begining", "ending", "middle", "single")
stat <- setDT(stat)[, paste0("id", 1:2) := tstrsplit(stat$id, "Gene_id.")]
stat_gene <- stat[,c("id2","begining", "middle", "ending", "single")]

colnames(stat_gene)<-c("Gene_id","begining", "middle", "ending", "single")

options(warn=-1)
join_gene <- left_join(join_gene, stat_gene, by = c('Gene_id'))
options(warn=0)

join_gene <- join_gene[,c("Transcript_id", "TE_gene", "Gene_id", "Gene_ref_id", "Total_number_transcripts", "Total_transcripts_with_this_TEs", "FPKM_transcript", "TPM_transcript", "Coverage_transcript", "Coverage_exon", "exon_chromosome", "exon_start", "exon_end", "TE_vs_exonlenght", "TE_id", "TE_chromosome", "TE_start", "TE_end" ,"TE_localisation", "method", "exon_strand", "exon_number", "total_exon_number", "TE_position_isoform" , "begining", "middle" , "ending" , "single" )]

#add a freq coulumn of TE isoform
join_gene$freq <- 0
join_gene$Freq_TE_isoform = ( join_gene$Total_transcripts_with_this_TEs/ join_gene$Total_number_transcripts )

#add a coulumn of the difference proportion size between TE length and exon length
join_gene$TE_vs_exon <- 0
join_gene$TE_vs_exonlenght <-  as.numeric(levels(join_gene$TE_vs_exonlenght))[join_gene$TE_vs_exonlenght]
join_gene$TE_end <-  as.numeric(levels(join_gene$TE_end))[join_gene$TE_end]
join_gene$TE_start <-  as.numeric(levels(join_gene$TE_start))[join_gene$TE_start]

join_gene$TE_vs_exon = ( join_gene$TE_vs_exonlenght / ( join_gene$TE_end - join_gene$TE_start ) *100 )

###### remove FP rules
###### genes with only a TE showing a single isoform
join_gene<-join_gene[!(join_gene$begining ==0 & join_gene$middle == 0 & join_gene$ending == 0) ,]

########################### pre-Detailed results

detailed_join_gene <- join_gene[,c("Transcript_id", "Gene_id", "Gene_ref_id", "Total_number_transcripts", "Total_transcripts_with_this_TEs", "Freq_TE_isoform", "FPKM_transcript", "TPM_transcript", "Coverage_transcript", "Coverage_exon", "exon_chromosome", "exon_start", "exon_end", "TE_id", "TE_chromosome", "TE_start", "TE_end", "TE_vs_exon","TE_localisation", "method", "exon_strand", "exon_number", "total_exon_number", "TE_position_isoform" )]

########################### pre-Overview results

stat_TE_gene <- do.call(data.frame, aggregate(TE_gene ~ TE_position_isoform, data=join_gene, FUN=table))
stat_TE_gene <- as.data.frame(t(stat_TE_gene))
stat_TE_gene <- setDT(stat_TE_gene, keep.rownames = TRUE )[]
stat_TE_gene = stat_TE_gene[-1,]
colnames(stat_TE_gene)<-c("id","begining", "ending", "middle", "single")
stat_TE_gene <- setDT(stat_TE_gene)[, paste0("id", 1:2) := tstrsplit(stat_TE_gene$id, "TE_gene.")]
stat_TE_gene <- stat_TE_gene[,c("id2","begining", "middle", "ending", "single")]
colnames(stat_TE_gene)<-c("TE_gene","begining", "middle", "ending", "single")

##################### Potential AS splicing sites check analyses
#check if the analyses TEs are overlaping a splicing site

splicing_sites <- read.table("../ParasiTE_output/splicing_sites.gff", header = FALSE, sep = '\t', check.names=FALSE, stringsAsFactors=F)
splicing_sites <- splicing_sites[1:3]
colnames(splicing_sites)<-c("chromosome","start", "end")

splicing_sites_start <- splicing_sites[,c("chromosome","start", "start")]
splicing_sites_end <- splicing_sites[,c("chromosome","end", "end")]

write.table(splicing_sites_start,file="../ParasiTE_output/splicing_sites_start.gff",sep = "\t",row.names=F, col.names=F)
write.table(splicing_sites_end,file="../ParasiTE_output/splicing_sites_end.gff",sep = "\t",row.names=F, col.names=F)
system(paste("cat ../ParasiTE_output/splicing_sites_start.gff ../ParasiTE_output/splicing_sites_end.gff > ../ParasiTE_output/splicing_sites_start_end.gff"))

ssites <- "awk 'BEGIN{FS=OFS=\"\t\"}{print $1,$2-1,$3}' ../ParasiTE_output/splicing_sites_start_end.gff > ../ParasiTE_output/pre-splicing_sites_start_end.bed"
system(paste(ssites))

splicing_sites_start_end <- read.table("../ParasiTE_output/pre-splicing_sites_start_end.bed", header = FALSE, sep = '\t', check.names=FALSE,  stringsAsFactors=F)
write.table(splicing_sites_start_end,file="../ParasiTE_output/splicing_sites_start_end.bed",sep = "\t",row.names=F, col.names=F)

candidates_TEs <- read.table("../ParasiTE_output/Results/Annotations_chimerick_TE_events/Candidate_TEs.bed", header = FALSE, sep = '\t', check.names=FALSE, stringsAsFactors=F)
candidates_TEs <- subset(candidates_TEs, select=c(1:4))
write.table(candidates_TEs,file="../ParasiTE_output/TEcandidates.bed",sep = "\t",row.names=F, col.names=F)

system(paste("sort ../ParasiTE_output/splicing_sites_start_end.bed | uniq > ../ParasiTE_output/splicing_sites_start_end.uniq.bed &&
bedtools sort -i ../ParasiTE_output/splicing_sites_start_end.uniq.bed > ../ParasiTE_output/splicing_sites_start_end.uniq.sorted.bed"))

system(paste("bedtools intersect -nonamecheck -wo -a ../ParasiTE_output/TEcandidates.bed -b ../ParasiTE_output/splicing_sites_start_end.uniq.sorted.bed > ../ParasiTE_output/TEcandidates_splicing_sites.bed"))
AS_TEcandidates_ss <- read.table("../ParasiTE_output/TEcandidates_splicing_sites.bed", header = FALSE, sep = '\t', check.names=FALSE)

AS_TEcandidates_ss <- AS_TEcandidates_ss[4]
colnames(AS_TEcandidates_ss)<-c("TE_id")
head(AS_TEcandidates_ss,10)
write.table(AS_TEcandidates_ss,file="../ParasiTE_output/Results/Annotations_chimerick_TE_events/Result_TE_candidates_splicing_sites.bed",sep = "\t",row.names=F, col.names=F)
head(AS_TEcandidates_ss,10)

########################### CATANA

system(paste("rm -Rf ../CATANA/Classification"))
system(paste("cd ../CATANA && perl CATANA.pl -i ../ParasiTE_output/Stringtie_annotation.gff3 -o Classification"))

system(paste("cd ../CATANA/Classification && cat A5SS.gff A3SS.gff SE.gff RI.gff MXE.gff MSE.gff > Alternative_splicing_events.gff &&
awk -F'\t' 'x$5' Alternative_splicing_events.gff > pre-Alternative_splicing_events.gff &&
bedtools sort -i  pre-Alternative_splicing_events.gff >  Alternative_splicing_events.sorted.gff && 
mkdir ../../ParasiTE_output/Isoform_classification &&
cp Alternative_splicing_events.sorted.gff ../../ParasiTE_output/Isoform_classification/Alternative_splicing_events.sorted.gff"))
#awk command to remove line that have a "." at column 5

system(paste("cd ../CATANA/Classification && cat ATSS.gff ATTS.gff AFE.gff ALE.gff > Alternative_transcription_products.gff &&
awk -F'\t' 'x$5' Alternative_transcription_products.gff > pre-Alternative_transcription_products.gff
bedtools sort -i  pre-Alternative_transcription_products.gff >  Alternative_transcription_products.sorted.gff && 
cp Alternative_transcription_products.sorted.gff ../../ParasiTE_output/Isoform_classification/Alternative_transcription_products.sorted.gff"))

########################### transcript alternative splicing events

system(paste("awk -F'\t' '$3~/exon/' ../ParasiTE_output/Isoform_classification/Alternative_splicing_events.sorted.gff > ../ParasiTE_output/Isoform_classification/Alternative_splicing_events.mRNA.sorted.gff"))

system(paste("convert2bed --input=GFF < ../ParasiTE_output/Isoform_classification/Alternative_splicing_events.mRNA.sorted.gff > ../ParasiTE_output/Isoform_classification/Alternative_splicing_events.mRNA.sorted.bed"))

CATANA <- read.table("../ParasiTE_output/Isoform_classification/Alternative_splicing_events.mRNA.sorted.bed", header = FALSE, sep = '\t', check.names=FALSE)

colnames(CATANA)<-c("chr","start", "end","info","point","strand","isoform","mRNA", "other","id")

#There is an error here, the number of column is different between long read "paste0("id", 1:5)" or short read "paste0("id", 1:6)", I used the fonction mtry to fix that
mtry <- try(setDT(CATANA)[, paste0("id", 1:4) := tstrsplit(CATANA$id, ";")])
if (inherits(mtry, "try-error")) {
CATANA <- setDT(CATANA)[, paste0("id", 1:5) := tstrsplit(CATANA$id, ";")]
CATANA <- setDT(CATANA)[, paste0("tid", 1:2) := tstrsplit(CATANA$id4, "=")]
CATANA <- subset( CATANA, select = -id5 )
} else {
CATANA <- setDT(CATANA)[, paste0("id", 1:4) := tstrsplit(CATANA$id, ";")]
CATANA <- setDT(CATANA)[, paste0("tid", 1:2) := tstrsplit(CATANA$id4, "=")]
}

head(CATANA,10)

CATANA <- cSplit(CATANA,"tid2",",")
CATANA <- subset( CATANA, select = -tid1 )
CATANA <- subset( CATANA, select = -id4 )
CATANA <- subset( CATANA, select = -id3 )
CATANA <- subset( CATANA, select = -id2 )
CATANA <- subset( CATANA, select = -id1 )
CATANA <- subset( CATANA, select = -id )
CATANA <- subset( CATANA, select = -other )
CATANA <- subset( CATANA, select = -mRNA )
CATANA <- subset( CATANA, select = -strand )
CATANA <- subset( CATANA, select = -point )
CATANA <- subset( CATANA, select = -info )
CATANA <- subset( CATANA, select = -chr )


CATANA$isoform <- paste(CATANA$isoform,CATANA$start,CATANA$end)
CATANA <- subset( CATANA, select = -start )
CATANA <- subset( CATANA, select = -end )
CATANA$isoform <- as.character(gsub(' ','.',CATANA$isoform))
CATANA = as_tibble(CATANA)  

CATANA = CATANA %>%
 pivot_longer(-isoform, 
              names_to = "income", 
              values_to = "count",
              values_drop_na = TRUE) %>% 
  as.data.frame
CATANA <- CATANA[,c("count", "isoform")]

CATANA <- cSplit(CATANA,"isoform",".")
CATANA$count<- paste(CATANA$count,CATANA$isoform_2,CATANA$isoform_3)
CATANA$count <- as.character(gsub(' ','_',CATANA$count))
CATANA <- subset( CATANA, select = -isoform_2 ) 
CATANA <- subset( CATANA, select = -isoform_3 )
CATANA<-CATANA[!duplicated(CATANA), ]
CATANA <- aggregate(data=CATANA,isoform_1~count,FUN=paste)
colnames(CATANA)<-c("isoform_id","Alternative_splicing")

Alternative_sp <- CATANA

CATANA$Alternative_splicing <- as.character(CATANA$Alternative_splicing)

detailed_join_gene$isoform_id<- paste(detailed_join_gene$Transcript_id,detailed_join_gene$exon_start,detailed_join_gene$exon_end)

detailed_join_gene$isoform_id <- as.character(gsub(' ','_',detailed_join_gene$isoform_id))

detailed_join_gene <- detailed_join_gene[,c("isoform_id", "Transcript_id", "Gene_id", "Gene_ref_id", "Total_number_transcripts", "Total_transcripts_with_this_TEs", "Freq_TE_isoform", "FPKM_transcript", "TPM_transcript", "Coverage_transcript", "Coverage_exon", "exon_chromosome", "exon_start", "exon_end", "TE_id", "TE_chromosome", "TE_start", "TE_end", "TE_localisation", "method", "TE_vs_exon" , "exon_strand", "exon_number", "total_exon_number", "TE_position_isoform")]

options(warn=-1)
detailed_join_gene <- left_join(detailed_join_gene, CATANA, by = c('isoform_id'))
options(warn=0)

detailed_join_gene$Alternative_splicing <- as.character(gsub('[()]','',detailed_join_gene$Alternative_splicing))
detailed_join_gene$Alternative_splicing <- as.character(gsub('c','',detailed_join_gene$Alternative_splicing))
detailed_join_gene$Alternative_splicing <- as.character(gsub('"','',detailed_join_gene$Alternative_splicing))
detailed_join_gene$Alternative_splicing <- as.character(gsub('<>','',detailed_join_gene$Alternative_splicing))

########################### transcript alternative transcription products


system(paste("awk -F'\t' '$3~/exon/' ../ParasiTE_output/Isoform_classification/Alternative_transcription_products.sorted.gff > ../ParasiTE_output/Isoform_classification/Alternative_transcription_products.mRNA.sorted.gff"))

system(paste("convert2bed --input=GFF < ../ParasiTE_output/Isoform_classification/Alternative_transcription_products.mRNA.sorted.gff > ../ParasiTE_output/Isoform_classification/Alternative_transcription_products.mRNA.sorted.bed"))

CATANA <- read.table("../ParasiTE_output/Isoform_classification/Alternative_transcription_products.mRNA.sorted.bed", header = FALSE, sep = '\t', check.names=FALSE)
colnames(CATANA)<-c("chr","start", "end","info","point","strand","isoform","mRNA", "other","id")

#There is an error here, the number of column is different between long read "paste0("id", 1:5)" or short read "paste0("id", 1:6)", I used the fonction mtry to fix that
mtry <- try(setDT(CATANA)[, paste0("id", 1:4) := tstrsplit(CATANA$id, ";")])
if (inherits(mtry, "try-error")) {
CATANA <- setDT(CATANA)[, paste0("id", 1:5) := tstrsplit(CATANA$id, ";")]
CATANA <- setDT(CATANA)[, paste0("tid", 1:2) := tstrsplit(CATANA$id4, "=")]
CATANA <- subset( CATANA, select = -id5 )
} else {
CATANA <- setDT(CATANA)[, paste0("id", 1:4) := tstrsplit(CATANA$id, ";")]
CATANA <- setDT(CATANA)[, paste0("tid", 1:2) := tstrsplit(CATANA$id4, "=")]
}

CATANA <- cSplit(CATANA,"tid2",",")
CATANA <- subset( CATANA, select = -tid1 )
CATANA <- subset( CATANA, select = -id4 )
CATANA <- subset( CATANA, select = -id3 )
CATANA <- subset( CATANA, select = -id2 )
CATANA <- subset( CATANA, select = -id1 )
CATANA <- subset( CATANA, select = -id )
CATANA <- subset( CATANA, select = -other )
CATANA <- subset( CATANA, select = -mRNA )
CATANA <- subset( CATANA, select = -strand )
CATANA <- subset( CATANA, select = -point )
CATANA <- subset( CATANA, select = -info )
CATANA <- subset( CATANA, select = -chr )

CATANA$isoform <- paste(CATANA$isoform,CATANA$start,CATANA$end)
CATANA <- subset( CATANA, select = -start )
CATANA <- subset( CATANA, select = -end )
CATANA$isoform <- as.character(gsub(' ','.',CATANA$isoform))
CATANA = as_tibble(CATANA)  

CATANA = CATANA %>%
 pivot_longer(-isoform, 
              names_to = "income", 
              values_to = "count",
              values_drop_na = TRUE) %>% 
  as.data.frame
CATANA <- CATANA[,c("count", "isoform")]

CATANA <- cSplit(CATANA,"isoform",".")
CATANA$count<- paste(CATANA$count,CATANA$isoform_2,CATANA$isoform_3)
CATANA$count <- as.character(gsub(' ','_',CATANA$count))
CATANA <- subset( CATANA, select = -isoform_2 ) 
CATANA <- subset( CATANA, select = -isoform_3 )
CATANA<-CATANA[!duplicated(CATANA), ]

CATANA <- aggregate(data=CATANA,isoform_1~count,FUN=paste)
colnames(CATANA)<-c("isoform_id","Alternative_transcription")
Alternative_tr <- CATANA

CATANA$Alternative_transcription <- as.character(CATANA$Alternative_transcription)

options(warn=-1)
detailed_join_gene <- left_join(detailed_join_gene, CATANA, by = c('isoform_id'))
options(warn=0)

detailed_join_gene$Alternative_transcription <- as.character(gsub('[()]','',detailed_join_gene$Alternative_transcription))
detailed_join_gene$Alternative_transcription <- as.character(gsub('c','',detailed_join_gene$Alternative_transcription))
detailed_join_gene$Alternative_transcription <- as.character(gsub('"','',detailed_join_gene$Alternative_transcription))
detailed_join_gene$Alternative_transcription <- as.character(gsub('<>','',detailed_join_gene$Alternative_transcription))

################################ splicing_site and exons_borders check analyses

########################### Adding Potential splicing_sites check analyses

AS_check <- detailed_join_gene$TE_id %in% AS_TEcandidates_ss$TE_id
detailed_join_gene <- cbind(detailed_join_gene,AS_check)

########################## Adding Potential exons_borders check analyses
detailed_join_gene_check <- detailed_join_gene

detailed_join_gene_check$exon_start <- as.numeric(levels(detailed_join_gene_check$exon_start))[detailed_join_gene_check$exon_start]
detailed_join_gene_check$exon_end <- as.numeric(levels(detailed_join_gene_check$exon_end))[detailed_join_gene_check$exon_end]

detailed_join_gene_check$Exon_check_start <- 0
detailed_join_gene_check$Exon_check_start <- ifelse((detailed_join_gene_check$TE_start<detailed_join_gene_check$exon_start) & (detailed_join_gene_check$TE_position_isoform=="begining" | detailed_join_gene_check$TE_position_isoform=="ending" | detailed_join_gene_check$TE_position_isoform=="single"),TRUE,FALSE)

detailed_join_gene_check$Exon_check_end <- 0
detailed_join_gene_check$Exon_check_end <- ifelse((detailed_join_gene_check$TE_end>detailed_join_gene_check$exon_end) & (detailed_join_gene_check$TE_position_isoform=="begining" | detailed_join_gene_check$TE_position_isoform=="ending" | detailed_join_gene_check$TE_position_isoform=="single"),TRUE,FALSE)

detailed_join_gene$Exon_check_start <- detailed_join_gene_check$Exon_check_start
detailed_join_gene$Exon_check_end <- detailed_join_gene_check$Exon_check_end

detailed_join_gene <- detailed_join_gene[,c("isoform_id", "Transcript_id", "Gene_id", "Gene_ref_id", "Total_number_transcripts", "Total_transcripts_with_this_TEs", "Freq_TE_isoform", "FPKM_transcript", "TPM_transcript", "Coverage_transcript", "Coverage_exon", "exon_chromosome", "exon_start", "exon_end", "TE_id", "TE_chromosome", "TE_start", "TE_end", "TE_localisation", "method", "AS_check", "Exon_check_start", "Exon_check_end","exon_strand", "exon_number", "total_exon_number", "TE_position_isoform", "Alternative_splicing", "Alternative_transcription")]

###### MXE or MSE are simplified to SE
#I prefere to simplify the AS, and to change MXE and MSE into SE
detailed_join_gene$Alternative_splicing <- as.character(gsub('MSE, SE','SE',detailed_join_gene$Alternative_splicing))
detailed_join_gene$Alternative_splicing <- as.character(gsub('SE, MSE','SE',detailed_join_gene$Alternative_splicing))
detailed_join_gene$Alternative_splicing <- as.character(gsub('MSE','SE',detailed_join_gene$Alternative_splicing))

detailed_join_gene$Alternative_splicing <- as.character(gsub('MXE, SE','SE',detailed_join_gene$Alternative_splicing))
detailed_join_gene$Alternative_splicing <- as.character(gsub('SE, MXE','SE',detailed_join_gene$Alternative_splicing))
detailed_join_gene$Alternative_splicing <- as.character(gsub('MXE','SE',detailed_join_gene$Alternative_splicing))

###### Rules to limit potential false positives predictions
###########################################################

##FOR AS

#Rule 1: If the Freq_TE_isoform is 1 and that the AS check is false (no overlap with a splicing site), do not give AS prediction
detailed_join_gene2 <- detailed_join_gene
detailed_join_gene$Alternative_splicing <- ifelse(detailed_join_gene$Freq_TE_isoform==1 & detailed_join_gene$AS_check==FALSE,NA, detailed_join_gene2$Alternative_splicing )

#Rule 2: If the intragenic TE has a TE freq of 1 and is located at the first or last exon, if it do not not overlap the correct (start or end) of the exon, there is no ATP prediction  
detailed_join_gene3 <- detailed_join_gene
detailed_join_gene$Alternative_transcription <- ifelse(detailed_join_gene$TE_localisation=="Intragenic" & detailed_join_gene$Freq_TE_isoform==1 & ((detailed_join_gene$TE_position_isoform=="ending" & detailed_join_gene$exon_strand=="+" & detailed_join_gene$Exon_check_end==FALSE) | (detailed_join_gene$TE_position_isoform=="begining" & detailed_join_gene$exon_strand=="+" & detailed_join_gene$Exon_check_start==FALSE) | (detailed_join_gene$TE_position_isoform=="ending" & detailed_join_gene$exon_strand=="-" & detailed_join_gene$Exon_check_start==FALSE) | (detailed_join_gene$TE_position_isoform=="begining" & detailed_join_gene$exon_strand=="-" & detailed_join_gene$Exon_check_end==FALSE) | ( detailed_join_gene$Exon_check_start==FALSE & detailed_join_gene$Exon_check_end==FALSE)),NA, detailed_join_gene3$Alternative_transcription )

#Rule 3
detailed_join_gene4 <- detailed_join_gene
detailed_join_gene$Alternative_transcription <- ifelse(detailed_join_gene$TE_localisation=="Intragenic" & detailed_join_gene$Exon_check_start==FALSE & detailed_join_gene$Exon_check_end==FALSE & detailed_join_gene$AS_check==FALSE & !is.na(detailed_join_gene$Alternative_splicing) ,NA, detailed_join_gene4$Alternative_transcription )

#Rule 4
detailed_join_gene5 <- detailed_join_gene
detailed_join_gene$Alternative_splicing <- ifelse(detailed_join_gene$TE_localisation=="Intragenic" & !is.na(detailed_join_gene$Alternative_transcription) & detailed_join_gene$AS_check==FALSE & (detailed_join_gene$TE_position_isoform=="begining" | detailed_join_gene$TE_position_isoform=="ending" | detailed_join_gene$TE_position_isoform=="single"), NA, detailed_join_gene5$Alternative_splicing )

#Rule 5
detailed_join_gene6 <- detailed_join_gene
detailed_join_gene$Alternative_splicing <- ifelse(detailed_join_gene$TE_localisation=="Intergenic" & ( detailed_join_gene$TE_position_isoform=="begining" | detailed_join_gene$TE_position_isoform=="ending" | detailed_join_gene$TE_position_isoform=="single") & detailed_join_gene$AS_check==FALSE,NA, detailed_join_gene6$Alternative_splicing )

#Rule c1: remove FP prediction for ATP prediction found for isoform composed of two exons
detailed_join_gene8 <- detailed_join_gene
detailed_join_gene$Alternative_transcription <- ifelse((detailed_join_gene$total_exon_number==2) & (detailed_join_gene$TE_position_isoform=="ending") & (grepl("ATTS",detailed_join_gene$Alternative_transcription)==TRUE) & (grepl("ATSS",detailed_join_gene$Alternative_transcription)==TRUE),as.character(gsub('ATSS','@',detailed_join_gene$Alternative_transcription)), detailed_join_gene8$Alternative_transcription )
detailed_join_gene9 <- detailed_join_gene
detailed_join_gene$Alternative_transcription <- ifelse((detailed_join_gene$total_exon_number==2) & (detailed_join_gene$TE_position_isoform=="begining") & (grepl("ATSS",detailed_join_gene$Alternative_transcription)==TRUE) & (grepl("ATTS",detailed_join_gene$Alternative_transcription)==TRUE),as.character(gsub('ATTS','@',detailed_join_gene$Alternative_transcription)), detailed_join_gene9$Alternative_transcription )
detailed_join_gene20 <- detailed_join_gene
detailed_join_gene$Alternative_transcription <- ifelse((detailed_join_gene$TE_position_isoform=="single") & (grepl("ATSS",detailed_join_gene$Alternative_transcription)==TRUE) & (grepl("ATTS",detailed_join_gene$Alternative_transcription)==TRUE) & ((detailed_join_gene$exon_strand=="+" & detailed_join_gene$Exon_check_start==FALSE) | (detailed_join_gene$exon_strand=="-" & detailed_join_gene$Exon_check_end==FALSE)),as.character(gsub('ATSS','@',detailed_join_gene$Alternative_transcription)), detailed_join_gene20$Alternative_transcription )
detailed_join_gene21 <- detailed_join_gene
detailed_join_gene$Alternative_transcription <- ifelse((detailed_join_gene$TE_position_isoform=="single") & (grepl("ATSS",detailed_join_gene$Alternative_transcription)==TRUE) & (grepl("ATTS",detailed_join_gene$Alternative_transcription)==TRUE) & ((detailed_join_gene$exon_strand=="+" & detailed_join_gene$Exon_check_end==FALSE) | (detailed_join_gene$exon_strand=="-" & detailed_join_gene$Exon_check_start==FALSE)),as.character(gsub('ATTS','@',detailed_join_gene$Alternative_transcription)), detailed_join_gene21$Alternative_transcription )
detailed_join_gene10 <- detailed_join_gene
detailed_join_gene$Alternative_transcription <- ifelse((detailed_join_gene$total_exon_number==2 | detailed_join_gene$total_exon_number==1) & (detailed_join_gene$TE_position_isoform=="ending") & (grepl("AFE",detailed_join_gene$Alternative_transcription)==TRUE) & (grepl("ALE",detailed_join_gene$Alternative_transcription)==TRUE),as.character(gsub('AFE','@',detailed_join_gene$Alternative_transcription)), detailed_join_gene10$Alternative_transcription )
detailed_join_gene11 <- detailed_join_gene
detailed_join_gene$Alternative_transcription <- ifelse((detailed_join_gene$total_exon_number==2 | detailed_join_gene$total_exon_number==1) & (detailed_join_gene$TE_position_isoform=="begining") & (grepl("ALE",detailed_join_gene$Alternative_transcription)==TRUE) & (grepl("AFE",detailed_join_gene$Alternative_transcription)==TRUE),as.character(gsub('ALE','@',detailed_join_gene$Alternative_transcription)), detailed_join_gene11$Alternative_transcription )
detailed_join_gene22 <- detailed_join_gene
detailed_join_gene$Alternative_transcription <- ifelse((detailed_join_gene$TE_position_isoform=="single") & (grepl("AFE",detailed_join_gene$Alternative_transcription)==TRUE) & (grepl("ALE",detailed_join_gene$Alternative_transcription)==TRUE) & ((detailed_join_gene$exon_strand=="+" & detailed_join_gene$Exon_check_start==FALSE) | (detailed_join_gene$exon_strand=="-" & detailed_join_gene$Exon_check_end==FALSE)),as.character(gsub('AFE','@',detailed_join_gene$Alternative_transcription)), detailed_join_gene22$Alternative_transcription )
detailed_join_gene23 <- detailed_join_gene
detailed_join_gene$Alternative_transcription <- ifelse((detailed_join_gene$TE_position_isoform=="single") & (grepl("AFE",detailed_join_gene$Alternative_transcription)==TRUE) & (grepl("ALE",detailed_join_gene$Alternative_transcription)==TRUE) & ((detailed_join_gene$exon_strand=="+" & detailed_join_gene$Exon_check_end==FALSE) | (detailed_join_gene$exon_strand=="-" & detailed_join_gene$Exon_check_start==FALSE)),as.character(gsub('ALE','@',detailed_join_gene$Alternative_transcription)), detailed_join_gene23$Alternative_transcription )
detailed_join_gene$Alternative_transcription <- as.character(gsub('\\b, @\\b','',detailed_join_gene$Alternative_transcription))
detailed_join_gene$Alternative_transcription <- as.character(gsub('\\b@, \\b','',detailed_join_gene$Alternative_transcription))
detailed_join_gene$Alternative_transcription <- as.character(gsub('\\b @, \\b',',',detailed_join_gene$Alternative_transcription))
detailed_join_gene$Alternative_transcription <- as.character(gsub('\\b, @, \\b',', ',detailed_join_gene$Alternative_transcription))
detailed_join_gene$Alternative_transcription <- as.character(gsub('\\b@\\b',NA,detailed_join_gene$Alternative_transcription))

#Rule c2: remove FP 3ASS and 5ASS
detailed_join_gene12 <- detailed_join_gene
detailed_join_gene$Alternative_splicing <- ifelse(( detailed_join_gene$TE_position_isoform=="begining" | detailed_join_gene$TE_position_isoform=="ending") & ((detailed_join_gene$exon_strand=="+" & grepl("A5SS",detailed_join_gene$Alternative_splicing)==TRUE & detailed_join_gene$Exon_check_end==FALSE) | (detailed_join_gene$exon_strand=="-" & grepl("A5SS",detailed_join_gene$Alternative_splicing)==TRUE & detailed_join_gene$Exon_check_start==FALSE)), as.character(gsub('A5SS','@',detailed_join_gene$Alternative_splicing)),detailed_join_gene12$Alternative_splicing)
detailed_join_gene13 <- detailed_join_gene
detailed_join_gene$Alternative_splicing <- ifelse(( detailed_join_gene$TE_position_isoform=="begining" | detailed_join_gene$TE_position_isoform=="ending") & ((detailed_join_gene$exon_strand=="+" & grepl("A3SS",detailed_join_gene$Alternative_splicing)==TRUE & detailed_join_gene$Exon_check_start==FALSE) | (detailed_join_gene$exon_strand=="-" & grepl("A3SS",detailed_join_gene$Alternative_splicing)==TRUE &
detailed_join_gene$Exon_check_end==FALSE)), as.character(gsub('A3SS','@',detailed_join_gene$Alternative_splicing)),detailed_join_gene13$Alternative_splicing)
detailed_join_gene$Alternative_splicing <- as.character(gsub('\\b, @\\b','',detailed_join_gene$Alternative_splicing))
detailed_join_gene$Alternative_splicing <- as.character(gsub('\\b@, \\b','',detailed_join_gene$Alternative_splicing))
detailed_join_gene$Alternative_splicing <- as.character(gsub('\\b @, \\b',',',detailed_join_gene$Alternative_splicing))
detailed_join_gene$Alternative_splicing <- as.character(gsub('\\b, @, \\b',', ',detailed_join_gene$Alternative_splicing))
detailed_join_gene$Alternative_splicing <- as.character(gsub('\\b@\\b',NA,detailed_join_gene$Alternative_splicing))

###### Check if exon candidates are overlaping CDS or 5UTR or 3UTR

TE_annotation_candidate <- detailed_join_gene[,c("TE_chromosome", "TE_start", "TE_end", "isoform_id")]
TE_annotation_candidate$TE_chromosome <- as.numeric(levels(TE_annotation_candidate$TE_chromosome))[TE_annotation_candidate$TE_chromosome]

write.table(TE_annotation_candidate,file="../ParasiTE_output/TE_annotation_candidate.txt",sep = "\t",row.names=F,col.names=F)
system(paste("bedtools sort -i  ../ParasiTE_output/TE_annotation_candidate.txt >  ../ParasiTE_output/TE_annotation_candidate.bed"))

##check CDS
system(paste("bedtools intersect -nonamecheck -wo -a ../ParasiTE_output/TE_annotation_candidate.bed -b ../ParasiTE_output/cds_annotation.bed  > ../ParasiTE_output/CDS_results.bed "))
CDS_results <-  read.table("../ParasiTE_output/CDS_results.bed", header = FALSE, sep = '\t', check.names=FALSE)
colnames(CDS_results)<-c("chromosome","start", "end", "isoform_id")
CDS_check <- detailed_join_gene$isoform_id %in% CDS_results$isoform_id
detailed_join_gene <- cbind(detailed_join_gene,CDS_check)

##check 5UTR
system(paste("bedtools intersect -nonamecheck -wo -a ../ParasiTE_output/TE_annotation_candidate.bed -b ../ParasiTE_output/five_annotation.bed  > ../ParasiTE_output/five_results.bed "))
five_results <-  read.table("../ParasiTE_output/five_results.bed", header = FALSE, sep = '\t', check.names=FALSE)
colnames(five_results)<-c("chromosome","start", "end", "isoform_id")
five_check <- detailed_join_gene$isoform_id %in% five_results$isoform_id
detailed_join_gene <- cbind(detailed_join_gene,five_check)

##check 3UTR
system(paste("bedtools intersect -nonamecheck -wo -a ../ParasiTE_output/TE_annotation_candidate.bed -b ../ParasiTE_output/three_annotation.bed  > ../ParasiTE_output/three_results.bed "))
three_results <-  read.table("../ParasiTE_output/three_results.bed", header = FALSE, sep = '\t', check.names=FALSE)
colnames(three_results)<-c("chromosome","start", "end", "isoform_id")
three_check <- detailed_join_gene$isoform_id %in% three_results$isoform_id
detailed_join_gene <- cbind(detailed_join_gene,three_check)


########### FINAL results
#########################
write.table(detailed_join_gene,file="../ParasiTE_output/Results/Annotations_chimerick_TE_events/Results_transcripts_TEs_all_events.tab",sep = "\t",row.names=F,col.names=T)

###############################################################
########################### GENE-TE alternative splicing events
###############################################################

Alternative_sp <- detailed_join_gene[,c("Gene_id", "TE_id","Alternative_splicing")]
Alternative_sp$Gene_TE <- paste(Alternative_sp$Gene_id,Alternative_sp$TE_id)
Alternative_sp$Gene_TE <- as.character(gsub(' ','_',Alternative_sp$Gene_TE ))
Alternative_sp <- Alternative_sp[,c("Gene_TE", "Alternative_splicing")]
Alternative_sp <- Alternative_sp[!duplicated(Alternative_sp), ]

Alternative_sp <- aggregate(Alternative_splicing ~ Gene_TE, data = Alternative_sp, c)
Alternative_sp$Alternative_splicing <- as.character(gsub('[()]','',Alternative_sp$Alternative_splicing))
Alternative_sp$Alternative_splicing <- as.character(gsub('c','',Alternative_sp$Alternative_splicing))
Alternative_sp$Alternative_splicing <- as.character(gsub('"','',Alternative_sp$Alternative_splicing))
Alternative_sp$Alternative_splicing <- as.character(gsub('<>','',Alternative_sp$Alternative_splicing))

Alternative_sp$Alternative_splicing <- sapply(Alternative_sp$Alternative_splicing, function(x) paste(unique(unlist(str_split(x,", "))), collapse = ", "))

overview_gene_TE <- join_gene[,c("Gene_id", "Total_number_transcripts", "Total_transcripts_with_this_TEs", "Freq_TE_isoform", "TE_id", "TE_chromosome", "TE_start", "TE_end", "TE_localisation", "method")]
overview_gene_TE <- overview_gene_TE[!duplicated(overview_gene_TE), ]

overview_gene_TE$Gene_TE <- paste(overview_gene_TE$Gene_id,overview_gene_TE$TE_id)
overview_gene_TE$Gene_TE <- as.character(gsub(' ','_',overview_gene_TE$Gene_TE ))

overview_gene_TE <- overview_gene_TE[,c("Gene_TE","Gene_id", "Total_number_transcripts", "Total_transcripts_with_this_TEs", "Freq_TE_isoform", "TE_id", "TE_chromosome", "TE_start", "TE_end", "TE_localisation" )]

options(warn=-1)
overview_gene_TE <- left_join(overview_gene_TE, Alternative_sp, by = c('Gene_TE'))
options(warn=0)

########################### GENE-TE alternative product events

Alternative_pr <- detailed_join_gene[,c("Gene_id", "TE_id","Alternative_transcription")]
Alternative_pr$Gene_TE <- paste(Alternative_pr$Gene_id,Alternative_pr$TE_id)
Alternative_pr$Gene_TE <- as.character(gsub(' ','_',Alternative_pr$Gene_TE ))
Alternative_pr <- Alternative_pr[,c("Gene_TE", "Alternative_transcription")]
Alternative_pr <- Alternative_pr[!duplicated(Alternative_pr), ]

Alternative_pr <- aggregate(Alternative_transcription ~ Gene_TE, data = Alternative_pr, c)
Alternative_pr$Alternative_transcription <- as.character(gsub('[()]','',Alternative_pr$Alternative_transcription))
Alternative_pr$Alternative_transcription <- as.character(gsub('c','',Alternative_pr$Alternative_transcription))
Alternative_pr$Alternative_transcription <- as.character(gsub('"','',Alternative_pr$Alternative_transcription))
Alternative_pr$Alternative_transcription <- as.character(gsub('<>','',Alternative_pr$Alternative_transcription))
Alternative_pr$Alternative_transcription <- sapply(Alternative_pr$Alternative_transcription, function(x) paste(unique(unlist(str_split(x,", "))), collapse = ", "))

options(warn=-1)
overview_gene_TE <- left_join(overview_gene_TE, Alternative_pr, by = c('Gene_TE'))
options(warn=0)

########################### Adding method information

method_info <- detailed_join_gene[,c("Gene_id", "TE_id","method")]
method_info$Gene_TE <- paste(method_info$Gene_id,method_info$TE_id)
method_info$Gene_TE <- as.character(gsub(' ','_',method_info$Gene_TE ))
method_info <- method_info[,c("Gene_TE", "method")]
method_info <- method_info[!duplicated(method_info), ]

method_info <- aggregate(method ~ Gene_TE, data = method_info, c)
head(method_info)

method_info$method <- as.character(gsub('[()]','',method_info$method))
method_info$method <- as.character(gsub('c','',method_info$method))
method_info$method <- as.character(gsub('"','',method_info$method))
method_info$method <- as.character(gsub('<>','',method_info$method))
method_info$method <- as.character(gsub(':','&',method_info$method))

head(method_info)

method_info$method <- sapply(method_info$method, function(x) paste(unique(unlist(str_split(x,", "))), collapse = ", "))

head(method_info)

method_info$method <- ifelse(method_info$method=="1&3","1&2&3", method_info$method )
method_info$method <- ifelse(method_info$method=="1, 3","1&3", method_info$method )

head(method_info)

options(warn=-1)
overview_gene_TE <- left_join(overview_gene_TE, method_info, by = c('Gene_TE'))
options(warn=0)

########################### Adding exon information

exon_info <- detailed_join_gene[,c("Gene_id", "TE_id","total_exon_number")]
exon_info$Gene_TE <- paste(exon_info$Gene_id,exon_info$TE_id)
exon_info$Gene_TE <- as.character(gsub(' ','_',exon_info$Gene_TE ))
exon_info <- exon_info[,c("Gene_TE", "total_exon_number")]
exon_info <- exon_info[!duplicated(exon_info), ]

exon_info <- aggregate(total_exon_number ~ Gene_TE, data = exon_info, c)
exon_info$total_exon_number <- as.character(gsub('[()]','',exon_info$total_exon_number))
exon_info$total_exon_number <- as.character(gsub('c','',exon_info$total_exon_number))
exon_info$total_exon_number <- as.character(gsub('"','',exon_info$total_exon_number))
exon_info$total_exon_number <- as.character(gsub('<>','',exon_info$total_exon_number))
exon_info$total_exon_number <- as.character(gsub(':',', ',exon_info$total_exon_number))

options(warn=-1)
overview_gene_TE <- left_join(overview_gene_TE, exon_info, by = c('Gene_TE'))
options(warn=0)

overview_gene_TE <- subset( overview_gene_TE, select = -Gene_TE )

########################### Adding position of TEs for every transcripts

overview_gene_TE$TE_gene <- paste(overview_gene_TE$TE_id,overview_gene_TE$Gene_id)
overview_gene_TE$TE_gene <- as.character(gsub(' ','.',overview_gene_TE$TE_gene ))
overview_gene_TE$TE_gene <- as.character(gsub("-",".",overview_gene_TE$TE_gene))
overview_gene_TE <- overview_gene_TE[,c("TE_gene", "Gene_id", "Total_number_transcripts", "Total_transcripts_with_this_TEs", "total_exon_number", "Freq_TE_isoform", "TE_id", "TE_chromosome", "TE_start", "TE_end", "TE_localisation", "method", "Alternative_splicing", "Alternative_transcription" )]

options(warn=-1)
overview_gene_TE <- left_join(overview_gene_TE, stat_TE_gene, by = c('TE_gene'))
options(warn=0)


########################### Adding information about CSD and UTR


##################


overview_gene_TE <- overview_gene_TE[!duplicated(overview_gene_TE), ]

overview_gene_TE <- overview_gene_TE[,c("TE_gene", "Gene_id", "Total_number_transcripts", "Total_transcripts_with_this_TEs", "total_exon_number", "Freq_TE_isoform", "TE_id", "TE_chromosome", "TE_start", "TE_end", "TE_localisation", "method","Alternative_splicing", "Alternative_transcription", "begining", "middle", "ending", "single" )]

write.table(overview_gene_TE,file="../ParasiTE_output/Results/Annotations_chimerick_TE_events/Result_genes_TEs_all_events.tab",sep = "\t",row.names=F,col.names=T)


########################### Separate Gene-TE alternative & no alternative transcript events
################## alternative transcript events

overview_gene_TE_alternative <- overview_gene_TE[!(is.na(overview_gene_TE$Alternative_transcription) & is.na(overview_gene_TE$Alternative_splicing)) ,]

write.table(overview_gene_TE_alternative,file="../ParasiTE_output/Results/Annotations_chimerick_TE_events/Result_genes_TEs_alternative_events.tab",sep = "\t",row.names=F,col.names=T)

########################### Separate Transcripts-TEs alternative & no alternative transcript events
detailed_join_gene_alternative <- detailed_join_gene[!(is.na(detailed_join_gene$Alternative_transcription) & is.na(detailed_join_gene$Alternative_splicing)) ,]

write.table(detailed_join_gene_alternative,file="../ParasiTE_output/Results/Annotations_chimerick_TE_events/Result_transcripts_TEs_alternative_events.tab",sep = "\t",row.names=F,col.names=T)

print("Finished")

