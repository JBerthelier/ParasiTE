
library(stringr)
library(data.table) ## v 1.9.6+ 
library(dplyr)
library(splitstackshape)
library(tidyr)


####### 1/ Loading of input data
## 1.1/ Stringtie transcript annotation is loaded
print("loading input data")
All_transcript <-  read.table("../ParasiTE_output/transcript_annotation.bed", header = FALSE, sep = '\t', check.names=FALSE)
colnames(All_transcript)<-c("chromosome","start", "end", "poin", "other", "strand","tool", "type", "point","id")

## 1.2/  ParasiTE annotation of exonic TE candidates is loaded
transposons <- read.table("../ParasiTE_output/Results/Annotations_chimerick_TE_events/Candidate_TEs.bed", header = FALSE, sep = '\t', check.names=FALSE)
colnames(transposons)<-c("techromosome","testart", "teend", "tename", "other", "strand","localisation","method", "point", "idte","exon_chromosome","exon_start", "exon_end", "exonpoint", "exonother", "exonstrand","exontool","exonmatch", "exonpoint2", "idexon", "TE_vs_exonlenght")

## 1.3/ Stringtie exon annotation is loaded
All_exons <- read.table("../ParasiTE_output/exon_transcript_annotation.bed", header = FALSE, sep = '\t', check.names=FALSE)
colnames(All_exons)<-c("chromosome","start", "end", "poin", "other", "strand","tool", "type", "point","id")

####### 2/ Important impormation are extracted
##2.1/ Extraction of important information from the Stringtie transcript annotation

All_transcript <- setDT(All_transcript)[, paste0("id", 1:2) := tstrsplit(All_transcript$id, "; cov")]
All_transcript <- setDT(All_transcript)[, paste0("indice", 1:3) := tstrsplit(All_transcript$id2, ";")]
All_transcript <- setDT(All_transcript)[, paste0("COV", 1:2) := tstrsplit(All_transcript$indice1, " ")]
All_transcript <- setDT(All_transcript)[, paste0("FPKM", 1:3) := tstrsplit(All_transcript$indice2, " ")]
All_transcript <- setDT(All_transcript)[, paste0("TPM", 1:3) := tstrsplit(All_transcript$indice3, " ")]
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

##2.3/ transcript with exon coverage < 1x can be removed here - OFF
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
join$newcolumn3 <- with(join, ifelse(newcolumn2 == "0" & exon_strand == "-" & exon_number == total_exon_number, "begining", newcolumn2 ))
join$newcolumn4 <- with(join, ifelse(newcolumn3 == "0" & exon_strand == "+" & exon_number == "1" & total_exon_number > 1, "begining", newcolumn3 ))
join$newcolumn5 <- with(join, ifelse(newcolumn4 == "0" & exon_strand == "-" & exon_number == "1" & total_exon_number > 1, "ending", newcolumn4 ))
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

########################### CATANA

system(paste("rm -Rf ../CATANA/Classification"))
system(paste("cd ../CATANA && perl CATANA.pl -i ../ParasiTE_output/Stringtie_annotation.gff3 -o Classification"))

system(paste("cd ../CATANA/Classification && cat A5SS.gff A3SS.gff SE.gff RI.gff MXE.gff MSE.gff > Alternative_splicing_events.gff &&
bedtools sort -i  Alternative_splicing_events.gff >  Alternative_splicing_events.sorted.gff && 
mkdir ../../ParasiTE_output/Isoform_classification &&
cp Alternative_splicing_events.sorted.gff ../../ParasiTE_output/Isoform_classification/Alternative_splicing_events.sorted.gff"))

system(paste("cd ../CATANA/Classification && cat ATSS.gff ATTS.gff AFE.gff ALE.gff > Alternative_transcription_products.gff &&
bedtools sort -i  Alternative_transcription_products.gff >  Alternative_transcription_products.sorted.gff && 
cp Alternative_transcription_products.sorted.gff ../../ParasiTE_output/Isoform_classification/Alternative_transcription_products.sorted.gff"))

########################### transcript alternative spliicing events

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

###########################

###### Remove FP rules
# Intergenic TEs showing ATP & AS (RI, SE, A5SS, A3SS)
detailed_join_gene2 <- detailed_join_gene
detailed_join_gene$Alternative_splicing <- ifelse(detailed_join_gene$TE_localisation=="Intergenic",NA, detailed_join_gene2$Alternative_splicing )

# Remove TE frequences = 1
# detailed_join_gene<-detailed_join_gene[!(detailed_join_gene$Freq_TE_isoform == 1) ,]

write.table(detailed_join_gene,file="../ParasiTE_output/Results/Annotations_chimerick_TE_events/Results_transcripts_TEs_all_events.tab",sep = "\t",row.names=F)

########################### GENE-TE alternative splicing events


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

method_info$method <- as.character(gsub('[()]','',method_info$method))
method_info$method <- as.character(gsub('c','',method_info$method))
method_info$method <- as.character(gsub('"','',method_info$method))
method_info$method <- as.character(gsub('<>','',method_info$method))
method_info$method <- as.character(gsub(':','&',method_info$method))

method_info$method <- sapply(method_info$method, function(x) paste(unique(unlist(str_split(x,", "))), collapse = ", "))

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

overview_gene_TE <- overview_gene_TE[!duplicated(overview_gene_TE), ]

###### Remove FP rules

# Intergenic TEs showing ATP & AS (RI, SE, A5SS, A3SS)
overview_gene_TE2 <- overview_gene_TE
overview_gene_TE$Alternative_splicing <- ifelse(overview_gene_TE$TE_localisation=="Intergenic",NA, overview_gene_TE2$Alternative_splicing )

# Remove TE frequences = 1
# overview_gene_TE<-overview_gene_TE[!(overview_gene_TE$Freq_TE_isoform == 1) ,]

write.table(overview_gene_TE,file="../ParasiTE_output/Results/Annotations_chimerick_TE_events/Result_genes_TEs_all_events.tab",sep = "\t",row.names=F)

########################### Separate Gene-TE alternative & no alternative transcript events

################## alternative transcript events

overview_gene_TE_alternative <- overview_gene_TE[!(is.na(overview_gene_TE$Alternative_transcription) & is.na(overview_gene_TE$Alternative_splicing)) ,]

write.table(overview_gene_TE_alternative,file="../ParasiTE_output/Results/Annotations_chimerick_TE_events/Result_genes_TEs_alternative_events.tab",sep = "\t",row.names=F)

########################### Separate Transcripts-TEs alternative & no alternative transcript events

detailed_join_gene_alternative <- detailed_join_gene[!(is.na(detailed_join_gene$Alternative_transcription) & is.na(detailed_join_gene$Alternative_splicing)) ,]

write.table(detailed_join_gene_alternative,file="../ParasiTE_output/Results/Annotations_chimerick_TE_events/Result_transcripts_TEs_alternative_events.tab",sep = "\t",row.names=F)


print("finished")
