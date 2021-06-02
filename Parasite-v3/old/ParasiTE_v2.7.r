#!/usr/bin/env Rscript
library("optparse")
library(data.table) ## v 1.9.6+ 
library(splitstackshape)

##ParasiTE has been wrote by Jeremy Berthelier, 2021, Plant Epigenetics Unit, Okinawa Institute of Science and Technology

###ParasiTE main launcher script 

##1/options and help message

option_list = list(
    make_option(c("-T", "--transposons"), type="character", default=NULL, 
              help="Annotation of transposable elements (.gff3)", metavar="character"),
    make_option(c("-L", "--transposonsgenes"), type="character", default=NULL, 
              help="Annotation of gene-like transposable elements (.gff3)", metavar="character"),
    make_option(c("-G", "--genes"), type="character", default=NULL, 
              help="Gene model annotation contening gene and exon positions (.gff3)", metavar="character"),   
    make_option(c("-R", "--transcripts"), type="character", default=NULL, 
              help="Transcript annotation obtain by Stringtie, contening transcript and exon positions (.gff3)", metavar="character"),  
    make_option(c("-C", "--cds"), type="character", default=NULL, 
              help="annotation of CDC (.gff3)", metavar="character"),  
    make_option(c("-f", "--five"), type="character", default=NULL, 
              help="annotation of five primer utr (.gff3)", metavar="character"),    
    make_option(c("-t", "--three"), type="character", default=NULL, 
              help="annotation of three prime utr (.gff3)", metavar="character"),    
    make_option(c("-F", "--Tfalsepositive"), type="character", default=0.8, 
              help="% of a TEs in a genes to be considered as false positive, [default= %default] (80%)", metavar="number"), 
    make_option(c("-I", "--Tintragenic"), type="character", default=0.8, 
              help="Proportion of TEs that overlap a gene to be view as intragenic [default= %default]", metavar="number"),
    make_option(c("-i", "--Tintron"), type="character", default=0.1,
              help="Min proportion of an exons that a TEs partialy overlap to be view as partial exonic, othervise it is view as intronic , [default= %default] (0.1%)", metavar="number"),
    make_option(c("-e", "--Texon1"), type="character", default=0.8, 
              help="Minimum proportion of a TEs in an exon to be considered as exonic, [default= %default] (50%)", metavar="number"),       
    make_option(c("-E", "--Texon2"), type="character", default=0.8, 
              help="Minimum proportion of a TEs in an exon to be considered as fragmented exonic, [default= %default] (50%)", metavar="number"), 
    make_option(c("-n", "--Tneighbor"), type="character", default=2000, 
              help="Maximum distance between a gene and a neighbor TE, [default= %default] (2000 bp)", metavar="number"),
    make_option(c("-X", "--MinLexons"), type="character", default=10, 
              help="Min length of exons',  [default= %default] (0)", metavar="number"),
    make_option(c("-m", "--MinLtransposons"), type="character", default=0, 
              help="Min length of transposons',  [default= %default] (0)", metavar="number"),
    make_option(c("-M", "--MaxLtransposons"), type="character", default=50000, 
              help="Maximun length of transposons',  [default= %default] (0)", metavar="number")
              
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$transposons)){
  print_help(opt_parser)
  stop("An annotation file of transposable elements must be supplied (option -T) ", call.=FALSE)
}

if (is.null(opt$genes)){
  print_help(opt_parser)
  stop("A gene model annotation file must be supplied (option -G)", call.=FALSE)
}

if (is.null(opt$transcript)){
  print_help(opt_parser)
  stop("A transcript annotation file must be supplied (option -R)", call.=FALSE)
}

## before to start

initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.basename <- dirname(script.name)
#print(script.basename)
#other.name <- file.path(script.basename, "ParasiTE.R")
#print(paste("Sourcing",other.name,"from",script.name))
#source(other.name)

setwd(script.basename)

if (file.exists("../ParasiTE_output"))
system(paste("rm -Rf ../ParasiTE_output"))

system(paste("mkdir ../ParasiTE_output"))
system(paste("mkdir ../ParasiTE_output/Inputs"))
system(paste("mkdir ../ParasiTE_output/Results"))
system(paste("mkdir ../ParasiTE_output/Results/Plots"))
system(paste("mkdir ../ParasiTE_output/Results/Annotations_TEs"))
system(paste("mkdir ../ParasiTE_output/Results/Annotations_TEs_events"))
system(paste("mkdir ../ParasiTE_output/Results/Annotations_chimerick_TE_events"))


###STEP 0 loading the files

#######1.1 load of the gene annotation

#Grep the "gene annotation" from the gene model file, then convert and sort the file
system(paste("awk -F'\t' '$3~/gene/'", opt$genes, "> ../ParasiTE_output/gene_annotation.gff3"))
system(paste("convert2bed --input=GFF < ../ParasiTE_output/gene_annotation.gff3 > ../ParasiTE_output/pre-gene_annotation.bed && 
bedtools sort -i ../ParasiTE_output/pre-gene_annotation.bed > ../ParasiTE_output/gene_annotation.bed &&
rm ../ParasiTE_output/pre-gene_annotation.bed ../ParasiTE_output/gene_annotation.gff3"))

#######1.2 load of the TE annotation

#Convert and sort the TE annotation file
system(paste("cp ",opt$transposons, "../ParasiTE_output/original_TE_annotation.gff3"))

#Option to remove shorts and/or huge transposons that can be false positives
system(paste("awk -F'\t' '{$10 = ($5-$4)+1} 1' OFS='\t' ../ParasiTE_output/original_TE_annotation.gff3 > ../ParasiTE_output/pre-TE_annotation.gff3 && awk -F'\t' '( ($10 >= ",opt$MinLtransposons, ") && ( $10 <= ",opt$MaxLtransposons," ) ) {print $1,$2,$3,$4,$5,$6,$7,$8,$9}' OFS='\t' ../ParasiTE_output/pre-TE_annotation.gff3  > ../ParasiTE_output/TE_annotation.gff3"))

system(paste("convert2bed --input=GFF < ../ParasiTE_output/TE_annotation.gff3 > ../ParasiTE_output/pre-transposons_annotation.bed && 
bedtools sort -i ../ParasiTE_output/pre-transposons_annotation.bed > ../ParasiTE_output/transposons_annotation.bed"))

#######1.3 (optional) load of the TE gene-like annotation

#Convert and sort the TE annotation file
system(paste("cp ",opt$transposonsgenes, "../ParasiTE_output/TE-genes-like_annotation.gff3"))
system(paste("convert2bed --input=GFF < ../ParasiTE_output/TE-genes-like_annotation.gff3 > ../ParasiTE_output/pre-TE-genes-like_annotation.bed && 
bedtools sort -i ../ParasiTE_output/pre-TE-genes-like_annotation.bed > ../ParasiTE_output/TE-genes-like_annotation.bed"))

#######1.4 load of the de novo transcript annotation (generated by Stringtie)

#2.1.1) Grep the transcript annotation convert and sort in bed file
system(paste("cp ",opt$transcripts," ../ParasiTE_output/Stringtie_annotation.gff3"))
system(paste("awk -F'\t' '$3~/transcript/'", opt$transcripts, "> ../ParasiTE_output/transcript_annotation.gff3"))
system(paste("convert2bed --input=GFF < ../ParasiTE_output/transcript_annotation.gff3 > ../ParasiTE_output/pre-transcript_annotation.bed && 
bedtools sort -i ../ParasiTE_output/pre-transcript_annotation.bed > ../ParasiTE_output/transcript_annotation.bed"))

#2.1.2) Grep exon annotation from the transcript annotation file, convert and sort in bed file
system(paste("awk -F'\t' '$3~/exon/'", opt$transcripts, "> ../ParasiTE_output/exon_transcript_annotation.gff3")) 

#convertion
system(paste("convert2bed --input=GFF < ../ParasiTE_output/exon_transcript_annotation.gff3 > ../ParasiTE_output/pre-exon_transcript_annotation.bed &&
bedtools sort -i ../ParasiTE_output/pre-exon_transcript_annotation.bed > ../ParasiTE_output/exon_transcript_annotation.bed")) 

######1.5 (optional) load of the CDS, 5prime UTR, 3prime UTR annotations

system(paste("cp ",opt$cds," ../ParasiTE_output/cds_annotation.gff3"))
system(paste("convert2bed --input=GFF < ../ParasiTE_output/cds_annotation.gff3 > ../ParasiTE_output/pre-cds_annotation.bed && 
bedtools sort -i ../ParasiTE_output/pre-cds_annotation.bed > ../ParasiTE_output/cds_annotation.bed"))

system(paste("cp ",opt$five," ../ParasiTE_output/five_annotation.gff3"))
system(paste("convert2bed --input=GFF < ../ParasiTE_output/five_annotation.gff3 > ../ParasiTE_output/pre-five_annotation.bed && 
bedtools sort -i ../ParasiTE_output/pre-five_annotation.bed > ../ParasiTE_output/five_annotation.bed"))

system(paste("cp ",opt$three," ../ParasiTE_output/three_annotation.gff3"))
system(paste("convert2bed --input=GFF < ../ParasiTE_output/three_annotation.gff3 > ../ParasiTE_output/pre-three_annotation.bed && 
bedtools sort -i ../ParasiTE_output/pre-three_annotation.bed > ../ParasiTE_output/three_annotation.bed"))

###STEP 1 Remove potential false positives of TE gene-like and the corresponding transcripts

#Round1 remove gene-like transposable elements isoforms
#Search for isoforms that cover at least XX% of a gene-like TEs annotation and remove the genes because it is considered as false positive
#remove transcript that match with TE like genes
#Search for TEs that cover at least XX% of a genes and remove the genes because it is considered as false positive
system(paste("bedtools intersect -nonamecheck -f ", opt$Tfalsepositive," -nonamecheck -wa -a ../ParasiTE_output/gene_annotation.bed -b ../ParasiTE_output/transposons_annotation.bed  > ../ParasiTE_output/gene_FP.bed &&
cat ../ParasiTE_output/gene_FP.bed ../ParasiTE_output/gene_annotation.bed > ../ParasiTE_output/All_genes_and_fp_genes.bed &&
sort ../ParasiTE_output/All_genes_and_fp_genes.bed | uniq -u > ../ParasiTE_output/genes_wo_fp.uniq.bed &&
bedtools sort -i ../ParasiTE_output/genes_wo_fp.uniq.bed > ../ParasiTE_output/genes_wo_fp.uniq.sorted.bed")) 

#remove transcript that match with TE like genes
system(paste("bedtools intersect -nonamecheck -f ", opt$Tfalsepositive," -wo -a ../ParasiTE_output/gene_FP.bed -b ../ParasiTE_output/transcript_annotation.bed  > ../ParasiTE_output/gene_transcript_FP.bed "))
xxxx <- "awk 'BEGIN{FS=OFS=\"\t\"} {print $11,$12,$13,$14,$15,$16,$17,$18,$19,$20}' ../ParasiTE_output/gene_transcript_FP.bed > ../ParasiTE_output/transcript_FP1.bed"
system(paste(xxxx))

#Search for TEs that cover at least XX% of a TE-like genes and remove it
system(paste("bedtools intersect -nonamecheck -F ", opt$Tfalsepositive," -wo -a ../ParasiTE_output/TE-genes-like_annotation.bed -b ../ParasiTE_output/transcript_annotation.bed  > ../ParasiTE_output/gene_transcript_FP.bed "))
xxxx <- "awk 'BEGIN{FS=OFS=\"\t\"} {print $11,$12,$13,$14,$15,$16,$17,$18,$19,$20}' ../ParasiTE_output/gene_transcript_FP.bed > ../ParasiTE_output/transcript_FP2.bed"
system(paste(xxxx))

system(paste("cat ../ParasiTE_output/transcript_FP1.bed ../ParasiTE_output/transcript_FP2.bed > ../ParasiTE_output/transcript_FP.bed &&
sort ../ParasiTE_output/transcript_FP.bed | uniq > ../ParasiTE_output/transcript_FP.uniq.bed &&
bedtools sort -i ../ParasiTE_output/transcript_FP.uniq.bed > ../ParasiTE_output/transcript_FP.uniq.sorted.bed"))

#Found transcripts that match with te-like gene or that are consisderate as False positifs, i retrieve the "STRG.XX" number and I remove the corresponding transcript and exons 
transcript_FP <-  read.table("../ParasiTE_output/transcript_FP.uniq.sorted.bed", header = FALSE, sep = '\t', check.names=FALSE)
colnames(transcript_FP)<-c("chromosome","start", "end", "poin", "other", "strand","tool", "type", "point","id")
transcript_FP <- cSplit(transcript_FP,"id",";")
transcript_FP <- setDT(transcript_FP)[, paste0("id_2", 1:2) := tstrsplit(transcript_FP$id_2, " ")]
transcript_FP <- cSplit(transcript_FP,"id_22",".")
FP_list <- paste(transcript_FP$id_22_1,".",transcript_FP$id_22_2)
FP_list <- as.character(gsub(' ','',FP_list))
write.table(FP_list,file="../ParasiTE_output/FP_list.txt",sep = "\t",row.names=F,col.names=F)
transcript_annotation <-  read.table("../ParasiTE_output/transcript_annotation.bed", header = FALSE, sep = '\t', check.names=FALSE)
 
system(paste("grep -vFwf ../ParasiTE_output/FP_list.txt ../ParasiTE_output/transcript_annotation.bed > ../ParasiTE_output/transcript_annotation_wo_FP.bed"))
system(paste("grep -vFwf ../ParasiTE_output/FP_list.txt ../ParasiTE_output/exon_transcript_annotation.bed > ../ParasiTE_output/exon_transcript_annotation_wo_FP.bed"))

###STEP 2 Classification of TEs regarding to their location
######2.1.Intergenic x Intragenic

#1.1) Identify intragenic TEs: 
#1.1.1) Intersect detect TEs located 80% in reference genes annotation, remove redundance and sort
system(paste("bedtools intersect -f ", opt$Tintragenic," -nonamecheck -wa -a ../ParasiTE_output/transposons_annotation.bed -b ../ParasiTE_output/genes_wo_fp.uniq.sorted.bed > ../ParasiTE_output/intragenicTEs.bed &&
sort ../ParasiTE_output/intragenicTEs.bed | uniq > ../ParasiTE_output/intragenicTEs.uniq.bed &&
bedtools sort -i ../ParasiTE_output/intragenicTEs.uniq.bed > ../ParasiTE_output/intragenicTEs.uniq.sorted.bed "))
#1.1.1) Create a new annotation without intragenic TEs (remove intragenic TEs from the total TE annotation)
system(paste("cat ../ParasiTE_output/transposons_annotation.bed ../ParasiTE_output/intragenicTEs.uniq.sorted.bed > ../ParasiTE_output/TEannotation_and_intragenicTEs.bed && 
sort ../ParasiTE_output/TEannotation_and_intragenicTEs.bed | uniq -u > ../ParasiTE_output/pre-TEannotation_wo_intragenicTEs.bed &&
bedtools sort -i ../ParasiTE_output/pre-TEannotation_wo_intragenicTEs.bed > ../ParasiTE_output/TEannotation_wo_intragenicTEs.bed &&
rm ../ParasiTE_output/pre-TEannotation_wo_intragenicTEs.bed"))

#1.3) Identify intergenic TEs
#1.3.2)label Intragenic TEs
system(paste("sed -i 's/\\S\\+/'Intragenic'/7' ../ParasiTE_output/intragenicTEs.uniq.sorted.bed &&
mv ../ParasiTE_output/intragenicTEs.uniq.sorted.bed ../ParasiTE_output/IntragenicTEs.labeled.bed &&
cp ../ParasiTE_output/IntragenicTEs.labeled.bed ../ParasiTE_output/Results/Annotations_TEs/All_Intragenic_TEs.bed"))

#1.3.2)label intergenic TEs
system(paste("sed -i 's/\\S\\+/'Intergenic'/7' ../ParasiTE_output/TEannotation_wo_intragenicTEs.bed  &&
mv ../ParasiTE_output/TEannotation_wo_intragenicTEs.bed ../ParasiTE_output/IntergenicTEs.labeled.bed &&
cp ../ParasiTE_output/IntergenicTEs.labeled.bed ../ParasiTE_output/Results/Annotations_TEs/All_Intergenic_TEs.bed"))

#1.4) Create a labeled final annotation of total TEs
system(paste("cat ../ParasiTE_output/IntragenicTEs.labeled.bed ../ParasiTE_output/IntergenicTEs.labeled.bed > ../ParasiTE_output/Total_TEs_annotation.labeled.bed &&
cp ../ParasiTE_output/Total_TEs_annotation.labeled.bed ../ParasiTE_output/Results/Annotations_TEs/All_TEs.bed"))

####STEP3 Selection of neighbor intergenic TEs
###### 3.1 Find the three first intergenic TEs close to genes

#round1 with overlap
#upstream
system(paste("bedtools closest -nonamecheck -D ref -id -a ../ParasiTE_output/genes_wo_fp.uniq.sorted.bed -b ../ParasiTE_output/IntergenicTEs.labeled.bed > ../ParasiTE_output/neighbor_intergenic_TEs_up_round1.bed"))

#The tool bedtool add line with "-1 -1" when no result were found. Thus I need to remove those lines.
ii <- "awk -F'\t' '$11 !=\".\"' ../ParasiTE_output/neighbor_intergenic_TEs_up_round1.bed > ../ParasiTE_output/neighbor_intergenic_TEs_up_round1.clean.bed " 
system(paste(ii))

#downstream
system(paste("bedtools closest -nonamecheck -D ref -iu -a ../ParasiTE_output/genes_wo_fp.uniq.sorted.bed -b ../ParasiTE_output/IntergenicTEs.labeled.bed > ../ParasiTE_output/neighbor_intergenic_TEs_down_round1.bed"))

#The tool bedtool add line with "-1 -1" when no result were found. Thus I need to remove those lines.
kk <- "awk -F'\t' '$11 !=\".\"' ../ParasiTE_output/neighbor_intergenic_TEs_down_round1.bed > ../ParasiTE_output/neighbor_intergenic_TEs_down_round1.clean.bed " 
system(paste(kk))

system(paste("cat ../ParasiTE_output/neighbor_intergenic_TEs_up_round1.clean.bed ../ParasiTE_output/neighbor_intergenic_TEs_down_round1.clean.bed > ../ParasiTE_output/neighbor_intergenic_TEs_up_down_round1.bed &&

sort ../ParasiTE_output/neighbor_intergenic_TEs_up_down_round1.bed | uniq > ../ParasiTE_output/neighbor_intergenic_TEs_up_down_round1.uniq.bed &&
bedtools sort -i ../ParasiTE_output/neighbor_intergenic_TEs_up_down_round1.uniq.bed > ../ParasiTE_output/neighbor_intergenic_TEs_up_down_round1.uniq.sorted.bed &&
awk -F'\t' '( ($21 >= -",opt$Tneighbor, ") && ( $21 <= ",opt$Tneighbor," ) )' ../ParasiTE_output/neighbor_intergenic_TEs_up_down_round1.uniq.sorted.bed > ../ParasiTE_output/neighbor_intergenic_TEs_up_down_round1.uniq.sorted.lenght.bed"))

dd <- "awk 'BEGIN{FS=OFS=\"\t\"} {print $11,$12,$13,$14,$15,$16,$17,$18,$19,$20}' ../ParasiTE_output/neighbor_intergenic_TEs_up_down_round1.uniq.sorted.lenght.bed > ../ParasiTE_output/neighbor_intergenic_TEs_up_down_round1.final.bed"
system(paste(dd))
system(paste("sort ../ParasiTE_output/neighbor_intergenic_TEs_up_down_round1.final.bed | uniq > ../ParasiTE_output/neighbor_intergenic_TEs_up_down_round1.final.uniq.bed &&
bedtools sort -i ../ParasiTE_output/neighbor_intergenic_TEs_up_down_round1.final.uniq.bed > ../ParasiTE_output/neighbor_intergenic_TEs_up_down_round1.final.uniq.sorted.bed"))

#remove round1 TEs from intergenic all annotation
system(paste("cat ../ParasiTE_output/IntergenicTEs.labeled.bed ../ParasiTE_output/neighbor_intergenic_TEs_up_down_round1.final.uniq.sorted.bed > ../ParasiTE_output/IntergenicTEs.labeled.with_round1.bed && 
sort ../ParasiTE_output/IntergenicTEs.labeled.with_round1.bed | uniq -u > ../ParasiTE_output/IntergenicTEs.labeled.wo_round1.bed &&
bedtools sort -i ../ParasiTE_output/IntergenicTEs.labeled.wo_round1.bed > ../ParasiTE_output/IntergenicTEs.labeled.wo_round1.sorted.bed"))

#round2 without overlap
#upstream

system(paste("bedtools closest -nonamecheck -D ref -id -io -a ../ParasiTE_output/genes_wo_fp.uniq.sorted.bed -b ../ParasiTE_output/IntergenicTEs.labeled.wo_round1.sorted.bed > ../ParasiTE_output/neighbor_intergenic_TEs_up_round2.bed"))

#The tool bedtool add line with "-1 -1" when no result were found. Thus I need to remove those lines.
ll <- "awk -F'\t' '$11 !=\".\"' ../ParasiTE_output/neighbor_intergenic_TEs_up_round2.bed > ../ParasiTE_output/neighbor_intergenic_TEs_up_round2.clean.bed " 
system(paste(ll))

#downstream
system(paste("bedtools closest -nonamecheck -D ref -iu -io -a ../ParasiTE_output/genes_wo_fp.uniq.sorted.bed -b ../ParasiTE_output/IntergenicTEs.labeled.wo_round1.sorted.bed > ../ParasiTE_output/neighbor_intergenic_TEs_down_round2.bed"))

#The tool bedtool add line with "-1 -1" when no result were found. Thus I need to remove those lines.
ll <- "awk -F'\t' '$11 !=\".\"' ../ParasiTE_output/neighbor_intergenic_TEs_down_round2.bed > ../ParasiTE_output/neighbor_intergenic_TEs_down_round2.clean.bed " 
system(paste(ll))

system(paste("cat ../ParasiTE_output/neighbor_intergenic_TEs_up_round2.clean.bed ../ParasiTE_output/neighbor_intergenic_TEs_down_round2.clean.bed > ../ParasiTE_output/neighbor_intergenic_TEs_up_down_round2.bed &&
sort ../ParasiTE_output/neighbor_intergenic_TEs_up_down_round2.bed | uniq > ../ParasiTE_output/neighbor_intergenic_TEs_up_down_round2.uniq.bed &&
bedtools sort -i ../ParasiTE_output/neighbor_intergenic_TEs_up_down_round2.uniq.bed > ../ParasiTE_output/neighbor_intergenic_TEs_up_down_round2.uniq.sorted.bed &&
awk -F'\t' '( ($21 >= -",opt$Tneighbor, ") && ( $21 <= ",opt$Tneighbor," ) )' ../ParasiTE_output/neighbor_intergenic_TEs_up_down_round2.uniq.sorted.bed > ../ParasiTE_output/neighbor_intergenic_TEs_up_down_round2.uniq.sorted.lenght.bed"))

dd <- "awk 'BEGIN{FS=OFS=\"\t\"} {print $11,$12,$13,$14,$15,$16,$17,$18,$19,$20}' ../ParasiTE_output/neighbor_intergenic_TEs_up_down_round2.uniq.sorted.lenght.bed > ../ParasiTE_output/neighbor_intergenic_TEs_up_down_round2.final.bed"
system(paste(dd))
system(paste("sort ../ParasiTE_output/neighbor_intergenic_TEs_up_down_round2.final.bed | uniq > ../ParasiTE_output/neighbor_intergenic_TEs_up_down_round2.final.uniq.bed &&
bedtools sort -i ../ParasiTE_output/neighbor_intergenic_TEs_up_down_round2.final.uniq.bed > ../ParasiTE_output/neighbor_intergenic_TEs_up_down_round2.final.uniq.sorted.bed"))

#remove round2 TEs from intergenic all annotation
system(paste("cat ../ParasiTE_output/IntergenicTEs.labeled.wo_round1.sorted.bed ../ParasiTE_output/neighbor_intergenic_TEs_up_down_round2.final.uniq.sorted.bed > ../ParasiTE_output/IntergenicTEs.labeled.wo_round1.with_round2.bed && 
sort ../ParasiTE_output/IntergenicTEs.labeled.wo_round1.with_round2.bed | uniq -u > ../ParasiTE_output/IntergenicTEs.labeled.wo_round1.wo_round2.bed &&
bedtools sort -i ../ParasiTE_output/IntergenicTEs.labeled.wo_round1.wo_round2.bed > ../ParasiTE_output/IntergenicTEs.labeled.wo_round1.wo_round2.sorted.bed"))

#round3 wo overlap
#upstream

system(paste("bedtools closest -nonamecheck -D ref -id -io -a ../ParasiTE_output/genes_wo_fp.uniq.sorted.bed -b ../ParasiTE_output/IntergenicTEs.labeled.wo_round1.wo_round2.sorted.bed > ../ParasiTE_output/neighbor_intergenic_TEs_up_round3.bed"))

#The tool bedtool add line with "-1 -1" when no result were found. Thus I need to remove those lines.
ll <- "awk -F'\t' '$11 !=\".\"' ../ParasiTE_output/neighbor_intergenic_TEs_up_round3.bed > ../ParasiTE_output/neighbor_intergenic_TEs_up_round3.clean.bed " 
system(paste(ll))

#downstream
system(paste("bedtools closest -nonamecheck -D ref -iu -io -a ../ParasiTE_output/genes_wo_fp.uniq.sorted.bed -b ../ParasiTE_output/IntergenicTEs.labeled.wo_round1.wo_round2.sorted.bed > ../ParasiTE_output/neighbor_intergenic_TEs_down_round3.bed"))

#The tool bedtool add line with "-1 -1" when no result were found. Thus I need to remove those lines.
ll <- "awk -F'\t' '$11 !=\".\"' ../ParasiTE_output/neighbor_intergenic_TEs_down_round3.bed > ../ParasiTE_output/neighbor_intergenic_TEs_down_round3.clean.bed " 
system(paste(ll))

system(paste("cat ../ParasiTE_output/neighbor_intergenic_TEs_up_round3.clean.bed ../ParasiTE_output/neighbor_intergenic_TEs_down_round3.clean.bed > ../ParasiTE_output/neighbor_intergenic_TEs_up_down_round3.bed &&

sort ../ParasiTE_output/neighbor_intergenic_TEs_up_down_round3.bed | uniq > ../ParasiTE_output/neighbor_intergenic_TEs_up_down_round3.uniq.bed &&
bedtools sort -i ../ParasiTE_output/neighbor_intergenic_TEs_up_down_round3.uniq.bed > ../ParasiTE_output/neighbor_intergenic_TEs_up_down_round3.uniq.sorted.bed &&
awk -F'\t' '( ($21 >= -",opt$Tneighbor, ") && ( $21 <= ",opt$Tneighbor," ) )' ../ParasiTE_output/neighbor_intergenic_TEs_up_down_round3.uniq.sorted.bed > ../ParasiTE_output/neighbor_intergenic_TEs_up_down_round3.uniq.sorted.lenght.bed"))

dd <- "awk 'BEGIN{FS=OFS=\"\t\"} {print $11,$12,$13,$14,$15,$16,$17,$18,$19,$20}' ../ParasiTE_output/neighbor_intergenic_TEs_up_down_round3.uniq.sorted.lenght.bed > ../ParasiTE_output/neighbor_intergenic_TEs_up_down_round3.final.bed"
system(paste(dd))
system(paste("sort ../ParasiTE_output/neighbor_intergenic_TEs_up_down_round3.final.bed | uniq > ../ParasiTE_output/neighbor_intergenic_TEs_up_down_round3.final.uniq.bed &&
bedtools sort -i ../ParasiTE_output/neighbor_intergenic_TEs_up_down_round3.final.uniq.bed > ../ParasiTE_output/neighbor_intergenic_TEs_up_down_round3.final.uniq.sorted.bed"))

#cat neihbor TEs found at each round and remove duplicates
system(paste("cat ../ParasiTE_output/neighbor_intergenic_TEs_up_down_round1.final.uniq.sorted.bed ../ParasiTE_output/neighbor_intergenic_TEs_up_down_round2.final.uniq.sorted.bed ../ParasiTE_output/neighbor_intergenic_TEs_up_down_round3.final.uniq.sorted.bed > ../ParasiTE_output/neighbor_intergenic_TEs_up_down_round_1_round2_round3.bed &&

sort ../ParasiTE_output/neighbor_intergenic_TEs_up_down_round_1_round2_round3.bed | uniq > ../ParasiTE_output/neighbor_intergenic_TEs_up_down_round_1_round2_round3.uniq.bed &&
bedtools sort -i ../ParasiTE_output/neighbor_intergenic_TEs_up_down_round_1_round2_round3.uniq.bed > ../ParasiTE_output/neighbor_intergenic_TEs_up_down_round_1_round2_round3.uniq.sorted.bed"))

###STEP4 Detection of exonic TE candidates

###### 4.1 Prepare dataset

#Create a dataset for the analyses
system(paste("cat ../ParasiTE_output/IntragenicTEs.labeled.bed ../ParasiTE_output/neighbor_intergenic_TEs_up_down_round_1_round2_round3.uniq.sorted.bed > ../ParasiTE_output/Intragenic_TEs_Neighbor_intergenic.bed"))

#Options to remove short exons, default no min
system(paste("awk -F'\t' '{$11 = ($3-$2)} 1' OFS='\t' ../ParasiTE_output/exon_transcript_annotation_wo_FP.bed > ../ParasiTE_output/pre-exon_transcript_annotation.woshort.bed && 
awk -F'\t' '( $11 >=",opt$MinLexons, ") {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}' OFS='\t' ../ParasiTE_output/pre-exon_transcript_annotation.woshort.bed  > ../ParasiTE_output/exon_transcript_annotation.woshort.bed "))

####### 4.2 Identification Method1 - Unfragmented exonic TEs
#2.2.1) The method 1 search TEs that are covered by at least 80% by an exon, one TEs can be in several transcript isoforms, thus the same TEs can be found several time
system(paste("bedtools intersect -f ", opt$Texon1, " -nonamecheck -wo -a ../ParasiTE_output/Intragenic_TEs_Neighbor_intergenic.bed -b ../ParasiTE_output/exon_transcript_annotation_wo_FP.bed > ../ParasiTE_output/TE_transcripts_exonic_events.bed &&
sort ../ParasiTE_output/TE_transcripts_exonic_events.bed | uniq > ../ParasiTE_output/TE_transcripts_exonic_events.uniq.bed &&
bedtools sort -i ../ParasiTE_output/TE_transcripts_exonic_events.uniq.bed > ../ParasiTE_output/TE_transcripts_exonic_events.sorted.bed"))

#2.2.2) Take every TEs involved in exonic events, then remove gene annotations and remove redundance of TEs, keep TE length
ttt <- "awk 'BEGIN{FS=OFS=\"\t\"} {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}' ../ParasiTE_output/TE_transcripts_exonic_events.sorted.bed > ../ParasiTE_output/TE_transcript_exonic.bed "
system(paste(ttt))
system(paste("sort ../ParasiTE_output/TE_transcript_exonic.bed | uniq > ../ParasiTE_output/TE_transcript_exonic.uniq.bed &&	    
bedtools sort -i ../ParasiTE_output/TE_transcript_exonic.uniq.bed > ../ParasiTE_output/TE_transcript_exonic.uniq.sorted.bed"))

#### 4.3 Identification Method2 - Fragemented exonic TEs
#The method 2 search TEs that covered exon by more than 80% the total length of the exon, those are also potential intragenic exonic TEs
system(paste("bedtools intersect -F ", opt$Texon2, " -nonamecheck -wo -a ../ParasiTE_output/Intragenic_TEs_Neighbor_intergenic.bed  -b ../ParasiTE_output/exon_transcript_annotation.woshort.bed > ../ParasiTE_output/exonic_ambigous_TEs.bed && sort ../ParasiTE_output/exonic_ambigous_TEs.bed | uniq > ../ParasiTE_output/exonic_ambigous_TEs.uniq.bed &&
bedtools sort -i ../ParasiTE_output/exonic_ambigous_TEs.uniq.bed > ../ParasiTE_output/TE_transcripts_exonic_method2_TEs_events.uniq.sorted.bed"))

#Take every TEs involved in exonic events, then remove gene annotations and remove redundance of TEs
ccc <- "awk 'BEGIN{FS=OFS=\"\t\"} {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}' ../ParasiTE_output/TE_transcripts_exonic_method2_TEs_events.uniq.sorted.bed > ../ParasiTE_output/TEs_transcipt_exonic_ambigous.bed "
system(paste(ccc))
system(paste("sort ../ParasiTE_output/TEs_transcipt_exonic_ambigous.bed | uniq > ../ParasiTE_output/TEs_transcipt_exonic_ambigous.uniq.bed &&	    
bedtools sort -i ../ParasiTE_output/TEs_transcipt_exonic_ambigous.uniq.bed > ../ParasiTE_output/TEs_transcipt_exonic_ambigous.uniq.sorted.bed"))

#### 4.3 Identification Method3 - Partially exonic TEs
#2.2.3 create a new annotation without the full exonic TEs detected with method 1 and method 2
#system(paste("cat ../ParasiTE_output/Intragenic_TEs_Neighbor_intergenic.bed ../ParasiTE_output/TE_transcript_exonic.uniq.sorted.bed ../ParasiTE_output/TEs_transcipt_exonic_ambigous.uniq.sorted.bed > ../ParasiTE_output/#TEannotation_and_exonicTEs.bed && 
#sort ../ParasiTE_output/TEannotation_and_exonicTEs.bed | uniq -u > ../ParasiTE_output/pre-TEannotation_wo_exonicTEs.bed &&
#bedtools sort -i ../ParasiTE_output/pre-TEannotation_wo_exonicTEs.bed > ../ParasiTE_output/TEannotation_wo_exonicTEs.bed"))
#Search TEs that covered exon by more than 80% the total length of the exon, those are potential intragenic_exonic_ambigous TEs

system(paste("bedtools intersect -F ", opt$Tintron," -nonamecheck -wo -a ../ParasiTE_output/Intragenic_TEs_Neighbor_intergenic.bed  -b ../ParasiTE_output/exon_transcript_annotation.woshort.bed > ../ParasiTE_output/partialy_exonic_TEs.bed && 
sort ../ParasiTE_output/partialy_exonic_TEs.bed | uniq > ../ParasiTE_output/partialy_exonic_TEs.uniq.bed &&
bedtools sort -i ../ParasiTE_output/partialy_exonic_TEs.uniq.bed > ../ParasiTE_output/partialy_exonic_TEs.uniq.sorted.bed"))

##### 4.4 Annotation file: Method1 and Method2 - Fragmented and Unfragmented exonic TEs
system(paste("cat ../ParasiTE_output/TE_transcripts_exonic_events.sorted.bed ../ParasiTE_output/TE_transcripts_exonic_method2_TEs_events.uniq.sorted.bed > ../ParasiTE_output/exonic_TEs.bed"))
#All full intragenic exonic identified are retreived and classified as intragenic exonic TEs
system(paste("awk -F'\t' '$7~/Intragenic/' ../ParasiTE_output/exonic_TEs.bed > ../ParasiTE_output/pre-Annotation_intragenic_exonic_TEs.bed"))
oo <- "awk 'BEGIN{FS=OFS=\"\t\"} {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}' ../ParasiTE_output/pre-Annotation_intragenic_exonic_TEs.bed > ../ParasiTE_output/Annotation_intragenic_exonic_TEs.bed"
system(paste(oo))
system(paste("sort ../ParasiTE_output/Annotation_intragenic_exonic_TEs.bed | uniq > ../ParasiTE_output/Annotation_intragenic_exonic_TEs.uniq.bed && 
bedtools sort -i ../ParasiTE_output/Annotation_intragenic_exonic_TEs.uniq.bed > ../ParasiTE_output/Annotation_intragenic_exonic_TEs.uniq.sorted.bed"))
system(paste("mv ../ParasiTE_output/Annotation_intragenic_exonic_TEs.uniq.sorted.bed ../ParasiTE_output/Results/Annotations_TEs/Intragenic_exonic_TEs.bed"))

#####4.5 Annotation file: Method3 for intragenic TEs - Partially exonic
#All partial intragenic exonic identified are retreived and classified as intragenic exonic TEs
system(paste("awk -F'\t' '$7~/Intragenic/' ../ParasiTE_output/partialy_exonic_TEs.uniq.sorted.bed > ../ParasiTE_output/pre-Annotation_intragenic_partial_exonic_TEs.bed"))
oo <- "awk 'BEGIN{FS=OFS=\"\t\"} {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}' ../ParasiTE_output/pre-Annotation_intragenic_partial_exonic_TEs.bed > ../ParasiTE_output/Annotation_intragenic_partial_exonic_TEs.bed"
system(paste(oo))

system(paste("sort ../ParasiTE_output/Annotation_intragenic_partial_exonic_TEs.bed | uniq > ../ParasiTE_output/Annotation_intragenic_partial_exonic_TEs.uniq.bed && 
bedtools sort -i ../ParasiTE_output/Annotation_intragenic_partial_exonic_TEs.uniq.bed > ../ParasiTE_output/Annotation_intragenic_partial_exonic_TEs.uniq.sorted.bed"))
# I need to change the word ambigous to partial
system(paste("mv ../ParasiTE_output/Annotation_intragenic_partial_exonic_TEs.uniq.sorted.bed ../ParasiTE_output/Results/Annotations_TEs/Intragenic_ambigous_TEs.bed"))

#Remove exonic TEs in the intragenic annotation
system(paste("cat ../ParasiTE_output/Results/Annotations_TEs/All_Intragenic_TEs.bed ../ParasiTE_output/Results/Annotations_TEs/Intragenic_exonic_TEs.bed ../ParasiTE_output/Results/Annotations_TEs/Intragenic_ambigous_TEs.bed > ../ParasiTE_output/CAT_Intragenic_and_exonic_TEs.bed "))
system(paste("sort ../ParasiTE_output/CAT_Intragenic_and_exonic_TEs.bed | uniq -u > ../ParasiTE_output/Intragenic_wo_exonic_TEs.bed &&
bedtools sort -i ../ParasiTE_output/Intragenic_wo_exonic_TEs.bed > ../ParasiTE_output/Intragenic_wo_partial_and_full_exonic_TEs.sorted.bed"))

#####4.5 Annotation file: Intronic TEs

#The TEs that overlap more 10% the length of an exon are excuded from the total intragenic TE annotation, thus the TEs left are potential intronic TEs  

system(paste("sort ../ParasiTE_output/Intragenic_wo_partial_and_full_exonic_TEs.sorted.bed | uniq -u > ../ParasiTE_output/TE_intronic.uniq.bed &&
bedtools sort -i ../ParasiTE_output/TE_intronic.uniq.bed > ../ParasiTE_output/TE_intronic.uniq.sorted.bed &&
cp ../ParasiTE_output/TE_intronic.uniq.sorted.bed ../ParasiTE_output/Results/Annotations_TEs/Intragenic_intronic_TEs.bed" ))

#Found every events of intronic TEs, one TEs can be in several gene, thus the same TEs can be found several time

system(paste("bedtools intersect -nonamecheck -wo -a ../ParasiTE_output/Results/Annotations_TEs/Intragenic_intronic_TEs.bed -b ../ParasiTE_output/transcript_annotation.bed > ../ParasiTE_output/TE_intronic_events.bed && 
sort ../ParasiTE_output/TE_intronic_events.bed | uniq > ../ParasiTE_output/TE_intronic_events.uniq.bed &&
bedtools sort -i ../ParasiTE_output/TE_intronic_events.uniq.bed > ../ParasiTE_output/TE_intronic_events.uniq.sorted.bed && 
cp ../ParasiTE_output/TE_intronic_events.uniq.sorted.bed ../ParasiTE_output/Results/Annotations_TEs_events/Intragenic_intronic_TEs_genes.bed"))

#####4.6 Remove redundance for transcripts, one transcript can be etheir detectected by method1 > method2 > method3

#remove the last column of each, because it is slighly differents according to the method
qw <- "awk 'BEGIN{FS=OFS=\"\t\"} {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,0}' ../ParasiTE_output/TE_transcripts_exonic_events.sorted.bed > ../ParasiTE_output/TE_transcripts_exonic_events2.sorted.bed"
system(paste(qw))
qe <- "awk 'BEGIN{FS=OFS=\"\t\"} {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,0}' ../ParasiTE_output/TE_transcripts_exonic_method2_TEs_events.uniq.sorted.bed > ../ParasiTE_output/TE_transcripts_exonic_method2_TEs_events2.uniq.sorted.bed"
system(paste(qe))
qr <- "awk 'BEGIN{FS=OFS=\"\t\"} {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,0}' ../ParasiTE_output/partialy_exonic_TEs.uniq.sorted.bed > ../ParasiTE_output/partialy_exonic_TEs2.uniq.sorted.bed"
system(paste(qr))

#sort files
system(paste("sort ../ParasiTE_output/TE_transcripts_exonic_events2.sorted.bed > ../ParasiTE_output/TE_transcripts_exonic_events.sort.bed &&
sort ../ParasiTE_output/TE_transcripts_exonic_method2_TEs_events2.uniq.sorted.bed > ../ParasiTE_output/TE_transcripts_exonic_method2_TEs_events.uniq.sort.bed &&
sort ../ParasiTE_output/partialy_exonic_TEs2.uniq.sorted.bed > ../ParasiTE_output/partialy_exonic_TEs.uniq.sort.bed"))

#remove method 2 candiates that are in methode 1
system(paste("comm -12 ../ParasiTE_output/TE_transcripts_exonic_events.sort.bed ../ParasiTE_output/TE_transcripts_exonic_method2_TEs_events.uniq.sort.bed > ../ParasiTE_output/comm1_2.bed &&
cat ../ParasiTE_output/TE_transcripts_exonic_method2_TEs_events.uniq.sort.bed ../ParasiTE_output/comm1_2.bed > ../ParasiTE_output/TE_transcripts_exonic_method2_TEs_events_comm1_2.bed &&
sort ../ParasiTE_output/TE_transcripts_exonic_method2_TEs_events_comm1_2.bed | uniq -u > ../ParasiTE_output/TE_transcripts_exonic_method2_TEs_events_wo_comm1.bed &&
bedtools sort -i ../ParasiTE_output/TE_transcripts_exonic_method2_TEs_events_wo_comm1.bed  > ../ParasiTE_output/TE_transcripts_exonic_method2_TEs_events_wo_comm1.sorted.bed "))

#remove method 1 and 2 candidates that are in methode 3
system(paste("cat ../ParasiTE_output/TE_transcripts_exonic_events.sort.bed ../ParasiTE_output/TE_transcripts_exonic_method2_TEs_events.uniq.sort.bed > ../ParasiTE_output/cat_TE_transcripts_exonic_method1_2_TEs_events.bed &&
sort ../ParasiTE_output/cat_TE_transcripts_exonic_method1_2_TEs_events.bed | uniq > ../ParasiTE_output/cat_TE_transcripts_exonic_method1_2_TEs_events.uniq.bed &&
comm -12 ../ParasiTE_output/cat_TE_transcripts_exonic_method1_2_TEs_events.uniq.bed ../ParasiTE_output/partialy_exonic_TEs.uniq.sort.bed > ../ParasiTE_output/comm1_2_3.bed &&
cat ../ParasiTE_output/partialy_exonic_TEs.uniq.sort.bed ../ParasiTE_output/comm1_2_3.bed > ../ParasiTE_output/partialy_exonic_TEs_comm1_2_3.bed &&
sort ../ParasiTE_output/partialy_exonic_TEs_comm1_2_3.bed | uniq -u > ../ParasiTE_output/partialy_exonic_TEs_wo_comm1_2_3.uniq.bed &&
bedtools sort -i ../ParasiTE_output/partialy_exonic_TEs_wo_comm1_2_3.uniq.bed > ../ParasiTE_output/partialy_exonic_TEs_wo_comm1_2_3.uniq.sorted.bed"))

#####4.7 Labelling TEs regarding their exonic proportion
#label colum 8 according to the detection method (1 or 2 or 3)

system(paste("sed -i 's/\\S\\+/'1'/8' ../ParasiTE_output/TE_transcripts_exonic_events2.sorted.bed &&
cp ../ParasiTE_output/TE_transcripts_exonic_events2.sorted.bed ../ParasiTE_output/TE_transcript_exonic.labeled.bed"))
system(paste("mv ../ParasiTE_output/TE_transcript_exonic.labeled.bed ../ParasiTE_output/Results/Annotations_chimerick_TE_events/TEs_exonic.bed"))

system(paste("sed -i 's/\\S\\+/'2'/8' ../ParasiTE_output/TE_transcripts_exonic_method2_TEs_events_wo_comm1.sorted.bed &&
cp ../ParasiTE_output/TE_transcripts_exonic_method2_TEs_events_wo_comm1.sorted.bed ../ParasiTE_output/TE_transcripts_exonic_method2_TEs_events_wo_comm1.sorted.labeled.bed"))
system(paste("mv ../ParasiTE_output/TE_transcripts_exonic_method2_TEs_events_wo_comm1.sorted.labeled.bed ../ParasiTE_output/Results/Annotations_chimerick_TE_events/TEs_exonic_method2.bed"))

system(paste("cat ../ParasiTE_output/Results/Annotations_chimerick_TE_events/TEs_exonic.bed ../ParasiTE_output/Results/Annotations_chimerick_TE_events/TEs_exonic_method2.bed > ../ParasiTE_output/Results/Annotations_chimerick_TE_events/Full_exonic_TEs.bed"))

system(paste("sed -i 's/\\S\\+/'3'/8' ../ParasiTE_output/partialy_exonic_TEs_wo_comm1_2_3.uniq.sorted.bed &&
cp ../ParasiTE_output/partialy_exonic_TEs_wo_comm1_2_3.uniq.sorted.bed ../ParasiTE_output/TE_transcript_partially_exonic.labeled.bed"))

system(paste("mv ../ParasiTE_output/TE_transcript_partially_exonic.labeled.bed ../ParasiTE_output/Results/Annotations_chimerick_TE_events/Partial_exonic_TEs.bed"))

system(paste("sed -i 's/\\S\\+/'0'/8' ../ParasiTE_output/Results/Annotations_TEs_events/Intragenic_intronic_TEs_genes.bed &&
cp ../ParasiTE_output/Results/Annotations_TEs_events/Intragenic_intronic_TEs_genes.bed ../ParasiTE_output/Intragenic_intronic_TEs_genes.labeled.bed"))
system(paste("mv  ../ParasiTE_output/Intragenic_intronic_TEs_genes.labeled.bed ../ParasiTE_output/Results/Annotations_chimerick_TE_events/TEs_intronic.bed"))

system(paste("cat ../ParasiTE_output/Results/Annotations_chimerick_TE_events/Full_exonic_TEs.bed ../ParasiTE_output/Results/Annotations_chimerick_TE_events/Partial_exonic_TEs.bed  > ../ParasiTE_output/Results/Annotations_chimerick_TE_events/Candidate_TEs.bed"))

## Launch the next scripts

#### Extraction of slicing_site
system(paste("../Splicing_sites/extract_splice_sites.py ", opt$transcripts," > ../ParasiTE_output/splicing_sites.gff"))
#### Extraction of classification information
system(paste("sh ../Extraction/Extract-intragenic-TE.sh"))
#### Generate plots
system(paste("sh ../Rscript/launcher-Rscript.sh"))

#### Final script to indentify chimeric TE-related transcripts and get outputs
system(paste("Rscript ./manage-exonic-TEs_v2.7.R"))

