#!/usr/bin/env Rscript
library("optparse")
library(data.table) ## v 1.9.6+ 
library(splitstackshape)

#ParasiTE has been wrote by Jeremy Berthelier, 2021, Plant Epigenetics Unit, Okinawa Institute of Science and Technology


###ParasiTE main launcher script 
##1/options and help message

option_list = list(
    make_option(c("-T", "--transposons"), type="character", default=NULL, 
              help="[required] Pathway to the annotation of transposable elements (.gff3)", metavar="character"),
    make_option(c("-G", "--genes"), type="character", default=NULL, 
              help="[required] Pathway to the gene model annotation (genes/exons) (.gff3)", metavar="character"),   
    make_option(c("-R", "--transcripts"), type="character", default=NULL, 
              help="[required] Pathway to the transcriptome annotation.
                    No need -P option if you use a Stringtie annotation obtained with short read, or long read with the option Assembly.
                    Use '-P SR' if you use a Stringtie annotation obtained with long reads and the option Raw.
                    Use '-P SM' if you want to use your own annotation that have been propery parsed. (.gff3)", metavar="character"),
    make_option(c("-P", "--Pmode"), type="character", default="0", 
              help="(optional) Add this option if: 
                    -the transcriptome annotation have been generated with long-reads by Stringtie using Raw mode (use -P SR).
                    -you are using a custom transcriptome annotation that have been properly formated (use -P SM)',  [default= %default] (0)", metavar="number")              
    make_option(c("-L", "--transposonsgenes"), type="character", default=NULL, 
              help="(optional) Pathway to the annotation of gene-like transposable elements (.gff3)", metavar="character"),
    make_option(c("-C", "--cds"), type="character", default=NULL, 
              help="(optional) Pathway to theannotation of CDS (.gff3)", metavar="character"),  
    make_option(c("-f", "--five"), type="character", default=NULL, 
              help="(optional) Pathway to theannotation of 5`UTR (.gff3)", metavar="character"),    
    make_option(c("-t", "--three"), type="character", default=NULL, 
              help="(optional) Pathway to theannotation of 3`UTR (.gff3)", metavar="character"),    
    make_option(c("-F", "--Tfalsepositive"), type="character", default=0.8, 
              help="Percentage of a TEs that overlap a gene annotation to be considered as false positive, [default= %default] (80%)", metavar="number"), 
    make_option(c("-I", "--Tintragenic"), type="character", default=0.8, 
              help="Proportion (%) of TEs that overlap a gene to be view as intragenic [default= %default]", metavar="number"),
    make_option(c("-i", "--Tintron"), type="character", default=0.01,
              help="Minimum proportion (%) of an exons that a TE have to overlap to be view as partial exonic, othervise it is view as intronic, [default= %default] (0.01%)", metavar="number"),
    make_option(c("-e", "--Texon1"), type="character", default=0.8, 
              help="Minimum proportion (%) of a TE in an exon to be considered as exonic, [default= %default] (80%)", metavar="number"),       
    make_option(c("-E", "--Texon2"), type="character", default=0.8, 
              help="Minimum proportion of a TEs in an exon to be considered as fragmented exonic, [default= %default] (80%)", metavar="number"), 
    make_option(c("-n", "--Tneighbor"), type="character", default=2000, 
              help="Maximum distance between a gene and a neighbor TE, [default= %default] (2000 bp)", metavar="number"),
    make_option(c("-X", "--MinLexons"), type="character", default=10, 
              help="Minimum length of exons',  [default= %default] (10 bp)", metavar="number"),
    make_option(c("-m", "--MinLtransposons"), type="character", default=0, 
              help="Minimum length of transposons ',  [default= %default] (0 bp)", metavar="number"),
    make_option(c("-M", "--MaxLtransposons"), type="character", default=50000, 
              help="Maximun length of transposons',  [default= %default] (50000 bp)", metavar="number"),
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

logo <- read.delim("../Rscript/logo.txt",header=F)
print(logo,row.names=F,col.names=F)

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
system(paste("awk -F'\t' '$3~/CDS/'", opt$genes, "> ../ParasiTE_output/cds_annotation.gff3"))
system(paste("convert2bed --input=GFF < ../ParasiTE_output/cds_annotation.gff3 > ../ParasiTE_output/pre-cds_annotation.bed && 
bedtools sort -i ../ParasiTE_output/pre-cds_annotation.bed > ../ParasiTE_output/cds_annotation.bed"))

system(paste("awk -F'\t' '$3~/five_prime_UTR/'", opt$genes, "> ../ParasiTE_output/five_annotation.gff3"))
system(paste("convert2bed --input=GFF < ../ParasiTE_output/five_annotation.gff3 > ../ParasiTE_output/pre-five_annotation.bed && 
bedtools sort -i ../ParasiTE_output/pre-five_annotation.bed > ../ParasiTE_output/five_annotation.bed"))

system(paste("awk -F'\t' '$3~/three_prime_UTR/'", opt$genes, "> ../ParasiTE_output/three_annotation.gff3"))
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

#######1.3 (optional) load of the TE gene-like annotation

if(is.null(opt$transposonsgenes)){
  print("No TE gene-like annotation provided", call.=FALSE)
cat(NULL,file="../ParasiTE_output/transcript_FP2.bed")
  } else{ 
  print("A TE gene-like annotation is provided", call.=FALSE)
#Convert and sort the TE annotation file
system(paste("cp ",opt$transposonsgenes, "../ParasiTE_output/TE-genes-like_annotation.gff3"))
system(paste("convert2bed --input=GFF < ../ParasiTE_output/TE-genes-like_annotation.gff3 > ../ParasiTE_output/pre-TE-genes-like_annotation.bed && 
bedtools sort -i ../ParasiTE_output/pre-TE-genes-like_annotation.bed > ../ParasiTE_output/TE-genes-like_annotation.bed"))

#Search for TEs that cover at least XX% of a TE-like genes and remove it
system(paste("bedtools intersect -nonamecheck -F ", opt$Tfalsepositive," -wo -a ../ParasiTE_output/TE-genes-like_annotation.bed -b ../ParasiTE_output/transcript_annotation.bed  > ../ParasiTE_output/gene_transcript_FP.bed "))
xxxx <- "awk 'BEGIN{FS=OFS=\"\t\"} {print $11,$12,$13,$14,$15,$16,$17,$18,$19,$20}' ../ParasiTE_output/gene_transcript_FP.bed > ../ParasiTE_output/transcript_FP2.bed"
system(paste(xxxx))
}

system(paste("cat ../ParasiTE_output/transcript_FP1.bed ../ParasiTE_output/transcript_FP2.bed > ../ParasiTE_output/transcript_FP.bed &&
sort ../ParasiTE_output/transcript_FP.bed | uniq > ../ParasiTE_output/transcript_FP.uniq.bed &&
bedtools sort -i ../ParasiTE_output/transcript_FP.uniq.bed > ../ParasiTE_output/transcript_FP.uniq.sorted.bed"))

#Load transcript annotation again
transcript_annotation <- read.table("../ParasiTE_output/transcript_annotation.bed", header = FALSE, sep = '\t', check.names=FALSE)

#Found isoforms that match with TE-like gene or that are consisderate as False positifs 
if (file.exists("../ParasiTE_output/transcript_FP.uniq.sorted.bed")&(file.size("../ParasiTE_output/transcript_FP.uniq.sorted.bed")>0)) {
print("TE-like genes False Positives found")
transcript_FP <-  read.table("../ParasiTE_output/transcript_FP.uniq.sorted.bed", header = FALSE, sep = '\t', check.names=FALSE)
colnames(transcript_FP)<-c("chromosome","start", "end", "poin", "other", "strand","tool", "type", "point","id")
transcript_FP <- cSplit(transcript_FP,"id",";")
transcript_FP <- setDT(transcript_FP)[, paste0("id_2", 1:2) := tstrsplit(transcript_FP$id_2, " ")]
FP_list <- transcript_FP$id_22
FP_list <- as.character(gsub(' ','',FP_list))
write.table(FP_list,file="../ParasiTE_output/FP_list.txt",sep = "\t",row.names=F,col.names=F)
#Retrieve the "STRG.XX.X number and remove the corresponding transcript and exons 
system(paste("grep -vFwf ../ParasiTE_output/FP_list.txt ../ParasiTE_output/transcript_annotation.bed > ../ParasiTE_output/transcript_annotation_wo_FP.bed"))
system(paste("grep -vFwf ../ParasiTE_output/FP_list.txt ../ParasiTE_output/exon_transcript_annotation.bed > ../ParasiTE_output/exon_transcript_annotation_wo_FP.bed"))
 } else{
print("NO TE-like genes false Positive found")
system(paste("cp ../ParasiTE_output/transcript_annotation.bed ../ParasiTE_output/transcript_annotation_wo_FP.bed"))
system(paste("cp ../ParasiTE_output/exon_transcript_annotation.bed ../ParasiTE_output/exon_transcript_annotation_wo_FP.bed"))
 }

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

### 4.4 Annotation file: Method1 and Method2 - Fragmented and Unfragmented exonic TEs
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
#system(paste("Rscript ./manage-exonic-TEs_v3.R"))


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
if (opt$Pmode=="SM"){
transposons$idexon <- ifelse(grepl("cov",transposons$idexon),transposons$idexon,paste0(transposons$idexon," cov 0.000000;"))
}
#FOR ANNOTATION OF STRINGTIE R 
if (opt$Pmode=="SR"){
transposons$idexon <- ifelse(grepl("cov",transposons$idexon),transposons$idexon,paste0(transposons$idexon,"cov0.000000;"))
}

## 1.3/ Stringtie exon annotation is loaded

All_exons <- read.table("../ParasiTE_output/exon_transcript_annotation_wo_FP.bed", header = FALSE, sep = '\t', check.names=FALSE)
colnames(All_exons)<-c("chromosome","start", "end", "poin", "other", "strand","tool", "type", "point","id")
#check if the cov is indicated for exon, beacause Stringtie do not give the coverage if we use the option -R in Stringtie. Therefore this cause problem for Parasite.
#Here, the script check if the coverage is indicated, if not it add a false one 0.000000, to make the script work.
#FOR ANNOTATION OF ARAPORT11 
if (opt$Pmode=="SM"){
All_exons$id <- ifelse(grepl("cov",All_exons$id),All_exons$id,paste0(All_exons$id," cov 0.000000;"))
}
#FOR ANNOTATION OF STRINGTIE R 
if (opt$Pmode=="SR"){
transposons$idexon <- ifelse(grepl("cov",transposons$idexon),transposons$idexon,paste0(transposons$idexon,"cov0.000000;"))
}

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
#join_gene<-join_gene[!(join_gene$begining ==0 & join_gene$middle == 0 & join_gene$ending == 0) ,]

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

print("coucou")

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



