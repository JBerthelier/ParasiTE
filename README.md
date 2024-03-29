
<p align="center">
  <img width="460" height="300" src="https://github.com/JBerthelier/ParasiTE/blob/master/logo.png">
</p>

# Inventory of transposon-related gene isoforms

In brief, ParasiTE is a tool aiming to: 

1) Identify TEs located in exonic or intronic regions of genes.
2) Detect TE sequences that are co-transcribed with gene mRNA (TE-Gene transcripts / TE-G transcripts).
3) Classify the ones contributing to alternative isoforms of genes (Alternative TE-Gene isoforms / ATE-G isoforms).

ParasiTE detects candidates for TE-AS and TE-ATP events as illustrated below using [CATANA](https://github.com/shiauck/CATANA) predictions:

![](https://github.com/JBerthelier/ParasiTE/blob/master/ParasiTE_altTE-Gi_illustration.png)


# Main steps of ParasiTE

ParasiTE is composed of five main steps:

1) Remove transcripts of gene-like TEs (transcripts of active TEs which are not involved in TE-gene transcripts)
2) Discrimination of intragenic and intergenic TEs
3) Discrimination of intronic and exonic TEs
4) Detection of TE-gene (TE-G) transcripts events
5) Detection of alternative TE-gene (ATE-G) isoform events

![](https://github.com/JBerthelier/ParasiTE/blob/master/ParasiTE_steps_illustration.png)

# Prerequisites to run ParasiTE in Linux 

An installation time of around 40 min if you need to install all dependencies

0) Download ParasiTE:

`git clone https://github.com/JBerthelier/ParasiTE.git`

1) R (versions 3.6.0 or 3.6.1) [https://cran.r-project.org/src/base/R-3/](https://cran.r-project.org/src/base/R-3/)
   or you can use conda to install R 3.6.1

  `conda create -n YourEnvironment -c conda-forge r-base=3.6.1`

3) R libraries: 
- [optparse](https://cran.r-project.org/web/packages/optparse/readme/README.html)  `install.packages("optparse")`
- [stringr](https://cran.r-project.org/web/packages/stringr/readme/README.html)  `install.packages("stringr")`
- [data.table](https://www.rdocumentation.org/packages/data.table/versions/1.14.2) `install.packages("data.table")`
- [dplyr](https://www.r-project.org/nosvn/pandoc/dplyr.html) `install.packages("dplyr")`
- [splitstackshape](https://www.r-project.org/nosvn/pandoc/splitstackshape.html) `install.packages("splitstackshape")`
- [tidyr](https://www.tidyverse.org/blog/2019/09/tidyr-1-0-0/) `install.packages("tidyr")`
 
3) [bedtools](https://bedtools.readthedocs.io/en/latest/index.html) (versions 2.27.1 and 2.29.2 were tested)
you can install it locally or with conda 
`conda install bedtools`

4) [BEDOPS](https://bedops.readthedocs.io/en/latest/index.html) (version 2.4.36 was tested)
you can install locally or with conda 
`conda install -c bioconda bedops`

# Input data

Three inputs are required and one is optional

1) TE annotation in .gff/gtf 
2) gene annotation in .gff/.gtf
3) de novo transcriptome annotation in .gtf file generated by Stringtie2 (version 2.1.4 was used) (https://ccb.jhu.edu/software/stringtie/)
or a gff/gtf file that has the same structure as files generated by Stringtie2 (see below Table 1).
4) [OPTIONAL] a gene-like transcript of TE annotation (transcripts of active TEs) such as the one of [Panda et al. 2020](https://academic.oup.com/plcell/article/32/9/2687/6115671?login=true) for A. thaliana

# Outputs

Results are in "ParasiTE_output" directory:

- Annotation of intergenic and intragenic TEs.
- Among Intragenic TEs, the annotation of intragenic (intronic and exonic) TEs.
- List of TE-genes candidates
- List of ATE-G isoforms candidates

# Run Parasite demo data to test your installation

The expected run time for demo data is 2 min (it was tested on a Linux Machine Ubuntu with 8Gb of memory)  

1. extract the folder
`cd /ParasiTE/`

`tar -xvf Demo_data_Araport11.tar.gz`

2. run the command

`Rscript /Fullpathway/ParasiTE/ParasiTE_v1/ParasiTE.R -T /Fullpathway/ParasiTE/Demo_data_Araport11/TEs_urgi_tair10.min200.gff3 -R /Fullpathway/ParasiTE/Demo_data_Araport11/Athaliana_447_Araport11.transcript_exons.for_ParasiTE_SC.gtf -G /Fullpathway/ParasiTE/Demo_data_Araport11/Athaliana_447_Araport11.gene.gff3 -L /Fullpathway/ParasiTE/Demo_data_Araport11/TAIR10-Panda_cat_TE_gene-like.gff3 -P SC`

**/Fullpathway/** has to be replaced by your data pathway

If the script finished without errors and you get files in five output directories (STEP1 to STEP5) in ParasiTE_output/Results your installation is well done.

# Run your data with ParasiTE (please read until the end)

**The basic command is:**

`Rscript /Fullpathway/ParasiTE/ParasiTE_v1/ParasiTE.R -T /Fullpathway/TE_annotation.gff3 -G /Fullpathway/gene_annotation.gtf 
-R /Fullpathway/transcripts_annotation.gtf -L /Fullpathway/gene-like_TE_annotation.gff3 -P {mode}`

**ParasiTE was built to work with Stringtie2 de novo transcriptome (but can use custom transcriptomes see how to use it below).** 

We choose [Stringtie2](https://ccb.jhu.edu/software/stringtie/), because it allows to identify chimeric transcripts such as ATE-G.
Moreover. Stringtie2 supports short reads (eg. Illumina) or long reads (eg. PacBio or Oxford Nanopore).

**For the -P {mode}**

1) Transcriptomes obtained with Stringtie2 from long reads alignment

- If the transcriptome was obtained with Stringtie2 with the -L mode you must use `-P SL`   

- If the transcriptome was obtained with Stringtie2 with the -R mode you must use `-P SR`  

2) Transcriptome annotation obtained with Stringtie2 from short reads alignment
 you must use `-P SL`  

3) Transcriptome obtained by stringtie2 Merge option or custom transcriptome following the same format as Stringtie2
you must use `-P SA`   

4) Custom transcriptome following the structure displayed below (such as "/Demo_data/Athaliana_447_Araport11.transcript_exons.for_ParasiTE.gtf" )
you must use `-P SC`

**Input data structure**

For gene annotation, transcriptome annotation, TE annotation and gene-like TE annotation, the seqname must be numeric (no "Chr1" but 1, no "Mit" but a number) (see Table 1) 

**Transcriptome annotation**

ParasiTE was developed to work with Stringtie2 structure transcriptome (see Table 1)
However, for some transcriptome annotations, the exon numbering differs from Stringtie2.

![](https://github.com/JBerthelier/ParasiTE/blob/master/Help_numbering.png)

For any custom transcriptome (SC or SA) the attribute must be as below:

Table 1:

|seqname|source|feature   |start |end  |score|strand|frame|attribute                                                  |
|------:|-----:|---------:|-----:|----:|----:|-----:|----:|----------------------------------------------------------:|
|1      |Ref   |transcript|6788  |9130 |1000 | +    | .   |gene_id "ID.1"; transcript_id "ID.1.1";                    |
|1      |Ref   |exon      |6788  |7069 |1000 | .    | .   |gene_id "ID.1"; transcript_id "ID.1.1"; exon_number_id "1";|
|1      |Ref   |exon      |8571  |9130 |1000 | .    | .   |gene_id "ID.1"; transcript_id "ID.1.1"; exon_number_id "2";|
|1      |Ref   |transcript|3631  |5899 |1000 | +    | .   |gene_id "ID.1"; transcript_id "ID.1.2";                    |
|1      |Ref   |exon      |6788  |7069 |1000 | .    | .   |gene_id "ID.1"; transcript_id "ID.1.2"; exon_number_id "1";|
|1      |Ref   |exon      |7157  |7450 |1000 | .    | .   |gene_id "ID.1"; transcript_id "ID.1.2"; exon_number_id "2";|
|1      |Ref   |exon      |8571  |8737 |1000 | .    | .   |gene_id "ID.1"; transcript_id "ID.1.2"; exon_number_id "3";|

etc...

*Be careful* 
- "seqname" must be a **number** shown in the table (no characters allowed, in example "1" corresponds to "chromosome 1").

-"Attribute" must follow the correct format as displayed.

-for transcript:  gene_id "ID.1"; transcript_id "ID.1.1"; ("ÏD" can be any other word, such as GENE.1, transcript_id "GENE.1.1")

-for exon: gene_id "ID.1"; transcript_id "ID.1.1"; exon_number_id "1";

Otherwise, ParasiTE will not work properly

**TE annotation**

We used the following structure for TE annotation, you need to provide a "ID=theTEid" in the attribute (see example in Table 2)

Table 2:
|seqname|source|feature   |start |end  |score|strand|frame|attribute                                                     |
|------:|-----:|---------:|-----:|----:|----:|-----:|----:|-------------------------------------------------------------:|
|1      |Ref   |TE        |17024 |18924|.    | .    | .   |ID=AT1TE00025;Number=594947;super_family=DHX;family=ATREP3    |
|1      |Ref   |TE        |18331 |18642|.    | .    | .   |ID=AT1TE00030;Number=597081;super_family=DTA;family=ATHATN7   |
|1      |Ref   |TE        |55676 |56576|.    | .    | .   |ID=AT1TE00150;Number=592885;super_family=DTA;family=SIMPLEHAT1|

etc...


# Detailled outputs

In /ParasiTE_output/Results/STEP4_TE-G_candidates/

- List_TE-G_exonic_level.tab ## List of Genes with exonic region overlaping with a TE annotation. 
- List_TE-G.tab ## List of exons region that overlaping with a TE annotation. 

In /ParasiTE_output/Results/STEP5_altTE-G_candidates/

- List_altTE-G.tab ## List of Genes that are involved in ATE-G isoforms. 
- List_altTE-G_exonic_level.tab ## List of exons that are involved in ATE-G isoforms. 

In "List_altTE-G.tab"

|Column name                     |description                                                                      |
|-------------------------------:|--------------------------------------------------------------------------------:|
|TE_gene                         |The id of the TE-gene event                                                      | 
|Gene_id                         |The id of the gene                                                               |
|Total_number_transcripts        |Number of isoforms transcribed by the gene                                       |
|Total_transcripts_with_TE       |Number of isoforms transcribed by the gene that are overlapped by the TE         |
|total_exon_number               |Number of exons                                                                  |
|Freq_TE_isoform                 |Frequence of isoforms overlapped by the TEs for the gene                         |
|TE_id                           |The id of the TE                                                                 |
|TE_chromosome                   |The Chromosome location of the TE annotation                                     |
|TE_start                        |The start location of the TE annotation                                          |
|TE_end                          |The end location of the TE annotation                                            |
|TE_localisation                 |intragenic or intergenic                                                         |
|method                          |ParasiTE method used to detect the exonic TE (M1 and/or M2 and/or M3)            |
|Alternative_splicing            |The predicted AS event caused by the TE                                          |
|Alternative_transcription       |The predicted ATP event caused by the TE                                         |
|first                           |Count of TE overlapping with exons of gene transcripts at the first exon (5'->3')|
|middle                          |Count of TE overlapping with exons of gene transcripts at middle exon (5'->3')   |
|last                            |Count of TE overlapping with exons of gene transcripts at the last exon (5'->3') |
|single                          |Count of TE overlapping with single-exon of gene transcripts (5'->3')            | 

# Parasite help

`Rscript /Fullpathway/ParasiTE.R -h`

![](https://github.com/JBerthelier/ParasiTE/blob/master/Help_illustration.PNG)


# Citation 

Please cite our work: 

Berthelier, J., Furci, L., Asai, S. et al. Long-read direct RNA sequencing reveals epigenetic regulation of chimeric gene-transposon transcripts in Arabidopsis thaliana. 
Nature Communications 14, 3248 (2023). https://doi.org/10.1038/s41467-023-38954-z

