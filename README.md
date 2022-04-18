![](https://github.com/JBerthelier/ParasiTE/blob/master/logo.png)

# Identification of alternative TE-genes isoforms events

ParasiTE is a tool that aim to 
1) identify TEs located in exonic or intronic region of genes 
2) detect TEs that are part of genes and that contribute in alternative isoforms of genes (altTE-Gi).

ParasiTE detect candidates for TE-AS and TE-ATP events as illustrate below:


![](https://github.com/JBerthelier/ParasiTE/blob/master/ParasiTE_altTE-Gi_illustration.png)


# Main steps of ParasiTE

ParasiTE is globally composed of five steps:

1) Removing of tramscripts produced by gene-like TEs (transcripts of free TEs)
2) Discrimination of intragenic and intergenic TEs
3) Discrimination of intronic and exonic TEs
4) Detection of TE-gene isoforms (TE-Gi) events
5) Detection of alternative TE-gene isoforms (altTE-Gi) events


![](https://github.com/JBerthelier/ParasiTE/blob/master/ParasiTE_steps_illustration.png)

# Prerequise to run ParasiTE in linux

R (version 3.6.1 was used, others were not tested)

R libraries: 
- optparse
- stringr 
- data.table 
- dplyr 
- splitstackshape 
- tidyr

bedtools

# Input data

Three input data are required

1) TE annotation in .gff/gtf 
2) gene annotation in .gff/.gtf
3) de novo transcriptome annotation in .gtf file generated by Stringtie2 (https://ccb.jhu.edu/software/stringtie/)
or a gff/gtf file that have the same structure than files generated by Stringtie2 (see below).

# Outputs

ParasiTE outputs are:
. Annotation of intergenic and intragenic TEs.
. Annotation of intragenic (intronic and exonic) TEs.
. List of TE-genes candidates
. List of altTE-gene candidates


# Run ParasiTE

The command is :

`Rscript /patway/to/ParasiTE_v1.r -T /Pathway/to/TE_annotation.gff3 -G /pathway/to/gene_annotation.gtf 
-R /pathway/to/transcripts_annotation.gtf -L /pathway/to/gene-like_TE_annotation.gff3 -P {mode}`

ParasiTE was built to work with Stringtie2 transcriptome annotation. 
Stringtie2 allows to generate a de novo transcriptome annotation and help to identify potential new isoform that have not been identificate in the reference genome annotation. It support short reads (eg. Illumina) or long reads (eg. PacBio or Oxford Nanopore).

For the -P mode

- Transcriptome obtained with Stringtie2 from long reads alignement
If the tramscript was obtained with  `Stringtie2 with the -L mode` you must use `-P SA`   
If the tramscript was obtained with  `Stringtie2 with the -R mode` you must use `-P SR`  

- Transcriptome obtained with Stringtie2 from short reads alignement
no `-P` must be supplied

- Transcriptome obtained with Stringtie2 'merged'
you must use `-P SM`   

- Home-made transciptome
The home made transcriptome must follow the stringtie2 format bellow







# Parasite help

`Rscript /patway/to/ParasiTE_v1.r -h`

![](https://github.com/JBerthelier/ParasiTE/blob/master/Help_illustration.PNG)



