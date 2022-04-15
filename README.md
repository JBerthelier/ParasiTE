# ParasiTE

![](https://github.com/JBerthelier/ParasiTE/blob/master/logo.png)

Identification of alternative TE-genes isoforms events from gtf/gff annotations

# Aim 

ParasiTE is a tool that aim to 
1) identify TEs located in exonic or intronic region of genes 
2) detect TEs that are part of genes and that contribute in alternative isoforms of genes.

#Steps

![](https://github.com/JBerthelier/ParasiTE/blob/master/ParasiTE_steps_illustration.png)

ParasiTE is globally composed of five steps:

1) Removing of tramscripts produced by gene-like TEs (transcripts of free TEs)
2) Discrimination of intragenic and intergenic TEs
3) Discrimination of intronic and exonic TEs
4) Detection of TE-gene isoforms (TE-Gi) events
5) Detection of alternative TE-gene isoforms (altTE-Gi) events

#Outputs
ParasiTE give in output:

. Annotation intergenic, intronic and exonic TEs.
. List of TE-related gene candidates
. List of TE-related alternative gene isoforms

______________________________________

#Prerequise to run ParasiTE

R (version XXX have been tested)
R packages: optparse , stringr , data.table , dplyr , splitstackshape , tidyr

#Input data

Three input data are required

1) TE annotation in .gff 
2) gene annotation in .gff
3) de novo transcriptome annotation in .gff file generated by Stringtie2 (https://ccb.jhu.edu/software/stringtie/)

_____________________________________

#How to generate a transcriptome annotation for ParasiTE ?

ParasiTE have been build to work with Stringtie2 transcriptome annotation. 
Stringtie2 allows to generate a de novo transcriptome annotation and help to identify potential new isoform that have not been identificate in the reference genome annotation.
It also support short reads (eg. Illumina) or long reads (eg. PacBio or Oxford Nanopore) based genome alignements. 

Here are the basic command lines to generate your transcriptome annotation for ParasiTE:

. For short-reads:

1) Hisat2

hisat2-build genome_assembly.fa index_name
hisat2 -x index_name --dta -1 forward.fq -2 reverse.fq -S short_reads_aln.sam
samtools view -Su short_reads_aln.sam > short_reads_aln.bam 
samtools sort short_reads_aln.bam -o short_reads_aln.sorted.bam

2) Stringtie2


. For long-reads:

1) minimap2

minimap2 -ax splice -uf -k14 -G 10000 genome_assembly.fa long_reads.fa > long_reads_aln.sam
samtools view -Su long_reads_aln.sam > long_reads_aln.bam
samtools sort -@ 10 long_reads_aln.bam > long_reads_aln.sorted.bam

2) Stringtie2

stringtie long_reads_aln.sorted.bam -G reference_genome_annotation.gff3 -L > Stringtie2_transcriptome_annotation.gff3
____________________________________________

#Parasite help#

Usage: ParasiTE.r [options]
Options:
        -T CHARACTER, --transposons=CHARACTER
                Annotation of transposable elements (.gff3)
        -G CHARACTER, --genes=CHARACTER
                Gene model annotation containing gene positions (.gff3)
        -R CHARACTER, --transcripts=CHARACTER
                Transcript annotation obtained by Stringtie, containing transcript and exon positions (.gff3)
        -F NUMBER, --Tfalsepositive=NUMBER
                Proportion of length of a TEs in a genes, to consider the gene as a TE-like gene, the gene annotation  will              
                be removed for the analyses [default= 0.8] (80%)
        -I NUMBER, --Tintragenic=NUMBER
                Proportion of length of a TEs that overlap a gene to be considerated as intragenic [default= 0.8] (80%)
        -i NUMBER, --Tambigous=NUMBER
                Proportion of length of an exons that an TEs can partialy overlap to be exonic ambigous, [default= 0.8]          
        -e NUMBER, --Texon=NUMBER
                Minimum proportion of a TEs in an exon to be considered as exonic, [default= 0.8] (80%)
        -n NUMBER, --Tneighbor=NUMBER
                Maximum distance between a TE to a gene to be considerate as a neighbor TE, [default= 1000] (1000 bp)
        -E NUMBER, --MinLexons=NUMBER
                Min length of exons',  [default= 200] (0)
        -m NUMBER, --MinLtransposons=NUMBER
                Min length of transposons',  [default= 0] (0)
        -M NUMBER, --MaxLtransposons=NUMBER
                Maximun length of transposons',  [default= 30000] (0)
        -h, --help                 Show this help message and exit




