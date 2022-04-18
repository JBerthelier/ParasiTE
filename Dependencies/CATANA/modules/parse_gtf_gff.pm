package parse_gtf_gff;

# Author: Cheng-Kai Shiau <shiauck at gmail dot com>

use strict;
use Exporter;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

$VERSION     = 0.1;
@ISA         = qw(Exporter);
@EXPORT      = ();
@EXPORT_OK   = qw(read_gtf_gff);
%EXPORT_TAGS = (DEFAULT => [qw(&read_gtf_gff)]);

sub read_gtf_gff($$) {
   my($f_name, $format) = @_;
   my $result_h;

   my $detect_format_flag = 0;
      $detect_format_flag = 1 if($format eq 'gff' || $format eq 'gtf');
   my $gff_flag = 0;
      $gff_flag = 1 if $format eq 'gff';
      $gff_flag = 0 if $format eq 'gtf';

   my $gff_gene_id = "";
   my @file_line_no = split(' ', `wc -l $f_name`);
   my $line_no = 0;
   print "\nLoading file $f_name...\n";
   open(FILE, $f_name);
   while(my $line = <FILE>) {
      printf("%.2f", ++$line_no / $file_line_no[0] * 100);
      print "\%\r";

      # ignore comment line
      next if $line =~ /^#/;
      chomp($line);

      my @line_arr = split("\t", $line);

      # ignore lines without enough columns
      next if $#line_arr != 8;

      # detect GTF/GFF
      ($detect_format_flag, $gff_flag) = detect_format($line_arr[8]) if $detect_format_flag == 0;

      # parse mRNA / transcript / exon
      if($line_arr[2] eq "gene" or $line_arr[2] eq "mRNA" or $line_arr[2] eq "transcript" or $line_arr[2] eq "exon") {
         my $attr_h;
         ($gff_flag == 1)?
            (%$attr_h = parse_gff_attr($line_arr[8])):
            (%$attr_h = parse_gtf_attr($line_arr[8]));

         if($line_arr[2] eq "gene") {
            # parse gene id for gff

            if($gff_flag == 1) {
               $gff_gene_id = $attr_h->{ID};
            }

         } elsif($line_arr[2] eq "mRNA" || $line_arr[2] eq "transcript") {
            # fill isoform structure hash

            my($gene_id, $tx_id) = ("", "");

            if($gff_flag == 1) {
               # GFF
               $gene_id = $attr_h->{Parent};
               $tx_id   = $attr_h->{ID};
            } else {
               #GTF
               $gene_id = $attr_h->{gene_id};
               $tx_id   = $attr_h->{transcript_id};
               $gene_id = $tx_id if $gene_id eq "";
            }

            $result_h->{$gene_id}->{chrom}  = $line_arr[0];
            $result_h->{$gene_id}->{strand} = $line_arr[6];
            $result_h->{$gene_id}->{transcripts}->{$tx_id} = 1;
            $result_h->{$gene_id}->{exons_boundaries}->{($line_arr[3] - 1)}->{$tx_id} = 1;
            $result_h->{$gene_id}->{exons_boundaries}->{($line_arr[4])}->{$tx_id} = 1;

         } elsif($line_arr[2] eq "exon") {

            my($gene_id, $tx_id, $exon_number, $ind_exon_number) = ("", "", "", "");

            if($gff_flag == 1) {
               # GFF
               $gene_id     = $gff_gene_id;
               $tx_id       = $attr_h->{Parent};
# calc exon no, deprecated
#               my @tmp_arr  = split(':', $attr_h->{ID});
#               $exon_number = $tmp_arr[2];
            } else {
               #GTF
               $gene_id     = $attr_h->{gene_id};
               $tx_id       = $attr_h->{transcript_id};
               $gene_id     = $tx_id if $gene_id eq "";
# calc exon no, deprecated
#               $exon_number = $attr_h->{exon_number};
            }

# calc exon no, deprecated
#            $ind_exon_number = $#{$result_h->{$gene_id}->{transcripts}->{$tx_id}->{exons}} + 1;

            $result_h->{$gene_id}->{chrom}  = $line_arr[0];
            $result_h->{$gene_id}->{strand} = $line_arr[6];
            $result_h->{$gene_id}->{transcripts}->{$tx_id} = 1;
            $result_h->{$gene_id}->{exons_boundaries}->{($line_arr[3] - 1)}->{$tx_id} = 1;
            $result_h->{$gene_id}->{exons_boundaries}->{($line_arr[4])}->{$tx_id} = 1;
         } else {
            # leave the rest tags: gene CDS start_codon stop_codon three_prime_UTR five_prime_UTR
         }
      }
   }
   close(FILE);
   print "\nLoading file $f_name Finished\n";

   return($result_h);
}

sub parse_gtf_attr($) {
   my $line = shift;

# example from ensembl
# gene_id "ENSG00000223972"; gene_name "DDX11L1"; gene_source "havana"; gene_biotype "transcribed_unprocessed_pseudogene"; 

   #remove the last white space if any
   $line =~ s/\s$//;
   #remove the last semi-colon if any
   $line =~ s/;$//;
   #remove last quot if any
   $line =~ s/\"$//;

   #remove white space in the begginning if any
   $line =~ s/^\s//;

   return(split(/\s\"|\";\s/, $line));

}

sub parse_gff_attr($) {
   my $line = shift;

# example from ensembl
# hid=trf; hstart=1; hend=21
# another example from ensembl
# ID=ENSMUSG00000005553;Name=Atp4a;Name=ENSMUSG00000005553

   #remove white space
   $line =~ s/\s//g;

   return(split(/=|;/, $line));
}

sub detect_format($) {
   my $line = shift;
   if($line =~ /\=/) {
      print "This file looks like GFF format\n\n";
      return(1, 1);
   } else {
      print "This file looks like GTF format\n\n";
      return(1, 0);
   }
}

1;
