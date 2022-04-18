package exon_usage;

# Author: Cheng-Kai Shiau <shiauck at gmail dot com>

use strict;
use Exporter;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

$VERSION     = 0.1;
@ISA         = qw(Exporter);
@EXPORT      = ();
@EXPORT_OK   = qw(generate_exon_usage_profile);
%EXPORT_TAGS = (DEFAULT => [qw(&generate_exon_usage_profile)]);

sub generate_exon_usage_profile($) {
   my $data_h = shift;

   print "\nComputing exon usage patterns...\n";
   my $gene_id_no = scalar(keys %$data_h);
   my $counter_no = 0;
   foreach my $gene_id(sort {$a cmp $b} keys %$data_h) {
      printf("%.2f", ++$counter_no / $gene_id_no * 100);
      print "\%\r";
      assign_exon_usage_flag(\%{$data_h->{$gene_id}});
   }
   print "\nComputing exon usage patterns finished\n";

   return;
}

sub assign_exon_usage_flag($) {
   my $gene_h = shift;

   my @tx_list = sort {$a cmp $b} keys %{$gene_h->{transcripts}};

   # initial exon usage_flag with 0 (no exon before 1st exon boundary)
   for(my $tx_ind = 0; $tx_ind <= $#tx_list; $tx_ind++) {
      $gene_h->{exon_usage_flag}->{$tx_list[$tx_ind]}->[0] = 0;
   }

   # also initial pairwise comparison array of exon usage
   for(my $tx_ind_1 = 0; $tx_ind_1 < $#tx_list; $tx_ind_1++) {
      for(my $tx_ind_2 = $tx_ind_1 + 1; $tx_ind_2 <= $#tx_list; $tx_ind_2++) {
         $gene_h->{exon_usage_comparison}->{$tx_list[$tx_ind_1]}->{$tx_list[$tx_ind_2]}->[0] = 0;
      }
   }

   # initiate flag 0 before the 1st exon boundary
   # change flag if there is an exon/intron boundary in the isoform
   # the exon usage flag array means the flags of exon containing
   # the exon boundaries array means the genomic position between adjacent 2 exon usage flag
   # so the length of exon usage flag array should be the length of exon boundaries array + 1

   # foreach exons boundaries from all isoforms
   foreach my $exon_pos(sort {$a <=> $b} keys %{$gene_h->{exons_boundaries}}) {
      my $exon_usage_flag_index = $#{$gene_h->{exon_usage_flag}->{$tx_list[0]}} + 1;

      # examine isoforms having this boundary
      for(my $tx_ind = 0; $tx_ind <= $#tx_list; $tx_ind++) {
         my $pre_exon_flag = $gene_h->{exon_usage_flag}->{$tx_list[$tx_ind]}->[$exon_usage_flag_index - 1];

         # change exon usage flag
         if($gene_h->{exons_boundaries}->{$exon_pos}->{$tx_list[$tx_ind]} == 1) {
            $gene_h->{exon_usage_flag}->{$tx_list[$tx_ind]}->[$exon_usage_flag_index] =
               reverse_flag($pre_exon_flag);
         } else {
         # for the rest isoforms do not have this boundary, keep exon usage flag
            $gene_h->{exon_usage_flag}->{$tx_list[$tx_ind]}->[$exon_usage_flag_index] =
               $pre_exon_flag;
         }
      }

      # pairwise comparison of exon usage
      for(my $tx_ind_1 = 0; $tx_ind_1 < $#tx_list; $tx_ind_1++) {
         for(my $tx_ind_2 = $tx_ind_1 + 1; $tx_ind_2 <= $#tx_list; $tx_ind_2++) {
            $gene_h->{exon_usage_comparison}->{$tx_list[$tx_ind_1]}->{$tx_list[$tx_ind_2]}->[$exon_usage_flag_index] =
               $gene_h->{exon_usage_flag}->{$tx_list[$tx_ind_1]}->[$exon_usage_flag_index] +
               $gene_h->{exon_usage_flag}->{$tx_list[$tx_ind_2]}->[$exon_usage_flag_index];
         }
      }

   }
}

sub reverse_flag($) {
   my $val = shift;
   ($val == 1)?(return 0):(return 1);
}

1;
