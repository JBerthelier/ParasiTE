package utils;

# Author: Cheng-Kai Shiau <shiauck at gmail dot com>

use strict;
use Exporter;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

$VERSION     = 0.1;
@ISA         = qw(Exporter);
@EXPORT      = ();
@EXPORT_OK   = qw(show_debug_msg show_exon_usage_profile_msg show_exon_usage_comparison_msg show_ae_msg dump_log);
%EXPORT_TAGS = (DEFAULT => [qw(&show_debug_msg &show_exon_usage_profile_msg &show_exon_usage_comparison_msg &show_ae_msg &dump_log)]);

sub show_debug_msg($$) {
   my($result_obj, $log_f) = @_;

   ($log_f eq "")?(open(LOG, ">&STDOUT")):(open(LOG, ">>$log_f"));

   print LOG "\$result_obj: ", $result_obj, "\n";

   foreach my $gene_id(sort {$a cmp $b} keys %$result_obj)
   {
      print LOG "$gene_id: $result_obj->{$gene_id}\n";
      foreach my $ele(sort {$a cmp $b} keys %{$result_obj->{$gene_id}})
      {
         print LOG "\t$ele: $result_obj->{$gene_id}->{$ele}\n";
      }

      foreach my $exons_pos(sort {$a cmp $b} keys %{$result_obj->{$gene_id}->{exons_boundaries}})
      {
         print LOG "\t\t$exons_pos: $result_obj->{$gene_id}->{exons_boundaries}\n";
         foreach my $tx_id(sort {$a cmp $b} keys %{$result_obj->{$gene_id}->{exons_boundaries}->{$exons_pos}})
         {
            print LOG "\t\t\t$tx_id: $result_obj->{$gene_id}->{exons_boundaries}->{$exons_pos}->{$tx_id}\n";
         }
      }
   }

   close(LOG);
}

sub show_exon_usage_profile_msg($$) {
   my($result_obj, $log_f) = @_;

   ($log_f eq "")?(open(LOG, ">&STDOUT")):(open(LOG, ">>$log_f"));

   print LOG "\$result_obj: ", $result_obj, "\n";

   foreach my $gene_id(sort {$a cmp $b} keys %$result_obj)
   {
      print LOG "$gene_id: $result_obj->{$gene_id}\n";

      foreach my $tx_id(sort {$a cmp $b} keys %{$result_obj->{$gene_id}->{exon_usage_flag}})
      {
         print LOG "\t\t$tx_id: ";
         print LOG join('', @{$result_obj->{$gene_id}->{exon_usage_flag}->{$tx_id}});
         print LOG "\n";
      }
   }

   close(LOG);
}

sub show_exon_usage_comparison_msg($$) {
   my($result_obj, $log_f) = @_;

   ($log_f eq "")?(open(LOG, ">&STDOUT")):(open(LOG, ">>$log_f"));

   print LOG "\$result_obj: ", $result_obj, "\n";

   foreach my $gene_id(sort {$a cmp $b} keys %$result_obj)
   {
      print LOG "$gene_id: $result_obj->{$gene_id}\n";

      my @tx_list = sort {$a cmp $b} keys %{$result_obj->{$gene_id}->{transcripts}};

      for(my $i = 0; $i < $#tx_list; $i++)
      {
         for(my $j = $i + 1; $j <= $#tx_list; $j++)
         {
            print LOG "\t\t$tx_list[$i]\n";
            print LOG "\t\t$tx_list[$j]: ";
            print LOG join('', @{$result_obj->{$gene_id}->{exon_usage_comparison}->{$tx_list[$i]}->{$tx_list[$j]}}), "\n";
         }
      }
   }

   close(LOG);
}

sub show_ae_msg($$) {
   my($result_obj, $log_f) = @_;

   ($log_f eq "")?(open(LOG, ">&STDOUT")):(open(LOG, ">>$log_f"));

   print LOG "\$result_obj: ", $result_obj, "\n";

   foreach my $gene_id(sort {$a cmp $b} keys %$result_obj)
   {
      print LOG "$gene_id: $result_obj->{$gene_id}\n";
      my @exon_boundaries_arr = sort {$a <=> $b} keys %{$result_obj->{$gene_id}->{exons_boundaries}};

      foreach my $ae_id(sort {$a cmp $b} keys %{$result_obj->{$gene_id}->{AE}})
      {
         if(defined $result_obj->{$gene_id}->{AE}->{$ae_id} && $#{$result_obj->{$gene_id}->{AE}->{$ae_id}} > -1)
         {
            print LOG "\t$ae_id:\n";
            for(my $i = 0; $i <= $#{$result_obj->{$gene_id}->{AE}->{$ae_id}}; $i++)
            {
               my $exon_start_ind = $result_obj->{$gene_id}->{AE}->{$ae_id}->[$i]->{start} - 1;
                  $exon_start_ind = 0 if $exon_start_ind < 0;
               my $exon_end_ind   = $result_obj->{$gene_id}->{AE}->{$ae_id}->[$i]->{end};
                  $exon_end_ind   = $#exon_boundaries_arr if $exon_end_ind > $#exon_boundaries_arr;

               print LOG "\t\t$i: $result_obj->{$gene_id}->{AE}->{$ae_id}->[$i]->{start}\t$result_obj->{$gene_id}->{AE}->{$ae_id}->[$i]->{end}\t";
               print LOG "$exon_boundaries_arr[$exon_start_ind]\t$exon_boundaries_arr[$exon_end_ind]\n";

               print LOG "\t\t\ttx1: ", join(",", @{$result_obj->{$gene_id}->{AE}->{$ae_id}->[$i]->{tx1}}), "\n";
               print LOG "\t\t\ttx2: ", join(",", @{$result_obj->{$gene_id}->{AE}->{$ae_id}->[$i]->{tx2}}), "\n";
            }
         }
      }
   }

   close(LOG);
}

sub dump_log($$$$$$$$$$$$) {
   my($pre_string, $mid_string, $suf_string, $seq_arr, $tx_1_arr, $tx_2_arr, $tx_1, $tx_2, $gene_id, $strand, $gene_h, $log_f) = @_;
   my(@tmp_arr, $pattern);
   my $ae_type_counter = 0;

   @tmp_arr = detect_ae::detect_se($pre_string, $mid_string, $suf_string,
                                   $tx_1_arr, $tx_2_arr, $tx_1, $tx_2);
   ++$ae_type_counter if $#tmp_arr > -1;
   draw_ae($gene_id, $tx_1, $tx_2, "SE", \@tmp_arr, $seq_arr, $strand, \%$gene_h, $log_f) if $#tmp_arr > -1;

   my($mse_arr, $mxe_arr) = detect_ae::detect_mse_mxe($pre_string, $mid_string, $suf_string,
                                           $tx_1_arr, $tx_2_arr, $tx_1, $tx_2);
   ++$ae_type_counter if $#$mse_arr > -1;
   draw_ae($gene_id, $tx_1, $tx_2, "MSE", $mse_arr, $seq_arr, $strand, \%$gene_h, $log_f) if $#$mse_arr > -1;

   ++$ae_type_counter if $#$mxe_arr > -1;
   draw_ae($gene_id, $tx_1, $tx_2, "MXE", $mxe_arr, $seq_arr, $strand, \%$gene_h, $log_f) if $#$mxe_arr > -1;

   @tmp_arr = detect_ae::detect_ri($pre_string,   $mid_string, $suf_string, $tx_1_arr, $tx_2_arr, $tx_1, $tx_2);
   ++$ae_type_counter if $#tmp_arr > -1;
   draw_ae($gene_id, $tx_1, $tx_2, "RI", \@tmp_arr, $seq_arr, $strand, \%$gene_h, $log_f) if $#tmp_arr > -1;

   @tmp_arr = detect_ae::detect_a5ss($pre_string, $mid_string, $suf_string, $tx_1_arr, $tx_2_arr, $tx_1, $tx_2);
   ++$ae_type_counter if $#tmp_arr > -1;
   $pattern = ($strand eq '+')?("A5SS"):("A3SS");
   draw_ae($gene_id, $tx_1, $tx_2, $pattern, \@tmp_arr, $seq_arr, $strand, \%$gene_h, $log_f) if $#tmp_arr > -1;

   @tmp_arr = detect_ae::detect_a3ss($pre_string, $mid_string, $suf_string, $tx_1_arr, $tx_2_arr, $tx_1, $tx_2);
   ++$ae_type_counter if $#tmp_arr > -1;
   $pattern = ($strand eq '+')?("A3SS"):("A5SS");
   draw_ae($gene_id, $tx_1, $tx_2, $pattern, \@tmp_arr, $seq_arr, $strand, \%$gene_h, $log_f) if $#tmp_arr > -1;

   @tmp_arr = detect_ae::detect_atss($pre_string, $mid_string, $suf_string, $tx_1, $tx_2);
   ++$ae_type_counter if $#tmp_arr > -1;
   $pattern = ($strand eq '+')?("ATSS"):("ATTS");
   draw_ae($gene_id, $tx_1, $tx_2, $pattern, \@tmp_arr, $seq_arr, $strand, \%$gene_h, $log_f) if $#tmp_arr > -1;

   @tmp_arr = detect_ae::detect_afe($pre_string,  $mid_string, $suf_string, $tx_1, $tx_2);
   ++$ae_type_counter if $#tmp_arr > -1;
   $pattern = ($strand eq '+')?("AFE"):("ALE");
   draw_ae($gene_id, $tx_1, $tx_2, $pattern, \@tmp_arr, $seq_arr, $strand, \%$gene_h, $log_f) if $#tmp_arr > -1;

   @tmp_arr = detect_ae::detect_atts($pre_string, $mid_string, $suf_string, $tx_1, $tx_2);
   ++$ae_type_counter if $#tmp_arr > -1;
   $pattern = ($strand eq '+')?("ATTS"):("ATSS");
   draw_ae($gene_id, $tx_1, $tx_2, $pattern, \@tmp_arr, $seq_arr, $strand, \%$gene_h, $log_f) if $#tmp_arr > -1;

   @tmp_arr = detect_ae::detect_ale($pre_string,  $mid_string, $suf_string, $tx_1, $tx_2);
   ++$ae_type_counter if $#tmp_arr > -1;
   $pattern = ($strand eq '+')?("ALE"):("AFE");
   draw_ae($gene_id, $tx_1, $tx_2, $pattern, \@tmp_arr, $seq_arr, $strand, \%$gene_h, $log_f) if $#tmp_arr > -1;

   open(LOG, ">>$log_f");
   print LOG "$gene_id\t$tx_1\t$tx_2\tNo_AS\t" if $ae_type_counter == 0;
   if($ae_type_counter == 0) {
      print LOG join('', @$seq_arr), "\t$strand\n";
   }
   close(LOG);

   return;
}

sub draw_ae($$$$$$$$$) {
   my($gene_id, $tx_1, $tx_2, $as_type, $pos_arr, $seq_arr, $strand, $ae_id, $log_f) = @_;

   open(LOG, ">>$log_f");

   for(my $i = 0; $i <= $#$pos_arr; $i++)
   {
      print LOG "$gene_id.$as_type.$ae_id\t$tx_1\t$tx_2\t$as_type\t";

      print LOG join('', @$seq_arr[$pos_arr->[$i]->{start} .. $pos_arr->[$i]->{end}]);
      print LOG "\t$strand\n";
   }
   close(LOG);

   return;
}

1;
