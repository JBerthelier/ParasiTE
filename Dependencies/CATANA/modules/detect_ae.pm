package detect_ae;

# Author: Cheng-Kai Shiau <shiauck at gmail dot com>

use strict;
use modules::utils;
use Exporter;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

$VERSION     = 0.8;
@ISA         = qw(Exporter);
@EXPORT      = ();
@EXPORT_OK   = qw(detect_ae);
%EXPORT_TAGS = (DEFAULT => [qw(&detect_ae)]);

sub detect_ae($$) {
   my($data_h, $log_f) = @_;

   my $gene_no = scalar(keys %$data_h);
   if($log_f ne "") {
      open(LOG, ">>$log_f");
      print LOG "Total gene number: $gene_no\n\n";
      close(LOG);
   }

   print "\nDetecting alternative events...\n";
   my $counter_no = 0;
   foreach my $gene_id(sort {$a cmp $b} keys %$data_h) {
      printf("%.2f", ++$counter_no / $gene_no * 100);
      print "\%\r";
      extract_exon_usage_comparison_string(\%{$data_h->{$gene_id}}, $gene_id, $log_f);
   }
   print "\nDetecting alternative events finished\n";

   return;
}

sub extract_exon_usage_comparison_string($$$) {
   my($gene_h, $gene_id, $log_f) = @_;
   my @tx_list = sort {$a cmp $b} keys %{$gene_h->{transcripts}};

   for(my $tx_ind_1 = 0; $tx_ind_1 < $#tx_list; $tx_ind_1++) {
      for(my $tx_ind_2 = $tx_ind_1 + 1; $tx_ind_2 <= $#tx_list; $tx_ind_2++) {
         my $exon_usage_comp_string = join('', @{$gene_h->{exon_usage_comparison}->{$tx_list[$tx_ind_1]}->{$tx_list[$tx_ind_2]}});

         #===separate segments by first and last common used exon===
         my($pre_string, $mid_string, $suf_string);
         if($exon_usage_comp_string =~ /^([^2]+)(2.+2)([^2]+)$/) {
            $pre_string = $1;
            $mid_string = $2;
            $suf_string = $3;

            push(@{$gene_h->{AE}->{SE}},    detect_se($pre_string,   $mid_string, $suf_string,
                                                      \@{$gene_h->{exon_usage_flag}->{$tx_list[$tx_ind_1]}},
                                                      \@{$gene_h->{exon_usage_flag}->{$tx_list[$tx_ind_2]}},
                                                      $tx_list[$tx_ind_1], $tx_list[$tx_ind_2]));
            push(@{$gene_h->{AE}->{RI}},    detect_ri($pre_string,   $mid_string, $suf_string,
                                                      \@{$gene_h->{exon_usage_flag}->{$tx_list[$tx_ind_1]}},
                                                      \@{$gene_h->{exon_usage_flag}->{$tx_list[$tx_ind_2]}},
                                                      $tx_list[$tx_ind_1], $tx_list[$tx_ind_2]));
            push(@{$gene_h->{AE}->{A5SS}},  detect_a5ss($pre_string, $mid_string, $suf_string,
                                                        \@{$gene_h->{exon_usage_flag}->{$tx_list[$tx_ind_1]}},
                                                        \@{$gene_h->{exon_usage_flag}->{$tx_list[$tx_ind_2]}},
                                                        $tx_list[$tx_ind_1], $tx_list[$tx_ind_2]));
            push(@{$gene_h->{AE}->{A3SS}},  detect_a3ss($pre_string, $mid_string, $suf_string,
                                                        \@{$gene_h->{exon_usage_flag}->{$tx_list[$tx_ind_1]}},
                                                        \@{$gene_h->{exon_usage_flag}->{$tx_list[$tx_ind_2]}},
                                                        $tx_list[$tx_ind_1], $tx_list[$tx_ind_2]));
            push(@{$gene_h->{AE}->{ATSS}},  detect_atss($pre_string, $mid_string, $suf_string, $tx_list[$tx_ind_1], $tx_list[$tx_ind_2]));
            push(@{$gene_h->{AE}->{AFE}},   detect_afe($pre_string,  $mid_string, $suf_string, $tx_list[$tx_ind_1], $tx_list[$tx_ind_2]));
            push(@{$gene_h->{AE}->{ATTS}},  detect_atts($pre_string, $mid_string, $suf_string, $tx_list[$tx_ind_1], $tx_list[$tx_ind_2]));
            push(@{$gene_h->{AE}->{ALE}},   detect_ale($pre_string,  $mid_string, $suf_string, $tx_list[$tx_ind_1], $tx_list[$tx_ind_2]));

            my($mse_arr, $mxe_arr) = detect_mse_mxe($pre_string, $mid_string, $suf_string,
                                                    \@{$gene_h->{exon_usage_flag}->{$tx_list[$tx_ind_1]}},
                                                    \@{$gene_h->{exon_usage_flag}->{$tx_list[$tx_ind_2]}},
                                                    $tx_list[$tx_ind_1], $tx_list[$tx_ind_2]);
            push(@{$gene_h->{AE}->{MSE}}, @$mse_arr);
            push(@{$gene_h->{AE}->{MXE}}, @$mxe_arr);
         } elsif($exon_usage_comp_string =~ /2/) {
            if($exon_usage_comp_string =~ /^([^2]+)(2+)/) {
               $pre_string = $1;
               $mid_string = $2;
               $suf_string = $';

               push(@{$gene_h->{AE}->{ATSS}},  detect_atss($pre_string, $mid_string, $suf_string, $tx_list[$tx_ind_1], $tx_list[$tx_ind_2]));
               push(@{$gene_h->{AE}->{AFE}},   detect_afe($pre_string,  $mid_string, $suf_string, $tx_list[$tx_ind_1], $tx_list[$tx_ind_2]));
            }

            if($exon_usage_comp_string =~ /(2+)([^2]+)$/) {
               $pre_string = $`;
               $mid_string = $1;
               $suf_string = $2;

               push(@{$gene_h->{AE}->{ATTS}},  detect_atts($pre_string, $mid_string, $suf_string, $tx_list[$tx_ind_1], $tx_list[$tx_ind_2]));
               push(@{$gene_h->{AE}->{ALE}},   detect_ale($pre_string,  $mid_string, $suf_string, $tx_list[$tx_ind_1], $tx_list[$tx_ind_2]));
            }
         }

         #===dump log===
#         utils::dump_log($pre_string, $mid_string, $suf_string,
#            \@{$gene_h->{exon_usage_comparison}->{$tx_list[$tx_ind_1]}->{$tx_list[$tx_ind_2]}},
#            \@{$gene_h->{exon_usage_flag}->{$tx_list[$tx_ind_1]}},
#            \@{$gene_h->{exon_usage_flag}->{$tx_list[$tx_ind_2]}},
#            $tx_list[$tx_ind_1], $tx_list[$tx_ind_2], $gene_id, $gene_h->{strand}, \%$gene_h, $log_f)
#            if $log_f ne "";
         #===dump log===

      }
   }

   reverse_strand(\%$gene_h);

   return;

}

sub detect_se($$$$$$$) {
   #===SE: skipped exon===
   #    [ A ]---[ B ]---[ C ]
   #    [ A ]-----------[ C ]
   #tx1:  1   0   1   0   1
   #tx2:  1   0   0   0   1
   #      2   0   1   0   2
   #     (                 )   $exon_usage_comp_string=$mid_string
   my($pre_string, $exon_usage_comp_string, $suf_string, $tx_1_arr, $tx_2_arr, $tx_1, $tx_2) = @_;
   my($pos_start,  $pos_end, @pos_arr);

   while($exon_usage_comp_string =~ /(2)(0+1+0+)(2)/g) {
      my @seq_arr = split('', $2);

      #===when exon usage comp. pattern having more than 1 exon, do MSE and MXE detection===
      $pos_start = length($pre_string) - 1 + length($`) + 1;
      $pos_end   = length($pre_string) - 1 + length($`) + length($1) + length($2) + length($3);

      if(sum(@$tx_1_arr[($pos_start + length($1)) .. ($pos_end - length($3))]) == 0 ||
         sum(@$tx_2_arr[($pos_start + length($1)) .. ($pos_end - length($3))]) == 0) {
      #===one of the two tx contains all 0 in the mid string
         $pos_arr[$#pos_arr + 1]->{start} = $pos_start;
         $pos_arr[$#pos_arr]->{end}       = $pos_end;
         $pos_arr[$#pos_arr]->{tx1}->[0]  = $tx_1;
         $pos_arr[$#pos_arr]->{tx2}->[0]  = $tx_2;
      }

      pos($exon_usage_comp_string) = $-[0] + 1;
   }

   return @pos_arr;
}

sub detect_mse_mxe($$$$$$$) {
   #===MSE: multiple skipped exons, MXE: mutually exclusive exons===
   #===using # of 0/1 transition, e.g. # of exon-intron junctions, for differentiation of MSE and MXE===
   #===MSE===
   #tx1: [ A ]---[ B ]---[ C ]---[ D ]   # of E-I junx: 6
   #tx2: [ A ]-------------------[ D ]   # of E-I junx: 2
   #tx1:   1   0   1   0   1   0   1
   #tx2:   1   0   0   0   0   0   1     <=== mid string contains only 0, hence, this pattern is MSE
   #       2   0   1   0   1   0   2
   #         (                   )      $exon_usage_comp_string=$mid_string
   #===MXE case 1===
   #tx1: [ A ]---[ B ]-----------[ D ]   # of E-I junx: 4
   #tx2: [ A ]-----------[ C ]---[ D ]   # of E-I junx: 4
   #tx1:   1   0   1   0   0   0   1
   #tx2:   1   0   0   0   1   0   1
   #       2   0   1   0   1   0   2
   #         (                   )      $exon_usage_comp_string=$mid_string

   #===MXE case 2===
   #tx1: [ A ]---[ B ]-------[ D ]   # of E-I junx: 4
   #tx2: [ A ]-------[ C ]---[ D ]   # of E-I junx: 4
   #tx1:   1   0   1   0   0   1
   #tx2:   1   0   0   1   0   1
   #       2   0   1   1   0   2
   #         (               )      $exon_usage_comp_string=$mid_string

   my($pre_string, $exon_usage_comp_string, $suf_string, $tx_1_arr, $tx_2_arr, $tx_1, $tx_2) = @_;
   my($pos_start,  $pos_end, @mse_arr, @mxe_arr);

#   while($exon_usage_comp_string =~ /(2)(0+1[^2]*10+)(2)/g) {
#   temporarily put MXE case 2 aside, rare case
   while($exon_usage_comp_string =~ /(2)(0+1+0[^2]*10+)(2)/g) {
      my @seq_arr = split('', $2);

      #===when exon usage comp. pattern having more than 1 exon, do MSE and MXE detection===
      $pos_start = length($pre_string) - 1 + length($`) + 1;
      $pos_end   = length($pre_string) - 1 + length($`) + length($1) + length($2) + length($3);

      my($E_I_junx_1, $E_I_junx_2) = (0, 0);
      for(my $i = $pos_start;  $i < $pos_end; $i++) {
         $E_I_junx_1 += 1 if $$tx_1_arr[$i] ne $$tx_1_arr[$i + 1];
         $E_I_junx_2 += 1 if $$tx_2_arr[$i] ne $$tx_2_arr[$i + 1];
      }

       if($E_I_junx_1 >= 4 || $E_I_junx_2 >= 4) {
         if($E_I_junx_1 == 2 || $E_I_junx_2 == 2) {
         #===one of the two tx contains all 0 in the mid string
            $mse_arr[$#mse_arr + 1]->{start} = $pos_start;
            $mse_arr[$#mse_arr]->{end}       = $pos_end;
            $mse_arr[$#mse_arr]->{tx1}->[0]  = $tx_1;
            $mse_arr[$#mse_arr]->{tx2}->[0]  = $tx_2;
         } elsif($E_I_junx_1 >= 4 && $E_I_junx_2 >= 4) {
            $mxe_arr[$#mxe_arr + 1]->{start} = $pos_start;
            $mxe_arr[$#mxe_arr]->{end}       = $pos_end;
            $mxe_arr[$#mxe_arr]->{tx1}->[0]  = $tx_1;
            $mxe_arr[$#mxe_arr]->{tx2}->[0]  = $tx_2;
         }
      }

      pos($exon_usage_comp_string) = $-[0] + 1;
   }

   return(\@mse_arr, \@mxe_arr);
}

sub sum(@) {
   my $res = 0;
   $res += $_ for @_;
   return $res;
}

sub detect_ri($$$$$$$) {
   #===RI: retained intron===
   #    [      A      ]
   #    [ A1 ]---[ A2 ]
   #tx1:  1    1   1
   #tx2:  1    0   1
   #      2    1   2
   #     (          )   $exon_usage_comp_string=$mid_string
   my($pre_string, $exon_usage_comp_string, $suf_string, $tx_1_arr, $tx_2_arr, $tx_1, $tx_2) = @_;
   my($pos_start, $pos_end, @pos_arr);

   while($exon_usage_comp_string =~ /(2)(1+)(2)/g) {
      $pos_start = length($pre_string) - 1 + length($`) + 1;
      $pos_end   = length($pre_string) - 1 + length($`) + length($1) + length($2) + length($3);

      if(sum(@$tx_1_arr[($pos_start + length($1)) .. ($pos_end - length($3))]) == 0 ||
         sum(@$tx_2_arr[($pos_start + length($1)) .. ($pos_end - length($3))]) == 0) {

         $pos_arr[$#pos_arr + 1]->{start} = $pos_start;
         $pos_arr[$#pos_arr]->{end}       = $pos_end;
         $pos_arr[$#pos_arr]->{tx1}->[0]  = $tx_1;
         $pos_arr[$#pos_arr]->{tx2}->[0]  = $tx_2;
      }

      pos($exon_usage_comp_string) = $-[0] + 1;
   }

   return @pos_arr;
}

sub detect_a5ss($$$$$$$) {
   #===A5SS: alternative 5' splice site===
   #    [ A1    ]---[ B ]
   #    [ A2 ]------[ B ]
   #tx1:  1    1  0   1
   #tx2:  1    0  0   1
   #      2    1  0   2
   #     (             )   $exon_usage_comp_string=$mid_string
   my($pre_string, $exon_usage_comp_string, $suf_string, $tx_1_arr, $tx_2_arr, $tx_1, $tx_2) = @_;
   my($pos_start, $pos_end, @pos_arr);

   while($exon_usage_comp_string =~ /(2)(1+0+)(2)/g) {
      $pos_start = length($pre_string) - 1 + length($`) + 1;
      $pos_end   = length($pre_string) - 1 + length($`) + length($1) + length($2) + length($3);

      if(sum(@$tx_1_arr[($pos_start + length($1)) .. ($pos_end - length($3))]) == 0 ||
         sum(@$tx_2_arr[($pos_start + length($1)) .. ($pos_end - length($3))]) == 0) {
         $pos_arr[$#pos_arr + 1]->{start} = $pos_start;
         $pos_arr[$#pos_arr]->{end}       = $pos_end;
         $pos_arr[$#pos_arr]->{tx1}->[0]  = $tx_1;
         $pos_arr[$#pos_arr]->{tx2}->[0]  = $tx_2;
      }

      pos($exon_usage_comp_string) = $-[0] + 1;
   }

   return @pos_arr;
}

sub detect_a3ss($$$$$$$) {
   #===A3SS: alternative 3' splice site===
   #    [ A ]---[    B1 ]
   #    [ A ]------[ B2 ]
   #tx1:  1   0  1   1
   #tx2:  1   0  0   1
   #      2   0  1   2
   #     (            )   $exon_usage_comp_string=$mid_string
   my($pre_string, $exon_usage_comp_string, $suf_string, $tx_1_arr, $tx_2_arr, $tx_1, $tx_2) = @_;
   my($pos_start, $pos_end, @pos_arr);

   while($exon_usage_comp_string =~ /(2)(0+1+)(2)/g) {
      $pos_start = length($pre_string) - 1 + length($`) + 1;
      $pos_end   = length($pre_string) - 1 + length($`) + length($1) + length($2) + length($3);

      if(sum(@$tx_1_arr[($pos_start + length($1)) .. ($pos_end - length($3))]) == 0 ||
         sum(@$tx_2_arr[($pos_start + length($1)) .. ($pos_end - length($3))]) == 0) {
         $pos_arr[$#pos_arr + 1]->{start} = $pos_start;
         $pos_arr[$#pos_arr]->{end}       = $pos_end;
         $pos_arr[$#pos_arr]->{tx1}->[0]  = $tx_1;
         $pos_arr[$#pos_arr]->{tx2}->[0]  = $tx_2;
      }

      pos($exon_usage_comp_string) = $-[0] + 1;
   }

   return @pos_arr;
}

sub detect_atss($$$$$) {
   #===ATSS: alternative transcription start site, including AFE===
   #        [    A1 ]
   #           [ A2 ]
   #tx1:  0  1   1
   #tx2:  0  0   1
   #      0  1   2
   #     (    )        $pre_string
   #         ^         Non-zero, meaning that an ATSS occured here
   #            (...   $exon_usage_comp_string=$mid_string
   #===AFE: alternative first exon===
   #        [ A ]-----------[ C ]
   #                [ B ]---[ C ]
   #tx1:  0   1   0   0   0   1
   #tx2:  0   0   0   1   0   1
   #      0   1   0   1   0   2
   #     (                 )        $pre_string
   #          ^       ^             Non-zero, meaning that an ATSS occured here
   #                         (...   $exon_usage_comp_string=$mid_string

   my($pre_string, $exon_usage_comp_string, $suf_string, $tx_1, $tx_2) = @_;
   my @pos_arr;

   if($pre_string =~ /[^0]/) {
      $pos_arr[0]->{start}    = 0;
      $pos_arr[0]->{end}      = length($pre_string) - 1 + 1;
      $pos_arr[0]->{tx1}->[0] = $tx_1;
      $pos_arr[0]->{tx2}->[0] = $tx_2;
   }

   return @pos_arr;
}

sub detect_afe($$$$$) {
   #===AFE: alternative first exon===
   #        [ A ]-----------[ C ]
   #                [ B ]---[ C ]
   #tx1:  0   1   0   0   0   1
   #tx2:  0   0   0   1   0   1
   #      0   1   0   1   0   2
   #     (                 )        $pre_string
   #          ^^^^^^^^^^^^^         pattern "10+1+0+$" occurred,
   #                                meaning that an exon appearred inside another tx's intron
   #                         (...   $exon_usage_comp_string=$mid_string

   my($pre_string, $exon_usage_comp_string, $suf_string, $tx_1, $tx_2) = @_;
   my @pos_arr;

   if($pre_string =~ /10+1+0+$/) {
      $pos_arr[0]->{start}    = 0;
      $pos_arr[0]->{end}      = length($pre_string) - 1 + 1;
      $pos_arr[0]->{tx1}->[0] = $tx_1;
      $pos_arr[0]->{tx2}->[0] = $tx_2;
   }

   return @pos_arr;
}

sub detect_atts($$$$$) {
   #===ATTS: alternative transcription termination site, including ALE===
   #    [ A1    ]
   #    [ A2 ]
   #tx1:  1    1    0
   #tx2:  1    0    0
   #      2    1    0
   #          (      )  $suf_string
   #           ^        Non-zero, meaning that an ATSS occured here
   #    ...)            $exon_usage_comp_string=$mid_string
   #===ALE: alternative last exon===
   #    [ A ]-----------[ C ]
   #    [ A ]---[ B ]
   #tx1:  1   0   0   0   1   0
   #tx2:  1   0   1   0   0   0
   #      2   0   1   0   1   0
   #         (                 )    $suf_string
   #              ^       ^         Non-zero, meaning that an ATSS occured here
   #    ...)                        $exon_usage_comp_string=$mid_string

   my($pre_string, $exon_usage_comp_string, $suf_string, $tx_1, $tx_2) = @_;
   my @pos_arr;

   if($suf_string =~ /[^0]/) {
      $pos_arr[0]->{start}    = length($pre_string) - 1 + length($exon_usage_comp_string);
      $pos_arr[0]->{end}      = length($pre_string) - 1 + length($exon_usage_comp_string) + length($suf_string);
      $pos_arr[0]->{tx1}->[0] = $tx_1;
      $pos_arr[0]->{tx2}->[0] = $tx_2;
   }

   return @pos_arr;
}

sub detect_ale($$$$$) {
   #===ALE: alternative last exon===
   #    [ A ]-----------[ C ]
   #    [ A ]---[ B ]
   #tx1:  1   0   0   0   1   0
   #tx2:  1   0   1   0   0   0
   #      2   0   1   0   1   0
   #         (                 )    $suf_string
   #          ^^^^^^^^^^^^^         pattern "^0+1+0+1" occurred,
   #    ...)                        $exon_usage_comp_string=$mid_string

   my($pre_string, $exon_usage_comp_string, $suf_string, $tx_1, $tx_2) = @_;
   my @pos_arr;

   if($suf_string =~ /^0+1+0+1/) {
      $pos_arr[0]->{start}    = length($pre_string) - 1 + length($exon_usage_comp_string);
      $pos_arr[0]->{end}      = length($pre_string) - 1 + length($exon_usage_comp_string) + length($suf_string);
      $pos_arr[0]->{tx1}->[0] = $tx_1;
      $pos_arr[0]->{tx2}->[0] = $tx_2;
   }

   return @pos_arr;
}

sub reverse_strand($) {
   my $gene_h = shift;
   my $tmp;

   if($gene_h->{'strand'} eq '-') {
      # reverse a5ss a3ss
      $tmp = $gene_h->{AE}->{A5SS};
      $gene_h->{AE}->{A5SS} = $gene_h->{AE}->{A3SS};
      $gene_h->{AE}->{A3SS} = $tmp;

      # reverse atss atts
      $tmp = $gene_h->{AE}->{ATSS};
      $gene_h->{AE}->{ATSS} = $gene_h->{AE}->{ATTS};
      $gene_h->{AE}->{ATTS} = $tmp;

      # reverse afe ale
      $tmp = $gene_h->{AE}->{AFE};
      $gene_h->{AE}->{AFE} = $gene_h->{AE}->{ALE};
      $gene_h->{AE}->{ALE} = $tmp;
   }

   return;
}

1;
