package report_result;

# Author: Cheng-Kai Shiau <shiauck at gmail dot com>

use strict;
use Exporter;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

$VERSION     = 0.1;
@ISA         = qw(Exporter);
@EXPORT      = ();
@EXPORT_OK   = qw(report_result);
%EXPORT_TAGS = (DEFAULT => [qw(&report_result)]);

sub report_result($$$) {
   my($data_h, $o_dir, $log_f) = @_;

   print "\nMerging identical alternative events...\n";
   my $gene_no = scalar(keys %$data_h);
   my $counter_no = 0;
   foreach my $gene_id(sort {$a cmp $b} keys %$data_h) {
      printf("%.2f", ++$counter_no / $gene_no * 100);
      print "\%\r";
      remove_redundant_ae(\%{$data_h->{$gene_id}});
   }
   print "\nMerging identical alternative events finished\n";

   output_gff_format(\%$data_h, $o_dir, $log_f);

   return;
}

sub remove_redundant_ae($) {
   my $gene_h = shift;

   foreach my $ae_type(keys %{$gene_h->{AE}}) {
      my $tmp_h;
      for(my $ae_id = 0; $ae_id <= $#{$gene_h->{AE}->{$ae_type}}; $ae_id++) {
         my $start_pos = sprintf("%s", $gene_h->{AE}->{$ae_type}->[$ae_id]->{start});
         my $end_pos   = sprintf("%s", $gene_h->{AE}->{$ae_type}->[$ae_id]->{end});
         #===guarantee exon_usage_1 is the transcript having the first left most exon inside
         my($exon_usage_1, $exon_usage_2) = retrieve_exon_usage_string(\%$gene_h, $ae_type, $ae_id);

         #===remove redundant===
         if(defined $tmp_h->{$start_pos}->{$end_pos}->{$exon_usage_1}->{$exon_usage_2}->{tx1}->[0]) {
            push(@{$tmp_h->{$start_pos}->{$end_pos}->{$exon_usage_1}->{$exon_usage_2}->{tx1}},
                 @{$gene_h->{AE}->{$ae_type}->[$ae_id]->{tx1}});
            push(@{$tmp_h->{$start_pos}->{$end_pos}->{$exon_usage_1}->{$exon_usage_2}->{tx2}},
                 @{$gene_h->{AE}->{$ae_type}->[$ae_id]->{tx2}});
         } else {
            $tmp_h->{$start_pos}->{$end_pos}->{$exon_usage_1}->{$exon_usage_2}->{tx1}->[0] =
               $gene_h->{AE}->{$ae_type}->[$ae_id]->{tx1}->[0];
            $tmp_h->{$start_pos}->{$end_pos}->{$exon_usage_1}->{$exon_usage_2}->{tx2}->[0] =
               $gene_h->{AE}->{$ae_type}->[$ae_id]->{tx2}->[0];
         }
         #===remove redundant===
      }

      #===clear data===
      $#{$gene_h->{AE}->{$ae_type}} = -1;

      #===reassign data===
      foreach my $start_pos(sort {$a <=> $b} keys %$tmp_h) {
         foreach my $end_pos(sort {$a <=> $b} keys %{$tmp_h->{$start_pos}}) {
            foreach my $exon_usage_1(keys %{$tmp_h->{$start_pos}->{$end_pos}}) {
               foreach my $exon_usage_2(keys %{$tmp_h->{$start_pos}->{$end_pos}->{$exon_usage_1}}) {
                  my $ae_type_ind = $#{$gene_h->{AE}->{$ae_type}} + 1;
                  $gene_h->{AE}->{$ae_type}->[$ae_type_ind]->{start} = $start_pos;
                  $gene_h->{AE}->{$ae_type}->[$ae_type_ind]->{end}   = $end_pos;
                  @{$gene_h->{AE}->{$ae_type}->[$ae_type_ind]->{tx1}} =
                     uniq(@{$tmp_h->{$start_pos}->{$end_pos}->{$exon_usage_1}->{$exon_usage_2}->{tx1}});
                  @{$gene_h->{AE}->{$ae_type}->[$ae_type_ind]->{tx2}} =
                     uniq(@{$tmp_h->{$start_pos}->{$end_pos}->{$exon_usage_1}->{$exon_usage_2}->{tx2}});
               }
            }
         }
      }
      #===reassign data===
   }

   return;
}

sub retrieve_exon_usage_string($$$) {
   my($gene_h, $ae_type, $ae_id) = @_;

   my $start_pos    = sprintf("%s", $gene_h->{AE}->{$ae_type}->[$ae_id]->{start});
   my $end_pos      = sprintf("%s", $gene_h->{AE}->{$ae_type}->[$ae_id]->{end});
   my $tx1          = $gene_h->{AE}->{$ae_type}->[$ae_id]->{tx1}->[0];
   my $tx2          = $gene_h->{AE}->{$ae_type}->[$ae_id]->{tx2}->[0];
   my $exon_usage_1 = join('', @{$gene_h->{exon_usage_flag}->{$tx1}}[$start_pos .. $end_pos]);
   my $exon_usage_2 = join('', @{$gene_h->{exon_usage_flag}->{$tx2}}[$start_pos .. $end_pos]);

   # sort tx1 tx2 by exon_usage_string
   if(comparing_exon_usage_string($exon_usage_1, $exon_usage_2)) {
      return($exon_usage_1, $exon_usage_2);
   } else {
      my @tmp = @{$gene_h->{AE}->{$ae_type}->[$ae_id]->{tx1}};
         @{$gene_h->{AE}->{$ae_type}->[$ae_id]->{tx1}} = @{$gene_h->{AE}->{$ae_type}->[$ae_id]->{tx2}};
         @{$gene_h->{AE}->{$ae_type}->[$ae_id]->{tx2}} = @tmp;
      return($exon_usage_2, $exon_usage_1);
   }
}

sub uniq {
    my %seen;
    grep !$seen{$_}++, @_;
}

sub comparing_exon_usage_string($$) {
   my($exon_usage_1, $exon_usage_2) = @_;

   my(@exon_usage_1) = split('', $exon_usage_1);
   my(@exon_usage_2) = split('', $exon_usage_2);

   for(my $i = 1; $i <= $#exon_usage_1; $i++) {
      if($exon_usage_1[$i] eq '1' && $exon_usage_2[$i] eq '0') {
         return 1;
      } elsif($exon_usage_1[$i] eq '0' && $exon_usage_2[$i] eq '1') {
         return 0;
      } else {}
   }
}

sub output_gff_format($$$) {
   my($data_h, $o_dir, $log_f) = @_;

   print "\nGenerating annotation file for alternative events...\n";

   my $gene_no = scalar(keys %$data_h);
   print "\nTotal gene number: $gene_no\n";

   my $gene_counter = 0;
   foreach my $gene_id(sort {$a cmp $b} keys %$data_h)
   {
      print ++$gene_counter, " of $gene_no (";
      printf("%.2f", $gene_counter / $gene_no * 100);
      print "%)\r";

      my @exon_boundaries_arr = sort {$a <=> $b} keys %{$data_h->{$gene_id}->{exons_boundaries}};

      foreach my $ae_type(sort {$a cmp $b} keys %{$data_h->{$gene_id}->{AE}})
      {
         if(defined $data_h->{$gene_id}->{AE}->{$ae_type} && $#{$data_h->{$gene_id}->{AE}->{$ae_type}} > -1)
         {
            my $o_name = "$o_dir/$ae_type.gff";

            for(my $ae_id = 0; $ae_id <= $#{$data_h->{$gene_id}->{AE}->{$ae_type}}; $ae_id++)
            {
               my $tx1_acc = $data_h->{$gene_id}->{AE}->{$ae_type}->[$ae_id]->{tx1}->[0];
               my $tx2_acc = $data_h->{$gene_id}->{AE}->{$ae_type}->[$ae_id]->{tx2}->[0];

               my $alt_event_id = "${gene_id}.${ae_type}.${ae_id}";

               #===calc gff gene left and right boundaries===
               my @exon_usage_comp_arr =
                  (defined $data_h->{$gene_id}->{exon_usage_comparison}->{$tx1_acc}->{$tx2_acc})?
                     (@{$data_h->{$gene_id}->{exon_usage_comparison}->{$tx1_acc}->{$tx2_acc}}):
                     (@{$data_h->{$gene_id}->{exon_usage_comparison}->{$tx2_acc}->{$tx1_acc}});
               my($AS_start_ind, $AS_end_ind) =
                  trace_back_left_right_boundary(\@exon_usage_comp_arr,
                     $data_h->{$gene_id}->{AE}->{$ae_type}->[$ae_id]->{start},
                     $data_h->{$gene_id}->{AE}->{$ae_type}->[$ae_id]->{end});
               $AS_start_ind--;
               $AS_start_ind = 0 if $AS_start_ind < 0;

               my $chrom_name = $data_h->{$gene_id}->{'chrom'};
               $chrom_name = "chr" . $chrom_name unless $chrom_name =~ /^chr/;

               open(OUT, ">>$o_name");

               #===event gff gene line===
               print OUT "$chrom_name\t$ae_type\tgene\t";
               print OUT ($exon_boundaries_arr[$AS_start_ind] + 1), "\t";
               print OUT $exon_boundaries_arr[$AS_end_ind], "\t";
               print OUT ".\t$data_h->{$gene_id}->{'strand'}\t.\t";
               print OUT "ID=", $alt_event_id, ";Name=", $gene_id, ";gid=", $gene_id, "\n";
               #===event gff gene line===

               print OUT generate_gff_mRNA_exon(\%$data_h, $gene_id, $tx1_acc, 'A', $ae_type, $ae_id, $alt_event_id);
               print OUT generate_gff_mRNA_exon(\%$data_h, $gene_id, $tx2_acc, 'B', $ae_type, $ae_id, $alt_event_id);

               close(OUT);

               #===dump log===
               if($log_f ne "") {

                  my $pos_h;
                  $pos_h->[0]->{start} = $data_h->{$gene_id}->{AE}->{$ae_type}->[$ae_id]->{start};
                  $pos_h->[0]->{end}   = $data_h->{$gene_id}->{AE}->{$ae_type}->[$ae_id]->{end};
                  my @seq_arr;
                  if(defined $data_h->{$gene_id}->{'exon_usage_comparison'}->{$tx1_acc}->{$tx2_acc}) {
                     @seq_arr = @{$data_h->{$gene_id}->{'exon_usage_comparison'}->{$tx1_acc}->{$tx2_acc}};
                  } else {
                     @seq_arr = @{$data_h->{$gene_id}->{'exon_usage_comparison'}->{$tx2_acc}->{$tx1_acc}};
                  }
                  utils::draw_ae($gene_id,
                                 join(',', @{$data_h->{$gene_id}->{AE}->{$ae_type}->[$ae_id]->{tx1}}),
                                 join(',', @{$data_h->{$gene_id}->{AE}->{$ae_type}->[$ae_id]->{tx2}}),
                                 $ae_type,
                                 \@$pos_h,
                                 \@seq_arr,
                                 $data_h->{$gene_id}->{'strand'}, $ae_id, $log_f);
               }
               #===dump log===
            }
         }
      }
   }
   print "\n";

   return;
}

sub generate_alt_event_id($$$$$$) {
   my($data_h, $gene_id, $tx1_acc, $tx2_acc, $ae_type, $ae_id) = @_;

   my $AS_start_ind = $data_h->{$gene_id}->{AE}->{$ae_type}->[$ae_id]->{start};
   my $AS_end_ind   = $data_h->{$gene_id}->{AE}->{$ae_type}->[$ae_id]->{end};
   my @exon_usage_comp_arr =
      (defined $data_h->{$gene_id}->{exon_usage_comparison}->{$tx1_acc}->{$tx2_acc})?
         (@{$data_h->{$gene_id}->{exon_usage_comparison}->{$tx1_acc}->{$tx2_acc}}):
         (@{$data_h->{$gene_id}->{exon_usage_comparison}->{$tx2_acc}->{$tx1_acc}});
   my @exon_bound_arr = sort {$a <=> $b} keys %{$data_h->{$gene_id}->{exons_boundaries}};

   #===trace AS event left boundary===
   ($AS_start_ind, $AS_end_ind) = trace_back_left_right_boundary(\@exon_usage_comp_arr, $AS_start_ind, $AS_end_ind);

   #===generate gff ID===
      #===initiate the left first exon boundary===
   my $alt_event_string = $exon_bound_arr[$AS_start_ind - 1];
   for(my $pos_ind = $AS_start_ind + 1; $pos_ind <= $AS_end_ind; $pos_ind++) {
      if($exon_usage_comp_arr[$pos_ind - 1] != $exon_usage_comp_arr[$pos_ind]) {
         #===2->1 transitions, A5SS/RI, using "|" symbol===
         if     ($exon_usage_comp_arr[$pos_ind - 1] == 2 && $exon_usage_comp_arr[$pos_ind] == 1) {
               $alt_event_string .= ':' . $exon_bound_arr[$pos_ind - 1];
         } elsif($exon_usage_comp_arr[$pos_ind - 1] == 1 && $exon_usage_comp_arr[$pos_ind] == 2) {
               $alt_event_string .= '|' . $exon_bound_arr[$pos_ind - 1];
         #===0 to 1/2, starting new exon, using "@" symbol===
         } elsif($exon_usage_comp_arr[$pos_ind - 1] == 0 &&  $exon_usage_comp_arr[$pos_ind] == 1) {
            $alt_event_string .= '@' . $exon_bound_arr[$pos_ind - 1];
         } elsif($exon_usage_comp_arr[$pos_ind - 1] == 0 &&  $exon_usage_comp_arr[$pos_ind] == 2) {
            $alt_event_string .= '@' . $exon_bound_arr[$pos_ind - 1];
         #===1/2 to 0, ending exon, using ":" symbol===
         } elsif($exon_usage_comp_arr[$pos_ind - 1] == 1 &&  $exon_usage_comp_arr[$pos_ind] == 0) {
               $alt_event_string .= '|' . $exon_bound_arr[$pos_ind - 1];
         } elsif($exon_usage_comp_arr[$pos_ind - 1] == 2 &&  $exon_usage_comp_arr[$pos_ind] == 0) {
               $alt_event_string .= ':' . $exon_bound_arr[$pos_ind - 1];
         } else {
            die "=======Fatal error in generating GFF ID!!======\n";
         }
      }
   }
      #===adding the right most exon boundary===
   if($exon_usage_comp_arr[$AS_end_ind] == 2) {
      $alt_event_string .= ':' . $exon_bound_arr[$AS_end_ind];
   } elsif($exon_usage_comp_arr[$AS_end_ind] == 1) {
      $alt_event_string .= '|' . $exon_bound_arr[$AS_end_ind];
   } else {
      #===do nothing for ATSS/ATTS===
   }
      #===polishing gff string===
   my @gff_string_arr = split('@', $alt_event_string);
   foreach(@gff_string_arr) {
      #===fix ATSS/ATTS 01+0 exons===
      $_ =~ s/^(\d+)\|(\d+)$/\1:\2/;
      #===adding chrom and strand===
      $_ = "$data_h->{$gene_id}->{chrom}:$_:$data_h->{$gene_id}->{strand}";
   }

   return(join('@', @gff_string_arr));
}

sub trace_back_left_right_boundary($$$) {
   my($exon_usage_arr, $AS_start_ind, $AS_end_ind) = @_;

   #===trace AS event left boundary===
   while(1) {
      #===move right if this position is at intron===
      if($exon_usage_arr->[$AS_start_ind] == '0') {
         $AS_start_ind++;
      #===move left if left position and this position are at exon===
      } elsif(defined $exon_usage_arr->[$AS_start_ind - 1] &&
                      $exon_usage_arr->[$AS_start_ind - 1] != '0') {
         $AS_start_ind--;
      } else {
         last;
      }
   }
   #===trace AS event right boundary===
   while(1) {
      #===move left if this position is at intron===
      if($exon_usage_arr->[$AS_end_ind] == '0') {
         $AS_end_ind--;
      #===move right if right position and this position are at exon===
      } elsif(defined $exon_usage_arr->[$AS_end_ind + 1] &&
                      $exon_usage_arr->[$AS_end_ind + 1] != '0') {
         $AS_end_ind++;
      } else {
         last;
      }
   }

   return($AS_start_ind, $AS_end_ind);
}

sub generate_gff_mRNA_exon($$$$$$$) {
   my($data_h, $gene_id, $tx_acc, $tx_type, $ae_type, $ae_id, $alt_event_id) = @_;

   my($AS_start_ind, $AS_end_ind) =
      trace_back_left_right_boundary(\@{$data_h->{$gene_id}->{exon_usage_flag}->{$tx_acc}},
         $data_h->{$gene_id}->{AE}->{$ae_type}->[$ae_id]->{start},
         $data_h->{$gene_id}->{AE}->{$ae_type}->[$ae_id]->{end});
   my @exon_boundaries_arr = sort {$a <=> $b} keys %{$data_h->{$gene_id}->{exons_boundaries}};

   my $t_id = ($tx_type eq 'A')?
              (join(',', @{$data_h->{$gene_id}->{AE}->{$ae_type}->[$ae_id]->{tx1}})):
              (join(',', @{$data_h->{$gene_id}->{AE}->{$ae_type}->[$ae_id]->{tx2}}));

   my $chrom_name = $data_h->{$gene_id}->{'chrom'};
   $chrom_name = "chr" . $chrom_name unless $chrom_name =~ /^chr/;

   #===gff mRNA===
   my $gff_string  = "$chrom_name\t$ae_type\tmRNA\t";
      $gff_string .= ($exon_boundaries_arr[$AS_start_ind - 1] + 1) . "\t";
      $gff_string .= "$exon_boundaries_arr[$AS_end_ind]\t";
      $gff_string .= ".\t$data_h->{$gene_id}->{'strand'}\t.\t";
      $gff_string .= "ID=$alt_event_id.$tx_type;Parent=$alt_event_id;gid=$gene_id;tid=$t_id\n";

   #===gff exon===
   my @exon_arr;
   for(my $pos_ind = $AS_start_ind; $pos_ind <= $AS_end_ind + 1; $pos_ind++) {
      #===exon start===
      if(     $data_h->{$gene_id}->{exon_usage_flag}->{$tx_acc}->[$pos_ind - 1] == 0 &&
              $data_h->{$gene_id}->{exon_usage_flag}->{$tx_acc}->[$pos_ind] == 1) {
         $exon_arr[$#exon_arr + 1]->{start} = $exon_boundaries_arr[$pos_ind - 1];
      } elsif($data_h->{$gene_id}->{exon_usage_flag}->{$tx_acc}->[$pos_ind - 1] == 1 &&
              $data_h->{$gene_id}->{exon_usage_flag}->{$tx_acc}->[$pos_ind] == 0) {
         $exon_arr[$#exon_arr]->{end} = $exon_boundaries_arr[$pos_ind - 1];
      }
   }

   for(my $exon_id = 0; $exon_id <= $#exon_arr; $exon_id++) {
      $gff_string .= "$chrom_name\t$ae_type\texon\t";
      $gff_string .= ($exon_arr[$exon_id]->{start} + 1) . "\t";
      $gff_string .= "$exon_arr[$exon_id]->{end}\t";
      $gff_string .= ".\t$data_h->{$gene_id}->{'strand'}\t.\t";
      $gff_string .= "ID=$alt_event_id.$tx_type";
      if($#exon_arr > 0) {
         if($exon_id == 0) {
            ($data_h->{$gene_id}->{'strand'} eq "+")?($gff_string .= ".up"):($gff_string .= ".dn");
         } elsif($exon_id == $#exon_arr) {
            ($data_h->{$gene_id}->{'strand'} eq "+")?($gff_string .= ".dn"):($gff_string .= ".up");
         } else {
            $gff_string .= ".${ae_type}${exon_id}";
         }
      }
      else {
         $gff_string .= ".${ae_type}";
      }
      $gff_string .= ";Parent=$alt_event_id.$tx_type;gid=$gene_id;tid=$t_id\n";
   }

   return $gff_string;
}

1;
