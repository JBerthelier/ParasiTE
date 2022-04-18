#! /usr/bin/env perl

# Author: Cheng-Kai Shiau <shiauck at gmail dot com>

use lib '.';

use warnings;
use strict;
use Getopt::Long 'HelpMessage';
use modules::parse_gtf_gff;
use modules::exon_usage;
use modules::detect_ae;
use modules::report_result;
use modules::utils;

#===cancel msg buffering===
$| = 1;

my $original_cmd = "$0 " . join(" ", @ARGV);

my $o_dir      = "alternative_event_annotation";
my $log_f      = "";
my $debug_mode = "";
my $format     = "auto";
GetOptions(
   'infile=s' => \ my $f_name,
   'outdir=s' => \    $o_dir,
   'format=s' => \    $format,
   'log=s'    => \    $log_f,
   'debug=s'  => \    $debug_mode,
   'help'     =>   sub { HelpMessage(0) },
) or HelpMessage(1);

#===die unless infile is given===
HelpMessage(1) unless $f_name;

die "\nError!!\nFile \"$f_name\" does not exist!!\n\n" unless -f $f_name;

#===format===
my $inp_format = "auto";
$inp_format = "auto" if($format =~ /auto/i || $format =~ /a/i);
$inp_format = "gff"  if($format =~ /gff/i  || $format =~ /gff3/i);
$inp_format = "gtf"  if($format =~ /gtf/i);

#===remove output folder===
print "\n\n";
unless(-d $o_dir) {
   print "Create folder: $o_dir\n";
   system("mkdir $o_dir");
} else {
   die "Folder \"$o_dir\" exists! Please give another one folder name to prevent accidentally re-write results!\n\n\n";
}

#===initialize log file===
if($log_f ne "") {
   open(LOG, ">>$log_f");
   print LOG "\n";
   print LOG "CMD: $original_cmd\n";
   my $datetime = localtime();
   print LOG "Timp stamp: $datetime\n";
   print LOG "\n";
   close(LOG);
}
#===initialize log file===

# parse gtf
my $data_h = parse_gtf_gff::read_gtf_gff($f_name, $inp_format);

# assumption 1: exons are non-overlapped regions
# assumption 2: all tx are on the same strand
# structure of $data_h
# $data_h->{$gene_id}->{'chrom'}
# $data_h->{$gene_id}->{'strand'}
# $data_h->{$gene_id}->{'transcripts'}->{$tx_id} = 1
# $data_h->{$gene_id}->{'exons_boundaries'}->{$exons_pos}->{$tx_id} = 1
# $data_h->{$gene_id}->{'exon_usage_flag'}->{$tx_id}->[$segments_by_exons_boundaries] = 0/1, 0 for intron, 1 for exon
# $data_h->{$gene_id}->{'exon_usage_comparison'}->{$tx_1}->{$tx_2}->[$segments_by_exons_boundaries] = 0/1/2, 0 for both intron, 1 for intron and exon, 2 for both exon
# $data_h->{$gene_id}->{AE}->{$ae_type}->[$ae_id]->{start} = # of starting segment by exons boundaries
# $data_h->{$gene_id}->{AE}->{$ae_type}->[$ae_id]->{end}   = # of ending segment by exons boundaries
# $data_h->{$gene_id}->{AE}->{$ae_type}->[$ae_id]->{tx1}   = acc of isoform 1
# $data_h->{$gene_id}->{AE}->{$ae_type}->[$ae_id]->{tx2}   = acc of isoform 2

# example of data structure:
# exon boundaries:
#   0th     1st     2nd     3rd     4th     5th
#   aaa     bbb     ccc     ddd     eee     fff
# exon usage array or exon usage comparison array:
# 0th	1st	2nd	3rd	4th	5th	6th
# 0	1	0	1	0	1	0
#^^^ initial value, always 0                   ^^^ always end with 0
# 0/1 flags are changed by having the exon boundary in the isoform

exon_usage::generate_exon_usage_profile(\%$data_h);

detect_ae::detect_ae(\%$data_h, $log_f);

report_result::report_result(\%$data_h, $o_dir, $log_f);

#==debug==
utils::show_exon_usage_profile_msg($data_h, $log_f) if $debug_mode ne "";
utils::show_exon_usage_comparison_msg($data_h, $log_f) if $debug_mode ne "";
utils::show_debug_msg($data_h, $log_f) if $debug_mode ne "";
utils::show_ae_msg(\%$data_h, $log_f) if $debug_mode ne "";
#==debug===





=head1 NAME

CATANA - 

=head1 SYNOPSIS

   --infile,-i   Input GTF/GFF file name
   --outdir,-o   (Optional) Output directory name, default: "alternative_event_annotation"
   --format,-f   Format of annotation, could be gff/gtf/auto, default: auto
   --log,   -l   Log file, default is on screen, e.g. standard output
   --debug, -d   Debug mode, default is off, to print debug message, set this parameter to "1"

=head1 VERSION

1.0

=cut
