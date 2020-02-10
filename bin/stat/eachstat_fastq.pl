#!/usr/bin/perl -w
use strict;

die "perl $0 <sample> <indir> <output> <genomesize>" unless @ARGV == 4;

my %stat;
my $fastqtable = "Raw reads,Raw bases(bp),Clean reads,Clean bases(bp),Q20,Q30,Total depth";
my $fragtable = "Total barcode type,Barcode number,Barcode type rate,Reads pair number,Reads pair number(after split),Barcode split rate,Split barcode(G)";

# raw fastq
print STDERR "processing $ARGV[0] raw fastq info ......\n";
open FQ,"$ARGV[1]/Basic_Statistics_of_Sequencing_Quality.txt";
while(<FQ>){
  chomp;
  next if /^#/ || !/\S/;
  my @b = split /\t/;
  if($b[0] =~ /Total number of reads/){
    $stat{"Raw reads"} += 2 * $1 if $b[1] =~ /^(\d+)\s+\(/;
    $stat{"Clean reads"} += 2 * $1 if $b[2] =~ /^(\d+)\s+\(/;
  }
  elsif($b[0] =~ /Read length/){
    $stat{"Raw read length"} = $b[2];
  }
  elsif($b[0] =~ /Total number of bases/){
    $stat{"Raw bases(bp)"} += 2 * $1 if $b[1] =~ /^(\d+)\s+\(/;
    $stat{"Clean bases(bp)"} += 2 * $1 if $b[2] =~ /^(\d+)\s+\(/;
  }
  elsif($b[0] =~ /Q20/){
    my $v1 = $1 if $b[2] =~ /.*\((.*)%\)/;
    my $v2 = $1 if $b[4] =~ /.*\((.*)%\)/;
    $stat{"Q20"} = int(($v1 + $v2) / 2 * 100 + 0.5) / 100;
  }
  elsif($b[0] =~ /Q30/){
    my $v1 = $1 if $b[2] =~ /.*\((.*)%\)/;
    my $v2 = $1 if $b[4] =~ /.*\((.*)%\)/;
    print "30\t$v1\t$v2\n";
    $stat{"Q30"} = int(($v1 + $v2) / 2 * 100 + 0.5) / 100;
  }
}
close FQ;

# barcode split
print STDERR "processing $ARGV[0] barcode info ......\n";
open SPLIT,"$ARGV[1]/split_stat_read1.log";
while(<SPLIT>){
  chomp;
  if(/^Barcode_types = .* = (\d+)$/){
    $stat{"Total barcode type"} = $1;
  }
  elsif(/^Real_Barcode_types = (\d+) \((\S+) \%\)/){
    $stat{"Barcode number"} = $1;
    $stat{"Barcode type rate"} = int($2 * 100 + 0.5) / 100;
    $stat{"Barcode type rate"} .= "%";
  }
  elsif(/^Reads_pair_num\s+= (\d+)/){
    $stat{"Reads pair number"} = $1;
  }
  elsif(/^Reads_pair_num\(after split\) = (\d+) \((\S+) \%\)/){
    $stat{"Reads pair number(after split)"} = $1;
    $stat{"Barcode split rate"} = int($2 * 100 + 0.5) / 100;
    $stat{"Barcode split rate"} .= "%";
  }
}
close SPLIT;

$stat{"Total depth"} = int($stat{"Raw bases(bp)"} / $ARGV[3] * 100 + 0.5) / 100;
$stat{"Split barcode(G)"} = $stat{"Reads pair number(after split)"} * $stat{"Raw read length"} * 2 / 10**9;
$stat{"Split barcode(G)"} = int($stat{"Split barcode(G)"} * 100 + 0.5) / 100;

open OT,">$ARGV[2]/$ARGV[0].fastqtable.xls";
print OT "Sample name\t$ARGV[0]\n";
foreach my $header (split /\,/, $fastqtable){
  print OT "$header\t$stat{$header}\n";
}
close OT;

open OT,">$ARGV[2]/$ARGV[0].fragtable.xls";
print OT "Sample name\t$ARGV[0]\n";
foreach my $header (split /\,/, $fragtable){
  print OT "$header\t$stat{$header}\n";
}
close OT;
