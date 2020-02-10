#!/usr/bin/perl -w
use strict;

die "perl $0 <in> <name>"unless @ARGV==2;

my $od = $ARGV[0];
my $flag = 0;

open IN,"$ARGV[0]/$ARGV[1].varianttable.xls";
open OT1,">$od/$ARGV[1].varianttable.xls.2";
while(<IN>){
  chomp;
  my @a = split /\t/;
  next if /^CNV|^SV/;
  if(/^Sample/){
    for(my $i = 1; $i < @a; $i++){
      $flag = $i if $a[$i] eq $ARGV[1];
      print OT1 "$a[0]\t$a[$flag]\n";
      last;
    }
    next;
  }
  print OT1 "$a[0]\t$a[$flag]\n";
}
close IN;

my @cnv = (0) x 2;
if(-e "$od/../../file/$ARGV[1]/CNV/$ARGV[1].CNV.result.xls"){
  open IN,"$od/../../file/$ARGV[1]/CNV/$ARGV[1].CNV.result.xls";
  while(<IN>){
    chomp;
    next if !/^>/;
    my @a = split;
    $a[3] eq "DEL" ? $cnv[0]++ : $cnv[1]++;
  }
  close IN;
}
print OT1 "CNV deletion\t$cnv[0]\n";
print OT1 "CNV duplication\t$cnv[1]\n";

my @sv = (0) x 5;
if(-e "$od/../../file/$ARGV[1]/SV/$ARGV[1].SV.result.xls"){
  open IN,"$od/../../file/$ARGV[1]/SV/$ARGV[1].SV.result.xls";
  while(<IN>){
    chomp;
    next if /^#/;
    my @a = split;
    next unless $a[11] =~ /PASS/ && $a[12] =~ /PASS/;
    $a[10] eq "DEL" ? $sv[0]++ :
    $a[10] eq "DUP" ? $sv[1]++ :
    $a[10] =~ /INV/ ? $sv[2]++ : 
    $a[10] =~ /TRA/ ? $sv[3]++ : $sv[4]++;
  }
  close IN;
}
print OT1 "SV DEL\t$sv[0]\n";
print OT1 "SV DUP\t$sv[1]\n";
print OT1 "SV INV\t$sv[2]\n";
print OT1 "SV TRA\t$sv[3]\n";

close OT1;

`mv $od/$ARGV[1].varianttable.xls.2 $od/$ARGV[1].varianttable.xls`;

