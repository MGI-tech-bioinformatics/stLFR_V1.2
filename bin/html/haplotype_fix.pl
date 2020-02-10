#!/usr/bin/perl -w
use strict;
die "perl $0 <in> "unless @ARGV==1;
my $num = 4;

open IN,$ARGV[0];
open OT,">$ARGV[0].tmp";
while(<IN>){
  chomp;
  if(/^chr\s+switch rate/){
    print OT "$_\n";
    next;
  }
  my @a = split;
  for(my $i = 1; $i < @a; $i++){
    if($i == 6){
      $a[$i] = int($a[$i]);
    }else{
      $a[$i] = int($a[$i] * 10**$num + 0.5) / 10**$num;
    }
  }
  my $out = join "\t", @a;
  print OT "$out\n";
}
close IN;
close OT;

`mv $ARGV[0].tmp $ARGV[0]`;
