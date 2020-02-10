#!/usr/bin/perl -w
use strict;
use FindBin qw($Bin);

die "perl $0 <name> <dir> <ref> <circos> <bcftools> "unless @ARGV==5;

my $name = $ARGV[0];
my $ref = $ARGV[2];
my $circos = $ARGV[3];
my $bcftools = $ARGV[4];
my $cd = "$ARGV[1]/tmp";
my %flag;
my $win = 1000000;
`mkdir -p $cd`;

my $header;
if($ref =~ /db\/reference\/hg19\/hg19.fa$|db\/reference\/hs37d5\/hs37d5.fa$/){
  $header = "
  karyotype = $circos/data/karyotype/karyotype.human.hg19.txt
  chromosomes = -hsY
  ";
}else{
  open IN,"$ref.fai";
  open OT,">$cd/karyotype.txt";
  while(<IN>){
    chomp;
    my @a = split;
    print OT "chr - $a[0] $a[0] 0 $a[1] $a[0]\n";
  }
  close IN;
  close OT;
  $header = "
  karyotype = $cd/karyotype.txt
  ";
}

  open MAIN,">$cd/circos.$name.conf";
  print MAIN "
  $header
  <<include $circos/etc/housekeeping.conf>>
  <<include $circos/etc/colors_fonts_patterns.conf>>
  <<include $Bin/circos.image.conf>>
  <<include $Bin/circos.ideogram.conf>>
  <<include $Bin/circos.ticks.conf>>
  chromosomes_units = 1000000
  ";

  # cnv and sv
  if($ref =~ /db\/reference\/hg19\/hg19.fa$|db\/reference\/hs37d5\/hs37d5.fa$/){
    if(-e "$ARGV[1]/../../file/$name/CNV/$name.CNV.result.xls"){
      open CNV,"$ARGV[1]/../../file/$name/CNV/$name.CNV.result.xls";
      open TXT,">$cd/cnv.$name.txt";
      while(<CNV>){
        chomp;
        next if /Start|^>/;
        my @a = split;
        if($a[0] =~ /^chr/){
          $a[0] =~ s/^chr/hs/;
        }else{
          $a[0] = "hs$a[0]";
        }
        print TXT "$a[0]\t$a[1]\t$a[2]\t$a[6]\n";
      }
      close CNV;
      close TXT;
    }else{
      `touch $cd/cnv.$name.txt`;
    }
    open CONF,">$cd/cnv.$name.conf";
    print CONF "
    <highlight>
    file = $cd/cnv.$name.txt
    <<include $Bin/circos.cnv.conf>>
    </highlight>
    ";
    close CONF;

    if(-e "$ARGV[1]/../../file/$name/SV/$name.SV.result.xls"){
      open SV,"$ARGV[1]/../../file/$name/SV/$name.SV.result.xls";
      open TXT,">$cd/sv.$name.txt";
      while(<SV>){
        chomp;
        my @a = split;
        next if /^#|^EventID/;
        next unless $a[11] =~ /PASS/ && $a[12] =~ /PASS/;
        if($a[4] =~ /^chr/){
          $a[4] =~ s/chr/hs/;
          $a[6] =~ s/chr/hs/;
        }else{
          $a[4] = "hs$a[4]";
          $a[6] = "hs$a[6]";
        }
        my $col = ($a[10] eq "DEL") ? 1 :
                  ($a[10] eq "DUP") ? 2 :
                  ($a[10] =~ /INV/) ? 3 :
                  ($a[10] =~ /TRA/) ? 4 :
                  5;
        print TXT "$a[4]\t$a[5]\t$a[5]\t$a[6]\t$a[7]\t$a[7]\tcolor=set2-5-qual-$col\n";
      }
      close TXT;
      close SV;
    }else{
      `touch $cd/sv.$name.txt`;
    }
    open CONF,">$cd/sv.$name.conf";
    print CONF "
    <link>
    file = $cd/sv.$name.txt
    <<include $Bin/circos.sv.conf>>
    </link>
    ";
    close CONF;

    print MAIN "
    <highlights>
    <<include $cd/cnv.$name.conf>>
    </highlights>
    <links>
    <<include $cd/sv.$name.conf>>
    </links>
    ";

  }

  # snp and indel
  treatvcf( "$ARGV[1]/../../file/$name/variant/$name.sortdup.bqsr.bam.HaplotypeCaller.vcf.gz", "snp", $name );
  treatvcf( "$ARGV[1]/../../file/$name/variant/$name.sortdup.bqsr.bam.HaplotypeCaller.vcf.gz", "indel", $name );
  print MAIN "
  <plots>
  <<include $cd/snp.$name.conf>>
  <<include $cd/indel.$name.conf>>
  </plots>
  ";

  # circos
  `$circos/bin/circos -conf $cd/circos.$name.conf -outputfile $name.circos.png -outputdir $cd`;

sub treatvcf{
  my ($vcf, $type, $name) = (@_);
  my %data;
  my $max = 0;
  open F,"$bcftools view -v ${type}s $vcf | grep -v ^# |";
  open TXT,">$cd/$type.$name.txt";
  while(<F>){
    chomp;
    my @a = split;
    my $key = int($a[1] / $win);
    my $chr = $a[0];
    $chr =~ s/^chr/hs/ if ($ref =~ /db\/reference\/hg19\/hg19.fa$/);
    $chr = "hs$chr" if ($ref =~ /db\/reference\/hs37d5\/hs37d5.fa$/);
    $data{$chr}{$key}++;
  }
  close F;
  foreach my $chr (sort keys %data){
    foreach my $pos (sort {$a <=> $b} keys %{$data{$chr}}){
      my $s = $pos * $win;
      my $e = $s + $win;
      my $r = $data{$chr}{$pos} / $win;
      $max = ($max > $r) ? $max : $r;
      print TXT "$chr\t$s\t$e\t$r\n";
    }
  }
  close TXT;
  open CONF,">$cd/$type.$name.conf";
  print CONF "
  <plot>
  file = $cd/$type.$name.txt
  max = $max
  <<include $Bin/circos.$type.conf>>
  <<include $Bin/circos.bg.conf>>
  </plot>
  ";
  close CONF;
};
