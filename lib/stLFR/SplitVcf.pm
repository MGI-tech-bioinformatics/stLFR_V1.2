package stLFR::SplitVcf;

use strict;
use warnings;
use FindBin;
use lib "$FindBin::RealBin";
use lib "$FindBin::RealBin/../";
use lib "$FindBin::RealBin/../lib/perl5";
use lib "$FindBin::RealBin/../perl5";
use lib "$FindBin::RealBin/perl5";
use MyModule::GlobalVar qw($REF_H_DBPATH $TOOL_PATH $QUEUE_PM);
use stLFR::ResultLink;

use base 'Exporter';
our @EXPORT = qw();
my (
  $bcftools,
);

$bcftools = "$TOOL_PATH/vcftools/bcftools";

sub phase_sv{
  my (
    $sample,
    $vcf,
    $reference,
    $outputdir,
    $shelldir,
    $prequeue,
    $monitortxt,
    $qsubsgetxt,
  ) = (@_);

  `mkdir -p $outputdir $shelldir`;

  open FAI,"$reference.fai";
  while(<FAI>){
    chomp;
    my @fai = split;
    next if $fai[0] =~ /^GL|NC|hs37d5|\_|MT|chrM/ 
         && $reference =~ /\/db\/reference\/hg19\/hg19.fa$|\/db\/reference\/hs37d5\/hs37d5.fa$/;

    open SUB,">$shelldir/sub031.$sample.splitvcf.$fai[0].sh";
    print SUB "
      $bcftools view \\
        -O v -g het -v snps \\
        $vcf \\
        $fai[0] \\
      > $outputdir/$sample.$fai[0].hetsnp.vcf
    ";
    close SUB;

    if($prequeue){
      print $monitortxt "$shelldir/sub021.$sample.sh:$prequeue\t$shelldir/sub031.$sample.splitvcf.$fai[0].sh:$$QUEUE_PM{defqid}:$$QUEUE_PM{s31mem}G:$$QUEUE_PM{s31cpu}cpu\n";
    }else{
      #print $monitortxt "$shelldir/sub031.$sample.splitvcf.$fai[0].sh:$$QUEUE_PM{defqid}:$$QUEUE_PM{s31mem}G:$$QUEUE_PM{s31cpu}cpu\n";
    }
    
    print $qsubsgetxt "sh $shelldir/sub031.$sample.splitvcf.$fai[0].sh \n";
  }
  close FAI;

};

1;
