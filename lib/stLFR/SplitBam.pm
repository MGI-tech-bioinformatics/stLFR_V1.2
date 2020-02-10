package stLFR::SplitBam;

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
  $samtools,
);

$samtools = "$TOOL_PATH/samtools/bin/samtools";

sub align_sb{
  my (
    $sample,
    $bam,
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

    open SUB,">$shelldir/sub022.$sample.splitbam.$fai[0].sh";
    print SUB "
      export LD_LIBRARY_PATH=$TOOL_PATH/cnv/lib:\$LD_LIBRARY_PATH

      $samtools view \\
        --threads $$QUEUE_PM{s31cpu} -h -F 0x400 \\
        $bam \\
        $fai[0] \\
      | awk -v OFS='\\t' '{if(\$1~/#/){split(\$1,a,\"#\"); if(a[2]!~/0_0_0/){print \$0,\"BX:Z:\"a[2]} }else{print}}' - \\
      | $samtools view \\
        --threads $$QUEUE_PM{s31cpu} -bhS - \\
      > $outputdir/$sample.$fai[0].bam

      $samtools index \\
        $outputdir/$sample.$fai[0].bam
    ";
    close SUB;

    if($prequeue){
      print $monitortxt "$shelldir/sub021.$sample.sh:$prequeue\t$shelldir/sub022.$sample.splitbam.$fai[0].sh:$$QUEUE_PM{defqid}:$$QUEUE_PM{s22mem}G:$$QUEUE_PM{s22cpu}cpu\n";
    }else{
      #print $monitortxt "$shelldir/sub022.$sample.splitbam.$fai[0].sh:$$QUEUE_PM{defqid}:$$QUEUE_PM{s22mem}G:$$QUEUE_PM{s22cpu}cpu\n";
    }

    print $qsubsgetxt "sh $shelldir/sub022.$sample.splitbam.$fai[0].sh \n";
  }
  close FAI;

};

1;
