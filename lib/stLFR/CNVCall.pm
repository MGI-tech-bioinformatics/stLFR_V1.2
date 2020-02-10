package stLFR::CNVCall;

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
  $python2, $bgzip, $tabix, $vcffilter, $vcffp, $lfrcnv, $lfrcnv_chr, $lfrcnv_ref, $minsize, $sp,
);

$bgzip     = "$TOOL_PATH/vcftools/bgzip";
$tabix     = "$TOOL_PATH/vcftools/tabix";
$lfrcnv    = "$TOOL_PATH/cnv/LFR-cnv";
$minsize   = 1000;
$sp        = 0.001;

sub cnvsv_cc{
  my (
    $sample,
    $alignprefix,
    $phaseprefix,
    $reference,
    $outputdir,
    $shelldir,
    $filedir,
    $prequeue,
    $monitortxt,
    $qsubsgetxt,
  ) = (@_);

  if($reference =~ /\/db\/reference\/hs37d5\/hs37d5.fa$/){
    $lfrcnv_chr = "N";
    $lfrcnv_ref = "GRCH37";
  }else{
    $lfrcnv_chr = "Y";
    $lfrcnv_ref = "GRCH37";
  }
  `mkdir -p $outputdir $shelldir $filedir`;

  open SUB,">$shelldir/sub041.$sample.cnv.sh";
  print SUB "
    rm -f $outputdir/cnv/$sample.CNV.result.xls
    export PYTHONPATH=\$PYTHONPATH:$TOOL_PATH/../lib/python2
    export LD_LIBRARY_PATH=$TOOL_PATH/cnv/lib:\$LD_LIBRARY_PATH
    export PATH=$TOOL_PATH/R/bin:/usr/local/bin:/usr/bin:\$PATH  

    mkdir -p $outputdir/tmp
    $lfrcnv \\
      -ncpu $$QUEUE_PM{s41cpu} \\
      -bam $alignprefix/$sample.sortdup.bqsr.bam \\
      -vcf $alignprefix/$sample.sortdup.bqsr.bam.HaplotypeCaller.vcf.gz \\
      -phase $phaseprefix/phasesplit \\
      -pname hapblock_$sample\_XXX \\
      -tmp $outputdir/tmp \\
      -sp $sp \\
      -out $outputdir/cnv \\
      -ref $lfrcnv_ref \\
      -chr $lfrcnv_chr \\
      -lcnv $minsize
      
    cp \\
      $outputdir/cnv/ALL.200.format.cnv.$minsize.highconfidence \\
      $outputdir/cnv/$sample.CNV.result.xls
  ";
  stLFR::ResultLink::filelink(
    "$sample.CNV.result.xls",
    "$outputdir/cnv",
    $filedir,
    \*SUB,
  );
  close SUB;

  if($prequeue){
    print $monitortxt "$shelldir/sub021.$sample.sh:$prequeue\t$shelldir/sub041.$sample.cnv.sh:$$QUEUE_PM{defqid}:$$QUEUE_PM{s41mem}G:$$QUEUE_PM{s41cpu}cpu\n";
    print $monitortxt "$shelldir/sub0322.$sample.phase.genome.sh:$$QUEUE_PM{defqid}:$$QUEUE_PM{s322mem}G:$$QUEUE_PM{s322cpu}cpu\t$shelldir/sub041.$sample.cnv.sh:$$QUEUE_PM{defqid}:$$QUEUE_PM{s41mem}G:$$QUEUE_PM{s41cpu}cpu\n";
  }else{
    #print $monitortxt "$shelldir/sub041.$sample.cnv.sh:$$QUEUE_PM{defqid}:$$QUEUE_PM{s41mem}G:$$QUEUE_PM{s41cpu}cpu\n";
  }
  print $qsubsgetxt "sh $shelldir/sub041.$sample.cnv.sh \n";

};

1;
