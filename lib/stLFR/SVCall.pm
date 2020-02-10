package stLFR::SVCall;

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
my ($svdir);

$svdir = "$TOOL_PATH/sv";

sub cnvsv_sc{
  my (
    $sample,
    $alignprefix,
    $phaseprefix,
    $reference,
    $blacklist,
    $controllist,
    $humanornot,
    $outputdir,
    $shelldir,
    $filedir,
    $prequeue,
    $monitortxt,
    $qsubsgetxt,
  ) = (@_);

  my $palist = "";
  $palist .= " -bl $blacklist" if $blacklist;
  $palist .= " -cl $controllist" if $controllist;
  $palist .= " -human Y" if $humanornot;

  `mkdir -p $outputdir $shelldir $filedir`;

  open SUB,">$shelldir/sub041.$sample.sv.sh";
  print SUB "
    rm -f $outputdir/$sample.SV.result.xls
    export LD_LIBRARY_PATH=$svdir/lib:$TOOL_PATH/R/bin:\$LD_LIBRARY_PATH

    $svdir/LFR-sv \\
      -bam $alignprefix/$sample.sortdup.bqsr.bam \\
      -out $outputdir \\
      -ncpu $$QUEUE_PM{s41cpu} $palist \\
      -phase $phaseprefix/svsplit

    cp \\
      $outputdir/$sample.sortdup.bqsr.bam.final \\
      $outputdir/$sample.SV.result.xls
  ";
  stLFR::ResultLink::filelink(
    "$sample.SV.result.xls",
    $outputdir,
    $filedir,
    \*SUB,
  );
  close SUB;

  if($prequeue){
    print $monitortxt "$shelldir/sub021.$sample.sh:$prequeue\t$shelldir/sub041.$sample.sv.sh:$$QUEUE_PM{defqid}:$$QUEUE_PM{s41mem}G:$$QUEUE_PM{s41cpu}cpu\n";
  }else{
    #print $monitortxt "$shelldir/sub041.$sample.sv.sh:$$QUEUE_PM{defqid}:$$QUEUE_PM{s41mem}G:$$QUEUE_PM{s41cpu}cpu\n";
  }
  print $qsubsgetxt "sh $shelldir/sub041.$sample.sv.sh \n";

};

1;
