package stLFR::LowQualityFilter;

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
  $soapnuke, $adapter1_default, $adapter2_default, $soapnukepa_default,
);

$soapnuke           = "$TOOL_PATH/SOAPnuke/SOAPnuke";
$adapter1_default   = "-f CTGTCTCTTATACACATCTTAGGAAGACAAGCACTGACGACATGA";
$adapter2_default   = "-r TCTGCTGAGTCGAGAACGTCTCTGTGAGCCAAGGAGTTGCTCTGG";
$soapnukepa_default = "-l 10 -q 0.1 -n 0.01 -Q 2 -G -T 4";

sub filter_lqf{
  my (
    $sample,
    $inputdir,
    $outputdir,
    $shelldir,
    $filedir,
    $soapnukepa,
    $adapter1,
    $adapter2,
    $monitortxt,
    $qsubsgetxt,
  ) = (@_);

  `mkdir -p $outputdir $shelldir $filedir`;

  open SUB,">$shelldir/sub013.$sample.sh";
  print SUB "
    $soapnuke filter \\
      $soapnukepa \\
      -f $adapter1 -r $adapter2 \\
      -1 $inputdir/split_read.1.fq.gz \\
      -2 $inputdir/split_read.2.fq.gz \\
      -o $outputdir \\
      -C $sample.clean_1.fq.gz \\
      -D $sample.clean_2.fq.gz
    rm $inputdir/split_read.*.fq.gz
  ";

  foreach my $lane (qw/1 2/){
    stLFR::ResultLink::filelink(
      "$sample.clean_$lane.fq.gz",
      $outputdir,
      $filedir,
      \*SUB,
    );
  }
  close SUB;

  print $monitortxt "$shelldir/sub012.$sample.sh:$$QUEUE_PM{defqid}:$$QUEUE_PM{s12mem}G:$$QUEUE_PM{s12cpu}cpu\t$shelldir/sub013.$sample.sh:$$QUEUE_PM{defqid}:$$QUEUE_PM{s13mem}G:$$QUEUE_PM{s13cpu}cpu\n";
  print $qsubsgetxt "sh $shelldir/sub013.$sample.sh \n";

};

1;
