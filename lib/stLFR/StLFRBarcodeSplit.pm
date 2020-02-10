package stLFR::StLFRBarcodeSplit;

use strict;
use warnings;
use FindBin;
use lib "$FindBin::RealBin/../lib/perl5";
use lib "$FindBin::RealBin/../perl5";
use lib "$FindBin::RealBin/perl5";
use lib "$FindBin::RealBin";
use MyModule::GlobalVar qw($REF_H_DBPATH $TOOL_PATH $QUEUE_PM);

use base 'Exporter';
our @EXPORT = qw();
my (
  $barcodesplitpl, $barcode_position_default,
);

$barcodesplitpl           = "$TOOL_PATH/../bin/split_barcode_stLFR.pl";
$barcode_position_default = "101_10,117_10,133_10";

sub filter_sbs{
  my (
    $sample,
    $inputdir,
    $outputdir,
    $shelldir,
    $barcode_position,
    $monitortxt,
    $qsubsgetxt,
  ) = (@_);
  
  `mkdir -p $outputdir $shelldir`;
  
  open SUB,">$shelldir/sub012.$sample.sh";
  print SUB "
    perl $barcodesplitpl \\
      -i1 $inputdir/$sample.raw.1.fq.gz \\
      -i2 $inputdir/$sample.raw.2.fq.gz \\
      -r $barcode_position \\
      -o $outputdir
    rm $inputdir/$sample.raw.*.fq.gz
  ";
  close SUB;

  print $monitortxt "$shelldir/sub011.$sample.1.sh:$$QUEUE_PM{defqid}:$$QUEUE_PM{s11mem}G:$$QUEUE_PM{s11cpu}cpu\t$shelldir/sub012.$sample.sh:$$QUEUE_PM{defqid}:$$QUEUE_PM{s12mem}G:$$QUEUE_PM{s12cpu}cpu\n";
  print $monitortxt "$shelldir/sub011.$sample.1.sh:$$QUEUE_PM{defqid}:$$QUEUE_PM{s11mem}G:$$QUEUE_PM{s11cpu}cpu\t$shelldir/sub012.$sample.sh:$$QUEUE_PM{defqid}:$$QUEUE_PM{s12mem}G:$$QUEUE_PM{s12cpu}cpu\n";
  print $qsubsgetxt "sh $shelldir/sub012.$sample.sh \n";

};

1;