package stLFR::FastqDistribution;

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
  $fqcheck, $fqcheckdis, $esfastqpl, $genomesize,
);

$fqcheck    = "$TOOL_PATH/fqcheck/fqcheck33";
$fqcheckdis = "$TOOL_PATH/fqcheck/fqcheck_distribute.pl";
$esfastqpl  = "$TOOL_PATH/../bin/stat/eachstat_fastq.pl";

sub filter_fd{
  my (
    $sample,
    $inputdir,
    $outputdir,
    $shelldir,
    $reportdir,
    $reference,
    $monitortxt,
    $qsubsgetxt1,
    $qsubsgetxt2,
  ) = (@_);

  if($reference =~ /\/db\/reference\/hg19\/hg19.fa$|\/db\/reference\/hs37d5\/hs37d5.fa$/){
    $genomesize = 3 * 10**9;
  }else{
    $genomesize = `less ${reference}.fai | awk '{a+=\$2}END{print a}'`;
  }
  chomp $genomesize;
  `mkdir -p $outputdir $shelldir $reportdir`;

  # figure of base/qual distribution
  foreach my $lane (qw/1 2/){
    open SUB1,">$shelldir/sub0141.$sample.$lane.sh";
    print SUB1 "
      $fqcheck \\
        -r $inputdir/$sample.clean_$lane.fq.gz \\
        -c $outputdir/$sample.clean_$lane.fqcheck
    ";
    close SUB1;
  }
  open SUB2,">$shelldir/sub0142.$sample.sh";
  print SUB2 "
    $fqcheckdis \\
      $outputdir/$sample.clean_1.fqcheck \\
      $outputdir/$sample.clean_2.fqcheck \\
      -o $outputdir/$sample.Cleanfq.
  ";
  foreach my $file (qw/base qual/){
    stLFR::ResultLink::reportlink(
      "$sample.Cleanfq.$file.png",
      $outputdir,
      $reportdir,
      \*SUB2,
    );
  }
  close SUB2;

  print $monitortxt "$shelldir/sub013.$sample.sh:$$QUEUE_PM{defqid}:$$QUEUE_PM{s13mem}G:$$QUEUE_PM{s13cpu}cpu\t$shelldir/sub0141.$sample.1.sh:$$QUEUE_PM{defqid}:$$QUEUE_PM{s141mem}G:$$QUEUE_PM{s141cpu}cpu\n";
  print $monitortxt "$shelldir/sub013.$sample.sh:$$QUEUE_PM{defqid}:$$QUEUE_PM{s13mem}G:$$QUEUE_PM{s13cpu}cpu\t$shelldir/sub0141.$sample.2.sh:$$QUEUE_PM{defqid}:$$QUEUE_PM{s141mem}G:$$QUEUE_PM{s141cpu}cpu\n";
  print $monitortxt "$shelldir/sub0141.$sample.1.sh:$$QUEUE_PM{defqid}:$$QUEUE_PM{s141mem}G:$$QUEUE_PM{s141cpu}cpu\t$shelldir/sub0142.$sample.sh:$$QUEUE_PM{defqid}:$$QUEUE_PM{s142mem}G:$$QUEUE_PM{s142cpu}cpu\n";
  print $monitortxt "$shelldir/sub0141.$sample.2.sh:$$QUEUE_PM{defqid}:$$QUEUE_PM{s141mem}G:$$QUEUE_PM{s141cpu}cpu\t$shelldir/sub0142.$sample.sh:$$QUEUE_PM{defqid}:$$QUEUE_PM{s142mem}G:$$QUEUE_PM{s142cpu}cpu\n";
  print $qsubsgetxt1 "sh $shelldir/sub0141.$sample.1.sh \n";
  print $qsubsgetxt1 "sh $shelldir/sub0141.$sample.2.sh \n";
  print $qsubsgetxt2 "sh $shelldir/sub0142.$sample.sh \n";

  # fasta and frag table
  open SUB1,">$shelldir/sub0141.$sample.table.sh";
  print SUB1 "
    perl $esfastqpl \\
    $sample \\
    $inputdir \\
    $outputdir \\
    $genomesize
  ";
  foreach my $file (qw/fastqtable fragtable/){
    stLFR::ResultLink::reportlink(
      "$sample.$file.xls",
      $outputdir,
      $reportdir,
      \*SUB1
    );
  }
  close SUB1;

  print $monitortxt "$shelldir/sub013.$sample.sh:$$QUEUE_PM{defqid}:$$QUEUE_PM{s13mem}G:$$QUEUE_PM{s13cpu}cpu\t$shelldir/sub0141.$sample.table.sh:$$QUEUE_PM{defqid}:$$QUEUE_PM{s141mem}G:$$QUEUE_PM{s141cpu}cpu\n";
  print $qsubsgetxt1 "sh $shelldir/sub0141.$sample.table.sh \n";

};

1;
