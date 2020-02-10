package stLFR::SampleBarcodeMerge;

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

sub filter_sbm{
  my (
    $sample, 
    $path, 
    $barcode, 
    $outputdir, 
    $shelldir, 
    $monitortxt, 
    $qsubsgetxt,
  ) = (@_);

  `mkdir -p $outputdir $shelldir`;

  my @lanepath = split /\:/, $path;
  my @samplebarcode = split /\:/, $barcode;
  my (@fq1, @fq2) = () x 2;
  my $fq;
  
  for(my $i = 0; $i < @lanepath; $i++){
    if($samplebarcode[$i] eq "0" && !-e "$lanepath[$i]/BarcodeStat.txt"){
      $fq = `ls $lanepath[$i]/*\_read_1.fq.gz`;
      chomp $fq;
      push @fq1, $fq;
      $fq =~ s/\_1.fq.gz/\_2.fq.gz/;
      push @fq2, $fq;
    }else{
      foreach my $barcodelist ( listGet( $lanepath[$i], $samplebarcode[$i] ) ){
        chomp $barcodelist;
        $fq = `ls $lanepath[$i]/*\_$barcodelist\_1.fq.gz`;
        chomp $fq;
        push @fq1, $fq;
        $fq =~ s/\_1.fq.gz/\_2.fq.gz/;
        push @fq2, $fq;
      }
    }
  }


  open  SUB1,">$shelldir/sub011.$sample.1.sh";
  print SUB1 "cat ".(join " ", @fq1)." > $outputdir/$sample.raw.1.fq.gz \n";
  close SUB1;
  open  SUB2,">$shelldir/sub011.$sample.2.sh";
  print SUB2 "cat ".(join " ", @fq2)." > $outputdir/$sample.raw.2.fq.gz \n";
  close SUB2;

  #print $monitortxt "$shelldir/sub011.$sample.1.sh:$$QUEUE_PM{defqid}:$$QUEUE_PM{s11mem}G:$$QUEUE_PM{s11cpu}cpu\n";
  #print $monitortxt "$shelldir/sub011.$sample.2.sh:$$QUEUE_PM{defqid}:$$QUEUE_PM{s11mem}G:$$QUEUE_PM{s11cpu}cpu\n";
  print $qsubsgetxt "sh $shelldir/sub011.$sample.1.sh \n";
  print $qsubsgetxt "sh $shelldir/sub011.$sample.2.sh \n";
  
};

sub listGet{
  my ($path, $barcode) = (@_);
  # A-B,C,D-E
  my @list = ();
  if($barcode eq "0"){
    if(-e "$path/BarcodeStat.txt"){
      foreach my $fq1 (`ls $path/*_1.fq.gz`){
        chomp $fq1;
        next if $fq1 =~ /unmap|undecoded|discarded/;
        push @list, $1 if $fq1 =~ /\_L0\d\_(\d+)\_1.fq.gz$/;
      }
    }else{
      push @list, "read";
    }
  }else{
    foreach my $x (split /\,/, $barcode){
      if($x =~ /\-/){
        my ($y, $z) = (split /\-/, $x)[0, 1];
        for(my $m = $y; $m <= $z; $m++){
          push @list, $m;
        }
      }else{
        push @list, $x;
      }
    }
  }

  return @list;
};

1;
