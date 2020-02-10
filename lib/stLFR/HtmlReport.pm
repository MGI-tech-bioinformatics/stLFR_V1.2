package stLFR::HtmlReport;

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
  $vfpl, $cpl, $convert, $size, $version, $htmlpy, $hfpl, $python3, $circos, $bcftools,
);

$vfpl     = "$TOOL_PATH/../bin/stat/variant_fix.pl";
$hfpl     = "$TOOL_PATH/../bin/html/haplotype_fix.pl";
$cpl      = "$TOOL_PATH/../bin/stat/circos.pl";
$htmlpy   = "$TOOL_PATH/../bin/html/stLFR/generate_DNApipe_report.py";
$convert  = "$TOOL_PATH/fqcheck/convert";
$python3  = "$TOOL_PATH/python3/python3";
$bcftools = "$TOOL_PATH/vcftools/bcftools";
$circos   = "$TOOL_PATH/circos";
$size     = "750x750";
$version  = "stLFR-reSeq V2.0.0.0";

sub report_hr{
  my (
    $sample,
    $fqprefix,
    $alignprefix,
    $phaseprefix,
    $cnvsvprefix,
    $reference,
    $reportdir,
    $shelldir,
    $monitortxt,
    $qsubsgetxt,
  ) = (@_);

  `mkdir -p $reportdir $shelldir`;

  open SUB,">$shelldir/sub051.$sample.sh";
  print SUB "export PATH=/ldfssz1/MGI_BIT/RUO/qiuwancen/local/bin:\$PATH \n";

  # fix variant table
  if($cnvsvprefix){
    print SUB "perl $vfpl $reportdir $sample \n";
  }

  # circos
  print SUB "perl $cpl $sample $reportdir $reference $circos $bcftools \n";
  print SUB "cp $TOOL_PATH/../bin/html/legend_circos.pdf $reportdir/$sample.legend_circos.pdf \n";
  foreach my $file (qw/svg png/){
    stLFR::ResultLink::reportlink(
      "$sample.circos.$file",
      "$reportdir/tmp",
      $reportdir,
      \*SUB,
    );
  }
  print SUB "rm -fr $reportdir/tmp \n";

  # haplotype fix
  if($phaseprefix){
    print SUB "perl $hfpl $reportdir/$sample.haplotype.xls \n";
  }

  # pdf to png
  my @pdf;
  if($fqprefix){
    push @pdf, "$sample.frag_cov.pdf";
    push @pdf, "$sample.fraglen_distribution_min5000.pdf";
    push @pdf, "$sample.frag_per_barcode.pdf";
  }
  if($alignprefix){
    push @pdf, "$sample.GCbias.pdf";
    push @pdf, "$sample.Insertsize.pdf";
    push @pdf, "$sample.Sequencing.depth.accumulation.pdf";
    push @pdf, "$sample.Sequencing.depth.pdf";
  }
  if($phaseprefix){
    push @pdf, "$sample.haplotype.pdf";
  }
  push @pdf, "$sample.legend_circos.pdf";
  foreach my $pdf (@pdf){
    chomp $pdf;
    next if $pdf !~ /.pdf/;
    my $png = $pdf;
    $png =~ s/.pdf/.png/;
    print SUB "$convert -resize $size $reportdir/$pdf $reportdir/$png \n";
  }

  # built html
  print SUB "
    $python3 \\
      $htmlpy \\
      \'$version\' \\
      $sample \\
      $reportdir/../ \\
      $reportdir/../
  ";

  close SUB;

  if($fqprefix){
    print $monitortxt "$shelldir/sub0142.$sample.sh:$$QUEUE_PM{defqid}:$$QUEUE_PM{s142mem}G:$$QUEUE_PM{s142cpu}cpu\t$shelldir/sub051.$sample.sh:$$QUEUE_PM{defqid}:$$QUEUE_PM{s51mem}G:$$QUEUE_PM{s51cpu}cpu\n";
  }
  if($alignprefix){
    print $monitortxt "$shelldir/sub0232.$sample.eachstat_aligncat.sh:$$QUEUE_PM{defqid}:$$QUEUE_PM{s232mem}G:$$QUEUE_PM{s232cpu}cpu\t$shelldir/sub051.$sample.sh:$$QUEUE_PM{defqid}:$$QUEUE_PM{s51mem}G:$$QUEUE_PM{s51cpu}cpu\n";
  }
  if($phaseprefix){
    print $monitortxt "$shelldir/sub0322.$sample.phase.genome.sh:$$QUEUE_PM{defqid}:$$QUEUE_PM{s322mem}G:$$QUEUE_PM{s322cpu}cpu\t$shelldir/sub051.$sample.sh:$$QUEUE_PM{defqid}:$$QUEUE_PM{s51mem}G:$$QUEUE_PM{s51cpu}cpu\n";
  }
  if($cnvsvprefix){
    print $monitortxt "$shelldir/sub041.$sample.sv.sh:$$QUEUE_PM{defqid}:$$QUEUE_PM{s41mem}G:$$QUEUE_PM{s41cpu}cpu\t$shelldir/sub051.$sample.sh:$$QUEUE_PM{defqid}:$$QUEUE_PM{s51mem}G:$$QUEUE_PM{s51cpu}cpu\n";
    print $monitortxt "$shelldir/sub041.$sample.cnv.sh:$$QUEUE_PM{defqid}:$$QUEUE_PM{s41mem}G:$$QUEUE_PM{s41cpu}cpu\t$shelldir/sub051.$sample.sh:$$QUEUE_PM{defqid}:$$QUEUE_PM{s51mem}G:$$QUEUE_PM{s51cpu}cpu\n";
  }
  print $qsubsgetxt "sh $shelldir/sub051.$sample.sh \n";

};

1;
