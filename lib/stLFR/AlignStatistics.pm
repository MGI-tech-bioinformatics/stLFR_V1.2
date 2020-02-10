package stLFR::AlignStatistics;

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
  $fileflag, $java, $picard, $bam2depth, $samtools, $bgzip, $tabix, $R, $bcftools, $circos, 
  $fragment1sh, $fragment2sh, $escovpl, $esdepthpl, $esaligncatpl, $esvcfpl,
);

$java        = "$TOOL_PATH/jre/bin/java";
$picard      = "$TOOL_PATH/picard/picard.jar";
$bam2depth   = "$TOOL_PATH/bam2depth/bam2depth";
$samtools    = "$TOOL_PATH/samtools/bin/samtools";
$bgzip       = "$TOOL_PATH/vcftools/bgzip";
$bcftools    = "$TOOL_PATH/vcftools/bcftools";
$circos      = "$TOOL_PATH/circos";
$tabix       = "$TOOL_PATH/vcftools/tabix";
$fragment1sh = "$TOOL_PATH/../bin/stat/eachstat_fragment_1.sh";
$fragment2sh = "$TOOL_PATH/../bin/stat/eachstat_fragment_2.sh";
$R           = "$TOOL_PATH/R/bin";
$escovpl     = "$TOOL_PATH/../bin/stat/eachstat_cov.pl";
$esdepthpl   = "$TOOL_PATH/../bin/stat/eachstat_depth.pl";
$esvcfpl     = "$TOOL_PATH/../bin/stat/eachstat_vcf.pl";
$esaligncatpl= "$TOOL_PATH/../bin/stat/eachstat_aligncat.pl";

sub align_as{
  my (
    $sample,
    $bam,
    $splitdir,
    $vcf,
    $reference,
    $outputdir,
    $shelldir,
    $reportdir,
    $prequeue,
    $monitortxt,
    $qsubsgetxt1,
    $qsubsgetxt2,
  ) = (@_);

  `mkdir -p $outputdir $shelldir $reportdir`;

  # 1. fragment statistics about stLFR-barcode in two-steps
  open FAI,"$reference.fai";
  while(<FAI>){
    chomp;
    my @fai = split;
    next if $fai[0] =~ /^GL|NC|hs37d5|\_|MT|chrM/ 
         && $reference =~ /\/db\/reference\/hg19\/hg19.fa$|\/db\/reference\/hs37d5\/hs37d5.fa$/;

    open SUB,">$shelldir/sub0231.$sample.fragment1.$fai[0].sh";
    print SUB "
      export LD_LIBRARY_PATH=$TOOL_PATH/cnv/lib:\$LD_LIBRARY_PATH
      $fragment1sh \\
        $splitdir/$sample.$fai[0].bam \\
        $outputdir \\
        $sample \\
        300000 5000 \\
        $samtools
    ";
    close SUB;

    print $qsubsgetxt1 "sh $shelldir/sub0231.$sample.fragment1.$fai[0].sh \n";
    print $monitortxt "$shelldir/sub022.$sample.splitbam.$fai[0].sh:$$QUEUE_PM{defqid}:$$QUEUE_PM{s22mem}G:$$QUEUE_PM{s22cpu}cpu\t$shelldir/sub0231.$sample.fragment1.$fai[0].sh:$$QUEUE_PM{defqid}:$$QUEUE_PM{s231mem}G:$$QUEUE_PM{s231cpu}cpu\n";
    # write step2 monitor because chromosome information: @fai
    print $monitortxt "$shelldir/sub0231.$sample.fragment1.$fai[0].sh:$$QUEUE_PM{defqid}:$$QUEUE_PM{s231mem}G:$$QUEUE_PM{s231cpu}cpu\t$shelldir/sub0232.$sample.fragment2.sh:$$QUEUE_PM{defqid}:$$QUEUE_PM{s232mem}G:$$QUEUE_PM{s232cpu}cpu\n";
  }
  close FAI;
  open SUB,">$shelldir/sub0232.$sample.fragment2.sh";
  print SUB "
    $fragment2sh \\
    $outputdir \\
    $sample \\
    5000 \\
    $R/Rscript
  ";
  foreach $fileflag (qw/frag_cov fraglen_distribution_min5000 frag_per_barcode/){
    stLFR::ResultLink::reportlink(
      "$sample.$fileflag.pdf",
      $outputdir,
      $reportdir,
      \*SUB,
    );
  }
  close SUB;
  print $qsubsgetxt2 "sh $shelldir/sub0232.$sample.fragment2.sh \n";

  # 2. samtools flagstat/stats, 
  open SUB,">$shelldir/sub0231.$sample.samtoolsflagstat.sh";
  print SUB "
  export LD_LIBRARY_PATH=$TOOL_PATH/cnv/lib:\$LD_LIBRARY_PATH
  $samtools flagstat $bam > $outputdir/$sample.sorted.bam.flagstat \n";
  close SUB;
  print $qsubsgetxt1 "sh $shelldir/sub0231.$sample.samtoolsflagstat.sh \n";

  open SUB,">$shelldir/sub0231.$sample.samtoolsstats.sh";
  print SUB "
  export LD_LIBRARY_PATH=$TOOL_PATH/cnv/lib:\$LD_LIBRARY_PATH
  $samtools stats $bam > $outputdir/$sample.sorted.bam.stats \n";
  close SUB;
  print $qsubsgetxt1 "sh $shelldir/sub0231.$sample.samtoolsstats.sh \n";

  # 2. mapping statistics
  open SUB,">$shelldir/sub0231.$sample.eachstat_cov.sh";
  print SUB "
  export LD_LIBRARY_PATH=$TOOL_PATH/cnv/lib:\$LD_LIBRARY_PATH
  perl $escovpl $bam $outputdir/$sample.sorted.bam.info_1.xls $samtools \n";
  close SUB;
  print $qsubsgetxt1 "sh $shelldir/sub0231.$sample.eachstat_cov.sh \n";

  # 2. coverage / depth
  open SUB,">$shelldir/sub0231.$sample.eachstat_depth.sh";
  if($reference =~ /\/db\/reference\/hg19\/hg19.fa$|\/db\/reference\/hs37d5\/hs37d5.fa$/){
    my $reftype = $reference =~ /hg19.fa$/ ? "hg19" : "hs37d5";
    my $nonfile = $$REF_H_DBPATH{"$reftype.non"};
    print SUB "perl $esdepthpl $bam $outputdir -hg $reftype -bd $bam2depth -n $nonfile -R $R -s $sample > $outputdir/$sample.sorted.bam.info_2.xls";
  }else{
    my $genomesize = `less ${reference}.fai | awk '{a+=\$2}END{print a}'`;chomp $genomesize;
    print SUB "perl $esdepthpl $bam $outputdir -l $genomesize -hg $sample -bd $bam2depth -R $R -s $sample > $outputdir/$sample.sorted.bam.info_2.xls\n";
  }
  foreach $fileflag (qw/Sequencing.depth.accumulation Sequencing.depth/){
    stLFR::ResultLink::reportlink(
      "$sample.$fileflag.pdf",
      $outputdir,
      $reportdir,
      \*SUB,
    );
  }
  close SUB;
  print $qsubsgetxt1 "sh $shelldir/sub0231.$sample.eachstat_depth.sh \n";

  # 3. Insertsize
  open SUB,">$shelldir/sub0231.$sample.Insertsize.sh";
  print SUB "export PATH=\$PATH:$R
    $java -Xms10g -Xmx10g -jar $picard CollectInsertSizeMetrics \\
      I=$bam \\
      O=$outputdir/$sample.Insertsize.metrics.txt \\
      H=$outputdir/$sample.Insertsize.pdf VALIDATION_STRINGENCY=SILENT
  ";
  foreach $fileflag (qw/Insertsize.metrics.txt Insertsize.pdf/){
    stLFR::ResultLink::reportlink(
      "$sample.$fileflag",
      $outputdir,
      $reportdir,
      \*SUB,
    );
  }
  close SUB;
  print $qsubsgetxt1 "sh $shelldir/sub0231.$sample.Insertsize.sh \n";

  # 3. GCbias
  open SUB,">$shelldir/sub0231.$sample.GCbias.sh";
  print SUB "export PATH=\$PATH:$R
    $java -Xms10g -Xmx10g -jar $picard CollectGcBiasMetrics \\
      R=$reference \\
      I=$bam \\
      O=$outputdir/$sample.GCbias.metrics.txt \\
      CHART=$outputdir/$sample.GCbias.pdf \\
      S=$outputdir/$sample.GCbias.summary_metrics.txt VALIDATION_STRINGENCY=SILENT
  ";
  foreach $fileflag (qw/GCbias.metrics.txt GCbias.pdf GCbias.summary_metrics.txt/){
    stLFR::ResultLink::reportlink(
      "$sample.$fileflag",
      $outputdir,
      $reportdir,
      \*SUB,
    );
  }
  close SUB;
  print $qsubsgetxt1 "sh $shelldir/sub0231.$sample.GCbias.sh \n";

  # 4. snp and indels statistics
  open SUB,">$shelldir/sub0231.$sample.eachstat_vcf.sh";
  print SUB "perl $esvcfpl $vcf $outputdir/$sample.varianttable.xls $bcftools $sample \n";
  stLFR::ResultLink::reportlink(
    "$sample.varianttable.xls",
    $outputdir,
    $reportdir,
    \*SUB,
  );  
  close SUB;
  print $qsubsgetxt1 "sh $shelldir/sub0231.$sample.eachstat_vcf.sh \n";

  # 2. aligntable fix
  open SUB,">$shelldir/sub0232.$sample.eachstat_aligncat.sh";
  print SUB " perl $esaligncatpl $outputdir $outputdir/$sample.aligntable.xls $sample \n";
  stLFR::ResultLink::reportlink(
    "$sample.aligntable.xls",
    $outputdir,
    $reportdir,
    \*SUB,
  );
  close SUB;
  print $qsubsgetxt2 "sh $shelldir/sub0232.$sample.eachstat_aligncat.sh \n";

  print $monitortxt "$shelldir/sub021.$sample.sh:$prequeue\t$shelldir/sub0231.$sample.samtoolsflagstat.sh:$$QUEUE_PM{defqid}:$$QUEUE_PM{s231mem}G:$$QUEUE_PM{s231cpu}cpu\n";
  print $monitortxt "$shelldir/sub021.$sample.sh:$prequeue\t$shelldir/sub0231.$sample.samtoolsstats.sh:$$QUEUE_PM{defqid}:$$QUEUE_PM{s231mem}G:$$QUEUE_PM{s231cpu}cpu\n";
  print $monitortxt "$shelldir/sub021.$sample.sh:$prequeue\t$shelldir/sub0231.$sample.eachstat_cov.sh:$$QUEUE_PM{defqid}:$$QUEUE_PM{s231mem}G:$$QUEUE_PM{s231cpu}cpu\n";
  print $monitortxt "$shelldir/sub021.$sample.sh:$prequeue\t$shelldir/sub0231.$sample.eachstat_depth.sh:$$QUEUE_PM{defqid}:$$QUEUE_PM{s231mem}G:$$QUEUE_PM{s231cpu}cpu\n";
  print $monitortxt "$shelldir/sub021.$sample.sh:$prequeue\t$shelldir/sub0231.$sample.Insertsize.sh:$$QUEUE_PM{defqid}:$$QUEUE_PM{s231mem}G:$$QUEUE_PM{s231cpu}cpu\n";
  print $monitortxt "$shelldir/sub021.$sample.sh:$prequeue\t$shelldir/sub0231.$sample.GCbias.sh:$$QUEUE_PM{defqid}:$$QUEUE_PM{s231mem}G:$$QUEUE_PM{s231cpu}cpu\n";
  print $monitortxt "$shelldir/sub021.$sample.sh:$prequeue\t$shelldir/sub0231.$sample.eachstat_vcf.sh:$$QUEUE_PM{defqid}:$$QUEUE_PM{s231mem}G:$$QUEUE_PM{s231cpu}cpu\n";
  print $monitortxt "$shelldir/sub0231.$sample.samtoolsflagstat.sh:$$QUEUE_PM{defqid}:$$QUEUE_PM{s231mem}G:$$QUEUE_PM{s231cpu}cpu\t$shelldir/sub0232.$sample.eachstat_aligncat.sh:$$QUEUE_PM{defqid}:$$QUEUE_PM{s232mem}G:$$QUEUE_PM{s232cpu}cpu\n";
  print $monitortxt "$shelldir/sub0231.$sample.samtoolsstats.sh:$$QUEUE_PM{defqid}:$$QUEUE_PM{s231mem}G:$$QUEUE_PM{s231cpu}cpu\t$shelldir/sub0232.$sample.eachstat_aligncat.sh:$$QUEUE_PM{defqid}:$$QUEUE_PM{s232mem}G:$$QUEUE_PM{s232cpu}cpu\n";
  print $monitortxt "$shelldir/sub0231.$sample.eachstat_cov.sh:$$QUEUE_PM{defqid}:$$QUEUE_PM{s231mem}G:$$QUEUE_PM{s231cpu}cpu\t$shelldir/sub0232.$sample.eachstat_aligncat.sh:$$QUEUE_PM{defqid}:$$QUEUE_PM{s232mem}G:$$QUEUE_PM{s232cpu}cpu\n";
  print $monitortxt "$shelldir/sub0231.$sample.eachstat_depth.sh:$$QUEUE_PM{defqid}:$$QUEUE_PM{s231mem}G:$$QUEUE_PM{s231cpu}cpu\t$shelldir/sub0232.$sample.eachstat_aligncat.sh:$$QUEUE_PM{defqid}:$$QUEUE_PM{s232mem}G:$$QUEUE_PM{s232cpu}cpu\n";
  print $monitortxt "$shelldir/sub0231.$sample.Insertsize.sh:$$QUEUE_PM{defqid}:$$QUEUE_PM{s231mem}G:$$QUEUE_PM{s231cpu}cpu\t$shelldir/sub0232.$sample.eachstat_aligncat.sh:$$QUEUE_PM{defqid}:$$QUEUE_PM{s232mem}G:$$QUEUE_PM{s232cpu}cpu\n";
  print $monitortxt "$shelldir/sub0231.$sample.GCbias.sh:$$QUEUE_PM{defqid}:$$QUEUE_PM{s231mem}G:$$QUEUE_PM{s231cpu}cpu\t$shelldir/sub0232.$sample.eachstat_aligncat.sh:$$QUEUE_PM{defqid}:$$QUEUE_PM{s232mem}G:$$QUEUE_PM{s232cpu}cpu\n";

}

1;
