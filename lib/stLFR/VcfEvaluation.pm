package stLFR::VcfEvaluation;

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
  $fileflag, $gatk, $rtg, $bcftools, $tabix, @type, @filterp,
);

$gatk     = "$TOOL_PATH/gatk4/gatk";
$rtg      = "$TOOL_PATH/rtg-tools/rtg";
$bcftools = "$TOOL_PATH/vcftools/bcftools";
$tabix    = "$TOOL_PATH/vcftools/tabix";
@type       = qw/snp indel/;
$filterp[0] = "QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0";
$filterp[1] = "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0";

sub align_ve{
  my (
    $sample,
    $vcf,
    $reference,
    $sdf,
    $baseline,
    $confbed,
    $outputdir,
    $shelldir,
    $reportdir,
    $prequeue,
    $deepvariant,
    $monitortxt,
    $qsubsgetxt1,
    $qsubsgetxt2,
  ) = (@_);

  `mkdir -p $outputdir $shelldir $reportdir`;

  # hard filter
  for(my $i = 0; $i < @type; $i++){
    open SUB1,">$shelldir/sub0231.$sample.$type[$i].sh";
    if($deepvariant){
      print SUB1 "
      $bcftools view -O z -f PASS --type $type[$i]s $vcf > $outputdir/$sample.$type[$i].filter.pass.vcf.gz
      $tabix -p vcf -f $outputdir/$sample.$type[$i].filter.pass.vcf.gz
      $bcftools view -O z --type $type[$i]s $baseline > $outputdir/$sample.tp.$type[$i].vcf.gz
      $tabix -p vcf -f $outputdir/$sample.tp.$type[$i].vcf.gz
      ";
    }else{
      print SUB1 "
      export PATH=$TOOL_PATH/jre/bin:\$PATH

      $bcftools view -O z --type $type[$i]s $vcf > $outputdir/$sample.$type[$i].filter.pass.vcf.gz
      $tabix -p vcf -f $outputdir/$sample.$type[$i].filter.pass.vcf.gz
      $bcftools view -O z --type $type[$i]s $baseline > $outputdir/$sample.tp.$type[$i].vcf.gz
      $tabix -p vcf -f $outputdir/$sample.tp.$type[$i].vcf.gz
      ";
    }

    print SUB1 "
    rm -fr $outputdir/$type[$i]
    $rtg RTG_MEM=$$QUEUE_PM{bwamem}G vcfeval \\
      -b $outputdir/$sample.tp.$type[$i].vcf.gz \\
      -c $outputdir/$sample.$type[$i].filter.pass.vcf.gz \\
      -e $confbed \\
      -t $sdf \\
      -o $outputdir/$type[$i] \\
      2> $outputdir/$type[$i].log

    rm -f $outputdir/$sample.tp.$type[$i].vcf.gz*
    ";
    close SUB1;

    print $monitortxt "$shelldir/sub021.$sample.sh:$prequeue\t$shelldir/sub0231.$sample.$type[$i].sh:$$QUEUE_PM{defqid}:$$QUEUE_PM{bwamem}G:$$QUEUE_PM{s231cpu}cpu\n";
    print $qsubsgetxt1 "sh $shelldir/sub0231.$sample.$type[$i].sh \n";
  }

  # statistics table
  open SUB2,">$shelldir/sub0232.$sample.vcfeval.sh";
  print SUB2 "
    head -1 $outputdir/snp/summary.txt \\
      | awk 'BEGIN{OFS=\"\\t\"}{print \$1,\$3,\$4,\$5,\$6,\$7,\$8}' \\
      > $outputdir/$sample.evaluation.xls

    sed -i \"s/None/SNP/g\" $outputdir/snp/summary.txt
    cat $outputdir/snp/summary.txt \\
      | tail -1 \\
      | awk 'BEGIN{OFS=\"\\t\"}{print \$1,\$3,\$4,\$5,\$6,\$7,\$8}' \\
      >> $outputdir/$sample.evaluation.xls

    sed -i \"s/None/Indel/g\" $outputdir/indel/summary.txt 
    cat $outputdir/indel/summary.txt \\
      | tail -1 \\
      | awk 'BEGIN{OFS=\"\\t\"}{print \$1,\$3,\$4,\$5,\$6,\$7,\$8}' \\
      >> $outputdir/$sample.evaluation.xls
  ";
  stLFR::ResultLink::reportlink(
    "$sample.evaluation.xls",
    $outputdir,
    $reportdir,
    \*SUB2,
  );
  close SUB2;

  print $monitortxt "$shelldir/sub0231.$sample.snp.sh:$$QUEUE_PM{defqid}:$$QUEUE_PM{bwamem}G:$$QUEUE_PM{s231cpu}cpu\t$shelldir/sub0232.$sample.vcfeval.sh:$$QUEUE_PM{defqid}:$$QUEUE_PM{s232mem}G:$$QUEUE_PM{s232cpu}cpu\n";
  print $monitortxt "$shelldir/sub0231.$sample.indel.sh:$$QUEUE_PM{defqid}:$$QUEUE_PM{bwamem}G:$$QUEUE_PM{s231cpu}cpu\t$shelldir/sub0232.$sample.vcfeval.sh:$$QUEUE_PM{defqid}:$$QUEUE_PM{s232mem}G:$$QUEUE_PM{s232cpu}cpu\n";
  print $qsubsgetxt2 "sh $shelldir/sub0232.$sample.vcfeval.sh \n";

};

1;
