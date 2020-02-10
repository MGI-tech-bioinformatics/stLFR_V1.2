package stLFR::AlignByBWA;

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
  $bwa, $gatk, $java, $samtools, $bgzip, $tabix,
);

$bwa            = "$TOOL_PATH/bwa/bwa";
$gatk           = "$TOOL_PATH/gatk4/gatk";
$java           = "$TOOL_PATH/jre/bin/java";
$samtools       = "$TOOL_PATH/samtools/bin/samtools";
$bgzip          = "$TOOL_PATH/vcftools/bgzip";
$tabix          = "$TOOL_PATH/vcftools/tabix";

sub align_abb{
  my (
    $sample,
    $fqprefix,
    $reference,
    $dbsnp,
    $kgindel,
    $kgmills,
    $outputdir,
    $shelldir,
    $filedir1,
    $filedir2,
    $monitortxt,
    $qsubsgetxt,
    $analysis,
  ) = (@_);

  `mkdir -p $outputdir $shelldir $filedir1 $filedir2`;

  open SUB,">$shelldir/sub021.$sample.sh";

  # bwa mem + sort + MarkDuplicates
  print SUB "
    export PATH=$TOOL_PATH/jre/bin:\$PATH
    export LD_LIBRARY_PATH=$TOOL_PATH/cnv/lib:\$LD_LIBRARY_PATH
    
    $bwa mem \\
      -t $$QUEUE_PM{bwacpu} \\
      -R \"\@RG\\tID:$sample\\tPL:COMPLETE\\tPU:COMPLETE\\tLB:COMPLETE\\tSM:$sample\" \\
      $reference \\
      $fqprefix/$sample.clean_1.fq.gz \\
      $fqprefix/$sample.clean_2.fq.gz \\
    | $samtools view \\
      -bhS --threads $$QUEUE_PM{bwacpu} \\
      -t $reference.fai \\
      -T $reference - \\
    | $samtools sort \\
      --threads $$QUEUE_PM{bwacpu} -m 1000000000 \\
      -T $outputdir/$sample.sort \\
      -o $outputdir/$sample.sort.bam -
    $samtools index \\
      $outputdir/$sample.sort.bam

    $gatk --java-options \"-Xmx$$QUEUE_PM{bwamem}G -Djava.io.tmpdir=$outputdir/tmpdir \" \\
      MarkDuplicates \\
      -I $outputdir/$sample.sort.bam \\
      -O $outputdir/$sample.sortdup.bam \\
      -M $outputdir/$sample.sortdup.mat
    $samtools index \\
      $outputdir/$sample.sortdup.bam
    rm $outputdir/$sample.sort.bam*
  ";

  # bqsr
  if($kgindel && $kgmills && $dbsnp){
    print SUB "
      $gatk --java-options \"-Xmx$$QUEUE_PM{bwamem}G -Djava.io.tmpdir=$outputdir/tmpdir\" \\
        BaseRecalibrator \\
        -R $reference \\
        -I $outputdir/$sample.sortdup.bam \\
        --known-sites $kgindel \\
        --known-sites $kgmills \\
        --known-sites $dbsnp \\
        -O $outputdir/$sample.sortdup.recal.table
      $gatk --java-options \"-Xmx$$QUEUE_PM{bwamem}G -Djava.io.tmpdir=$outputdir/tmpdir\" \\
        ApplyBQSR \\
        -R $reference \\
        -I    $outputdir/$sample.sortdup.bam \\
        -bqsr $outputdir/$sample.sortdup.recal.table \\
        -O    $outputdir/$sample.sortdup.bqsr.bam
      $samtools index \\
        $outputdir/$sample.sortdup.bqsr.bam
      rm $outputdir/$sample.sortdup.bam*
    ";
  }elsif($dbsnp){
    print SUB "
      $gatk --java-options \"-Xmx$$QUEUE_PM{bwamem}G -Djava.io.tmpdir=$outputdir/tmpdir\" \\
        BaseRecalibrator \\
        -R $reference \\
        -I $outputdir/$sample.sortdup.bam \\
        --known-sites $dbsnp \\
        -O $outputdir/$sample.sortdup.recal.table
      $gatk --java-options \"-Xmx$$QUEUE_PM{bwamem}G -Djava.io.tmpdir=$outputdir/tmpdir\" \\
        ApplyBQSR \\
        -R $reference \\
        -I    $outputdir/$sample.sortdup.bam \\
        -bqsr $outputdir/$sample.sortdup.recal.table \\
        -O    $outputdir/$sample.sortdup.bqsr.bam
      $samtools index \\
        $outputdir/$sample.sortdup.bqsr.bam
      rm $outputdir/$sample.sortdup.bam*
    ";
  }else{
    print SUB "
      mv \\
        $outputdir/$sample.sortdup.bam \\
        $outputdir/$sample.sortdup.bqsr.bam
      mv \\
        $outputdir/$sample.sortdup.bam.bai \\
        $outputdir/$sample.sortdup.bqsr.bam.bai
    ";  
  }

  # HaplotypeCaller
  if($dbsnp){
    print SUB "
      $gatk --java-options \"-Xmx$$QUEUE_PM{bwamem}G -Djava.io.tmpdir=$outputdir/tmpdir\" \\
        HaplotypeCaller \\
        -R $reference \\
        -I $outputdir/$sample.sortdup.bqsr.bam \\
        --dbsnp $dbsnp \\
        -O $outputdir/$sample.sortdup.bqsr.bam.HaplotypeCaller.vcf
    ";
  }else{
    print SUB "
      $gatk --java-options \"-Xmx$$QUEUE_PM{bwamem}G -Djava.io.tmpdir=$outputdir/tmpdir\" \\
        HaplotypeCaller \\
        -R $reference \\
        -I $outputdir/$sample.sortdup.bqsr.bam \\
        -O $outputdir/$sample.sortdup.bqsr.bam.HaplotypeCaller.vcf
    ";
  }
  print SUB "
    $bgzip \\
      $outputdir/$sample.sortdup.bqsr.bam.HaplotypeCaller.vcf
    $tabix -p vcf \\
      $outputdir/$sample.sortdup.bqsr.bam.HaplotypeCaller.vcf.gz
  ";

  foreach my $file (qw/sortdup.bqsr.bam
                  sortdup.bqsr.bam.bai/){
    stLFR::ResultLink::filelink(
      "$sample.$file",
      $outputdir,
      $filedir1,
      \*SUB,
    );
  }
  foreach my $file (qw/sortdup.bqsr.bam.HaplotypeCaller.vcf.gz 
                  sortdup.bqsr.bam.HaplotypeCaller.vcf.gz.tbi/){
    stLFR::ResultLink::filelink(
      "$sample.$file",
      $outputdir,
      $filedir2,
      \*SUB,
    );
  }  
  close SUB;

  if($analysis =~ /all|filter/){
    print $monitortxt "$shelldir/sub013.$sample.sh:$$QUEUE_PM{defqid}:$$QUEUE_PM{s13mem}G:$$QUEUE_PM{s13cpu}cpu\t$shelldir/sub021.$sample.sh:$$QUEUE_PM{defqid}:$$QUEUE_PM{bwamem}G:$$QUEUE_PM{bwacpu}cpu\n";
  }else{
    #print $monitortxt "$shelldir/sub021.$sample.sh:$$QUEUE_PM{defqid}:$$QUEUE_PM{bwamem}G:$$QUEUE_PM{bwacpu}cpu \n";
  }
  print $qsubsgetxt "sh $shelldir/sub021.$sample.sh \n";

};

1;
