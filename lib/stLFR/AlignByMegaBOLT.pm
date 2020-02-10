package stLFR::AlignByMegaBOLT;

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
  $bolt, $boltbwa, $megaboltparameter,
);

$boltbwa = "/ldfssz1/MGI_ALGORITHM/MegaBOLT/source/V1.5.6.01/MegaBOLT/bin/bwa";
$bolt    = "/ldfssz1/MGI_ALGORITHM/MegaBOLT/V1.5.6.01/MegaBOLT";

sub align_abm{
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
    $bwa,
    $deepvariant,
  ) = (@_);

  `mkdir -p $outputdir $shelldir $filedir1 $filedir2`;

  open SUB,">$shelldir/sub021.$sample.sh";
  print SUB "
    echo -e \\
      \"$sample,$sample,$sample,COMPLETE,$fqprefix/$sample.clean_1.fq.gz,$fqprefix/$sample.clean_2.fq.gz\" \\
      > $outputdir/$sample.boltlist

    source /etc/bashrc
    export LD_LIBRARY_PATH=/mnt/ssd/MegaBOLT/bin/:\$LD_LIBRARY_PATH
  ";
  if($bwa){
    print SUB "
    cd $outputdir
    ln -s $reference reference.fa
    $boltbwa index reference.fa
    ";
    $megaboltparameter = "-bwa -ref $outputdir/reference.fa -stand_call_conf 10";
  }elsif($deepvariant){
    print SUB "
    cd $outputdir
    ln -s $reference reference.fa
    $boltbwa index reference.fa
    ";
    $megaboltparameter = "-bwa -deepvariant -fast_model 0 -model $TOOL_PATH/MegaBOLT/stLFR_HG001_sample1_2_HG004-200000/model.ckpt-190168_frozen.xml -ref $outputdir/reference.fa -stand_call_conf 10";
  }else{
    $megaboltparameter = "--stLFR 1 -ref $reference -stand_call_conf 10";
  }
  if($dbsnp && $kgindel && $kgmills){
    print SUB "
    $bolt \\
      $megaboltparameter \\
      -list $outputdir/$sample.boltlist \\
      -vcf $dbsnp \\
      -outputprefix $sample \\
      -outputDir $outputdir \\
      -knownSites $kgindel \\
      -knownSites $kgmills
    ";
  }elsif($dbsnp){
    print SUB "
    $bolt \\
      $megaboltparameter \\
      -list $outputdir/$sample.boltlist \\
      -vcf $dbsnp \\
      -outputprefix $sample \\
      -outputDir $outputdir 
    ";
  }else{
    print SUB "
    $bolt \\
      $megaboltparameter \\
      -list $outputdir/$sample.boltlist \\
      -outputprefix $sample \\
      -outputDir $outputdir
    ";
  }

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
    print $monitortxt "$shelldir/sub013.$sample.sh:$$QUEUE_PM{defqid}:$$QUEUE_PM{s13mem}G:$$QUEUE_PM{s13cpu}cpu\t$shelldir/sub021.$sample.sh:$$QUEUE_PM{boltqid}:$$QUEUE_PM{boltmem}G:$$QUEUE_PM{boltcpu}cpu\n";
  }else{
    #print $monitortxt "$shelldir/sub021.$sample.sh:$$QUEUE_PM{boltqid}:$$QUEUE_PM{boltmem}G:$$QUEUE_PM{boltcpu}cpu \n";
  }
  print $qsubsgetxt "sh $shelldir/sub021.$sample.sh \n";

};

1;
