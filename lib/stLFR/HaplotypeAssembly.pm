package stLFR::HaplotypeAssembly;

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
  $extractHAIRS, $python3, $linkfragment, $hapcut2, $compareblock, $linkdist, $R,
  $pv, @h1, @v1, @f1, @pv, $esphasepl, $hkpl, $svbin,
);

$extractHAIRS  = "$TOOL_PATH/HapCUT2-master/build/extractHAIRS";
$hapcut2       = "$TOOL_PATH/HapCUT2-master/build/HAPCUT2";
$linkfragment  = "$TOOL_PATH/HapCUT2-master/utilities/LinkFragments.py";
$compareblock  = "$TOOL_PATH/HapCUT2-master/utilities/calculate_haplotype_statistics.py";
$python3       = "$TOOL_PATH/python3/python3";
$linkdist      = 1000000;
$esphasepl     = "$TOOL_PATH/../bin/stat/eachstat_phase.pl";
$hkpl          = "$TOOL_PATH/../bin/stat/haplotype.karyotype.pl";
$R             = "$TOOL_PATH/R/bin";
$svbin         = "$TOOL_PATH/sv/tools/gen_phase";

sub phase_ha{
  my (
    $sample,
    $reference,
    $inputdir,
    $outputdir,
    $shelldir,
    $reportdir,
    $filedir,
    $monitortxt,
    $qsubsgetxt1,
    $qsubsgetxt2,
  ) = (@_);

  `mkdir -p $outputdir $shelldir $reportdir $filedir $outputdir/../svsplit`;

  open SUB2,">$shelldir/sub0322.$sample.phase.genome.sh";
  print SUB2 "rm -f $outputdir/../$sample.* \n";

  (@h1, @v1, @f1, @pv) = () x 4;
  open FAI,"$reference.fai";
  while(<FAI>){
    chomp;
    my @fai = split;
    next if $fai[0] =~ /^GL|NC|hs37d5|\_|MT|chrM/ 
         && $reference =~ /\/db\/reference\/hg19\/hg19.fa$|\/db\/reference\/hs37d5\/hs37d5.fa$/;

    if($reference =~ /\/db\/reference\/hg19\/hg19.fa$|\/db\/reference\/hs37d5\/hs37d5.fa$/){
      my $type = `basename $reference`;
      chomp $type;
      $type =~ s/.fa$//;
      $type = "$type.pvcf";
      $pv = "$$REF_H_DBPATH{$type}.$fai[0].vcf";
    }else{
      $pv = "$inputdir/$sample.$fai[0].hetsnp.vcf";
    }

    # phase on chromosome
    open SUB,">$shelldir/sub0321.$sample.phase.$fai[0].sh";
    print SUB "
      export LD_LIBRARY_PATH=$TOOL_PATH/sv/lib:\$LD_LIBRARY_PATH

      $extractHAIRS --10X 1 \\
        --bam $inputdir/$sample.$fai[0].bam \\
        --VCF $inputdir/$sample.$fai[0].hetsnp.vcf \\
        --out $outputdir/1.$sample.$fai[0].unlinked_frag \n
      $python3 \\
        $linkfragment \\
        --bam $inputdir/$sample.$fai[0].bam \\
        --vcf $inputdir/$sample.$fai[0].hetsnp.vcf \\
        --fragments $outputdir/1.$sample.$fai[0].unlinked_frag \\
        --out $outputdir/linked_fragment.$sample.$fai[0] \\
        -d $linkdist \n
      $hapcut2 --nf 1 \\
        --fragments $outputdir/linked_fragment.$sample.$fai[0] \\
        --vcf $inputdir/$sample.$fai[0].hetsnp.vcf \\
        --output $outputdir/hapblock_$sample\_$fai[0] \n
      echo $fai[0] > $outputdir/2.$sample.$fai[0].hapcut_stat.txt \n
      $python3 \\
        $compareblock \\
        -h1 $outputdir/hapblock_$sample\_$fai[0] \\
        -v1 $inputdir/$sample.$fai[0].hetsnp.vcf \\
        -f1 $outputdir/linked_fragment.$sample.$fai[0] \\
        -pv $pv \\
        -c $reference.fai \\
        >> $outputdir/2.$sample.$fai[0].hapcut_stat.txt

      $svbin/format_phase \\
        $outputdir/hapblock_$sample\_$fai[0] \\
        $outputdir/../svsplit/$fai[0].region \\
        $outputdir/../svsplit/$fai[0].vcf
      $svbin/get_barcode_from_phase \\
        $reference \\
        $inputdir/$sample.$fai[0].bam \\
        $outputdir/../svsplit/$fai[0].vcf \\
        $outputdir/../svsplit/$fai[0].barcode.phase
    ";
    close SUB;

    print $monitortxt "$shelldir/sub031.$sample.splitvcf.$fai[0].sh:$$QUEUE_PM{defqid}:$$QUEUE_PM{s31mem}G:$$QUEUE_PM{s31cpu}cpu\t$shelldir/sub0321.$sample.phase.$fai[0].sh:$$QUEUE_PM{defqid}:$$QUEUE_PM{s321mem}G:$$QUEUE_PM{s321cpu}cpu\n";
    if(-e "$shelldir/sub031.$sample.splitbam.$fai[0].sh"){
      print $monitortxt "$shelldir/sub031.$sample.splitbam.$fai[0].sh:$$QUEUE_PM{defqid}:$$QUEUE_PM{s31mem}G:$$QUEUE_PM{s31cpu}cpu\t$shelldir/sub0321.$sample.phase.$fai[0].sh:$$QUEUE_PM{defqid}:$$QUEUE_PM{s321mem}G:$$QUEUE_PM{s321cpu}cpu\n";
    }else{
      print $monitortxt "$shelldir/sub022.$sample.splitbam.$fai[0].sh:$$QUEUE_PM{defqid}:$$QUEUE_PM{s22mem}G:$$QUEUE_PM{s22cpu}cpu\t$shelldir/sub0321.$sample.phase.$fai[0].sh:$$QUEUE_PM{defqid}:$$QUEUE_PM{s321mem}G:$$QUEUE_PM{s321cpu}cpu\n"
    }
    print $qsubsgetxt1 "sh $shelldir/sub0321.$sample.phase.$fai[0].sh \n";

    # get information for genome
    print SUB2 "
      cat \\
        $outputdir/linked_fragment.$sample.$fai[0] \\
        >> $outputdir/../$sample.linked_fragment
      cat \\
        $outputdir/hapblock_$sample\_$fai[0] \\
        >> $outputdir/../$sample.hapblock
      cat \\
        $outputdir/2.$sample.$fai[0].hapcut_stat.txt \\
        >> $outputdir/../$sample.hapcut_stat.txt
    ";
    push @h1, "$outputdir/hapblock_$sample\_$fai[0]";
    push @f1, "$outputdir/linked_fragment.$sample.$fai[0]";
    push @v1, "$inputdir/$sample.$fai[0].hetsnp.vcf";
    push @pv, "$pv";

    # write monitortxt as chromosome
    print $monitortxt "$shelldir/sub0321.$sample.phase.$fai[0].sh:$$QUEUE_PM{defqid}:$$QUEUE_PM{s321mem}G:$$QUEUE_PM{s321cpu}cpu\t$shelldir/sub0322.$sample.phase.genome.sh:$$QUEUE_PM{defqid}:$$QUEUE_PM{s322mem}G:$$QUEUE_PM{s322cpu}cpu\n";

  }
  close FAI;

  print SUB2 "
    echo \"combine all chrs\" \\
      >> $outputdir/../$sample.hapcut_stat.txt
    $python3 \\
      $compareblock \\
      -h1 ".(join " ", @h1)." -v1 ".(join " ", @v1)." -f1 ".(join " ",@f1)." -pv ".(join " ",@pv)." -c $reference.fai \\
      >> $outputdir/../$sample.hapcut_stat.txt

    perl $esphasepl \\
      $sample \\
      $outputdir/../ \\
      $inputdir \\
      $outputdir/../
      
    perl $hkpl \\
      $sample \\
      $outputdir/../$sample.hapblock \\
      $outputdir/../ \\
      $reference \\
      $R
  ";
  foreach my $file (qw/hapblock
                       hapcut_stat.txt
                       linked_fragment/){
    stLFR::ResultLink::filelink(
      "$sample.$file",
      "$outputdir/../",
      $filedir,
      \*SUB2,
    );
  }
  foreach my $file (qw/haplotype.pdf haplotype.xls/){
    stLFR::ResultLink::reportlink(
      "$sample.$file",
      "$outputdir/../",
      $reportdir,
      \*SUB2,
    );
  }
  close SUB2;
  print $qsubsgetxt2 "sh $shelldir/sub0322.$sample.phase.genome.sh \n";

};

1;
