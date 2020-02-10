package stLFR;

our $VERSION       = "0.4.2";
our $VERSIOIN_DATE = "18 Nov 2019";

use strict;
use warnings;
use FindBin;
use lib "$FindBin::RealBin";
use lib "$FindBin::RealBin/../lib";
use lib "$FindBin::RealBin/../lib/perl5";
use lib "$FindBin::RealBin/perl5";
use MyModule::GlobalVar qw($REF_H_DBPATH $TOOL_PATH $QUEUE_PM);

use stLFR::ParameterCheck;
use stLFR::Directory;
use stLFR::ResultLink;
use stLFR::SampleBarcodeMerge;
use stLFR::StLFRBarcodeSplit;
use stLFR::LowQualityFilter;
use stLFR::FastqDistribution;
use stLFR::AlignByBWA;
use stLFR::AlignByMegaBOLT;
use stLFR::AlignStatistics;
use stLFR::SplitBam;
use stLFR::VcfEvaluation;
use stLFR::SplitVcf;
use stLFR::HaplotypeAssembly;
use stLFR::CNVCall;
use stLFR::SVCall;
use stLFR::HtmlReport;
use stLFR::CleanShell;

use base 'Exporter';
our @EXPORT = qw();
our (
  %OPT, %sampleflag, $qsubsgepl, $watchdogpl,
);

$qsubsgepl  = "$TOOL_PATH/monitor/qsub-sge.pl";
$watchdogpl = "$TOOL_PATH/monitor/watchDog.pl";

sub run{
  %OPT = ref $_[0] eq "HASH" ? %{$_[0]} : @_;

  stLFR::ParameterCheck::run_pc(%OPT);

  stLFR::Directory::run_d(%OPT);

  #=============================================================================#
  # open all file handles
  #=============================================================================#
  open MAIN1 ,">$OPT{outputdir}/00.shell/pipeline.sh";
  open MAIN2 ,">$OPT{outputdir}/00.shell/pipeline.txt";
  open SUB11 ,">$OPT{outputdir}/00.shell/t011.filter_samplebarcodemerge.sh";
  open SUB12 ,">$OPT{outputdir}/00.shell/t012.filter_stLFRbarcodesplit.sh";
  open SUB13 ,">$OPT{outputdir}/00.shell/t013.filter_lowqualityfilter.sh";
  open SUB141,">$OPT{outputdir}/00.shell/t0141.filter_fastqdistribution.sh";
  open SUB142,">$OPT{outputdir}/00.shell/t0142.filter_fastqdistribution.sh";
  open SUB21 ,">$OPT{outputdir}/00.shell/t021.align_fastqalign.sh";
  open SUB22 ,">$OPT{outputdir}/00.shell/t022.align_bamsplit.sh";
  open SUB231,">$OPT{outputdir}/00.shell/t0231.align_alignstat.sh";
  open SUB232,">$OPT{outputdir}/00.shell/t0232.align_alignstat.sh";
  open SUB31 ,">$OPT{outputdir}/00.shell/t031.phase_chrsplit.sh";
  open SUB321,">$OPT{outputdir}/00.shell/t0321.phase_phase.sh";
  open SUB322,">$OPT{outputdir}/00.shell/t0322.phase_phasestat.sh";
  open SUB41 ,">$OPT{outputdir}/00.shell/t041.cnvsv_call.sh";
  open SUB51 ,">$OPT{outputdir}/00.shell/t051.report_html.sh";

  open LIST,$OPT{'list'};
  while(<LIST>){
    # sam1  lane1:lane2 0:1-4,5,6-7
    # sam2  lane1       0              reference   dbsnp
    chomp;
    next if /^#|^$|^sample.*path/;
    my @tmp = split;
    $sampleflag{$tmp[0]}++;
    next if $sampleflag{$tmp[0]} > 1;

    #=============================================================================#
    # Description: 
    #    process stLFR fastq from MGI DNBSEQ platforms in 4 steps.
    #      SampleBarcodeMerge + 
    #      StLFRBarcodeSplit +
    #      LowQualityFilter +
    #      FastqDistribution
    # Input:
    #    sample name +
    #    fastq path +
    #    barcode information + 
    #    output directory +
    #    reference for coverage
    # Output:
    #    clean fastq files in FILE + 
    #    clean fastq figures in REPORT
    #=============================================================================#
    if($OPT{'analysis'} =~ /all|base|filter/){

      #=============================================================================#
      # Description: 
      #    merge all fastq files to one according to its sample-barcode
      # Barcode format:
      #    A-B,C : barcodes from A to B, and C belong to THIS sample
      #    0     : all barcodes belong to THIS sample, no matter split or not
      # Input:
      #    1. sample name
      #    2. fastq path(s), split by colon
      #    3. sample barcode(s), split by colon
      #    4. output directory
      #    5. shell directory
      #    6. monitor txt
      #    7. qsubsge txt
      # Output :
      #    raw fastq files in output directory named 'SAMPLE.raw.{1,2}.fq.gz'
      #=============================================================================#
      stLFR::SampleBarcodeMerge::filter_sbm(
        $tmp[0],
        $tmp[1],
        $tmp[2],
        "$OPT{outputdir}/01.filter/$tmp[0]",
        "$OPT{outputdir}/00.shell/fragment",
        \*MAIN2,
        \*SUB11,
      );

      #=============================================================================#
      # Description:
      #    split stLFR barocde for each raw fastq files
      # Barcode location format:
      #    the stLFR barcode location in read2: 101_10,117_10,133_10
      #    101_10 : barocde start on 101bp in read2 with 10bp length
      # Input:
      #    1. sample name
      #    2. input directory with raw fastq files in pre-step
      #    3. output directory
      #    4. shell directory
      #    5. barcode position
      #    6. monitor txt
      #    7. qsubsge txt 
      # Output:
      #    split fastq files in output directory named 'split_read.{1,2}.fq.gz'
      #=============================================================================#
      stLFR::StLFRBarcodeSplit::filter_sbs(
        $tmp[0],
        "$OPT{outputdir}/01.filter/$tmp[0]",
        "$OPT{outputdir}/01.filter/$tmp[0]",
        "$OPT{outputdir}/00.shell/fragment",
        $OPT{'position'},
        \*MAIN2,
        \*SUB12,
      );

      #=============================================================================#
      # Description:
      #    filter low quality read by SOAPnuke filter
      # Paramter in SOAPnuke filter:
      #    -l 10 -q 0.1 -n 0.01 -Q 2 -G -T 4
      #    -f CTGTCTCTTATACACATCTTAGGAAGACAAGCACTGACGACATGA 
      #    -r TCTGCTGAGTCGAGAACGTCTCTGTGAGCCAAGGAGTTGCTCTGG
      # Annotations of SOAPnuke filter:
      #    -l : low quality threshold
      #    -q : low quality rate
      #    -n : N rate threshold
      #    -Q : quality system 1:illumina, 2:sanger
      #    -G : set clean data quality system to sanger
      #    -T : threads
      #    -f : 3' adapter sequence of fq1 file
      #    -r : 5' adapter sequence of fq2 file
      # Input:
      #    1. sample name
      #    2. input directory with split fastq files in pre-step
      #    3. output directory
      #    4. shell directory
      #    5. FILE directory
      #    6. soapnuke parameter
      #    7. adapter1
      #    8. adapter2
      #    9. monitor txt
      #    10. qsubsge txt 
      # Output:
      #    clean fastq files in FILE for result-uploading
      #    named 'SAMPLE.clean_{1,2}.fq.gz'
      #=============================================================================#
      stLFR::LowQualityFilter::filter_lqf(
        $tmp[0],
        "$OPT{outputdir}/01.filter/$tmp[0]",
        "$OPT{outputdir}/01.filter/$tmp[0]",
        "$OPT{outputdir}/00.shell/fragment",
        "$OPT{outputdir}/file/$tmp[0]/sequence",
        '-l 10 -q 0.1 -n 0.01 -Q 2 -G -T 4',
        'CTGTCTCTTATACACATCTTAGGAAGACAAGCACTGACGACATGA',
        'TCTGCTGAGTCGAGAACGTCTCTGTGAGCCAAGGAGTTGCTCTGG',
        \*MAIN2,
        \*SUB13,
      );

      #=============================================================================#
      # Description:
      #    calculate fastq base and qual distribution within two steps
      # Note:
      #    These steps are a time-consuming process, so it will be done parallel
      #    The priority is: 
      #      align with bwa > phase > cnv > report
      # Input:
      #    1. sample name
      #    2. input directory with clean fastq files in pre-step
      #    3. output directory
      #    4. shell directory
      #    5. REPORT directory
      #    6. reference
      #    6. monitor txt
      #    7. qsubsge txt 1
      #    8. qsubsge txt 2
      # Output:
      #    base distribution figure in REPORT named 'SAMPLE.Cleanfq.base.png'
      #    quality distribution figure in REPORT named 'SAMPLE.Cleanfq.qual.png'
      #    fatq statistics table in REPORT named 'SAMPLE.fastqtable.xls'
      #    fragment statistics table in REPORT named 'SAMPLE.fragtable.xls'
      #=============================================================================#
      my $reference = "";
      if(@tmp == 3 || $tmp[3] =~ /^hg19$|^hs37d5$|^$/){
        $reference = @tmp == 3        ? $$REF_H_DBPATH{"hs37d5.fa"} :
                     $tmp[3] =~ /\S+/ ? $$REF_H_DBPATH{"$tmp[3].fa"} :
                                        $$REF_H_DBPATH{"hs37d5.fa"};
      }else{
        $reference = $tmp[3] if @tmp > 3 && $tmp[3] =~ /\S+/;
      }
      stLFR::FastqDistribution::filter_fd(
        $tmp[0],
        "$OPT{outputdir}/01.filter/$tmp[0]",
        "$OPT{outputdir}/01.filter/$tmp[0]",
        "$OPT{outputdir}/00.shell/fragment",
        "$OPT{outputdir}/report/$tmp[0]",
        $reference,
        \*MAIN2,
        \*SUB141,
        \*SUB142,
      );
    }

    #=============================================================================#
    # Description: 
    #    align the clean fastq to reference and follow with variant call
    # Input:
    #    clean fastq from pre-step or INPUTDIR +
    #    reference and dbsnp information from LIST table +
    #    output directory
    # Output:
    #    BAM and VCF in FILE + 
    #    statistics reports in REPORT
    #=============================================================================#
    if($OPT{'analysis'} =~ /all|base|align/){

      #=============================================================================#
      # Two methods are ready for alignment: BWA + GATK and MegaBOLT, which is 
      #  appointed by --software.
      # The reference and dbsnp are defined by 4th and 5th column in LIST:
      #    no 4th and 5th in LIST: use inner hs37d5 reference and dbsnp
      #    4th is hg19/hs37d5: use inner hg19/hs37d5 reference and dbsnp
      #    4th is fasta, 5th is vcf: use 4th as reference and 5th as dbsnp
      #    4th/5th is not hg19 or hs37d5: only phase but no CNV and SV
      # The input fastq file is detected according to --analysis
      #    1. --analysis contains all/filter : get fastq from previous
      #    2. --analysis not contains filter : get fastq from input --inputdir
      # The module about split BAM/VCF will be built in following modules, 
      #    because --inputdir contains BAM/VCF will not be split if built here
      #=============================================================================#
      my ($reference, $dbsnp, $kgindel, $kgmills, $fqprefix, $type, $prequeue, $sdf, 
         $baseline, $confbed) = 
         ("") x 10;
      $baseline = $OPT{'baseline'} if $OPT{'baseline'};
      $confbed  = $OPT{'confbed'} if $OPT{'confbed'};
      if(@tmp == 3 || $tmp[3] =~ /^hg19$|^hs37d5$|^$/){
        $type = @tmp == 3        ? "hs37d5" :
                $tmp[3] =~ /\S+/ ? $tmp[3] :
                                   "hs37d5";
        $reference = $$REF_H_DBPATH{"$type.fa"};
        $sdf       = $$REF_H_DBPATH{"$type.sdf"};
        $dbsnp     = $$REF_H_DBPATH{"$type.dbsnp"};
        $kgindel   = $$REF_H_DBPATH{"$type.1kgindel"};
        $kgmills   = $$REF_H_DBPATH{"$type.1kgmills"};
        $baseline  = $$REF_H_DBPATH{"$type.baseline"};
        $confbed   = $$REF_H_DBPATH{"$type.confbed"};
      }else{
        $reference = $tmp[3] if @tmp > 3 && $tmp[3] =~ /\S+/;
        $dbsnp     = $tmp[4] if @tmp > 4 && $tmp[4] =~ /\S+/;
      }
      $fqprefix =    $OPT{'inputdir'} 
                  && $OPT{'analysis'} !~ /all|base|filter/ 
                  ?  "$OPT{inputdir}/01.filter/$tmp[0]" 
                  :  "$OPT{outputdir}/01.filter/$tmp[0]";
      $prequeue = $OPT{'software'} eq "bwa" ?
                  "$$QUEUE_PM{defqid}:$$QUEUE_PM{bwamem}G:$$QUEUE_PM{bwacpu}cpu" :
                  "$$QUEUE_PM{boltqid}:$$QUEUE_PM{boltmem}G:$$QUEUE_PM{boltcpu}cpu";

      #=============================================================================#
      # Description:
      #    Porcess fastq with MegaBOLT, get alignment and variant vcf
      # Input:
      #    1. sample name
      #    2. fastq directory
      #    3. reference
      #    4. dbsnp
      #    ++. 1KG indels for BQSR if hg19/hs37d5
      #    ++. 1KG mills indels for BQSR if hg19/hs37d5
      #    5. output directory
      #    6. shell directory
      #    7. FILE directory for alignment
      #    8. FILE directory for variant
      #    9. monitor txt
      #    10. qsubsge txt
      #    11. analysis modules, used for monitor
      # Output:
      #    BAM file in FILE named 'SAMPLE.sortdup.bqsr.bam'
      #    VCF file in FILE named 'SAMPLE.sortdup.bqsr.bam.HaplotypeCaller.vcf.gz'
      #=============================================================================#
      stLFR::AlignByMegaBOLT::align_abm(
        $tmp[0],
        $fqprefix,
        $reference,
        $dbsnp,
        $kgindel,
        $kgmills,
        "$OPT{outputdir}/02.align/$tmp[0]",
        "$OPT{outputdir}/00.shell/fragment",
        "$OPT{outputdir}/file/$tmp[0]/alignment",
        "$OPT{outputdir}/file/$tmp[0]/variant",
        \*MAIN2,
        \*SUB21,
        $OPT{'analysis'},
        $OPT{'BWAinMegaBOLT'},
        $OPT{'DeepVariantinMegaBOLT'},
      ) if $OPT{'software'} =~ /MegaBOLT/i;

      #=============================================================================#
      # Description:
      #    Porcess fastq with BWA and GATK, BQSR only work if hg19/hs37d5
      # Input:
      #    1. sample name
      #    2. fastq directory
      #    3. reference
      #    4. dbsnp
      #    5. 1KG indels for BQSR if hg19/hs37d5
      #    6. 1KG mills indels for BQSR if hg19/hs37d5
      #    7. output directory
      #    8. shell directory
      #    9. FILE directory for alignment
      #    10. FILE directory for variant
      #    11. monitor txt
      #    12. qsubsge txt
      #    13. analysis modules, used for monitor
      # Output:
      #    BAM file in FILE named 'SAMPLE.sortdup.bqsr.bam'
      #    VCF file in FILE named 'SAMPLE.sortdup.bqsr.bam.HaplotypeCaller.vcf.gz'
      #=============================================================================#
      stLFR::AlignByBWA::align_abb(
        $tmp[0],
        $fqprefix,
        $reference,
        $dbsnp,
        $kgindel,
        $kgmills,
        "$OPT{outputdir}/02.align/$tmp[0]",
        "$OPT{outputdir}/00.shell/fragment",
        "$OPT{outputdir}/file/$tmp[0]/alignment",
        "$OPT{outputdir}/file/$tmp[0]/variant",
        \*MAIN2,
        \*SUB21,
        $OPT{'analysis'},
      ) if $OPT{'software'} =~ /BWA/i;

      #=============================================================================#
      # Description:
      #    split BAM of each chromosome for fragment statistics
      # Input:
      #    1. sample name
      #    2. BAM
      #    3. reference
      #    4. output directory
      #    5. shell directory
      #    6. pre-step queue parameter, used in monitor
      #    7. monitor txt
      #    8. qsubsge txt
      # Output:
      #    split bam in output directory for phase
      #=============================================================================#
      stLFR::SplitBam::align_sb(
        $tmp[0],
        "$OPT{outputdir}/02.align/$tmp[0]/$tmp[0].sortdup.bqsr.bam",
        $reference,
        "$OPT{outputdir}/02.align/$tmp[0]/alignsplit",
        "$OPT{outputdir}/00.shell/fragment",
        $prequeue,
        \*MAIN2,
        \*SUB22,
      );

      #=============================================================================#
      # Description:
      #    get statistics of BAM and VCF within several modules.
      # Modules:
      #    1. fragment statistics about stLFR-barcode in two-steps
      #    2. alignment statistics in two-steps and the first step contains 4 parts
      #    3. Insert size and GC bias distribution by Picard
      #    4. snp and indels statistics
      # Input:
      #    1. sample name
      #    2. BAM from pre-step
      #    3. BAM split directory
      #    4. VCF from pre-step
      #    5. reference
      #    6. output directory
      #    7. shell directory
      #    8. REPORT directory
      #    9. pre-step queue parameter, used in monitor
      #    10. monitor txt
      #    11. qsubsge txt 1
      #    12. qsubsge txt 2
      # Output:
      #    fragment figures in REPORT named 'SAMPLE.frag_cov.pdf', 
      #                                     'SAMPLE.fraglen_distribution_min5000.pdf',
      #                                     'SAMPLE.frag_per_barcode.pdf'
      #    coverage figures in REPORT named 'SAMPLE.Sequencing.depth.accumulation.pdf'
      #                                     'SAMPLE.Sequencing.depth.pdf'
      #    Insert size figure in REPORT named 'SAMPLE.Insertsize.pdf'
      #                txt in REPORT named    'SAMPLE.Insertsize.metrics.txt'
      #    GC bias figure in REPORT named 'SAMPLE.GCbias.pdf'
      #            txt in REPORT named    'SAMPLE.GCbias.metrics.txt'
      #                                   'SAMPLE.GCbias.summary_metrics.txt'
      #    alignment statistics table in REPORT named 'SAMPLE.aligntable.xls'
      #    variant statistics table in REPORT named 'SAMPLE.varianttable.xls'
      #=============================================================================#
      stLFR::AlignStatistics::align_as(
        $tmp[0],
        "$OPT{outputdir}/02.align/$tmp[0]/$tmp[0].sortdup.bqsr.bam",
        "$OPT{outputdir}/02.align/$tmp[0]/alignsplit",
        "$OPT{outputdir}/02.align/$tmp[0]/$tmp[0].sortdup.bqsr.bam.HaplotypeCaller.vcf.gz",
        $reference,
        "$OPT{outputdir}/02.align/$tmp[0]/stat",
        "$OPT{outputdir}/00.shell/fragment",
        "$OPT{outputdir}/report/$tmp[0]",
        $prequeue,
        \*MAIN2,
        \*SUB231,
        \*SUB232,
      );

      #=============================================================================#
      # Description:
      #    SNP/Indel evaluation based on baseline VCF
      # Input:
      #    1. sample name
      #    2. BAM from pre-step
      #    3. VCF from pre-step
      #    4. reference
      #    5. fastq directory, used for split_barcode
      #    6. output directory
      #    7. shell directory
      #    8. REPORT directory
      #    9. pre-step queue parameter, used in monitor
      #    ++. DeepVariantinMegaBOLT for PASS filter
      #    10. monitor txt
      #    11. qsubsge txt 1
      #    12. qsubsge txt 2
      # Output:
      #    evaluation table in REPORT named 'SAMPLE.evaluation.xls'
      #=============================================================================#
      stLFR::VcfEvaluation::align_ve(
        $tmp[0],
        "$OPT{outputdir}/02.align/$tmp[0]/$tmp[0].sortdup.bqsr.bam.HaplotypeCaller.vcf.gz",
        $reference,
        $sdf,
        $baseline,
        $confbed,
        "$OPT{outputdir}/02.align/$tmp[0]/stat",
        "$OPT{outputdir}/00.shell/fragment",
        "$OPT{outputdir}/report/$tmp[0]",
        $prequeue,
        $OPT{'DeepVariantinMegaBOLT'},
        \*MAIN2,
        \*SUB231,
        \*SUB232,
      ) if @tmp == 3 || $tmp[3] =~ /^hg19$|^hs37d5$|^$/;
    }

    #=============================================================================#
    # Description: 
    #    haplotype assembly based on BAM + het SNP of each chromosome by Hapcut2
    #    statistics will be done with baseline if reference is hg19/hs37d5
    # Input:
    #    BAM + VCF + 
    #    baseline VCF if hg19/hs37d5 +
    #    output directory
    # Output:
    #    haplotype statistics table in REPORT +
    #    haplotype figure in REPORT +
    #    haplotype files in FILE
    #=============================================================================#
    if($OPT{'analysis'} =~ /all|base|phase/){

      #=============================================================================#
      # For phase, reference is need for split BAM and VCF
      #   More, the align queue information is need for monitor
      # The BAM and VCF could get from INPUTDIR if -analysis is not contain all|align
      #=============================================================================#
      my ($reference, $type, $prequeue, $alignprefix) = ("") x 4;
      if(@tmp == 3 || $tmp[3] =~ /^hg19$|^hs37d5$|^$/){
        $type = @tmp == 3        ? "hs37d5" :
                $tmp[3] =~ /\S+/ ? $tmp[3] :
                                   "hs37d5";
        $reference = $$REF_H_DBPATH{"$type.fa"};
      }else{
        $reference = $tmp[3] if @tmp > 3 && $tmp[3] =~ /\S+/;
      }      
      $prequeue = $OPT{'analysis'} !~ /all|base|align/ ? ""
                : $OPT{'software'} eq "bwa" 
                ? "$$QUEUE_PM{defqid}:$$QUEUE_PM{bwamem}G:$$QUEUE_PM{bwacpu}cpu" 
                : "$$QUEUE_PM{boltqid}:$$QUEUE_PM{boltmem}G:$$QUEUE_PM{boltcpu}cpu";
      $alignprefix = $OPT{'inputdir'} 
                  && $OPT{'analysis'} !~ /all|base|align/ 
                  ?  "$OPT{inputdir}/02.align/$tmp[0]" 
                  :  "$OPT{outputdir}/02.align/$tmp[0]";

      #=============================================================================#
      # Description:
      #    split BAM of each chromosome or link split BAM from previous module
      # Input:
      #    1. sample name
      #    2. BAM
      #    3. reference
      #    4. output directory
      #    5. shell directory
      #    6. pre-step queue parameter, used in monitor
      #    7. monitor txt
      #    8. qsubsge txt
      # Output:
      #    split bam in output directory for phase
      #=============================================================================#                  
      if($OPT{'analysis'} !~ /all|base|align/){
        stLFR::SplitBam::align_sb(
          $tmp[0],
          "$alignprefix/$tmp[0].sortdup.bqsr.bam",
          $reference,
          "$OPT{outputdir}/03.phase/$tmp[0]/alignsplit",
          "$OPT{outputdir}/00.shell/fragment",
          $prequeue,
          \*MAIN2,
          \*SUB31,
        );
      }else{
        `mkdir -p $OPT{outputdir}/03.phase/$tmp[0]`;
        `ln -s ../../02.align/$tmp[0]/alignsplit $OPT{outputdir}/03.phase/$tmp[0]/`
         if !-e "$OPT{outputdir}/03.phase/$tmp[0]/alignsplit";
      }

      #=============================================================================#
      # Description:
      #    split heterozygous snp VCF of each chromosome as only het-snp is useful
      # Input:
      #    1. sample name
      #    2. VCF
      #    3. reference
      #    4. output directory
      #    5. shell directory
      #    6. pre-step queue parameter, used in monitor
      #    7. monitor txt
      #    8. qsubsge txt
      # Output:
      #    split vcf in output directory for phase      
      #=============================================================================#
      stLFR::SplitVcf::phase_sv(
        $tmp[0],
        "$alignprefix/$tmp[0].sortdup.bqsr.bam.HaplotypeCaller.vcf.gz",
        $reference,
        "$OPT{outputdir}/03.phase/$tmp[0]/alignsplit",
        "$OPT{outputdir}/00.shell/fragment",
        $prequeue,
        \*MAIN2,
        \*SUB31,
      );

      #=============================================================================#
      # Description:
      #    haplotype assembly and statistics on chromosome in step one,
      #    genome statistics and figure/table in step two
      # Note:
      #    statistics will be done with baseline if reference is hg19/hs37d5
      # Input:
      #    1. sample name
      #    2. reference
      #    3. input directory
      #    4. output directory
      #    5. shell directory
      #    6. REPORT directory
      #    7. FILE directory
      #    8. monitor txt
      #    9. qsubsge txt 1
      #    10. qsubsge txt 2
      # Output:
      #    haplotype block in FILE named 'SAMPLE.hapblock'
      #    statistics file in FILE named 'SAMPLE.hapcut_stat.txt'
      #    fragment file in FILE named 'SAMPLE.linked_fragment'
      #    statistics table in REPORT named 'SAMPLE.haplotype.xls'
      #    haplotype block figure in REPORT named 'SAMPLE.haplotype.pdf'
      #=============================================================================#
      stLFR::HaplotypeAssembly::phase_ha(
        $tmp[0],
        $reference,
        "$OPT{outputdir}/03.phase/$tmp[0]/alignsplit",
        "$OPT{outputdir}/03.phase/$tmp[0]/phasesplit",
        "$OPT{outputdir}/00.shell/fragment",
        "$OPT{outputdir}/report/$tmp[0]",
        "$OPT{outputdir}/file/$tmp[0]/haplotype",
        \*MAIN2,
        \*SUB321,
        \*SUB322,
      );
    }

    #=============================================================================#
    # Description: 
    #    CNV/SV calling based on BAM + Filtered VCF + haplotype block
    # Input:
    #    BAM + VCF + phasesplit +
    #    reference + output directory
    # Output:
    #    CNV xls result in FILE +
    #    SV xls result in FILE
    #=============================================================================#
    if($OPT{'analysis'} =~ /all|cnvsv/){
      
      #=============================================================================#
      # If no all|align|phase in -analysis, then -inputdir INPUTDIR shows
      #   the previous results
      #=============================================================================#
      my ($reference, $type, $prequeue, $alignprefix, $phaseprefix,
          $svblacklist, $svcontrollist, $svhumanornot) = ("") x 8;
      if(@tmp == 3 || $tmp[3] =~ /^hg19$|^hs37d5$|^\-$|^$/){
        $type = @tmp == 3        ? "hs37d5" :
                $tmp[3] =~ /\S+/ ? $tmp[3] :
                                   "hs37d5";
        $reference = $$REF_H_DBPATH{"$type.fa"};
        $svblacklist = $type eq "hs37d5" 
                     ? "$TOOL_PATH/sv/data/bl_region_hg19_nochr"
                     : "$TOOL_PATH/sv/data/bl_region_hg19_withchr";
        $svcontrollist = $type eq "hs37d5" 
                       ? "$TOOL_PATH/sv/data/con_hg19_nochr" 
                       : "$TOOL_PATH/sv/data/con_hg19_withchr";
        $svhumanornot = "Y";
      }else{
        $reference = $tmp[3] if @tmp > 3 && $tmp[3] =~ /\S+/ && $tmp[3] ne "-";
      }
      $prequeue = $OPT{'analysis'} !~ /all|base|align/ ? ""
                : $OPT{'software'} eq "bwa" 
                ? "$$QUEUE_PM{defqid}:$$QUEUE_PM{bwamem}G:$$QUEUE_PM{bwacpu}cpu" 
                : "$$QUEUE_PM{boltqid}:$$QUEUE_PM{boltmem}G:$$QUEUE_PM{boltcpu}cpu";
      $alignprefix = $OPT{'inputdir'} 
                  && $OPT{'analysis'} !~ /all|base|align/ 
                  ?  "$OPT{inputdir}/02.align/$tmp[0]" 
                  :  "$OPT{outputdir}/02.align/$tmp[0]";
      $phaseprefix = $OPT{'inputdir'} 
                  && $OPT{'analysis'} !~ /all|base|phase/
                  ?  "$OPT{inputdir}/03.phase/$tmp[0]"
                  :  "$OPT{outputdir}/03.phase/$tmp[0]";
      $svblacklist   = $tmp[5] if @tmp >= 6;
      $svcontrollist = $tmp[6] if @tmp == 7;

      #=============================================================================#
      # Description:
      #    CNV calling based on BAM, Filtered VCF and haplotype block
      # Input:
      #    1. sample name
      #    2. align directory
      #    3. phase directory
      #    4. reference      
      #    5. ouptut directory
      #    6. shell directory
      #    7. FILE directory
      #    8. pre-step queue parameter, used in monitor
      #    9. monitor txt
      #    10. qsubsge txt
      # Output:
      #    CNV file in FILE directory named 'SAMPLE.CNV.result.xls'
      #=============================================================================#
      stLFR::CNVCall::cnvsv_cc(
        $tmp[0],
        $alignprefix,
        $phaseprefix,
        $reference,
        "$OPT{outputdir}/04.cnv/$tmp[0]",
        "$OPT{outputdir}/00.shell/fragment",
        "$OPT{outputdir}/file/$tmp[0]/CNV",
        $prequeue,
        \*MAIN2,
        \*SUB41,
      );

      #=============================================================================#
      # Description:
      #    SV calling based on BAM, Filtered VCF and haplotype block
      # Input:
      #    1. sample name
      #    2. align directory
      #    ++. phase directory
      #    ++. reference      
      #    ++. black list
      #    ++. control list
      #    ++. human or not
      #    3. ouptut directory
      #    4. shell directory
      #    5. FILE directory
      #    6. pre-step queue parameter, used in monitor
      #    7. monitor txt
      #    8. qsubsge txt
      # Output:
      #    SV files in FILE directory named 'SAMPLE.SV.result.xls'
      #=============================================================================#
      stLFR::SVCall::cnvsv_sc(
        $tmp[0],
        $alignprefix,
        $phaseprefix,
        $reference,        
        $svblacklist,
        $svcontrollist,
        $svhumanornot,
        "$OPT{outputdir}/04.sv/$tmp[0]",
        "$OPT{outputdir}/00.shell/fragment",
        "$OPT{outputdir}/file/$tmp[0]/SV",
        $prequeue,
        \*MAIN2,
        \*SUB41,
      );
    }

    #=============================================================================#
    # Description: built html report
    #=============================================================================#
    if($OPT{'analysis'} =~ /all|base|report/){
      #=============================================================================#
      # Description:
      #    built html report
      # Input:
      #    1. sample name
      #    2. filter flag
      #    3. align flag
      #    4. phase flag
      #    5. cnvsv flag
      #    6. reference
      #    7. REPORT directory
      #    8. shell directory
      #    9. monitor txt
      #    10. qsubsge txt
      # Output:
      #    CIRCOS in REPORT named 'SAMPLE.circos.svg'
      #                           'SAMPLE.circos.png'
      #                           'SAMPLE.legend_circos.pdf'
      #    html in REPORT/../ named 'SAMPLE_cn.html'
      #                             'SAMPLE_en.html'
      #=============================================================================#
      my ($reference, $fqprefix, $alignprefix, $phaseprefix, $cnvsvprefix) = ("") x 5;
      $fqprefix = 1 if $OPT{'analysis'} =~ /all|base|filter/;
      $alignprefix = 1 if $OPT{'analysis'} =~ /all|base|align/;
      $phaseprefix = 1 if $OPT{'analysis'} =~ /all|base|phase/;
      $cnvsvprefix = 1 if $OPT{'analysis'} =~ /all|cnvsv/;  
      if(@tmp == 3 || $tmp[3] =~ /^hg19$|^hs37d5$|^$/){
        my $type = @tmp == 3        ? "hs37d5" :
                   $tmp[3] =~ /\S+/ ? $tmp[3] :
                                     "hs37d5";
        $reference = $$REF_H_DBPATH{"$type.fa"};
      }else{
        $reference = $tmp[3] if @tmp > 3 && $tmp[3] =~ /\S+/;
      }    
      stLFR::HtmlReport::report_hr(
        $tmp[0],
        $fqprefix,
        $alignprefix,
        $phaseprefix,
        $cnvsvprefix,
        $reference,
        "$OPT{outputdir}/report/$tmp[0]",
        "$OPT{outputdir}/00.shell/fragment",
        \*MAIN2,
        \*SUB51,
      );
    }

  }
  close LIST;

  #=============================================================================#
  # output qubsge txt
  #=============================================================================#
  if($OPT{'type'} =~ /blc/i){
    if($OPT{'task'} == 1 && $OPT{'analysis'} =~ /all|base|filter/){
      print MAIN1 "\necho ========== 1.filter start at : `date` ==========\n";
      print MAIN1 "perl $qsubsgepl --convert no --jobprefix s11 \\
          -resource=\"vf=$$QUEUE_PM{s11mem}G,num_proc=$$QUEUE_PM{s11cpu} \\
          -binding linear:$$QUEUE_PM{s11cpu} -P $OPT{Project} -q $OPT{queue}\" \\
          $OPT{outputdir}/00.shell/t011.filter_samplebarcodemerge.sh \n";
      print MAIN1 "perl $qsubsgepl --convert no --jobprefix s12 \\
          -resource=\"vf=$$QUEUE_PM{s12mem}G,num_proc=$$QUEUE_PM{s12cpu} \\
          -binding linear:$$QUEUE_PM{s12cpu} -P $OPT{Project} -q $OPT{queue}\" \\
          $OPT{outputdir}/00.shell/t012.filter_stLFRbarcodesplit.sh \n";
      print MAIN1 "perl $qsubsgepl --convert no --jobprefix s13 \\
          -resource=\"vf=$$QUEUE_PM{s13mem}G,num_proc=$$QUEUE_PM{s13cpu} \\
          -binding linear:$$QUEUE_PM{s13cpu} -P $OPT{Project} -q $OPT{queue}\" \\
          $OPT{outputdir}/00.shell/t013.filter_lowqualityfilter.sh \n";
      print MAIN1 "echo ========== 1.filter end   at : `date` ==========\n";
    }
    if($OPT{'task'} == 1 && $OPT{'analysis'} =~ /all|base|align/){
      print MAIN1 "\necho ========== 2.align start at : `date` ==========\n";
      print MAIN1 "perl $qsubsgepl --convert no --jobprefix s21 \\
          -resource=\"vf=$$QUEUE_PM{bwamem}G,num_proc=$$QUEUE_PM{bwacpu} \\
          -binding linear:$$QUEUE_PM{bwacpu} -P $OPT{Project} -q $OPT{queue}\" \\
          $OPT{outputdir}/00.shell/t021.align_fastqalign.sh \n" if $OPT{'software'} =~ /BWA/i;
      print MAIN1 "perl $qsubsgepl --convert no --jobprefix s21 \\
          -resource=\"vf=$$QUEUE_PM{boltmem}G,num_proc=$$QUEUE_PM{boltcpu} \\
          -binding linear:$$QUEUE_PM{boltcpu} -P $OPT{Project} -q $OPT{boltqueue}\" \\
          $OPT{outputdir}/00.shell/t021.align_fastqalign.sh \n" if $OPT{'software'} =~ /MegaBOLT/i;
      print MAIN1 "perl $qsubsgepl --convert no --jobprefix s22 \\
          -resource=\"vf=$$QUEUE_PM{s22mem}G,num_proc=$$QUEUE_PM{s22cpu} \\
          -binding linear:$$QUEUE_PM{s22cpu} -P $OPT{Project} -q $OPT{queue}\" \\
          $OPT{outputdir}/00.shell/t022.align_bamsplit.sh \n";        
      print MAIN1 "echo ========== 2.align end   at : `date` ==========\n";
    }
    if($OPT{'task'} == 1 && $OPT{'analysis'} =~ /all|base|phase/){
      print MAIN1 "\necho ========== 3.phase start at : `date` ==========\n";
      print MAIN1 "perl $qsubsgepl --convert no --jobprefix s31 \\
          -resource=\"vf=$$QUEUE_PM{s31mem}G,num_proc=$$QUEUE_PM{s31cpu} \\
          -binding linear:$$QUEUE_PM{s31cpu} -P $OPT{Project} -q $OPT{queue}\" \\
          $OPT{outputdir}/00.shell/t031.phase_chrsplit.sh \n";
      print MAIN1 "echo ========== 3.phase end   at : `date` ==========\n";
    }
    if($OPT{'task'} == 1 && $OPT{'analysis'} =~ /all|base|filter|align|phase/){
      print MAIN1 "\necho ========== 3.stat start at : `date` ==========\n";
      print MAIN1 "perl $qsubsgepl --convert no --jobprefix s321 \\
          -resource=\"vf=$$QUEUE_PM{s321mem}G,num_proc=$$QUEUE_PM{s321cpu} \\
          -binding linear:$$QUEUE_PM{s321cpu} -P $OPT{Project} -q $OPT{queue}\" \\
          $OPT{outputdir}/00.shell/t0321.phase_phase.sh \n";
      print MAIN1 "perl $qsubsgepl --convert no --jobprefix s322 \\
          -resource=\"vf=$$QUEUE_PM{s322mem}G,num_proc=$$QUEUE_PM{s322cpu} \\
          -binding linear:$$QUEUE_PM{s322cpu} -P $OPT{Project} -q $OPT{queue}\" \\
          $OPT{outputdir}/00.shell/t0322.phase_phasestat.sh \n";         
      print MAIN1 "echo ========== 3.stat end   at : `date` ==========\n";
    }
    if($OPT{'task'} == 1 && $OPT{'analysis'} =~ /all|cnvsv/){
      print MAIN1 "\necho ========== 4.cnvsv start at : `date` ==========\n";
      print MAIN1 "perl $qsubsgepl --convert no --jobprefix s41 \\
          -resource=\"vf=$$QUEUE_PM{s41mem}G,num_proc=$$QUEUE_PM{s41cpu} \\
          -binding linear:$$QUEUE_PM{s41cpu} -P $OPT{Project} -q $OPT{queue}\" \\
          $OPT{outputdir}/00.shell/t041.cnvsv_call.sh \n";
      print MAIN1 "echo ========== 4.cnvsv end   at : `date` ==========\n";    
    }
    if($OPT{'task'} == 1 && $OPT{'analysis'} =~ /all|base|report/){
      print MAIN1 "\necho ========== 5.html start at : `date` ==========\n";
      print MAIN1 "perl $qsubsgepl --convert no --jobprefix s51 \\
          -resource=\"vf=$$QUEUE_PM{s51mem}G,num_proc=$$QUEUE_PM{s51cpu} \\
          -binding linear:$$QUEUE_PM{s51cpu} -P $OPT{Project} -q $OPT{queue}\" \\
          $OPT{outputdir}/00.shell/t051.report_html.sh \n";
      print MAIN1 "echo ========== 5.html end   at : `date` ==========\n";  
    }
  }elsif($OPT{'type'} =~ /fpga|local/i){
    if($OPT{'task'} == 1 && $OPT{'analysis'} =~ /all|base|filter/){
      print MAIN1 "\necho ========== 1.filter start at : `date` ==========\n";
      print MAIN1 "perl $watchdogpl --mem $$QUEUE_PM{s11mem}G --num_paral $OPT{'cpu'} \\
          $OPT{outputdir}/00.shell/t011.filter_samplebarcodemerge.sh \n";
      print MAIN1 "perl $watchdogpl --mem $$QUEUE_PM{s12mem}G --num_paral $OPT{'cpu'} \\
          $OPT{outputdir}/00.shell/t012.filter_stLFRbarcodesplit.sh \n";
      print MAIN1 "perl $watchdogpl --mem $$QUEUE_PM{s13mem}G --num_paral $OPT{'cpu'} \\
          $OPT{outputdir}/00.shell/t013.filter_lowqualityfilter.sh \n";
      print MAIN1 "echo ========== 1.filter end   at : `date` ==========\n";
    }
    if($OPT{'task'} == 1 && $OPT{'analysis'} =~ /all|base|align/){
      print MAIN1 "\necho ========== 2.align start at : `date` ==========\n";
      print MAIN1 "perl $watchdogpl --mem $$QUEUE_PM{bwamem}G --num_paral $OPT{'cpu'} \\
          $OPT{outputdir}/00.shell/t021.align_fastqalign.sh \n" if $OPT{'software'} =~ /BWA/i;
      print MAIN1 "perl $watchdogpl --mem $$QUEUE_PM{boltmem}G --num_paral $OPT{'cpu'} \\
          $OPT{outputdir}/00.shell/t021.align_fastqalign.sh \n" if $OPT{'software'} =~ /MegaBOLT/i;
      print MAIN1 "perl $watchdogpl --mem $$QUEUE_PM{s22mem}G --num_paral $OPT{'cpu'} \\
          $OPT{outputdir}/00.shell/t022.align_bamsplit.sh \n";        
      print MAIN1 "echo ========== 2.align end   at : `date` ==========\n";
    }
    if($OPT{'task'} == 1 && $OPT{'analysis'} =~ /all|base|phase/){
      print MAIN1 "\necho ========== 3.phase start at : `date` ==========\n";
      print MAIN1 "perl $watchdogpl --mem $$QUEUE_PM{s31mem}G --num_paral $OPT{'cpu'} \\
          $OPT{outputdir}/00.shell/t031.phase_chrsplit.sh \n";
      print MAIN1 "echo ========== 3.phase end   at : `date` ==========\n";
    }
    if($OPT{'task'} == 1 && $OPT{'analysis'} =~ /all|base|filter|align|phase/){
      print MAIN1 "\necho ========== 3.stat start at : `date` ==========\n";
      print MAIN1 "perl $watchdogpl --mem $$QUEUE_PM{s321mem}G --num_paral $OPT{'cpu'} \\
          $OPT{outputdir}/00.shell/t0321.phase_phase.sh \n";
      print MAIN1 "perl $watchdogpl --mem $$QUEUE_PM{s322mem}G --num_paral $OPT{'cpu'} \\
          $OPT{outputdir}/00.shell/t0322.phase_phasestat.sh \n";         
      print MAIN1 "echo ========== 3.stat end   at : `date` ==========\n";
    }
    if($OPT{'task'} == 1 && $OPT{'analysis'} =~ /all|cnvsv/){
      print MAIN1 "\necho ========== 4.cnvsv start at : `date` ==========\n";
      print MAIN1 "perl $watchdogpl --mem $$QUEUE_PM{s41mem}G --num_paral $OPT{'cpu'} \\
          $OPT{outputdir}/00.shell/t041.cnvsv_call.sh \n";
      print MAIN1 "echo ========== 4.cnvsv end   at : `date` ==========\n";    
    }
    if($OPT{'task'} == 1 && $OPT{'analysis'} =~ /all|base|report/){
      print MAIN1 "\necho ========== 5.html start at : `date` ==========\n";
      print MAIN1 "perl $watchdogpl --mem $$QUEUE_PM{s51mem}G --num_paral $OPT{'cpu'} \\
          $OPT{outputdir}/00.shell/t051.report_html.sh \n";
      print MAIN1 "echo ========== 5.html end   at : `date` ==========\n";  
    }
  }

  #=============================================================================#
  # close file handles
  #=============================================================================#
  close MAIN1;
  close MAIN2;
  close SUB11;
  close SUB12;
  close SUB13;
  close SUB141;
  close SUB142;
  close SUB21;
  close SUB22;
  close SUB231;
  close SUB232;
  close SUB31;
  close SUB321;
  close SUB322;
  close SUB41;
  close SUB51;

  #=============================================================================#
  # Description:
  #    fix qsubsge txt such as:
  #      gather some task for parallel
  #=============================================================================#
  my ($fqprefix, $alignprefix, $phaseprefix, $cnvsvprefix) = ("") x 4;
  $fqprefix = 1 if $OPT{'analysis'} =~ /all|base|filter/;
  $alignprefix = ($OPT{'software'} =~ /bwa/i) ? $OPT{'queue'} : "fpga.q"
    if $OPT{'analysis'} =~ /all|base|align/;
  $phaseprefix = 1 if $OPT{'analysis'} =~ /all|base|phase/;
  $cnvsvprefix = 1 if $OPT{'analysis'} =~ /all|cnvsv/;  
  stLFR::CleanShell::report_cs(
    $fqprefix,
    $alignprefix,
    $phaseprefix,
    $cnvsvprefix,
    "$OPT{outputdir}/00.shell",
    $OPT{'task'},
    $OPT{'Project'},
    $OPT{'queue'},
    $OPT{'name'},
  );

};

1;
