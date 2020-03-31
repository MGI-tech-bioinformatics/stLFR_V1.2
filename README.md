# stLFR_V1.2

Introduction
----------------
Tools of stLFR(Single Tube Long Fragment Reads) data analysis.

stLFR FAQs is directed to bgi-MGITech_Bioinfor@genomics.cn.

Download source code package from https://github.com/MGI-tech-bioinformatics/stLFR_V1.2.

Updates 
----------------
Mar 31, 2020
1. add a Docker version of stLFR V1.2.

Jan 10, 2020
1. Updated SV module (SV2.1, https://github.com/MGI-tech-bioinformatics/stLFR_SV2.1_module)
2. Novel summary report of HMTL format.
3. Redesigned pipeline structure and changed some parameters.

May 6, 2019
There are several updates in stLFR_v1.1 comparing with v1:
1. Users could use an alternative reference type (hg19 or hs37d5) in stLFR_v1.1 by --ref option instead of only hg19.
2. Updated CNV and SV detection tools are implied in stLFR_v1.1 for decreasing false discovery rate.
3. Three figures used for illustrating stLFR fragment distribution and coverage are added.
4. NA12878 benchmark VCF by GIAB is used for haplotype phasing error calculation.

Run in Docker
----------------
  1. Install docker follow the official website: https://www.docker.com/
  2. Then do the following for the workflow:
      docker pull rjunhua/stlfr_reseq_v1.2:v1
  3. Download and unzip the database from https://pan.genomics.cn/ucdisk/s/jQJFVn,
     and MD5 from https://pan.genomics.cn/ucdisk/s/rAva6z
  4. Run the command:
      docker run -d -P --name STLFRNAME \
      -v /USER/DB:/stLFR/db -v /USER/DATA:/USER/DATA -v /USER/RESULT:/USER/RESULT \
      rjunhua/stlfr_reseq_v1.2:v1 /bin/bash /stLFR/bin/stLFR_SGE \
      /USER/DATA/SAMPLELIST /USER/RSEULT
  5. Close docker
      docker rm STLFRNAME
  
  Notes:
  1. Please make sure that you run the docker container with at least 4GB memory and 15 CPU.
  2. The input is sample list and output directory which descripted below (Main progarm arguments).

Run in local server
----------------

Preinstallation
----------------
More than 20 softwares/tools are used in this pipeline, and some of them are difficult to build static binary.
The tools directory descriped below is just an example. So, please make sure these softwares/tools are installed firstly.

      1. software/tool list
         the directory/path is required for pipeline, * means the software/tool need pre-install
         tools
         |--- bam2depth
              |___ bam2depth
         |--- bwa
              |___ bwa
         |--- circos            * pre-install, some perl-packages are required
              |___ bin
                   |___ circos
         |--- cnv
         |--- fqcheck
              |--- convert      * pre-install
              |--- fqcheck33
              |--- fqcheck_distribute.pl
              |--- gnuplot      * pre-install
              |___ PLOT.pm
         |---gatk4              * pre-install
             |___ gatk
         |---HapCUT2-master
             |___ utilitie
                  |--- LinkFragments.py
                  |___ calculate_haplotype_statistics.py
         |--- jre
              |___ bin
                   |___ java
         |--- monitor
              |___ watchDog.pl
         |--- picard             * pre-install
              |___ picard.jar
         |--- Python2            * pre-install
              |___ python
         |--- python3            * pre-install
              |___ python3
         |--- R                  * pre-install 
              |___ bin
                   |___ R
         |--- rtg-tools
              |___ rtg
         |--- samtools
              |___ bin
                   |___ samtools
         |--- SOAPnuke
              |___ SOAPnuke
         |--- sv
         |___ vcftools
              |--- bcftools
              |--- bgzip
              |___ tabix
              
      2. some packages are required for software/tool, such as:
            R:          ggplot2, scales, regioneR, karyoploteR(https://bioconductor.org/packages/release/bioc/html/karyoploteR.html)
            Python2:    vcf, pysam, numpy
            Python3:    pysam

Download/Install
----------------
Due to the size limitation of GitHub repository, the database directory ('stLFR_V1.2/db') and tools directory ('stLFR_V1.2/tools') are provided below:

      1. tools:
         
         For China mainland users, please using BGI Cloud Drive link:
         https://pan.genomics.cn/ucdisk/s/ryUvuq
         
         For other region users, please using OneDrive link:
         https://dwz.cn/DWEiVTD6
         
      2. database:
      
         For China mainland users, please using BGI Cloud Drive link:
         https://pan.genomics.cn/ucdisk/s/3mmUzy

         For other region users, please using OneDrive link:
         https://dwz.cn/UKb9SRWU

Tool list in directory ('stLFR_V1.2/tools'):

      1. download tools from BGI Cloud or OneDrive and check up MD5 values in the related MD5.txt files.
      2. install some softwares/tools by yourself.

Database:

      1. download database from BGI Cloud or OneDrive.
      2. check and prepare database follow the 'readme.txt' in BGI Cloud.
      3. place database follow the 'db.tree.list' in GitHub.

Meanwhile, two demo stLFR libraries are provided for testing, and every library consists two lanes:

      1. T0001-2:
            ftp://ftp.cngb.org/pub/CNSA/CNP0000387/CNS0057111/
      2. T0001-4:
            ftp://ftp.cngb.org/pub/CNSA/CNP0000387/CNS0094773/

Usage
----------------
1. Make sure 'SAMPLELIST' file is in a right format.
2. Run script with default parameters:

         perl bin/stLFR_SGE -l SAMPLELIST -outputdir OUTPUTDIR
         
3. Or analyze one data in more than one steps:

         # only run fastq filter module first
         perl bin/stLFR_SGE -l SAMPLELIST -outputdir ANALYSIS_FQ -analysis filter
         # do other modules later
         perl bin/stLFR_SGE -l SAMPLELIST -outputdir ANALYSIS_OTHER -analysis align,phase,cnvsv,report -inputdir ANALYSIS_FQ
         

Main progarm arguments:
----------------

   Sample List
   
       -list FILE
                     Name of input file. This is required.

                     Five columns in the list file, which the front three columns are required:
                     1. name    : unique sample ID in this analysis
                     2. path    : fastq path(s) for this sample split with colon(:)
                     3. barcode : sample-barcode for each path split with colon(:), 0 means all used
                     4. reffile : reference with index, two inner options are 'hg19' and 'hs37d5', NULL or '-' means 'hs37d5'
                     5. vcffile : dbsnp file, default is NULL or '-'
                     6. blacklist : black list file(BED format) for SV
                     7. controllist : sorted control list file(BEDPE format) for SV

                     eg:  SAM1   /DATA/slide1/L01                    1-4,5,7-9
                          SAM2   /DATA/slide1/L01:/DATA/slide2/L01   0:1-8         hg19
                          SAM3   /DATA/slide2/L02                    0             REFERENCE/ref.fa      DBSNP/dbsnp.vcf
                          SAM4   /DATA/slide2/L03                    0             hs37d5                -                   BLACKLIST    CONTROLLIST

   Output Directory
   
       -outputdir [ ./ ]
                     Output directory path.

                     The Format of Output directory (also input directory of this workflow):
                     Inputdir/
                     |-- 01.filter   // for align
                     |   |__ SAMPLE
                     |       |__ SAMPLE.clean_1.fq.gz
                     |       |__ SAMPLE.clean_2.fq.gz
                     |-- 02.align    // for phase, cnv and sv
                     |   |__ SAMPLE
                     |       |__ SAMPLE.sortdup.bqsr.bam
                     |       |__ SAMPLE.sortdup.bqsr.bam.bai
                     |       |__ SAMPLE.sortdup.bqsr.bam.HaplotypeCaller.vcf.gz
                     |       |__ SAMPLE.sortdup.bqsr.bam.HaplotypeCaller.vcf.gz.tbi
                     |-- 03.phase    // for cnv
                     |   |__ SAMPLE
                     |       |__ phasesplit
                     |           |__ hapblock_SAMPLE_CHROMOSOME
                     |-- 04.cnv
                     |   |__ SAMPLE
                     |       |__ SAMPLE.CNV.result.xls
                     |-- 05.sv
                     |   |__ SAMPLE
                     |       |__ SAMPLE.SV.simple.result.xls
                     |       |__ SAMPLE.SV.complex.result.xls
                     |__ file
                         |__ SAMPLE
                             |-- alignment
                             |   |__ SAMPLE.sortdup.bqsr.bam
                             |   |__ SAMPLE.sortdup.bqsr.bam.bai
                             |-- CNV
                             |   |__ SAMPLE.CNV.result.xls
                             |-- haplotype
                             |   |__ SAMPLE.hapblock
                             |   |__ SAMPLE.hapcut_stat.txt
                             |   |__ SAMPLE.inked_fragment
                             |-- sequence
                             |   |__ SAMPLE.clean_1.fq.gz
                             |   |__ SAMPLE.clean_2.fq.gz
                             |-- SV
                             |   |__ SAMPLE.SV.result.xls
                             |__ variant
                                 |__ SAMPLE.sortdup.bqsr.bam.HaplotypeCaller.vcf.gz
                                 |__ SAMPLE.sortdup.bqsr.bam.HaplotypeCaller.vcf.gz.tbi

   Analysis Modules
   
       -analysis  [ all ]
                     There are 5 modules in this program: filter -> align -> phase -> cnvsv -> report.
                       -analysis all  == -analysis filter,align,phase,cnvsv,report.
                       -analysis base == -analysis filter,align,phase,report.

                     eg: all                : filter + align + phase + cnvsv + report
                         base               : filter + align + phase + report
                         align              : align + report
                         filter,align,phase : filter + align + phase + report

       -inputdir
                     Input directory path with results in previous process. [ ]
                     1. align:  need filter result
                     2. phase:  need align result
                     3. cnvsv:  need align and phase result

   stLFR barcode position
   
       -position [ 101_10,117_10,133_10 ]
                     Position of stLFR barcodes on read2.

   Task Monitor Type
   
       -cpu      [ 60 ]
                     CPU number on server

   Baseline for SNP/INDEL evaluation
   
       -baseline   Baseline VCF for NA12878 on hg19 or hs37d5.
       -confbed    High confidence BED for NA12878 on hg19 or hs37d5.

   Usage
   
       -help       Show brief usage synopsis.
       -man        Show man page.
       -run        Run workflow after main shells built.

Result
----------------
After all analysis processes ending, you will get these files below:

      1.  HTML report:                              *_cn.html, *_en.html
      2.  raw data summary:                         *.fastqtable.xls
      3.  stLFR barcode summary:                    *.fragtable.xls
      4.  alignment summary:                        *.aligntable.xls
      5.  variant summary:                          *.varianttable.xls
      6.  haplotype phasing summary:                *.haplotype.xls, *.haplotype.pdf
      7.  evaluation summary (NA12878):             *.evaluation.xls
      8.  quality distrubution in cleanfq:          *.Cleanfq.qual.png
      9.  base distribution in cleanfq:             *.Cleanfq.base.png
      10. depth distribution in alignment:          *.Sequencing.depth.pdf
      11. accumulated depth distribution:           *.Sequencing.depth.accumulation.pdf
      12. insert size distrubition:                 *.Insertsize.metrics.txt, *.Insertsize.pdf
      13. GC bias distrubution:                     *.GCbias.metrics.txt, *.GCbias.pdf
      14. Fragment coverage figure:                 *.frag_cov.pdf
      15. Fragment length distribution figure:      *.fraglen_distribution_min5000.pdf
      16. Fragment per barcode distribution figure: *.frag_per_barcode.pdf
      17. variant CIRCOS:                           *.circos.svg, *.circos.png, *.legend_circos.pdf

Additional Information
----------------
1. If user has "Permission denied" problem in the process of running，you can use the command "chmod +x -R stLFR_v2.1/tools" to get executable permission of tools.


License
----------------
Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions： 
  
The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
  
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
