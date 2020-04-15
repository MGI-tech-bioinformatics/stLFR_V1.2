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
    2. fix one BUG occur in SV detection.

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

With Docker
----------------
    Setup
    1. Install docker follow the official website
        https://www.docker.com/
    2. Then do the following for the workflow:
        docker pull rjunhua/stlfr_reseq_v1.2:v1
    3. Download and unzip the database
        From BGI Cloud Drive:
          https://pan.genomics.cn/ucdisk/s/Jvmuii
        Or from OneDrive:
          https://dwz.cn/ZPlGA0eJ
    Notes:
        1. Please make sure that you run the docker container with at least 45GB memory and 30 CPU.
        2. The input is sample list and output directory which descripted below (Main progarm arguments).    
    
    Running
    1. Please set the following variables on your machine:
      (a) $DB_LOCAL: directory on your local machine that has the database files.
      (b) $DATA_LOCAL: directory on your local machine that has the sequence data and "samplelist" file.
          "samplelist" must follow the format descripted bellow,
          and the *PATH* in "samplelist" must be absolute dicrtory of $DATA_LOCAL.
      (c) $RESULT_LOCAL: directory for result.
    2. Run the command:
        docker run -d -P \
        --name $STLFRNAME \
        -v $DB_LOCAL:/stLFR/db \
        -v $DATA_LOCAL:$DATA_LOCAL \
        -v $RESULT_LOCAL:$RESULT_LOCAL \
        rjunhua/stlfr_reseq_v1.2:v1 \
        /bin/bash \
        /stLFR/bin/stLFR_SGE \
        $DATA_LOCAL/samplelist \
        $RESULT_LOCAL
    3. After report is generated:
        docker rm $STLFRNAME

Without Docker but run in local server
----------------

    Download
    ----------------
    1.  download main scripts:
        git clone git@github.com:MGI-tech-bioinformatics/stLFR_V1.2.git
    2.  download database and unzip into 'db':
        For China mainland users, please using BGI Cloud Drive link:
             https://pan.genomics.cn/ucdisk/s/3mmUzy
        For other region users, please using OneDrive link:
             https://dwz.cn/UKb9SRWU
    3.  download tools and unzip into 'tools':
        For China mainland users, please using BGI Cloud Drive link:
             https://pan.genomics.cn/ucdisk/s/ErYf6v
        For other region users, please using OneDrive link:
             https://dwz.cn/dVkXg7No
    4.  preinstall tools
        More than 20 softwares/tools are used in this pipeline, and some of them are difficult to build static binary.
        The tools directory descriped below is just an example. So, please make sure these softwares/tools are installed firstly.
        (a). software/tool list
            bam2depth, bwa, circos, convert(ImageMagick), gnuplot,
            gatk4, java, picard, python2, python3, R, rtg-tools, samtools
        (b). some packages are required for software/tool, such as:
            R:          ggplot2,
                        karyoploteR(https://bioconductor.org/packages/release/bioc/html/karyoploteR.html)
            Python2:    vcf, pysam, numpy
            Python3:    pysam

    Running
    ----------------
    1. Make sure 'SAMPLELIST' file is in a right format.
    2. Run script with default parameters:
         perl bin/stLFR_SGE -l SAMPLELIST -outputdir OUTPUTDIR

Demo data
----------------
    two demo stLFR libraries are provided for testing, and every library consists two lanes:
      1. T0001-2:
            ftp://ftp.cngb.org/pub/CNSA/CNP0000387/CNS0057111/
      2. T0001-4:
            ftp://ftp.cngb.org/pub/CNSA/CNP0000387/CNS0094773/

Input: Sample List
----------------
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
