#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import base64
import re
import os
import sys
import platform
import logging
import hashlib
import html_util
from io import open

logging.basicConfig(format='%(asctime)s - %(name)s - %(levelname)s - %(message)s', level=logging.INFO)

def upload(host, username, password, cwd, file):
    try:
        name = file.split('/')[-1]
        logging.info('Upload the file: ' + name)
        with FTP(host, timeout=30) as ftp:
            ftp.login(username, password)
            ftp.cwd(cwd)
            with open(file, 'rb') as f:
                ftp.storbinary('STOR ' + name[:-4] + hashlib.md5(open(file, 'rb').read()).hexdigest(), f)
            ftp.dir()
    except:
        logging.error('Error in using FTP to upload report.')

def getTableContent(filePath, isWarpTitle = True, lang = 'cn'):
    try:
        content = ""
        with open(filePath, 'r', encoding='utf-8') as file:
            index = 1
            for line in file.readlines():
                if index == 1:
                    content += "<tr>"
                    for field in line.split("\t"):
                        field = getContentLang(field, lang)
                        if isWarpTitle == True:
                            content += "<th scope='col'>" + setWarpTitle(field) + "</th>"           
                        else:
                            content += "<th scope='col'>" + field + "</th>"
                    content += "</tr>"
                if index != 1 :
                    content += "<tr>"
                    tdindex = 0
                    for field in line.split("\t"):
                        content += "<td scope='row'>"
                        # remvoe the str of line break. 
                        field = field.replace('\r', '').replace('\n', '')
                        field = getContentLang(field, lang)
                        #check the str is int and format it.
                        if field.isdigit():
                            content += format(int(field),',')
                        else:
                            content += field
                        content += "</td>"
                        tdindex += 1      
                    content += "</tr>"
                index = index + 1    
    except:
        content = '<div class="noDataTitle">无数据文件:' + filePath + '</div>'  
    return content

def getSpecialTableContent(path, directoryName, fimeName, isWarpTitle = 'false', lang = 'cn', italicIndx = 1):
    try:
        content = ""
        # print(path + directoryName + fimeName)
        with open(path + directoryName + fimeName, 'r', encoding='utf-8') as file:
            index = 1
            for line in file.readlines():
                if index == 1:
                    thIndx = 0
                    content += "<tr>"
                    style="";
                    for field in line.split("\t"):
                        field = getFiledNameByLang(field, lang)
                        if thIndx == 0:
                            style = "style='min-width: 45px;max-width: 45px;'"
                        elif thIndx == 1:
                            style = "style='min-width: 200px;max-width: 200px;'"
                        else:
                            style = "style='min-width: 100px;max-width: 100px;'"  
                        if isWarpTitle == 'false':
                            content += "<th scope='col'" + style + ">" + field + "</th>"
                        else:
                            content += "<th scope='col'" + style + ">" + setWarpTitle(field) + "</th>"
                        thIndx += 1         
                    content += "</tr>"
                if index != 1:
                    content += "<tr>"
                    # content += "<td>" + str(index - 1) + "</td>"
                    colIndex = 0
                    for field in line.split("\t"):
                        if colIndex == italicIndx: 
                            content += "<td scope='row' style='font-style:italic;'>" + field + "</td>"
                        else:
                            field = field.replace('\r', '').replace('\n', '')
                            content += "<td scope='row'>" 
                            if field.isdigit():
                                content += format(int(field),',')
                            else:
                                content += field
                            content += "</td>"
                        colIndex += 1       
                    content += "</tr>"
                index = index + 1
                if index == 12:
                    break
    except:
        content = '<div class="noDataTitle">无数据文件:' + path + directoryName + fimeName + '</div>'  
    return content 

def setWarpTitle(title):
    # print("filed:" +title)
    if "_" in title:
        titleName = ""
        for name in title.split("_"):
            titleName += name + "<br>"
        return titleName
    if "/" in title:
        titleName = ""
        index = 1
        for name in title.split("/"):
            if len(title.split("/")) == index:
                titleName += "/" + name + "<br>"
            else:
                titleName += name + "<br>"
            index += 1
        return titleName
    return title

def isInList(name, files):
    #print("name is " + name)
    #print("files is ")
    #print(files)
    if name in files:
        return {"hasFile": True, "name": name}
    return {"hasFile": False, "name": name}

def replaceTxtTOXls(name):
    return name.replace(".txt",".xls")

def getContentLang(contentName, lang):
    contentDict = {
        "SAMPLE"       : ["样品", "Sample"],
        "DOWNLOAD"     : ["下载", "Download"],
        "stLFR"        : ["stLFR重测序分析报告",  "stLFR resequencing analysis report"],
        "FASTQTITLE"   : ["数据产出",       "Fastq product"],
        "FASTQTDES"    : ["对高分子量gDNA样品使用华大智造的MGIEasy stLFR文库构建试剂盒进行文库制备，随后使用华大智造DNBSEQ测序仪对stLFR文库测序最终产出stLFR数据。stLFR数据将使用SOAPnuke进行低质量序列过滤得到有效数据，对有效数据进行碱基分布、质量分布、Q20、Q30等指标析以评估数据质量。", 
                          "Generation of stLFR resequencing data is started from the HMW gDNA sample. A library is constructed by MGIEasy stLFR Library Prep Kit by using HMW gDNA, then the library is sequenced on DNBSEQ platforms which  manufactured by MGI. The stLFR fastq files are got at the end. After removing low-quality reads by SOAPnuke, we got clean data of stLFR sample. Base distribution, quality distribution, Q20 and Q30 of clean stLFR data are showed below."],
        "FRAGTITLE"    : ["stLFR长片段产出",     "stLFR long-fragment product"],
        "FRAGTDES"     : ["使用虚拟隔离共标记技术实现的长片段序列是stLFR产品的特色。通过对长片段与barcode的组合关系、长片段长度及覆盖分析来展示stLFR产品的长片段性能。", 
                          "Long-fragment read is a feature of stLFR technology, built by virtual co-barcoding method. The performance of Long-fragment read is showed by analyzing relationship between long-fragment and barcode, long-fragment length distribution and long-fragment coverage."],
        "ALIGNTITLE"   : ["比对信息",       "Alignment Statistics"],
        "ALIGNTDES"    : ["stLFR有效序列通过BWA或MegaBOLT软件与参考基因组进行比对定位，并通过比对率、深度覆盖、插入片段等指标进行性能评估。", 
                          "Clean stLFR data is aligned to reference genome by BWA or MegaBOLT. The performance of alignment is evaluated by mapping rate, coverage, insert size distribution, etc."],
        "VARIANTTITLE" : ["变异信息",       "Variant Statistics"],
        "VARIANTTDES"  : ["通过GATK或MegaBLOT以及华大自主开发的CNV及SV软件对stLFR数据进行多种变异检测，得到SNP、INDEL、CNV、SV突变信息，并通过CIRCOS进行可视化展示。", 
                          "The variants of stLFR data are called by GATK or MegaBOLT and BGI self-developed CNV/SV calling softwares, and visualized by CIRCOS."],
        "PHASETITLE"   : ["单倍体组装信息", "Phasing Statistics"],
        "PHASETDES"    : ["使用HapCUT2软件结合stLFR数据的比对及突变信息进行单倍体组装，获取高质量的单倍体组装结果。", 
                          "High quality haplotype assembly results are achieved by HapCUT2 from stLFR data."],

        "BASEDIS"     : ["碱基分布", "Base distribution"],
        "BASEDISDES"  : ["碱基比例分布图（x轴表示序列上的位置，y轴显示序列具体位置上的五种碱基比例）。", 
                         "Base distribution (x-axis is the sequencing cycle, y-axis is the proportion of five bases at a certain cycle)."],
        "BASEQUAL"    : ["碱基质量分布", "Base quality distribution"],
        "BASEQUALDES" : ["碱基质量分布图（x轴表示序列上的位置，y轴显示序列具体位置上的碱基质量分布热图，从白到绿到红依次表示碱基比例从低到高）。", 
                         "Quality distribution (x-axis is the cycle along reads, y-axis is the heatmap of quality at given cycle)."],

        "FRAGBAR"     : ["文库片段barcode分布图", "Fragment barcode distribution"],
        "FRAGBARDES"  : ["stLFR文库片段与barcode关系示意图（x轴表示长片段中包含barcode的个数，y轴表示包含特定barcode个数的长片段数目）。", 
                         "Long-fragment vs. barcode in stLFR data (x-axis is average long-fragment DNA molecular number captured by barcode beads, y-axis is the density of barcode number)."],
        "FRAGLEN"     : ["文库片段长度分布图", "Fragment length distribution"],
        "FRAGLENDES"  : ["stLFR文库片段长度分布（x轴表示长度，y轴表示特定长度片段的比例）。", 
                         "Length distribution of long-fragment in stLFR data (x-axis is the length, y-axis is the count of long-fragment with given length)."],
        "FRAGCOV"     : ["文库片段覆盖分布图", "Fragment coverage distribution"],
        "FRAGCOVDES"  : ["stLFR文库片段覆盖度分布（x轴表示覆盖度，y轴表示特定覆盖度的片段数目）。", 
                         "Coverage of long-fragment of stLFR data (x-axis is the coverage of long-fragment, y-axis is the density of long-fragment at a certain coverage)."],

        "INSERT"      : ["插入片段", "Insert size distribution"],
        "INSERTDES"   : ["插入片段分布图（x轴表示插入片段大小，y轴表示特定插入片段的比例）。", 
                         "Insert size distribution (x-axis is insert size, y-axis is the count of paired reads at a certain insert size)."],
        "DEEP"        : ["深度累积分析", "Depth accumulative distribution"],
        "DEEPDES"     : ["累积深度分布图（x轴表示测序覆盖深度，y轴表示全基因组范围不低于特定覆盖深度的比例）。", 
                         "Cumulatived sequencing depth distribution (x-axis is the depth, y-axis is the proportion of reference genome that achieves at or above certain depth)."],
        "DEEP2"       : ["深度分布", "Depth distribution"],
        "DEEP2DES"    : ["深度分布图（x轴表示测序覆盖深度，y轴表示特定覆盖深度在全基因组范围中所占比例）。",
                         "Sequencing depth distribution (x-axis is the depth, y-axis is the proportion of reference genome at a certain depth)."],
        "GCBIAS"      : ["GC偏差", "GC bias distribution"],
        "GCBIASDES"   : ["GC-bias分布图（x轴表示GC比例，y轴蓝点表示特定GC比例的归一化覆盖深度，y轴绿线表示特定GC区域的平均质量值，y轴红线表示参考基因组中特定GC所占比例）。", 
                         "GC-bias distribution (x-aixs shows GC content, y-axis plots the NORMALIZED_COVERAGE, the distribution of WINDOWs corresponding to GC percentages, and base qualities corresponding to each GC content bin)."],

        "FASTQTAB"    : ["序列统计", "Fastq report"],
        "FASTQTABDES" : ["<p>&emsp;a. 原始序列经过低质量过滤后得到有效序列，低质量序列为包含以下任何一条：（1）含有接头序列；（2）序列中N比例超过1%；（3）序列中质量值低于10比例超过10%。</p><p>&emsp;b. Q20、Q30为有效序列中质量高于20或30的碱基比例。</p><p>&emsp;c. 测序深度根据原始序列数统计。</p>", 
                         "<p>&emsp;a. Clean fastq is built by removing low-quality read from raw fastq. One pair of read is defined as low-quality when it:(1) contains adapter sequences, (2) N ratio > 1%, (3) the ratio of the base (base quality is lower than 10) > 10%.</p><p>&emsp;b. Q20 and Q30 means the ratio of base whose quality > 20 or 30 of the clean fastq.</p><p>&emsp;c. Total depth is calculated by the raw fastq.</p>"],
        "FRAGTAB"     : ["文库片段Barcode统计", "Fragment & Barcode report"],
        "FRAGTABDES"  : ["<p>&emsp;a. stLFR数据是通过barcode进行长片段序列构建，因此需要首先对stLFR数据进行barcode拆分及过滤，并对每条序列进行barcode标记。表格显示了stLFR数据barcode拆分及过滤结果。</p>",
                         "<p>&emsp;a. The barcodes of stLFR data are related to the Long-fragment read. We should do and filter and split the barcodes before we do next analysis. After this step, each read will get a barcode ID. The table above shows the barcode split statistics results.</p>"],
        "ALIGNTAB"    : ["比对统计", "Alignment report"],
        "ALIGNTABDES" : ["<p>&emsp;a. 比对结果基于低质量过滤及barcode拆分后的有效数据进行统计。</p><p>&emsp;b. 比对率和成对比对率分别表示正确和成对正确的比对到参考基因组的序列比例。</p><p>&emsp;c. 错配表示序列中与参考基因组不同的碱基信息。</p><p>&emsp;d. 重复序列表示来自PCR或其他因素导致的重复序列。</p><p>&emsp;e. 有效平均深度表示去除重复序列后基因组平均覆盖深度。</p><p>&emsp;f. 覆盖度分别统计参考基因组中覆盖超过1层、4层、10层、20层的比例。</p><p>&emsp;g. 平均插入片段表示成对比对序列的平均片段长度。</p>", 
                         "<p>&emsp;a. The statistics results shown here are calculated based on clean data after low-quality filter and barcode-split filter.</p><p>&emsp;b. The mapping rate and paired mapping rate are the proportion of read which is mapped or paired-mapped to reference genome.</p><p>&emsp;c. The mismatch means the unmapped bases between read and reference genome.</p><p>&emsp;d. The duplication means duplicated reads introduced by PCR or other processes.</p><p>&emsp;e. The average sequencing depth is calculated by mapped reads without duplicated reads.</p><p>&emsp;f. The coverage is calculated as proportion of reference genome covered by more than 1, 4, 10 or 20 folds reads.</p><p>&emsp;g. The mean insert size is the mean insert size of paired mapped reads.</p>"],

        "Sample name"                    : ["样本", "Sample"],
        "Raw reads"                      : ["原始序列数", "Raw reads"],
        "Raw bases(bp)"                  : ["原始碱基数(bp)", "Raw bases(bp)"],
        "Clean reads"                    : ["有效序列数", "Clean reads"],
        "Clean bases(bp)"                : ["有效碱基数(bp)", "Clean bases(bp)"],
        "Q20"                            : ["Q20(%)", "Q20(%)"],
        "Q30"                            : ["Q30(%)", "Q30(%)"],
        "Total barcode type"             : ["理论Barcode种类", "Total barcode type"],
        "Barcode number"                 : ["实际Barcode种类", "Barcode number"],
        "Barcode type rate"              : ["实际Barcode比例", "Barcode type rate"],
        "Reads pair number"              : ["序列对数", "Reads pair number"],
        "Reads pair number(after split)" : ["拆分后有效序列对数", "Reads pair number(after split)"],
        "Barcode split rate"             : ["拆分率", "Barcode split rate"],
        "Mapped reads"                   : ["比对序列数", "Mapped reads"],
        "Mapped bases(bp)"               : ["比对碱基数(bp)", "Mapped bases(bp)"],
        "Mapping rate"                   : ["比对率", "Mapping rate"],
        "Paired mapping rate"            : ["成对比对率", "Paired mapping rate"],
        "Mismatch bases(bp)"             : ["错配碱基数", "Mismatch bases(bp)"],
        "Mismatch rate"                  : ["错配率", "Mismatch rate"],
        "Duplicate reads"                : ["重复序列数", "Duplicate reads"],
        "Duplicate rate"                 : ["重复率", "Duplicate rate"],
        "Total depth"                    : ["测序深度(X)", "Total depth(X)"],
        "Split barcode(G)"               : ["拆分后数据量(G)", "Split barcode(G)"],
        "Dup depth"                      : ["平均测序深度(X)", "Dup depth(X)"],
        "Average sequencing depth"       : ["有效平均深度(X)", "Average sequencing depth(X)"],
        "Coverage"                       : ["覆盖度(≥1X)",  "Coverage(≥1X)"],
        "Coverage at least 4X"           : ["覆盖度(≥4X)",  "Coverage(≥4X)"],
        "Coverage at least 10X"          : ["覆盖度(≥10X)", "Coverage(≥10X)"],
        "Coverage at least 20X"          : ["覆盖度(≥20X)", "Coverage(≥20X)"],
        "Mean insert size"               : ["平均插入片段(bp)", "Mean insert size(bp)"],

        "CIRCOS"           : ["变异结果CIRCOS示意图", "Circos of variants"],
        "CIRCOSDES"        : ["全基因组变异示意图及图例（图形由6个圆环组成，从外到内依次是(i)染色体、(ii)SNP密度曲线、(iii)INDEL密度曲线、(iv)CNV中deletion分布图、(v)CNV中duplication分布图、(vi)SV分布图）。", 
                              "CIRCOS of variant from whole genome (From the outside to the inside are: (i) chromosomes, (ii) SNP densitity, (iii) INDEL densitity, (iv) CNV deletions, (v) CNV duplications and (vi) structure variants)."],
        "VARIANTTAB"       : ["变异统计", "Variant statistics"],
        "VARIANTTABDES"    : ["<p>&emsp;a. SNP： 单核苷酸多态性。</p><p>&emsp;b. Ti/TV：SNP的转换颠换比例。</p><p>&emsp;c. INDEL：插入与缺失。</p><p>&emsp;d. CNV：拷贝数变异。</p><p>&emsp;e. SV：结构变异。</p><p>&emsp;f. DEL： 缺失突变。</p><p>&emsp;g. DUP：重复突变。</p><p>&emsp;h. INV： 颠倒突变。</p><p>&emsp;i. TRA：移位突变。</p>", 
                              "<p>&emsp;a. SNP: Insertion& Deletion.</p><p>&emsp;b. Ti/TV：the ratio of transition (Ti) to transversion (Tv) SNPs.</p><p>&emsp;c. INDEL: Insertion & Deletion.</p><p>&emsp;d. CNV：copy number variant.</p><p>&emsp;e. SV: Structure Variant.</p><p>&emsp;f. DEL: deletion.</p><p>&emsp;g. DUP: duplication.</p><p>&emsp;h. INV: inversion.</p><p>&emsp;i. TRA: translocation.</p>"],
        "VCFEVALTAB"       : ["变异评估", "Variant evaluation"],
        "VCFEVALTABDES"    : ["<p>&emsp;a. PPV： 阳性预测值，准确度，TP/(TP+FP)。</p><p>&emsp;b. Sensitivity： 灵敏度, TP/(TP+FN).</p><p>&emsp;c. F-measure： 灵敏度与准确度的调和平均数, 2*TP/(2*TP+FP+FN).</p>",
                              "<p>&emsp;a. PPV: positive predictive value, precision, TP/(TP+FP).</p><p>&emsp;b. Sensitivity: true positive rate, TP/(TP+FN).</p><p>&emsp;c. F-measure: the harmonic mean of precision and sensitivity, 2*TP/(2*TP+FP+FN).</p>"],
        "Sample"           : ["样本", "Sample"],
        "Total_SNP"        : ["SNP个数", "Total_SNP"],
        "dbSNP_rate"       : ["dbSNP比例", "dbSNP_rate"],
        "Novel_SNP"        : ["未知SNP", "Novel_SNP"],
        "Novel_SNP_Rate"   : ["未知SNP比例", "Novel_SNP_Rate"],
        "Ti/Tv"            : ["Ti/Tv", "Ti/Tv"],
        "Total_INDEL"      : ["INDEL个数", "Total_INDEL"],
        "dbINDEL_Rate"     : ["dbINDEL比例", "dbINDEL_Rate"],
        "CNV deletion"     : ["CNV deletion", "CNV deletion"],
        "CNV duplication"  : ["CNV duplication", "CNV duplication"],
        "SV DEL"           : ["SV DEL", "SV DEL"],
        "SV DUP"           : ["SV DUP", "SV DUP"],
        "SV INV"           : ["SV INV", "SV INV"],
        "SV TRA1"          : ["SV TRA1", "SV TRA1"],
        "SV TRA2"          : ["SV TRA2", "SV TRA2"],

        "SNP"              : ["SNP", "SNP"],
        "INDEL"            : ["INDEL", "INDEL"],
        "Threshold"        : ["&emsp;", "&emsp;"],
        "True-pos-call"    : ["TP", "TP"],
        "False-pos"        : ["FP", "FP"],
        "False-neg"        : ["FN", "FN"],
        "Precision"        : ["PPV", "PPV"],
        "Sensitivity"      : ["Sensitivity", "Sensitivity"],
        "F-measure"        : ["F-score", "F-score"],

        "PHASE"              : ["单倍体组装示意图", "Phasing distribution"],
        "PHASEDES"           : ["单倍体组装示意图（每一行表示一个染色体，白色表示没有组装区域，灰色、深蓝色表示独立组装好的单倍体block）。", 
                                "The haplotype phasing plot (Each line shows one chrmosome. White shows no phasing block, where grey and dark-bule show phasing blocks)."],
        "PHASETAB"           : ["单倍体组装统计", "Phasing statistics"],
        "PHASETABDES"        : ["<p>&emsp;a. 转换率: 发生转换错误位点的比例。</p><p>&emsp;b. 错配率：发生不一致错配位点的比例。</p><p>&emsp;c. 一致率： 单倍体与参考序列的最小汉明距离。</p><p>&emsp;d. 缺失率：所有覆盖位点上发生错配的比例。</p><p>&emsp;e. SNP组装数： 用于组装成单倍体的SNV数据。</p><p>&emsp;f. AN50： 单倍体结果中AN50。</p><p>&emsp;g. N50： 单倍体结果中N50。</p><p>&emsp;h. block组装率： 最大组装单倍体block中SNP比例。</p><p>&emsp;i. SNP组装率： 用于单倍体组装的SNP比例。</p>", 
                                "<p>&emsp;a. switch rate: the fraction of switch errors.</p><p>&emsp;b. mismatch rate: the fraction of mismatch errors.</p><p>&emsp;c. flat rate: the fraction of flat errors.</p><p>&emsp;d. missing rate: the fraction missing errors.</p><p>&emsp;e. phased count: counts of total SNVs phased in the test haplotype.</p><p>&emsp;f. AN50: the phasing block AN50 length of haplotype completeness.</p><p>&emsp;g. N50: the phasing block N50 length of haplotype completeness.</p><p>&emsp;h. max block snp frac: the fraction of SNVs in the largest (most variants phased) block.</p><p>&emsp;i. phasing rate: the fraction of SNVs in all blocks.</p>"],
        "chr"                : ["染色体", "chr"],
        "chrAll"             : ["全基因组", "Genome"],
        "switch rate"        : ["转换率", "switch rate"],
        "mismatch rate"      : ["错配率", "mismatch rate"],
        "flat rate"          : ["一致率", "flat rate"],
        "missing rate"       : ["缺失率", "missing rate"],
        "phased count"       : ["SNP组装数", "phased count"],
        "AN50"               : ["AN50", "AN50"],
        "N50"                : ["N50", "N50"],
        "max block snp frac" : ["block组装率", "max block snp frac"],
        "phasing rate"       : ["SNP组装率", "phasing rate"],
    }

    if contentName in contentDict.keys():
        index = 0
        if lang == 'en':
            index = 1
        return contentDict[contentName][index]
    return contentName

def usePlatform( ):
    sysstr = platform.system()
    if(sysstr =="Windows"):
        return '\\\\'    
    else:
        return '\\/'

def create_hide_box(title,  box_name, content, descript):
    html = '''
            <div class="secOne">
               <h1 class="headText headTextOne">
                 %(title)s
                 <div class="headTextIconDIV">
                   <img id="iDiv%(box_name)s" src="%(up)s" onclick="showAndHidden_%(box_name)s();">
                 </div>
               </h1>
               <div class = "content">
                 <p>
                   %(descript)s
                 </p>
               </div>
               <div id="%(box_name)s">
                   %(content)s                            
               </div>             
               <script type="text/javascript">
                 var %(box_name)s = document.getElementById('%(box_name)s');
                 %(box_name)s.style.display = 'block';
                 function showAndHidden_%(box_name)s() {
                   if (%(box_name)s.style.display == 'block') {
                      %(box_name)s.style.display = 'none';
                      document.getElementById('iDiv%(box_name)s').src = "%(down)s";
                   } else {
                      %(box_name)s.style.display = 'block';
                      document.getElementById('iDiv%(box_name)s').src = "%(up)s";
                   }
                 }
		           </script>
            </div>
    ''' 
    html = html % {
        "title"    : title,
        "box_name" : box_name,
        "content"  : content,
        "descript" : descript,
        "down"     : html_util.HtmlUtil.getPNGBinary('arrow-down.png'),
        "up"       : html_util.HtmlUtil.getPNGBinary('arrow-up.png') 
    }

    return html

def get_png_div(path):
    if os.path.exists(path):
        return '''
                <div class="secThree">
                    <img src="%s">
                </div>
        ''' % html_util.HtmlUtil.getPNGBinary(path, False)
    else:
        return '<div class="noDataTitle">无图片文件:'+ path + '</div>'

def create_sec_box(filePath, title, context_name, lang, type, linkPath):
    html = '''
           <div class="secBox">
               <h1>%s</h1>
           </div>
           <div class = "dataLink">
                <a href = "%s" download>%s</a>
           </div>
    ''' % (getContentLang(title, lang), 
           linkPath, 
           getContentLang("DOWNLOAD", lang)
    )
    if type == 'table':
        html += '''
            <div class="secThree">
                <table cellspacing="0" style="width: 1050px;">
                    <tbody>%s</tbody>
                </table>
            </div>''' % getTableContent(filePath, False ,lang)
    else:
        html += get_png_div(filePath)
    if getContentLang(context_name, lang):
        html += '''
                <div class="content">
                     <p>%s</p>
                </div>
             ''' % getContentLang(context_name, lang)
    return html

def create_sample(name, lang, type):
    html = '''
        <div class = "secBox">
            <h1>%s</h1>
        </div>
        <div class = "content">
            <p>%s</p>
        </div>
    ''' % (getContentLang("SAMPLE", lang), name)
    return html

def generate_html_report(path, version, name='rem', output_path='E:/codezlims/rem/Result/sample1', lang="cn"):
    base_png        = path + "/" + name + "/" + name + ".Cleanfq.base.png"
    qual_png        = path + "/" + name + "/" + name + ".Cleanfq.qual.png"
    fragbarcode_png = path + "/" + name + "/" + name + ".frag_per_barcode.png"
    fraglen1_png    = path + "/" + name + "/" + name + ".fraglen_distribution_min5000-0.png"
    fragcov_png     = path + "/" + name + "/" + name + ".frag_cov.png"
    insert_png      = path + "/" + name + "/" + name + ".Insertsize.png"
    deep_png1       = path + "/" + name + "/" + name + ".Sequencing.depth.accumulation.png"
    deep_png2       = path + "/" + name + "/" + name + ".Sequencing.depth.png"
    gc_png          = path + "/" + name + "/" + name + ".GCbias.png"
    circos_png      = path + "/" + name + "/" + name + ".circos.png"
    circosleg_png   = path + "/" + name + "/" + name + ".legend_circos.png"
    phase_png       = path + "/" + name + "/" + name + ".haplotype.png"

    fragbarcode_pdf = path + "/" + name + "/" + name + ".frag_per_barcode.pdf"
    fraglen_pdf     = path + "/" + name + "/" + name + ".fraglen_distribution_min5000.pdf"
    fragcov_pdf     = path + "/" + name + "/" + name + ".frag_cov.pdf"
    insert_pdf      = path + "/" + name + "/" + name + ".Insertsize.pdf"
    deep_pdf1       = path + "/" + name + "/" + name + ".Sequencing.depth.accumulation.pdf"
    deep_pdf2       = path + "/" + name + "/" + name + ".Sequencing.depth.pdf"
    gc_pdf          = path + "/" + name + "/" + name + ".GCbias.pdf"
    circos_svg      = path + "/" + name + "/" + name + ".circos.svg"
    phase_pdf       = path + "/" + name + "/" + name + ".haplotype.pdf"
    
    fastq_table     = path + "/" + name + "/" + name + ".fastqtable.xls"
    frag_table      = path + "/" + name + "/" + name + ".fragtable.xls"
    align_table     = path + "/" + name + "/" + name + ".aligntable.xls"
    variant_table   = path + "/" + name + "/" + name + ".varianttable.xls"
    phase_table     = path + "/" + name + "/" + name + ".haplotype.xls"
    vcfeval_table   = path + "/" + name + "/" + name + ".evaluation.xls"

    fastqBoxContent = (
                        create_sample(name, lang, "png")
                      + create_sec_box(fastq_table, "FASTQTAB", "FASTQTABDES", lang, "table", fastq_table)
                      + create_sec_box(base_png,    "BASEDIS",  "BASEDISDES",  lang, "png", base_png) 
                      + create_sec_box(qual_png,    "BASEQUAL", "BASEQUALDES", lang, "png", qual_png)
                      )

    fragmentBoxContent = (
                           create_sec_box(frag_table,      "FRAGTAB", "FRAGTABDES", lang, "table", frag_table)
                         + create_sec_box(fragbarcode_png, "FRAGBAR", "FRAGBARDES", lang, "png", fragbarcode_pdf)
                         + create_sec_box(fraglen1_png,    "FRAGLEN", "FRAGLENDES", lang, "png", fraglen_pdf)
                         + create_sec_box(fragcov_png,     "FRAGCOV", "FRAGCOVDES", lang, "png", fragcov_pdf)
                         )

    alignmentBoxContent = (
                            create_sec_box(align_table, "ALIGNTAB", "ALIGNTABDES", lang, "table", align_table)
                          + create_sec_box(insert_png,  "INSERT",   "INSERTDES",   lang, "png", insert_pdf)
                          + create_sec_box(deep_png1,   "DEEP",     "DEEPDES",     lang, "png", deep_pdf1)
                          + create_sec_box(deep_png2,   "DEEP2",    "DEEP2DES",    lang, "png", deep_pdf2)
                          + create_sec_box(gc_png,      "GCBIAS",   "GCBIASDES",   lang, "png", gc_pdf)
                          )

    variantBoxContent = ( 
                          create_sec_box(variant_table, "VARIANTTAB", "VARIANTTABDES", lang, "table", variant_table)
                        + create_sec_box(vcfeval_table, "VCFEVALTAB", "VCFEVALTABDES", lang, "table", vcfeval_table)
                        + create_sec_box(circos_png, "CIRCOS", "CIRCOSDES", lang, "png", circos_svg)
                        + get_png_div(circosleg_png)
                        )

    phaseBoxContent = (
                        create_sec_box(phase_table, "PHASETAB", "PHASETABDES", lang, "table", phase_table)
                      + create_sec_box(phase_png, "PHASE", "PHASEDES", lang, "png", phase_pdf)
                      )

    fastqBox     = create_hide_box("1." + getContentLang("FASTQTITLE", lang),   "boxOne",   fastqBoxContent, getContentLang("FASTQTDES", lang))
    fragmentBox  = create_hide_box("2." + getContentLang("FRAGTITLE", lang),    "boxTow",   fragmentBoxContent, getContentLang("FRAGTDES", lang))
    alignmentBox = create_hide_box("3." + getContentLang("ALIGNTITLE", lang),   "boxThree", alignmentBoxContent, getContentLang("ALIGNTDES", lang))
    variantBox   = create_hide_box("4." + getContentLang("VARIANTTITLE", lang), "boxFour",  variantBoxContent, getContentLang("VARIANTTDES", lang))
    phaseBox     = create_hide_box("5." + getContentLang("PHASETITLE", lang),   "boxFive",  phaseBoxContent, getContentLang("PHASETDES", lang))

    html=   '''
            <!DOCTYPE html>
            <html>
            <head>
                <meta http-equiv="Content-Type" content="text/html; charset=UTF-8">
                <title>'''+ getContentLang(product_type, lang) + '''</title>
                <meta http-equiv="X-UA-Compatible" content="IE=edge,chrome=1">
                <meta name="viewport" content="width=device-width, initial-scale=1, minimum-scale=1, maximum-scale=1, user-scalable=no">
                <link href="https://fonts.googleapis.com/css?family=Roboto" rel="stylesheet">
                <!-- load css file -->
                <style type="text/css">
            ''' + html_util.HtmlUtil.getFileContent('base.css') + '''
            ''' + html_util.HtmlUtil.getFileContent('common.css') + '''
            ''' + html_util.HtmlUtil.getFileContent('table.css') + '''
                </style>
            </head>
            <body>
            <!--header starts-->
            <div class="header">
                <div class="repeater"></div>
                <div class="wrapper">
                    <div class="headerBox clearfloat">
                        <div class="headLeft fl">
                        <h1>'''+ getContentLang(product_type, lang) + '''</h1>
                        <div style="float:right; margin:-40px -130px 10px 160px; background-color: #f7faeb; width: 100px; height: 30px; border-radius: 8px;text-align:center;"><a href="./'''+ name 
                        
    if lang == "cn":
        html+= "_en"
    else:
        html+= "_cn"
    html+=               '''.html" style="color: #1c567f; font-size: 20px; padding-top: 2px;">'''
    if lang == "cn":
        html+= "English"
    else:
        html+= "中文"
    html+=          '''</a></div>
                        <h2>%(version)s</h2>
                        </div>
                        <!-- headLeft -->
                        <div class="headRight fr">
                        <div class="logo">
                        <img src="%(logo)s">
                        </div>
                        </div>
                    <!-- headRight -->
                    </div>
                </div>
            </div>
            <!--header ends-->

            <!--container starts-->
            <div class="container">
            <div class="wrapper">

            %(fastqBox)s
            %(fragmentBox)s
            %(alignmentBox)s
            %(variantBox)s
            %(phaseBox)s

        <div class="repeater"></div>
        <div class="secBottom"></div>
        </div>
        </body>
        </html>
    ''' %{
        "version": version,
        "logo": html_util.HtmlUtil.getPNGBinary('logo.png'),
        "fastqBox": fastqBox,
        "fragmentBox" : fragmentBox,
        "alignmentBox" : alignmentBox,
        "variantBox" : variantBox,
        "phaseBox" : phaseBox
    }

    with open(output_path + '/' + name + '_' + lang +'.html', 'w+', encoding='utf-8') as report:
        print(output_path + '/' + name + '_' + lang +'.html') 
        report.write(html)


if __name__ == '__main__':

    if len(sys.argv) < 4:
        logging.error("Usage: python3 %s [version] [name] [result_path] [outdir|default:result_path]\n" %sys.argv[0])
        sys.exit(-1)
    version = sys.argv[1]
    name    = sys.argv[2]
    path    = sys.argv[3]
    outdir  = sys.argv[4]
    product_type = "stLFR"

    logging.info('\nResult folder is: %s\nSample name is: %s' % (path, name))
    # Local test
    # html_util.HtmlUtil.getPNGBinary('link.png')
    generate_html_report(path, version, name, outdir, "en")
    generate_html_report(path, version, name, outdir, "cn")
