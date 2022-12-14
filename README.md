# Length variation in short tandem repeats affects gene expression in natural populations of Arabidopsis thaliana


## 1、寻找STR

## 1.1 数据下载
```bash
mkdir sequence
cd sequence

prefetch SRR1945487 SRR1945478 SRR1945477 SRR1945464

parallel -j 4 "
    mv {1} ./
" ::: $(find ./ -maxdepth 2 -type f -name "*.sra")

parallel -j 4 "
    rm -r {1}
" ::: $(find ./ -maxdepth 1 -type d -name "SRR*")
```

+ 格式转换（sra ===> fastq.gz）
```bash
parallel -j 4 "
    fastq-dump --split-3 --gzip {1}
" ::: $(ls *.sra)

rm *.sra

parallel -j 4 "
    gzip -d {1}
" ::: $(ls *.gz)
```

## 1.2.1 利用 BWA 将 reads 匹配到 TAIR10 的参考基因组（Col-0）上
+ BWA (Burrows-Wheeler-Alignment Tool) 简介

BWA 是一种能够将差异度较小的序列比对到一个较大的参考基因组上的软件包。它由三个不同的算法：

```
BWA-backtrack: 是用来比对 Illumina 的序列的，reads 长度最长能到 100bp。

BWA-SW: 用于比对 long-read ，支持的长度为 70bp-1Mbp；同时支持剪接性比对。

BWA-MEM: 推荐使用的算法，支持较长的read长度，同时支持剪接性比对（split alignments)，但是BWA-MEM是更新的算法，也更快，更准确，且 BWA-MEM 对于 70bp-100bp 的 Illumina 数据来说，效果也更好些。
```

在运用这三种算法之前，需要先利用 BWA 的 index 命令，对参考基因组构建索引
```
bwa index [-p prefix] [-a algoType] <in.db.fasta>
-p 输出数据库的前缀
-a 构建index的算法
    -a is 是默认的算法，虽然相对较快，但需要较大的内存，并且不能构建大于 2GB 的数据库；
    -a bwtsw 对与短的参考序列是不工作的，必须要大于等于 10MB ，能用于较大的基因组数据。
```

三种算法对应于三种不同的命令
```
aln/samse/sampe ----> BWA-backtrack (samse 中的 se 是 single-end 的简写，而 sampe 中的 pe 是 paired-end 的简写）。
bwasw ----> BWA-SW
mem ----> BWA-MEM
```

+ TAIR10 参考基因组（Col-0）下载
```bash
mkdir sequence/reference
cd sequence/reference

wget http://ftp.ensemblgenomes.org/pub/plants/release-54/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
gzip -d Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
```

+ BWA 的使用
```bash
brew install bwa-mem2

mkdir BWA
cd BWA

bwa index ../sequence/reference/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa -p ref
# for i in SRR1945464 SRR1945477 SRR1945478 SRR1945487;do
#     bwa mem -t 4 ref ../sequence/$i'_1'.fastq  ../sequence/$i'_2'.fastq > $i.sam
# done


for i in SRR1945464 SRR1945477 SRR1945478 SRR1945487;do
    bwa mem -t 4 -M -R "@RG\tID:$i\tLB:$i\tPL:Illumina\tPU:$i\tSM:$i"  ref ../sequence/$i'_1'.fastq  ../sequence/$i'_2'.fastq > $i.sam 2> $i.log
done

parallel -k -j 2 "
    samtools sort -@ 2 {1}.sam > {1}.sort.bam
    samtools index {1}.sort.bam
" ::: $(ls *.sam | perl -p -e 's/\.sam$//')
```

## 1.2.2 利用hisat2进行序列比对
```bash
mkdir hisat2
cd hisat2

# 建立索引
hisat2-build  -p 6 ../sequence/reference/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa Arabidopsis_thaliana

# 比对
mkdir output
cd output

parallel -k -j 4 "
    hisat2 -t -x ../Arabidopsis_thaliana \
      -1 ../../sequence/{1}_1.fastq.gz -2 ../../sequence/{1}_2.fastq.gz -S ./{1}.sam \
      2>./{1}.log
" ::: $(ls ../../sequence/*.gz | perl -p -e 's/..\/..\/sequence\///;
s/_.+//' | uniq)

# 格式转化与排序
parallel -k -j 2 "
    samtools sort -@ 2 {1}.sam > {1}.sort.bam
    samtools index {1}.sort.bam
" ::: $(ls *.sam | perl -p -e 's/\.sam$//')
```

## 1.3 使用 Tandem Repeats Finder 查找参考基因组的 STRs 并且构建 bed 文件（HipSTR 参考文件）
HipSTR 通过分析测序 reads 与参考基因组（此处为登录号 Col-0）中检测到的 STR 的比对情况来查找STR

首先需要明确参考基因组中有哪些区域是 STR，并进行注释（bed 文件）

bed 文件格式如下，前1-5列是必须的，分别是 STR 所在染色体，起始位置，终止位置，最小重复单元碱基数，重复次数。后边列为可选，分别是 STR 名称，以及序列。
```
# 以人中的 STR 为例
chr11   2171088 2171115 4   7   TH01    N/A
chr12   5983969 5984044 4   19  vWa N/A
```

+ Tandem Repeats Finder 查找参考基因组的 STRs
```bash
mkdir TRF
cd TRF

# 从官网下载 trf409.linux64 文件
mv trf409.linux64 trf
trf # 运行
```
+ TRF使用说明
```
Usage：trf File Match Mismatch Delta PM PI Minscore MaxPeriod [options]

Match, Mismatch, and Delta: 匹配上，没匹配上和插入的权重，建议2 7 7
PM and PI ：PM是指比上的概率，可选择数值为80 和75，PI 是插入的概率，可选择数值为10 和20，最好效果的参数是PM=80 和PI=10
Minscore: 被匹配上的串联重复序列的最小分值。比如，我们设定了Match=2，Minscore=50， 那么就要求最少有25bp 被完全比上（比如，5bp 的重复单元，重复 5 次）
Maxperiod: 最大的重复单元 bp 数

参数
-m: 该参数将输入文件中trf序列屏蔽为N输出
-f: 该参数将输出每一串联重复序列两侧200bp 的侧翼序列，输出到比对文件中
-d: 该参数将产生一个屏蔽文件，记录了与列表文件一样的信息，及比对信息，可用于后续程序的处理
-h：不输出网页版结果
```

+ 注释
```bash
cd reference
faops split-name Arabidopsis_thaliana.TAIR10.dna.toplevel.fa .
# 参考基因组被拆分成了 7 个文件（ 5条染色体 + 线粒体Mt + 叶绿体Pt ）

cd TRF
for i in $(seq 5);do
    trf ../sequence/reference/$i.fa 2 7 7 80 10 50 500 -f -d -h
done
```

+ .dat结果解读
```
1.Indices of the repeat relative to the start of the sequence.
相对于序列开头的距离（两列）

2.Period size of the repeat.
每一个重复单元的大小

Number of copies aligned with the consensus pattern.
重复单元的数量

Size of consensus pattern (may differ slightly from the period size).
因为重复的次数可能不是整数倍，所有重复单元共有的序列，可能略小于重复单元大小

Percent of matches between adjacent copies overall.


Percent of indels between adjacent copies overall.


Alignment score.
比对得分

Percent composition for each of the four nucleotides.
四种核苷酸中每一种的百分比组成（四列）

Entropy measure based on percent composition.
基于百分比组成的熵的测量
```

+ bed文件制作
```bash
# 预处理（整理成矩阵）
JOB=$(ls *.dat)
for J in $JOB;do
    sed -i '1,15d' $J
    export CHR=$(echo $J | cut -d "." -f 1)
    cat $J | perl -ne'
        print "chr" . "$ENV{'CHR'}\t" . "$_"
    ' > tem&&
    mv tem $J
    sed -i 's/\s/\t/g' $J
    mv $J $CHR.dat
done

# 筛选出其中的短串联重复位点（ 1-6bp ）
JOB=$(ls *.dat)
for J in $JOB;do
    tsv-filter --le 4:6 $J > tem&&
        mv tem $J
done

# 构建bed文件
for J in $JOB;do
    cat $J >> Arabidopsis_thaliana_STR.tsv
done 

cat Arabidopsis_thaliana_STR.tsv | perl -a -F"\t" -ne'
    $i = $i + 1;
    print "@F[0]\t@F[1]\t@F[2]\t@F[3]\t@F[4]\t" . "Arabidopsis_thaliana_STR_$i\t" . "@F[14]\n";
' > Arabidopsis_thaliana_STR.bed

sed -i 's/chr//1' Arabidopsis_thaliana_STR.bed
```

## 1.4 使用 HipSTR 查找不同 accession 中的 STR
+ HipSTR下载
```bash
mkdir biosoft
cd biosoft

git clone https://github.com/HipSTR-Tool/HipSTR
cd HipSTR
make
export PATH="$(pwd):$PATH"
source $HOME/.bashrc

./HipSTR --help
```

+ 查找STR
```bash
mkdir STR
cd STR

HipSTR --bams ../BWA/SRR1945464.sort.bam,../BWA/SRR1945477.sort.bam,../BWA/SRR1945478.sort.bam,../BWA/SRR1945487.sort.bam --fasta ../sequence/reference/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa --region ../TRF/Arabidopsis_thaliana_STR.bed --str-vcf AT_vcf.gz 

# 出现了错误
# ERROR: Provided BAM/CRAM files don't contain read groups in the header and the --bam-samps flag was not specified
# 利用 BWA 的 -R 参数或者 hisat2 的 --rg-id 参数再次进行比对或者
# 在 HipSTR 后面加上--bam-samps 与 --bam-libs 参数

samtools faidx ../sequence/reference/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa

HipSTR --bams ../hisat2/output/SRR1945464.sort.bam,../hisat2/output/SRR1945477.sort.bam,../hisat2/output/SRR1945478.sort.bam,../hisat2/output/SRR1945487.sort.bam --fasta ../sequence/reference/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa --region ../TRF/Arabidopsis_thaliana_STR.bed --str-vcf AT_vcf.gz --bam-samps SAMPLE1,SAMPLE2,SAMPLE3,SAMPLE4 --bam-libs SRR1945464,SRR1945477,SRR1945478,SRR1945487

# 出现了错误
# ERROR: Invalid CIGAR option encountered in trimAlignment
```

+ BAM文件中的CIGAR参数

>CIGAR 是 Compact Idiosyncratic Gapped Alignment Report的首字母缩写，称为“雪茄”字符串。

>作为一个字符串，它用数字和几个字符的组合形象记录了read比对到参考序列上的细节情况，读起来要比FLAG直观友好许多，只是记录的是不同的信息。比如，一条150bp长的read比对到基因组之后，假如看到它的CIGAR字符串为：33S117M，其意思是说在比对的时候这条read开头的33bp在被跳过了（S），紧接其后的117bp则比对上了参考序列（M）。这里的S代表软跳过（Soft clip），M代表匹配（Match）。CIGAR的标记字符有“MIDNSHP=XB”这10个，分别代表read比对时的不同情况.











## 1.5 查找 R 基因
+ 直接使用文章附录中的数据（STR文件：210217.SuppDataSet2.DiploidUnitNumberCalls.tsv 以及基因表达量文件：210217.ExtraMaterial.logX.tsv）
### 1.6 样本 R 基因表达量统计
>**R基因**
>
>植物在漫长的生命历程中往往面临着多种多样逆境的威胁，所以植物不得不进化出一系列复杂的防御机制以抵抗环境中的真菌、细菌、病毒等病原体。抗性基因（Resistance gene，R gene）作为一类重要的植物防御基因，可以直接或间接地通过基因和基因之间的互作作用特异性识别病原体的产物来产生信号并传导信号，引起强烈的防御作用。
>
>大部分植物的 R 蛋白的特征结构域包括核酸结合位点（nucleotide-binding site，NBS）,富含亮氨酸（GAC --> CUG）的重复序列（leucine-rich repeat，LRR）,Toll-白细胞介素-1受体结构域（TIR）等，NBS 结构域在植物产生和传递防御信号中有着非常重要的作用。
>
>植物先天免疫系统由两个主要的免疫反应组成，即病原相关分子模式激发的免疫反应（PTI）和效应蛋白激发的免疫反应（ETI）。其中，PTI主要由病原微生物表面的病原相关分子模式（如多糖、鞭毛蛋白等）刺激诱导，可导致植物产生非特异性的防卫反应（基础防卫反应）；ETI则由植物的抗病蛋白（R蛋白）识别病原微生物产生的效应蛋白引发，可使植物产生特异性的防卫反应。除此而外，植物小RNA途径则通过RNA沉默方式，参与了对病毒和细菌等病原的抗性反应。

+ R基因数据库[PRGDB](http://prgdb.crg.eu/wiki/Download)
```bash
wget http://prgdb.crg.eu/data/annotation/all_annotation_all.tsv.gz
gzip -d all_annotation_all.tsv.gz
# 这个文件里面包含了以及确定的R基因（reference，25个）以及推断的R基因（putative）
# 格式：PRGID\tName\tType\tSpecies\tClass\tGenBank ID\tGenBank Locus\tDescription

tsv-filter --str-eq 4:'Arabidopsis thaliana' all_annotation_all.tsv | wc -l
# 2871
tsv-filter --str-eq 4:'Arabidopsis thaliana' all_annotation_all.tsv > ATRgene.tsv

# 210217.ExtraMaterial.logX.tsv 中的基因为 Gene symbol 格式而 ATRgene.tsv 中只含有 GenBank ID 和GenBank Locus编号，进行ID转换
# 在线网站[bioDBnet](https://biodbnet-abcc.ncifcrf.gov/)进行 ID 转换
```
+ [TAIR 数据库](https://www.arabidopsis.org/index.jsp)

输入关键词“NBS-LRR”进行查询

```bash
sed '1d' NBS-LRR.csv | cut -f 2 | sort | uniq > NBS-LRR.lst

cat ../210217.ExtraMaterial.logX.tsv | datamash transpose | sed '1d' | cut -f 1 | grep -f NBS-LRR.lst | wc -l
# 86 个基因有表达量

LIST=$(cat ../210217.ExtraMaterial.logX.tsv | datamash transpose | sed '1d' | cut -f 1 | grep -f NBS-LRR.lst | perl -ne'
    chomp($_);
    print "$_,";
' | sed 's/,$//')

tsv-select -H --fields 1,$LIST ../210217.ExtraMaterial.logX.tsv > Rexp.tsv

# 构建表型文件
sed '1d' Rexp.tsv | perl -a -F"\t" -ne'
    print "@F[0]\t$_";
' > pheno.txt
```

## 1.7 关联分析
+ STR 信息整理
> 210217.SuppDataSet2.DiploidUnitNumberCalls.tsv 中统计了 STR 重复单元的数量

```bash
# 构建 .map 文件
cat ../210217.SuppDataSet2.DiploidUnitNumberCalls.tsv | datamash transpose | 
    sed '1d' | cut -f 1 | 
    perl -a -F"_" -ne'
    print "@F[0]\t.\t0\t@F[1]";
    ' > STR.map
sed -i 's/^chr//' STR.map

# 构建 .ped 文件，从第七列开始，每两列代表一个基因型
sed '1d' ../210217.SuppDataSet2.DiploidUnitNumberCalls.tsv | cut -f 1 | perl -ne'
    chomp($_);
    print "$_\t$_\t0\t0\t0\t-9\n";
' > STR.ped

REPEAT_TIMES=$(cat ../210217.SuppDataSet2.DiploidUnitNumberCalls.tsv | datamash transpose | sed '1d' | wc -l)
for i in $(seq $REPEAT_TIMES);do
    for j in 1 2;do
    START_LINE=$(( $i + 1 ))
#    echo $START_LINE
    tsv-join --filter-file ../210217.SuppDataSet2.DiploidUnitNumberCalls.tsv --key-fields 1 --append-fields $START_LINE STR.ped > tem&&
        mv tem STR.ped
    done
done
# 计算时间较长，是否可以用 parallel 并行
# 运算结束发现 .ped 中间有些空缺值，用 0 填补
sed -i 's/\t\t/\t0\t/g' STR.ped  # 重复至少两次
cat STR.ped | perl -ne'
    s/\t\n/\t0\n/g;
    print "$_";
' > tem&&
mv tem STR.ped 
# 由于 sed 无法匹配行末尾的换行符，所以用 perl 对行末尾的空缺值进行填补


# 关联分析
parallel -k -j 4 "
#    echo {1}
    plink2 -file ./STR --linear --pheno ../TAIR/pheno.txt --mpheno {1} -noweb --allow-no-sex --out mydata{1}
" ::: $(seq 86)
# Error: More than 4 different alleles at variant 303 (post-sort/filter).
```
> STR 的变异不像 SNP，只有四种形式，STR 重复次数的多样性更高，如何解决？

```bash


```


## 1.8 可视化

```bash


```

## 参考
[1、BWA 使用详解](https://www.jianshu.com/p/3b86615d647b?utm_campaign=maleskine&utm_content=note&utm_medium=seo_notes&utm_source=recommendation)

[2、通过HipSTR 对二代测序（NGS）数据进行 短串联重复序列（short tandem repeats，STR）检测](https://blog.csdn.net/qq_42962326/article/details/105058661)

[3、TRF（tandem repeats finder ）安装与使用（二）](https://www.jianshu.com/p/66d562fb47de)

[4、基因组注释--重复序列注释（一）：Trf软件安装与使用](https://www.jianshu.com/p/efe2d66f8969)

[5、HipSTR](https://github.com/HipSTR-Tool/HipSTR)

[6、Length variation in short tandem repeats affects gene expression in natural populations of Arabidopsis thaliana](https://academic.oup.com/plcell/article/33/7/2221/6225030?login=true)

[7、BEM文件的：flags、CIGAR、MAPQ](https://www.jianshu.com/p/ede7735eda20)

