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
for i in SRR1945464 SRR1945477 SRR1945478 SRR1945487;do
    bwa mem -t 4 ref ../sequence/$i'_1'.fastq  ../sequence/$i'_2'.fastq > $i.sam
done

parallel -k -j 4 "
    samtools sort -@ 4 {1}.sam > {1}.sort.bam
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
parallel -k -j 4 "
    samtools sort -@ 4 {1}.sam > {1}.sort.bam
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
# 尝试利用hisat2再次进行比对


HipSTR --bams ../hisat2/output/SRR1945464.sort.bam,../hisat2/output/SRR1945477.sort.bam,../hisat2/output/SRR1945478.sort.bam,../hisat2/output/SRR1945487.sort.bam --fasta ../sequence/reference/genome.fa --region ../TRF/Arabidopsis_thaliana_STR.bed --str-vcf AT_vcf.gz 
```




## 参考
[1、BWA 使用详解](https://www.jianshu.com/p/3b86615d647b?utm_campaign=maleskine&utm_content=note&utm_medium=seo_notes&utm_source=recommendation)

[2、通过HipSTR 对二代测序（NGS）数据进行 短串联重复序列（short tandem repeats，STR）检测](https://blog.csdn.net/qq_42962326/article/details/105058661)

[3、TRF（tandem repeats finder ）安装与使用（二）](https://www.jianshu.com/p/66d562fb47de)

[4、基因组注释--重复序列注释（一）：Trf软件安装与使用](https://www.jianshu.com/p/efe2d66f8969)

[5、HipSTR](https://github.com/HipSTR-Tool/HipSTR)

[6、Length variation in short tandem repeats affects gene expression in natural populations of Arabidopsis thaliana](https://academic.oup.com/plcell/article/33/7/2221/6225030?login=true)



