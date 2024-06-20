########################################################
#                安装软件                          #
########################################################
#安装bioconda
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh  
sh Miniconda3-latest-Linux-x86_64.sh  
source ~/.bashrc

#添加软件源
conda config --add channels bioconda 
conda config --add channels conda-forge

#使用bioconda进行安装
conda install -y mamba
mamba install -y spades

mamba install -y soapdenovo2
#SOAPdenovo2 installation error 

#/usr/bin/ld: failed to set dynamic section sizes: bad value
#collect2: error: ld returned 1 exit status
#make: *** [Makefile:58: SOAPdenovo-63mer] Error 1)

#make clean
#Edit makefile. Add -no-pie in Line 58 and Line 62.
#@$(CC) sparsePregraph/*.o standardPregraph/*.o -no-pie $(LIBPATH) $(LIBS) $(EXTRA_FLAGS) -o SOAPdenovo-63mer
#@$(CC) sparsePregraph/*.o standardPregraph/*.o -no-pie $(LIBPATH) $(LIBS) $(EXTRA_FLAGS) -o SOAPdenovo-127mer
#create a new conda env and install gcc9 by conda install  -c conda-forge cxx-compiler=1.2.0 -y
#make

mamba install -y flye
mamba install -y canu
mamba install -y wgsim
mamba install -y megahit
mamba install -y gapfiller
mamba install -y jellyfish
mamba install -y gapfiller
mamba install -y soapdenovo2-gapcloser
mamba create -n medaka -y medaka
mamba create -n  -y quast
mamba install -y genomescope2
mamba create -n unicycler -y unicycler

#安装busco
#软件官网：https://gitlab.com/ezlab/busco
#数据库下载：https://www.orthodb.org/?page=filelist
#数据库下载：https://busco-data.ezlab.org/v4/data/lineages/
mamba create -n busco -y busco=5.2.2

#列出数据库
busco --list-datasets
#下载数据
cd miniconda3/envs/busco/busco_downloads
grep "bacteria" file_versions.tsv
busco --download bacteria_odb10

########################################################
#                 准备数据                              #
########################################################
mkdir data;cd data;
#文章地址
https://www.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000595

#1 下载参考序列
#基因组下载
#Klebsiella pneumoniae MGH78578
#基因组： NC_009648.1 
https://www.ncbi.nlm.nih.gov/nuccore/NC_009648.1/
#质粒： NC_009649.1
https://www.ncbi.nlm.nih.gov/nuccore/NC_009649

efetch -db nuccore -format fasta -id NC_009648.1 > MGH78578.fasta
efetch -db nuccore -format fasta -id NC_009649 >>MGH78578.fasta

#Escherichia coli CFT073
efetch -db nuccore -format fasta -id NC_004431.1 >CFT073.fasta


#下载测序数据
#https://www.ncbi.nlm.nih.gov/bioproject/PRJNA422511 

#下载illumina数据
prefetch SRR8482586 -O ./ -o illumina.sra
fastq-dump --split-files --gzip illumina.sra

#下载pacbio数据
prefetch SRR8494912 -O ./ -o pacbio.sra
fasterq-dump pacbio.sra
gzip pacbio.sra.fastq
#下载nanopore数据
prefetch SRR8494939 -O ./ -o nanopore.sra
fasterq-dump nonopore.sra
gzip nanopore.sra.fastq

########################################################
#             35 二代数据处理                           #
########################################################
#fastqc质控
mkdir 1.qc
fastqc -f fastq -o qc illumina.sra_1.fastq.gz illumina.sra_2.fastq.gz

#利用fastp进行数据过滤  
fastp -i illumina.sra_1.fastq.gz -I illumina.sra_2.fastq.gz -o clean.1.fq.gz -O clean.2.fq.gz -D -z 4 -q 20 -u 30 -n 10 -f 20 -t 10 -F 20 -T 10 -h clean.html

#过滤完在质控一遍
fastqc -f fastq -o 1.qc clean.1.fq.gz clean.2.fq.gz

#  kmer估计基因组大小
mkdir 34.kmer;cd 34.kmer/
ln -s ../data/illumina.sra_1.fastq.gz .
ln -s ../data/illumina.sra_2.fastq.gz .
#jellyfish统计kmer
#不支持压缩格式
jellyfish count -m 15 -o kmer15.count  -s 2G -C <(zcat ../data/clean.1.fq.gz )
jellyfish count -m 15 -o kmer15.count  -s 2G -C ../data/clean.fq
jellyfish stats -o kmer15.stats  kmer15.count
jellyfish histo -t 3  kmer15.count > kmer15.histo

#genomescope2估计基因组大小
genomescope2 -i kmer15.histo -o gscope -p 1 -k 15

########################################################
#            36 二代测序数据拼接                           #
########################################################
mkdir 36.illumina;cd 36.illumina
#spaces拼接
/ifs1/Software/biosoft/SPAdes-3.14.0-Linux/bin/spades.py -o spades_result -t 12 -1 ../data/clean.1.fq.gz -2 ../data/clean.2.fq.gz  1>spades.log 2>spades.err

#SOAPdenovo配置文件
max_rd_len=150
[LIB]
avg_ins=439
reverse_seq=0
asm_flags=3
rank=1
pair_num_cutoff=3
q1=#reads1文件路径 (),这个需要替换，不要那么教条直接复制这段
q2=#reads2文件路径 (),这个需要替换，不要那么教条直接复制这段

#soapdenovo2运行
mkdir kmer45  
SOAPdenovo-63mer all -s lib.list -K 45 -o kmer45/kmer45 -D 1 -d 1 -u 2 -p 12 >kmer45.log 

# 补洞 
GapCloser -a kmer45/kmer45.scafSeq -b lib.list -o kmer45.fill.fa -l 100 -p 25 -t 12 >gapcloser.log

########################################################
#           37 拼接结果评估                              #
########################################################
mkdir 37.evaluate;cd 37.evaluate
# fasta文件格式处理案例
#案例一：统计
seqkit stats kmer45.scafSeq 
#分别统计每一条序列长度
seqkit fx2tab kmer45.scafSeq |awk '{print $1"\t"length($2)}'

#案例二：格式化
seqkit seq -w 0 kmer45.scafSeq  #seqtk  seq -l 0 kmer45.scafSeq
#每行显示50个碱基
seqtk  seq -l 50 kmer45.scafSeq

#案例三：逐条统计
seqtk  seq -l 0 kmer45.scafSeq  | grep -v ">" | awk '{print length($0)}' | head
#统计长度并按照长度计算频数
seqtk  seq -l 0 kmer45.scafSeq  |grep -v ">" | awk '{print length($0)}' | sort | uniq -c

#案例四：成分分析
seqtk comp kmer45.scafSeq | head

#案例五：提取序列
seqkit grep -r -p "C2877" kmer45.scafSeq 

#案例六：截取序列
seqkit subseq -r 1000:3000 kmer45.scafSeq 
seqkit subseq -r 1000:3000 kmer45.scafSeq --chr C2689

#案例七：排序
seqkit sort -l -r kmer45.scafSeq | less -S

#案例八：按照长度过滤
seqkit seq -m 1000 kmer45.scafSeq
#过滤长度大于1000bp序列
seqkit seq -M 1000 kmer45.scafSeq

# 案例九：反向互补
#seqkit取反向序列
seqkit seq -r test.fasta
#seqkit seq 加-r -p同时取反向互补序列
seqkit seq -r -p test.fasta

#案例十：转换大小写
seqkit seq -l kmer45.scafSeq| head
seqkit seq -u kmer45.scafSeq| head

#2 quast评估
mkdir 1.quast
#quast
conda activate quast
quast.py -r MGH78578.fasta spades.fa soapdenovo.fa -o quast
quast -o quast -r GCF_000240185.1_ASM24018v2_genomic.fna -t 12 -g GCF_000240185.1_ASM24018v2_genomic.gff soapdenovo.fa spades.fa  --glimmer
#3 busco
mkdir 2.busco
#激活虚拟环境
conda activate busco
#列出数据库
busco --list-datasets

#运行busco
busco -i spades.fa -o busco_result -m geno -l busco_downloads/lineages/bacteria_odb10/ -c 20 --offline 
#利用busco结果绘图
generate_plot.py -wd busco_result

#tablet
#与拼接结果比对
bwa-mem2 index spades.fa
bwa-mem2 mem -t 12 spades.fa ../data/clean.1.fq.gz ../data/clean.2.fq.gz >spades.sam

#samtools处理比对结果
samtools sort -@ 12 -o spades.sorted.bam spades.sam 
samtools index spades.sorted.bam

#对spades.fa建立索引
samtools faidx spades.fa
#tablet可视化
#将文件拷贝至windows下使用tablet可视化
mkdir tablet
mv spades.fa spades.fa.fai spades.sorted.bam spades.sorted.bam.bai tablet

#Bandage可视化
https://rrwick.github.io/Bandage/

########################################################
#                    38 基因组拼接探索                  #
########################################################
mkdir 38.Test;cd 38.test;
#不同kmer大小
mkdir 1.kmer
SOAPdenovo-31mer all -s lib.list -K 31 -o kmer31 -D 1 -d 1 -u 2 -p 12 >kmer31.log 
SOAPdenovo-63mer all -s lib.list -K 63 -o kmer63 -D 1 -d 1 -u 2 -p 12 >kmer63.log 
SOAPdenovo-127mer all -s lib.list -K 127 -o kmer127 -D 1 -d 1 -u 2 -p 12 >kmer127.log 

#大片段文库
max_rd_len=90
[LIB]
avg_ins=500
reverse_seq=0
asm_flags=3
rank=1
pair_num_cutoff=3
q1=/ifs1/Sequencing/cleandata/500_clean.1.fq.gz
q2=/ifs1/Sequencing/cleandata/500_clean.2.fq.gz

[LIB]
avg_ins=2000
reverse_seq=1
asm_flags=2
rank=3
pair_num_cutoff=3
q1=/ifs1/Sequencing/cleandata/2000_clean.1.fq.gz
q2=/ifs1/Sequencing/cleandata/2000_clean.2.fq.gz

#数据质量
/ifs1/Software/biosoft/SPAdes-3.14.0-Linux/bin/spades.py -o raw_result -t 12 -1 ../data/illumina.sra_1.fastq.gz -2 ../data/illumina.sra_1.fastq.gz  1>spades.log 2>spades.err

#数据量大小
#分别抽取10%，30%，50%，80%进行比较
seqkit sample -p 0.1 -s 1234 ../data/illumina.sra_1.fastq.gz | gzip >reads.0.1_1.fq.gz
seqkit sample -p 0.1 -s 1234 ../data/illumina.sra_2.fastq.gz | gzip >reads.0.1_2.fq.gz

for i in {0.1,0.3,0.5,0.8};do 
    seqkit sample -p ${i} -s 1234 ../data/illumina.sra_1.fastq.gz | gzip >reads.${i}_1.fq.gz
    seqkit sample -p ${i} -s 1234 ../data/illumina.sra_2.fastq.gz | gzip >reads.${i}_2.fq.gz;
done;

ls -1 *.fq.gz | xargs -n 2 | while read {i,j};do echo /ifs1/Software/biosoft/SPAdes-3.14.0-Linux/bin/spades.py -o spades_${i} -t 12 -1 ${i} -2 ${j};done;


#不同reads长度
#利用wgsim分别模拟长度70与150bp长度reads
wgsim  -1 70 -2 70 MGH78578.fasta read.70_1.fq read.70_2.fq 
wgsim  -1 150 -2 150 MGH78578.fasta read.150_1.fq read.150_2.fq 

#不同选项参数


#不同错误率的影响
wgsim -e 0.01 -N 1000000 -1 150 -2 150 MGH78578.fasta read.0.01_1.fq read.0.01_2.fq 
wgsim -e 0.1 -N 1000000 -1 150 -2 150 MGH78578.fasta read.0.1_1.fq read.0.1_2.fq


########################################################
#                         39 混合拼接                     #
########################################################
mkdir 39.hybrid;cd 39.hybrid;
#单独使用illumina数据拼接，写到脚本里
/ifs1/Software/biosoft/SPAdes-3.14.0-Linux/bin/spades.py -o illumina -t 24 -1 ../data/clean_1.fastq.gz -2 ../data/clean_2.fastq.gz

#利用illumina数据+pacbio数据拼接
/ifs1/Software/biosoft/SPAdes-3.14.0-Linux/bin/spades.py -o ill_ont -t 24 -1 ../data/clean_1.fastq.gz -2 ../data/clean_2.fastq.gz --pacbio ../data/pacbio.sra.fastq.gz
#利用illumina数据+nanopore数据拼接
/ifs1/Software/biosoft/SPAdes-3.14.0-Linux/bin/spades.py -o ill_pac -t 24 -1 ../data/clean_1.fastq.gz -2 ../data/clean_2.fastq.gz --nanopore ../data/nanopore.sra.fastq.gz 

#unicycler混合拼接
conda activate unicycler
unicycler -1 ../data/clean_1.fastq.gz -2 ../data/clean_2.fastq.gz --long ../data/nanopore.sra.fastq.gz -o unicycler -t 12

#quast比较不同拼接方案
quast.py -r MGH78578.fasta illumina.fa ill_pac.fa ill_ont.fa unicycler.fa -o quast

########################################################
#                  40 三代测序数据拼接                    #
########################################################
mkdir 40.nanopore;cd 40.nanopore
#pacbio测序拼接
flye --pacbio-raw ../data/pacbio.sra.fastq.gz --genome-size 5.4m --threads 12 --out-dir flye

#nanopore数据处理
#激活nanoplot环境
conda activate nanoplot

#NanoPlot质控
NanoPlot --fastq ../data/nanopore.fastq.gz -o nanoplot
#过滤数据
filtlong --min_length 2000 --min_mean_q 90 ../data/nanopore.fastq.gz |  gzip >clean.filtlong.fq.gz
#过滤完质控
NanoPlot --fastq clean.filtlong.fq.gz -o clean
conda deactivate


#nannopre序列拼接
#1 flye案例
time flye --nano-raw ../data/clean.filtlong.fq.gz --genome-size 5.4m --threads 12 --out-dir flye

#2 canu案例
time canu -d canu -p canu genomeSize=5.4m maxThreads=24 -nanopore-raw ../data/clean.filtlong.fq.gz

#3 wtdbg2案例
mkdir wtdbg2;
perl /ifs1/Software/biosoft/wtdbg2/wtdbg2.pl -t 12 -x ont -g 5.4m -o wtdbg2/wtdbg2 ../data/clean.filtlong.fq.gz
