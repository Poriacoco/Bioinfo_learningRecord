#############################
#      下载数据             #
#############################
mkdir 03.sequencing
cd 03.sequencing
cp /ifs1/VipData/03.sequencing/data/*.sra ./
#1 安装sratools
wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.10.8/sratoolkit.current-centos_linux64.tar.gz
tar -zxvf sratoolkit.current-centos_linux64.tar.gz 

#2 安装aspera
wget https://download.asperasoft.com/download/sw/connect/3.9.9/ibm-aspera-connect-3.9.9.177872-linux-g2.12-64.tar.gz
tar -zxvf ibm-aspera-connect-3.9.9.177872-linux-g2.12-64.tar.gz
sh ibm-aspera-connect-3.9.9.177872-linux-g2.12-64.sh

#下载illumina数据
prefetch SRR8482586 -O ./ -o illumina.sra
fastq-dump --split-files --gzip illumina.sra
#下载pacbio数据
prefetch SRR8494912 -O ./ -o pacbio.sra
fasterq-dump pacbio.sra
#下载nanopore数据
prefetch SRR8494939 -O ./ -o nanopore.sra
fasterq-dump nonopore.sra
gzip nanopore.fastq

#基因组下载
#Klebsiella pneumoniae MGH78578
#基因组： NC_009648.1 
https://www.ncbi.nlm.nih.gov/nuccore/NC_009648.1/
#质粒： NC_009649.1
https://www.ncbi.nlm.nih.gov/nuccore/NC_009649

#批量下载
prefetch --list SRR_Acc_List.txt -O ./


#############################
#     fastq文件处理          #
#############################
#1 压缩与解压缩
#解压缩
gunzip illumina_1.fastq.gz
gzip -d  illumina_2.fastq.gz
#压缩
gzip illumina_1.fastq
gzip illumina_2.fastq

#2 fastq文件统计
seqkit stats illumina_1.fastq.gz illumina_2.fastq.gz

#3 统计fastq文件每条序列ATCG四种碱基组成以及质量值分布
seqtk comp illumina_1.fastq.gz illumina_2.fastq.gz

#4 ATCG以及质量值分布
seqtk fqchk illumina_1.fastq.gz
seqtk fqchk illumina_2.fastq.gz

#57 交叉合并pairend文件
seqtk mergepe illumina_1.fastq.gz illumina_2.fastq.gz >merge.fastq

#6 过滤短的序列
#过滤小于150bp序列，并压缩输出 
seqkit seq -m 150 nanopore.fastq.gz | gzip -  >filter_150.fq.gz
#保留小于150bp序列
seqkit seq -M 150 nanopore.fastq.gz

#7 转换为列表格式
seqkit fx2tab nanopore.fastq.gz

#8 分别统计每一条序列长度
seqkit fx2tab nanopore.fastq.gz |awk -F "\t" '{print length($2) }'

#9 质量值转换
#将illumina 1.8转换为1.5
seqkit convert --to Illumina-1.5+ illumina_1.fastq.gz |head -4
#将illumina 1.5转换为1.8，什么都不加就是转换为1.8
seqkit convert  illumina_illmina1.5.gz

#10 排序，按照长度
seqkit sort -l -r nanopore.fastq.gz

#11 #seqkit抽样，按照百分比
seqkit sample -p 0.1 illumina_1.fastq.gz

#12 seqkit抽样，按照条数
seqkit sample -n 1000 illumina_1.fastq.gz 

#13 拆分数据
seqkit split2 -1 illumina_1.fastq.gz -2 illumina_2.fastq.gz -p 2 -f

#14 转换为fasta
#seqkit工具
seqkit fq2fa nanopore.fastq.gz >nanopore.fasta

#15 只输出20行ID
seqkit seq -n -i nanopore.fastq.gz |head -20 >id.list

#16 提取序列
seqkit grep -f id.list nanopore.fastq.gz

#17 截取头尾
seqtk trimfq -b 15 -e 15 -Q illumina_1.fastq.gz

#17 修改reads ID
seqkit replace -p "SRR8494939\.sra" -r 'reads'  nanopore.fastq.gz 

#18 长度分布直方图
seqkit watch -L -f ReadLen hairpin.fa

#19 平均质量直方图
seqkit watch -p 500 -O qhist.pdf -f MeanQual nanopore.fastq.gz 

#20 选取固定范围
seqkit range -r 200:300 nanopore.fastq.gz

#21 移除重名序列
seqkit rmdup -n nanopore.fastq.gz  -o clean.fa.gz

#22 将小于Q20的替换为小写字母
seqtk seq -q 20 illumina_1.fastq.gz


#############################
#     illumina         #
#############################
#fastqc质控
mkdir illumina_qc
fastqc -f fastq -o illumina_qc/ illumina_1.fastq.gz illumina_2.fastq.gz

#fastp 数据过滤
fastp -i illumina_1.fastq.gz  -I illumina_2.fastq.gz -o clean.1.fq.gz -O clean.2.fq.gz -z 4 -q 20 -u 40 -n 10 -f 15 -t 15 -F 15 -T 15 -h fastp.html

#过滤完质控
mkdir illumina_clean
fastqc -f fastq -o illumina_clean/ clean.1.fq.gz  clean.2.fq.gz


#############################
#     pacbio          #
#############################
#fastqc质控
mkdir pacbio_qc/
fastqc -f fastq -o pacbio_qc/ pacbio.fastq.gz

#过滤数据
filtlong --min_length 300 --min_mean_q 90 pacbio.fastq.gz |  gzip >pacbio.filtlong.fq.gz

#质控完过滤
mkdir pacbio_clean
fastqc -f fastq -o pacbio_clean/ pacbio.filtlong.fq.gz

#############################
#     Basecalling           #
#############################
#链接数据到当前目录下
ln -s /ifs1/VipData/03.sequencing/data/fast5_files/ ./
ll
#1 Basecalling
guppy_basecaller -i fast5_files/ -s fastq_files --config dna_r9.4.1_450bps_fast.cfg -r 
#合并全部fastq并压缩
cat fastq_files/*.fastq | gzip >lambda.fastq.gz

#############################
#        barcoding          #
#############################
#basecalling同时拆分barcode
guppy_basecaller -i /ifs1/VipData/03.sequencing/data/16s_fast5/ -s fastq --config dna_r9.4.1_450bps_hac.cfg -r  -x cuda:all --trim_barcodes --barcode_kits SQK-RBK004 --min_score 10

#basecalling与拆分barcode分布完成
#首先进行basecalling
guppy_basecaller -i /ifs1/VipData/03.sequencing/data/16s_fast5/ -s fastq --config dna_r9.4.1_450bps_hac.cfg -r  -x cuda:all
#拆分barcode
guppy_barcoder -i  fastq -s barcode --trim_barcodes --barcode_kits SQK-RBK004 --min_score 10

#############################
#          数据质控与过滤    #
#############################
#激活nanoplot环境
conda activate nanoplot

#NanoPlot质控
NanoPlot --fastq nanopore.fastq.gz -o nanoplot
#过滤数据
filtlong --min_length 2000 --min_mean_q 90 nanopore.fastq.gz |  gzip >nanopore.filtlong.fq.gz
#过滤完质控
NanoPlot --fastq nanopore.filtlong.fq.gz -o clean
conda deactivate
