#############################
#      准备工作：安装软件     #
#############################
mkdir 04.ncov;cd 04.ncov
conda install -y blast
conda install -y entrez-direct
conda install -y tablet
conda install -y megahit
conda install -y blat
conda install -y lastz
conda install -y seqkit
conda install -y bwa
conda install -y bwa-mem2
conda install -y minimap2
conda install -y samtools
conda install -y bcftools
conda install -y freebayes
conda install -y ivar
conda install -y treebest

#artic network
#https://github.com/artic-network/fieldbioinformatics
#1 下载
git clone https://github.com/artic-network/fieldbioinformatics
cd fieldbioinformatics
#2 安装依赖
conda env create -f environment.yml
conda activate artic

#3 安装流程
python setup.py install
#4 测试流程
artic -v

#centrifuge
wget https://codeload.github.com/DaehwanKimLab/centrifuge/tar.gz/refs/tags/v1.0.4 -O centrifuge-1.0.4.tar.gz
tar -zxvf centrifuge-1.0.4.tar.gz
cd centrifuge-1.0.4
make

#pangolin
git clone https://github.com/cov-lineages/pangolin.git
cd pangolin
conda env create -f environment.yml
conda activate pangolin
pip install ./

#############################
#         28下载参考数据     #
#############################
mkdir 28.data;cd 28.data;
mkdir reference;cd reference;
#1 下载新冠参考序列
efetch -db nuccore -format fasta -id NC_045512 > NC_045512.fa
#下载SARS序列
efetch -db nuccore -format fasta -id KT444582 > SARS_ref.fa
#下载batSARS序列
efetch -db nuccore -format fasta -id MG772933 > batSARS_ref.fa
#下载新冠突变株
efetch -db nuccore -format fasta -id MZ310552 >UK_B117.fa;sed -i "s/>.*/>UK_B117	MZ310552/g" UK_B117.fa
efetch -db nuccore -format fasta -id MZ202314 >SouthAfrica_B1351.fa;sed -i "s/>.*/>SouthAfrica_B1351	MZ202314/g" SouthAfrica_B1351.fa
efetch -db nuccore -format fasta -id MZ169911 >Brazil_P1.fa;sed -i "s/>.*/>Brazil_P1	MZ169911/g" Brazil_P1.fa
efetch -db nuccore -format fasta -id MZ318159 >India_B16172.fa;sed -i "s/>.*/>India_B16172	MZ318159/g" India_B16172.fa
efetch -db nuccore -format fasta -id MZ373479 >USA_B1427.fa;sed -i "s/>.*/>USA_B1427	MZ373479/g" USA_B1427.fa
efetch -db nuccore -format fasta -id MZ169912 >Brazil_P2.fa;sed -i "s/>.*/>Brazil_P2	MZ169912/g" Brazil_P2.fa
efetch -db nuccore -format fasta -id MW852494 >USA_B1525.fa;sed -i "s/>.*/>USA_B1525	MW852494/g" USA_B1525.fa
efetch -db nuccore -format fasta -id MZ257684 >Philippines_P3.fa;sed -i "s/>.*/>Philippines_P3	MZ257684/g" Philippines_P3.fa
efetch -db nuccore -format fasta -id MZ310903 >USA_B1526.fa;sed -i "s/>.*/>USA_B1526	MZ310903/g" USA_B1526.fa
efetch -db nuccore -format fasta -id MZ310580 >india_B16171.fa;sed -i "s/>.*/>india_B16171	MZ310580/g" india_B16171.fa

#############################
#  准备工作：下载测序数据     #
#############################
#参考基因组测序数据
mkdir sra;cd sra
#prefetch SRR10971381 -O ./ 
wget https://storage.googleapis.com/nih-sequence-read-archive/run/SRR10971381/SRR10971381 -O SRR10971381.sra

#API 下载
ID=$1
wget https://storage.googleapis.com/nih-sequence-read-archive/run/$1/$1 -O ${1}.sra

#1 宏基因组测序
#illumina：
prefetch SRR11092056 -O ./ 
prefetch SRR11092057 -O ./
prefetch SRR11092058 -O ./
prefetch SRR11092059 -O ./
prefetch SRR11092060 -O ./
prefetch SRR11092061 -O ./
prefetch SRR11092062 -O ./
prefetch SRR11092063 -O ./
prefetch SRR11092064 -O ./

#nanopore：SRR11178050 -O ./
prefetch SRR11178051 -O ./
prefetch SRR11178052 -O ./
prefetch SRR11178053 -O ./
prefetch SRR11178054 -O ./
prefetch SRR11178055 -O ./
prefetch SRR11178056 -O ./
prefetch SRR11178057 -O ./

#2 PCR扩增产物
#illumina：
prefetch SRR14867660 -O ./	
#nanopore：
prefetch SRR14859525 -O ./
prefetch ERR6055891 -O ./	

#############################
#     google api下载数据     #
#############################
vim sra.sh 
#/usr/bin/bash

if [ $# == 0 ];then 
	echo "This shell is used to download sra data from Google apis
    Usage: sh sra.sh SRA1 SRA2 SRA3...
    ";	

fi

for i in $*;
	do echo "Download SRA Data ${i}\n";
	wget -c https://storage.googleapis.com/nih-sequence-read-archive/run/${i}/${i} -O ${i}.sra;
done;

sh sra.sh SRR11092056 SRR11092057 SRR11092058

#############################
#  准备工作：下载鉴定数据库   #
#############################
#blast数据库
~/.aspera/connect/bin/ascp -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh --overwrite=diff -QTr -l6000m anonftp@ftp.ncbi.nlm.nih.gov:blast/db/Betacoronavirus.00.tar.gz ./
~/.aspera/connect/bin/ascp -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh --overwrite=diff -QTr -l6000m anonftp@ftp.ncbi.nlm.nih.gov:blast/db/Betacoronavirus.01.tar.gz ./
~/.aspera/connect/bin/ascp -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh --overwrite=diff -QTr -l6000m anonftp@ftp.ncbi.nlm.nih.gov:blast/db/Betacoronavirus.02.tar.gz ./
~/.aspera/connect/bin/ascp -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh --overwrite=diff -QTr -l6000m anonftp@ftp.ncbi.nlm.nih.gov:blast/db/Betacoronavirus.03.tar.gz ./

#centrifuge数据库
wget https://zenodo.org/record/3732127/files/h+v+c.tar.gz?download=1
mkdir hpv
tar -zxvf h+v+c.tar.gz -C hpv

#############################
#     29 快速鉴定            #
#############################
mkdir 29.identify;cd 29.identify;
#SRA数据转换
fastq-dump --split-files --gzip SRR10971381.sra

#centrifuge进行物种分类鉴定
centrifuge -x /ifs1/MetaDatabase/centrifuge_h+v_v2_20200327/hv2 -1 SRR10971381.sra_1.fastq.gz -2 SRR10971381.sra_2.fastq.gz  --report-file report.tsv -S result.tsv -p 12 >centrifuge.log
centrifuge-kreport -x /ifs1/MetaDatabase/centrifuge_h+v_v2_20200327/hv2 result.tsv >kreport.txt

#############################
#      快速鉴定              #
#############################
#在线工具：https://fbreitwieser.shinyapps.io/pavian/
#将centrifuge-kreport的结果上传到pavian网站上

#############################
#    30 拼接新冠参考基因组    #
#############################
mkdir 30.ncov;cd 30.ncov;
#数据质控
mkdir qc
fastqc -f fastq  -o qc SRR10971381_1.fastq SRR10971381_2.fastq
#数据过滤
fastp -i SRR10971381.sra_1.fastq.gz -I SRR10971381.sra_2.fastq.gz -o clean.1.fq.gz -O clean.2.fq.gz -z 4 -q 20 -u 30 -n 10 -h clean.html
#过滤之后质控
mkdir clean
fastqc -f fastq -o clean clean.1.fq.gz clean.2.fq.gz 

#下载人基因组序列建立索引
~/.aspera/connect/bin/ascp -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh --overwrite=diff -QTr -l6000m anonftp@ftp.ncbi.nlm.nih.gov:1000genomes/ftp/technical/reference/human_g1k_v37.fasta.gz ./
gunzip -zxvf human_g1k_v37.fasta.gz
nohup time bwa-mem2 index human_g1k_v37.fasta &  #较长时间

#bwa-mem2 加快比对速度
bwa-mem2 mem -t 24 /ifs1/VipData/04.ncov/database/human_g1k_v37_bwamem2/human_g1k_v37.fasta clean.1.fq.gz clean.2.fq.gz >ncov.sam

#过滤宿主序列
samtools fastq -G 2 ncov.sam -1 ncov.1.fq -2 ncov.2.fq
gzip ncov.1.fq ncov.2.fq
#统计过滤前后数变化
seqkit stat clean.1.fq.gz clean.2.fq.gz
seqkit stat ncov.1.fq.gz ncov.2.fq.gz


#拼接基因组
vim megahit.sh
megahit -t 24 -o megahit/ -1 ncov.1.fq.gz -2 ncov.2.fq.gz 1>megahit.log 2>megahit.err

#统计每一条长度
seqkit fx2tab final.contigs.fa  |awk -F"\t" '{print $1,length($2)}'

#NCBI Blast比对
blastn -query final.contigs.fa -out blast.out -db /ifs1/MetaDatabase/Betacoronavirus/Betacoronavirus -outfmt 6 -evalue 1e-5 -num_threads 24

#筛选最长序列
#每次拼接出来最长序列ID不一样
seqkit grep -p "ID" final.contigs.fa #将ID替换为最长序列ID
seqkit grep -p "ID" final.contigs.fa >longest.fa  #将ID替换为最长序列ID

#修改名字
mv longest.fa ncov.fa

#与已发表NC_045512进行比对
seqkit stat NC_045512.fa ncov.fa
lastz NC_045512.fa ncov.fa  --out=lastz.axt

#PCR产物验证
#下载PCR扩增测序数据
prefetch SRR14867660 -O ./
locate SRR14867660.sra
fasterq-dump SRR14867660.sra 

#与拼接结果比对
bwa-mem2 index ncov.fa
bwa-mem2 mem -t 12 ncov.fa SRR14867660_1.fastq SRR14867660_2.fastq >pcr.sam
#samtools处理比对结果
samtools sort -@ 12 -o pcr.sorted.bam pcr.sam 
samtools index pcr.sorted.bam

#对ncov.fa建立索引
samtools faidx ncov.fa
#tablet可视化
#将文件拷贝至windows下使用tablet可视化
mkdir tablet
mv ncov.fa ncov.fa.fai pcr.sorted.bam pcr.sorted.bam.bai tablet

#############################
#    31 拼接新冠病毒基因组    #
#############################
mkdir 31.assembly;cd 31.assembly;
#下载参考序列
efetch -db nuccore -format fasta -id MN908947.3 > MN908947.3.fa
#https://www.ncbi.nlm.nih.gov/sars-cov-2/
#下载illumina测序数据
/ifs1/Software/bin/prefetch SRR14867660 -O ./
fasterq-dump SRR14867660.sra 

#下载引物bed文件
git clone https://github.com/artic-network/artic-ncov2019.git

#比对
#illumina测序数据比对，排序一步完成
bwa-mem2 index MN908947.3.fa 
bwa-mem2 mem -t 12 MN908947.3.fa SRR14867660.sra_1.fastq SRR14867660.sra_2.fastq | samtools sort  | samtools view -F 4 -o ncov.sorted.bam

#切出引物部分
ivar trim -e -i ncov.sorted.bam -b nCoV-2019.bed -p ncov.primertrim
#排序建立索引
samtools sort -@ 12 -o ncov.primertrim.sorted.bam  ncov.primertrim.bam
samtools index ncov.primertrim.sorted.bam
#统计覆盖度
samtools coverage ncov.primertrim.sorted.bam -o ncov.samcov.txt -m

#ivar
samtools mpileup -A -d 1000 -B -Q 0 --reference MN908947.3.fa ncov.primertrim.sorted.bam | ivar consensus -p ivar_consensus -n N

#统计长度
seqkit stat MN908947.3.fa ivar_consensus.fa 

#比较二者差别
blat MN908947.3.fa ivar_consensus.fa -out=axt blat.out
blat MN908947.3.fa ivar_consensus.fa -out=blast blat.out

#############################
#      Artic Network        #
#############################
# 创建工作目录
mkdir artic;cd artic;
#激活工作环境
conda activate /ifs1/Software/miniconda3/envs/artic
# 下载数据
#下载nanopore测序数据
/ifs1/Software/bin/prefetch SRR14800265 -O ./
fasterq-dump SRR14800265.sra

# 数据过滤
artic guppyplex --min-length 400 --max-length 700 --directory ../artic/ --prefix clean 
# 运行流程
#medaka模式
cp -r /ifs1/Software/biosoft/fieldbioinformatics/test-data/primer-schemes/ ./
artic minion --normalise 200 --threads 12 --scheme-directory primer_schemes/ --read-file clean_test.fastq nCoV-2019/V3 out  --medaka --medaka-model r941_min_high_g351

#############################
#      纳米孔测序拼接基因组   #
#############################
#1 下载数据
mkdir ebov;cd ebov;
wget http://artic.s3.climb.ac.uk/run-folders/EBOV_Amplicons_flongle.tar.gz
tar -zxvf EBOV_Amplicons_flongle.tar.gz 

#拷贝配置文件
cp -r /ifs1/Software/biosoft/fieldbioinformatics/test-data/primer-schemes/ ./

#2 basecalling以及demultiplex
#激活环境
conda activate /ifs1/Software/miniconda3/envs/artic
#basecalling
artic gather --min-length 400 --max-length 800 --prefix ebov-mayinga --directory ./20190830_1509_MN22126_AAQ411_9efc5448 --fast5-directory ./20190830_1509_MN22126_AAQ411_9efc5448/fast5_pass
#拆分barcode
artic demultiplex --threads 2 ebov-mayinga_fastq_pass.fastq
#3 获得基因组
artic minion --normalise 200 --threads 12 --scheme-directory primer-schemes/ --read-file ebov-mayinga_fastq_pass-NB03.fastq --fast5-directory ./20190830_1509_MN22126_AAQ411_9efc5448/fast5_pass --sequencing-summary ./20190830_1509_MN22126_AAQ411_9efc5448/lab-on-an-ssd_20190830_160932_AAQ411_minion_sequencing_run_EBOV_Amplicons_flongle_sequencing_summary.txt IturiEBOV/V1 ebov

#############################
#     32 构建系统发育树       #
#############################
#下载测序数据
mkdir 32.tree;cd 32.tree;
efetch -db nuccore -format fasta -id "MN908947 MG772933 KT444582 MZ310552 MZ202314 MZ169911 MZ318159 MZ373479 MZ169912 MW852494 MZ257684 MZ310903 MZ310580">ncov.fasta
muscle -in ncov.fasta -out ncov.clw -clw
treebest nj ncov.clw >ncov.nwk

#利用mega构建系统发育树
#megacc比对
megacc -a muscle_align_nucleotide.mao -d ncov.fasta -o mega/ncov -f MEGA   
#mega画树
megacc -a infer_UPGMA_nucleotide.mao -d mega/ncov.meg -o mega/ncov

#上传至itol可视化
#https://itol.embl.de/


#############################
#      33 变种病毒识别       #
#############################
mkdir 33.pangolin;cd 33.pangolin
#下载参考序列
efetch -db nuccore -format fasta -id MN908947.3 > MN908947.3.fa
#下载待分析样本测序数据
/ifs1/Software/bin/prefetch ERR6057099 -O ./
fasterq-dump ERR6057099.sra 

#获得基因组序列
#比对
#illumina测序数据比对，排序一步完成
bwa-mem2 index MN908947.3.fa 
bwa-mem2 mem -t 12 MN908947.3.fa ERR6057099.sra_1.fastq ERR6057099.sra_2.fastq | samtools sort | samtools view -F 4 -o ncov.sorted.bam
#切出引物部分
ivar trim -e -i ncov.sorted.bam -b nCoV-2019.bed -p ncov.primertrim
#排序建立索引
samtools sort -@ 12 -o ncov.primertrim.sorted.bam  ncov.primertrim.bam
samtools index ncov.primertrim.sorted.bam
#ivar
samtools mpileup -A -d 1000 -B -Q 0 --reference MN908947.3.fa ncov.primertrim.sorted.bam | ivar consensus -p ERR6057099 -n N

#序列比对，blat比对，输出blast列表结果
blat /ifs1/VipData/04.ncov/nCov.fasta ERR6057099.fa -out=blast8 blat.out
#对结果进行过滤，筛选，分组计算平均值
cat blat.out |awk '{if ($4 >= 10000) print $2"\t"$3"\t"$4}' | datamash --sort groupby 1 mean 2 |sort -nr -k2


#==================Pango Lineages===============
#运行软件
conda activate /ifs1/Software/miniconda3/envs/pangolin
pangolin ERR6057099.fa --alignment --outfile pangolin.csv

#在线分析
#https://pangolin.cog-uk.io/