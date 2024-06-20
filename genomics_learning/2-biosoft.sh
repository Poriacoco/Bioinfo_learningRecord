
#创建文件夹
mkdir bin biosoft src

#拷贝数据
cp -r /ifs1/VipData/02.biosoft/biosoft/* ~/src

#############################
#       编译软件             #
#############################
#R软件
#安装R语言依赖
yum install -y --skip-broken zlib java gcc-gfortran gcc gcc-c++ readline-devel libXt-devel bzip2-devel.x86_64 bzip2-libs.x86_64 xz-devel.x86_64 pcre-devel.x86_64 libcurl-devel.x86_64

#下载
wget https://cloud.r-project.org/src/base/R-4/R-4.1.1.tar.gz
tar -zxvf R-4.1.1.tar.gz -C ~/biosoft
cd ~/biosoft/ R-4.1.1
#检测配置
./configure 
#编译
make
make install


#############################
#        安装预编译软件      #
#############################

#1 blast+
wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.12.0+-x64-linux.tar.gz
tar -zxvf ncbi-blast-2.12.0+-x64-linux.tar.gz
cd ~/biosoft/ncbi-blast-2.10.1+/bin
ls -1  | while read i;do ln -s $PWD/$i ~/bin/;done;

#2 edirect
wget  https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/edirect.tar.gz
tar -zxvf edirect.tar.gz 
cd ~/bin/
ln -s ~/biosoft/edirect/efetch .
ln -s ~/biosoft/edirect/edirect.pl .

#3 sratookit
#下载指定版本
wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.10.8/sratoolkit.current-centos_linux64.tar.gz
tar -zxvf sratoolkit.current-centos_linux64.tar.gz
cd ~/bin
ln -s  ~/biosoft/sratoolkit.2.11.1-centos_linux64/bin/prefetch ./
ln -s  ~/biosoft/sratoolkit.2.11.1-centos_linux64/bin/fasterq-dump ./
ln -s  ~/biosoft/sratoolkit.2.11.1-centos_linux64/bin/fastq-dump ./

#4 aspera
wget https://download.asperasoft.com/download/sw/connect/3.9.9/ibm-aspera-connect-3.9.9.177872-linux-g2.12-64.tar.gz
tar -zxvf ibm-aspera-connect-3.9.9.177872-linux-g2.12-64.tar.gz
sh ibm-aspera-connect-3.9.9.177872-linux-g2.12-64.sh

#5 bowtie2
wget https://jaist.dl.sourceforge.net/project/bowtie-bio/bowtie2/2.3.5.1/bowtie2-2.3.5.1-linux-x86_64.zip
unzip bowtie2-2.3.5.1-linux-x86_64.zip 
ln -s ~/biosoft/bowtie2-2.3.5.1-linux-x86_64/bowtie2-build ~/bin
ln -s ~/biosoft/bowtie2-2.3.5.1-linux-x86_64/bowtie2 ~/bin

#6 fastqc
wget http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip
unzip fastqc_v0.11.9.zip 
chmod u+x  ~/biosoft/FastQC/fastqc 
ln -s ~/biosoft/FastQC/fastqc ~/bin

#7 fastp
#https://github.com/OpenGene/fastp#get-fastp
cd ~/bin
​wget http://opengene.org/fastp/fastp
chmod a+x ./fastp

#8 spades
#https://github.com/ablab/spades#sec2.1
wget http://cab.spbu.ru/files/release3.14.1/SPAdes-3.14.1-Linux.tar.gz
tar -zxvf SPAdes-3.14.1-Linux.tar.gz
ln -s ~/biosoft/SPAdes-3.14.1-Linux/bin/spades.py ~/bin

#9 haist2
#https://daehwankimlab.github.io/hisat2/
curl https://cloud.biohpc.swmed.edu/index.php/s/hisat2-220-Linux_x86_64/download -o hisat2-2.2.0-Linux_x86_64.zip
unzip hisat2-2.2.0-Linux_x86_64.zip -d ~/biosoft
ln -s ~/biosoft/hisat2-2.2.0/hisat2-build ~/bin/
ln -s ~/biosoft/hisat2-2.2.0/hisat2 ~/bin/

#10 stringtie
#http://ccb.jhu.edu/software/stringtie/index.shtml
wget http://ccb.jhu.edu/software/stringtie/dl/stringtie-2.1.3b.Linux_x86_64.tar.gz
tar -zxvf stringtie-2.1.3b.Linux_x86_64.tar.gz/
ln -s ~/biosoft/stringtie-2.1.3b.Linux_x86_64/stringtie ~/bin/

#11 diamond
#http://www.diamondsearch.org/index.php
wget http://github.com/bbuchfink/diamond/releases/download/v0.9.36/diamond-linux64.tar.gz
tar -zxvf diamond-linux64.tar.gz  

#12 vsearch
#https://github.com/torognes/vsearch/
wget https://github.com/torognes/vsearch/releases/download/v2.15.0/vsearch-2.15.0-linux-x86_64.tar.gz
tar -zxvf vsearch-2.15.0-linux-x86_64.tar.gz/
ln -s ~/biosoft/vsearch-2.15.0-linux-x86_64/bin/vsearch ~/bin

#13 mothur
# https://github.com/mothur/mothur
wget https://github.com/mothur/mothur/releases/download/v.1.44.1/Mothur.linux.zip
unzip Mothur.linux.zip -d ~/biosoft/
ln -s ~/biosoft/mothur/mothur ~/bin/

#14 muscle
#http://www.drive5.com/muscle/
wget http://www.drive5.com/muscle/downloads3.8.31/muscle3.8.31_i86linux64.tar.gz
tar -zxvf muscle3.8.31_i86linux64.tar.gz -C ~/bin

#15 gatk4
#https://github.com/broadinstitute/gatk/
wget https://github.com/broadinstitute/gatk/releases/download/4.1.8.0/gatk-4.1.8.0.zip
unzip gatk-4.1.8.0.zip -d ~/biosoft/
ln -s ~/biosoft/gatk-4.1.8.0/gatk ~/bin/

#16 安装guppy，需要官网购买仪器注册用户下载
wget https://mirror.oxfordnanoportal.com/software/analysis/ont-guppy-cpu_5.0.11_linux64.tar.gz
tar -zxvf ont-guppy-cpu_5.0.11_linux64.tar.gz
#自行编译软件

#17 安装centrifuge
wget ftp://ftp.ccb.jhu.edu/pub/infphilo/centrifuge/downloads/centrifuge-1.0.3-beta-Linux_x86_64.zip
unzip centrifuge-1.0.3-beta-Linux_x86_64.zip 


#############################
#      自行编译             #
#############################
#拷贝数据
cp -r /ifs1/VipData/02.biosoft/git/* ~/biosoft
#1 bwa
cd ~/biosoft
git clone https://github.com/lh3/bwa.git
cd bwa; make

#2 minimap2
git clone https://github.com/lh3/minimap2
cd minimap2 && make

#3 prodigal
git clone https://github.com/hyattpd/Prodigal.git
cd Prodigal
make install 

#14 安装filtlong
git clone https://github.com/rrwick/Filtlong.git
cd Filtlong
make -j

#5 flye
git clone https://github.com/fenderglass/Flye
cd Flye
make

#6 mummer4
wget https://github.com/mummer4/mummer/releases/download/v4.0.0rc1/mummer-4.0.0rc1.tar.gz
tar -zxvf mummer-4.0.0rc1.tar.gz
cd tar -zxvf mummer-4.0.0rc1
./configure
make
make install


#7 htslib
#包括samtools和bcftools在内的很多软件都需要依赖htslib，所以需要先安装htslib
cd ~biosoft
git clone https://github.com/samtools/htslib.git
cd htslib
autoreconf -i
git submodule update --init --recursive 
./configure   
make
make install

#8 安装samtools
cd ~/biosoft
git clone https://github.com/samtools/samtools.git
cd samtools
autoheader     
autoconf -Wno-syntax
./configure          
make

#9 安装bcftools
cd ~/biosoft
git clone https://github.com/samtools/bcftools.git
cd bcftools
autoheader 
autoconf 
./configure --enable-libgsl --enable-perl-filters
make

#10 bedtools
wget https://github.com/arq5x/bedtools2/archive/v2.29.2.tar.gz
tar -zxvf v2.29.2.tar.gz 
cd bedtools2-2.29.2/
make

#11 deeptools
git clone https://github.com/deeptools/deepTools.git
python setup.py install --prefix ~/biosoft/deepTools2.0

#12 Augustus
# https://github.com/Gaius-Augustus/Augustus
#安装依赖
yum install -y libboost-iostreams-dev zlib1g-dev libgsl-dev  libboost-all-dev libsuitesparse-dev liblpsolve55-dev libsqlite3-dev libmysql++-dev libbamtools-dev libboost-all-dev  libboost-all-dev
git clone https://github.com/Gaius-Augustus/Augustus.git
cd Augustus
make

#13  sniffles
wget https://github.com/fritzsedlazeck/Sniffles/archive/master.tar.gz -O Sniffles.tar.gz
tar xzvf Sniffles.tar.gz
cd Sniffles-master/
mkdir -p build/
cd build/
cmake ..
make

#14 rocon网址
git clone --recursive https://github.com/lbcb-sci/racon.git 
cd racon
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
#############################
#       环境设置             #
#############################
#首先备份一下
cp ~/.bashrc ~/.bashrc.bak
#打开vim修改
vim ~/.bashrc
#保存退出，刷新设置
source ~/.bashrc

#设置快捷命令
alias rm='rm -i'
alias cp='cp -i'
alias mv='mv -i'
alias df='df -h'
alias du='du -skh'
alias grep='grep --color'
alias ls='ls -hF --color=tty'                 # classify files in colour
alias dir='ls --color=auto --format=vertical'
alias ll='ls -lh -rt --file-type'                              # long list
alias l='ls -CF'                              #
alias lla='ls -a -l'
alias e='exit'
alias le='less -S'
alias d='display'
alias t='top -u $USER'

#2 修改命令行显示
if [ -f /etc/bashrc ]; then
    . /etc/bashrc
fi

if [ "$TERM" = "xterm" ]
then
 export PS1="\[\e[31;1m\]\u\[\e[0m\] \[\e[32;1m\]\t \[\e[0m\]\[\e[34;1m\]\w\[\e[0m\]\n\[\e[31;1m\]$ \[\e[0m\]"
else
 export PS2="\[\e[31;1m\]\u\[\e[0m\] \[\e[32;1m\]\t \[\e[0m\]\[\e[34;1m\]\w\[\e[0m\]\n\[\e[31;1m\]$ \[\e[0m\]"
fi


#3 export环境变量
export PATH="$PATH:./:/usr/bin:$PATH"
export PATH="$PATH:/ifs1/Software/bin/:$PATH"
export PERL5LIB="/ifs1/Software/biosoft/tRNAscan-SE-1.3.1/"

# added by Miniconda3 installer
export LD_LIBRARY_PATH="/ifs1/Software/boost-1.60.0-py27_3/lib/"

#4 恢复初始配置
#恢复备份
cp ~/.bashrc.bak ~/.bashrc

#也可以恢复到系统初始化
cp /etc/skel/.bashrc ~/.bashrc
#选择覆盖原始文件


#############################
#      Bioconda             #
#############################
#软件官网
http://bioconda.github.io/
#软件搜索
https://anaconda.org/bioconda/repo

#1下载biconda
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh  
#2安装
sh Miniconda3-latest-Linux-x86_64.sh  
source ~/.bashrc
#3添加软件源
conda config --add channels bioconda 
conda config --add channels conda-forge

#安装软件(bwa，samtools，bcftools软件为例)
#搜索软件    
conda search bwa    
#安装软件    
conda install -y samtools=1.9 bcftools=1.9    
#查看已安装软件  
conda list  
#升级软件  
conda update bwa  
#移除软件  
conda remove bwa 

#############################
# 使用bioconda安装常用软件   #
#############################
#安装mamba
conda install -c conda-forge -y mamba
#利用mamba安装软件
mamba install -y bwa 
mamba install -y samtools
mamba install -y bcftools
mamba install -y blast 
mamba install -y blat 
mamba install -y mummer 
mamba install -y mafft 
mamba install -y muscle 
mamba install -y lastz
mamba install -y sra-tools
mamba install -y seqkit
mamba install -y seqtk
mamba install -y bedtools
mamba install -y bedops
mamba install -y gfatools
mamba install -y circos
mamba install -y entrez-direct
mamba install -y emboss

#安装数据质控软件
mamba install -y fastqc multiqc 
mamba install -y trimmomatic
mamba install -y fastp

#安装基因组拼接相关工具
mamba install -y velvet
mamba install -y flye
mamba install -y miniasm
mamba install -y canu
mamba install -y megahit
mamba install -y spades
mamba install -y quast
mamba install -y racon
mamba install -y miniasm
mamba install -y nanopolish

#安装基因功能分析软件
mamba install -y prodigal
mamba install -y glimmer
mamba install -y augustus
mamba install -y trf

#############################
#      虚拟环境             #
#############################

#1 安装不同版本软件
#查看虚拟环境
conda env list
#安装指定版本软件blast 2.7.1，samtools 1.7
conda create -y -n test 
#激活虚拟环境
conda activate test
#安装软件
conda install -y blast=2.7.1 samtools=1.7
#退出虚拟环境
conda deactivate

#2 创建python 2.7环境
conda create -n py27 -y python=2.7
#查看现有虚拟环境
conda env list
#激活python2.7环境
conda activate py27
#查看python版本
python -V

#2 安装nanoplot
conda create -n nanoplot -y nanoplot

#3 拼接结果纠错
#medaka网址：https://github.com/nanoporetech/medaka
conda create -y -n medaka -c conda-forge -c bioconda medaka


#artic network
git clone https://github.com/artic-network/artic-ncov2019.git
cd artic-ncov2019/
conda env create -f environment.yml

#pangolin
git clone https://github.com/cov-lineages/pangolin.git
cd pangolin
conda env create -f environment.yml
conda activate pangolin
pip install .


#安装prokka，建立虚拟环境
conda create -n prokka  -y
conda activate prokka
conda install -y prokka

#############################
#      下载数据             #
#############################
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

