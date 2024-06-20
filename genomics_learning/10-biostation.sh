########################################################
#                 六十五 环境配置                       #
########################################################
#=========================
#       1 服务器设置     #
#=========================
#1 默认文本界面
systemctl get-default
graphical.target #显示的是graphical.target，也就是图形模式，多用户模式，multi-user.target
systemctl set-default multi-user.target

#升级到centos-stream
dnf swap centos-linux-repos centos-stream-repos 
dnf distro-sync

cat /etc/centos-release

#3 安装epel源
yum install -y epel-release.noarch
yum clean metadata
yum makecache

#4 安装ntfs支持ntfs格式磁盘
yum install -y ntfs-3g.x86_64 ntfs-3g-devel.x86_64 ntfsprogs.x86_64 ntfs-3g-system-compression.x86_64

#5 安装openssh-server远程登录
yum install -y openssh-server.x86_64
systemctl start sshd.service

#=========================
#       2 安装基本配置     #
#=========================
yum install -y gcc* 
yum install -y zlib* 
yum install -y glibc* 
yum install -y boost-* --skip-broken
yum install -y lib* --skip-broken
yum install -y compat-* 
yum install -y java-*  --skip-broken
yum install -y pip* 
yum install -y python2-pip.noarch python34-pip.noarch
yum install -y cmake3.x86_64 cmake.x86_64
yum install -y ncurses-devel.i686 ncurses-devel.x86_64 
yum install -y root.x86_64
yum install -y libcurl-devel.i686 libcurl-devel.x86_64 libcurl.i686 libcurl.x86_64
yum install -y libXScrnSaver.i686 libXScrnSaver.x86_64
yum  install -y argtable.x86_64  argtable-devel.x86_64
yum install -y openssl-devel
yum install -y libstdc++.i686 libstdc++.x86_64 libstdc++-devel.i686 libstdc++-devel.x86_64 libstdc++-static.i686 libstdc++-static.x86_64 compat-libstdc++-33.i686 compat-libstdc++-33.x86_64
yum install -y tbb-devel.x86_64
yum install -y gsl-devel.i686   gsl-devel.x86_64 gsl.i686 gsl.x86_64
yum install -y perl-Sys-SigAction.noarch
yum install -y build-essentials
yum install -y swig-doc.noarch
yum install -y swig.x86_64
yum install -y libcurl-devel.i686
yum install -y libcurl-devel.x86_64
yum install -y libcurl.i686
yum install -y libcurl.x86_64
yum install -y glibc-devel.i686
yum install -y glibc-devel.x86_64
yum install -y libXtst-devel.i686
yum install -y libXtst-devel.x86_64
yum install -y openssl-devel.i686
yum install -y openssl-devel.x86_64
yum install -y xmlsec1.x86_64
yum install -y xmlsec1.i686
yum install -y xmlsec1-nss.x86_64 
yum install -y xmlsec1-nss.i686
yum install -y xmlsec1-openssl.x86_64
yum install -y xmlsec1-openssl.i686
yum install -y xmlsec1-openssl-devel.x86_64
yum install -y libquadmath.i686
yum install -y libquadmath.x86_64
yum install -y libquadmath-devel.i686
yum install -y libquadmath-devel.x86_64
yum install -y libstdc++-static
yum install -y pandoc.x86_64
yum install -y pcre2.x86_64
yum install -y pcre.x86_64
yum install -y pcre2.x86_64
yum install -y ImageMagick* --skip-broken
#=========================
#       安装perl模块  #
#=========================
yum install -y perl-* --skip-broken
#=========================
#       安装python模块   #
#=========================
yum install -y python3-* --skip-broken
#=========================
#       其他应用   #
#=========================
yum install -y dnf
yum install -y git 
yum install -y tree 
yum install -y htop 
yum install -y screen
yum install -y lftp
yum install -y iftop
yum install -y tmux
yum install -y podman-docker.noarch
yum install -y pbzip2 nethogs

########################################################
#                   六十六 用户管理                   #
########################################################
#创建名为tests1的用户
useradd test1 -d /ifs1/User/test1 
echo "Pass1234" | passwd --stdin test1
su - test1
chown -R test1:test1 /ifs1/User/test1;
chmod -R 700 /ifs1/User/test1

#创建用户组
groupadd vip

#将用户添加到docker组中
usermod -aG docker ${Yonghu}

#修改属主权限
chown -R test1:vip /ifs1/User/test1
#修改文件权限
chmod -R 700 /ifs1/User/test1

########################################################
#                六十七 安装生物软件和配置数据库          #
########################################################
#=========================
#       1 安装bioconda   #
#=========================
#下载安装：
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
sh Miniconda3-latest-Linux-x86_64.sh  
#刷新环境
source ~/.bashrc

#运行以下命令，添加软件源：
conda config --add channels bioconda 
conda config --add channels conda-forge

#4 安装mamba
conda install -y mamba

#=========================
#       2 安装生物软件    #
#=========================
#安装mamba
conda install -y mamba 
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
mamba create -n nanoplot -y nanoplot

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

#=========================
#       3 切换bioconda    #
#=========================
conda activate /ifs1/Software/miniconda3/
#查看bwa路径
(/ifs1/Software/miniconda3) meta 07:44:36 ~
$ which bwa
/ifs1/Software/miniconda3/bin/bwa
#激活nanoplot环境
(/ifs1/Software/miniconda3) meta 07:45:09 ~
$ conda activate /ifs1/Software/miniconda3/envs/nanoplot
(nanoplot) meta 07:45:19 ~
$ which NanoPlot 
/ifs1/Software/miniconda3/envs/nanoplot/bin/NanoPlot
(nanoplot) meta 07:45:22 ~

#=========================
#       4 管理生物数据库   #
#=========================
mkdir -p /ifs1/Database

########################################################
#                   六十八 配置R环境                           #
########################################################
#=========================
#       1 安装R          #
#=========================
#yum搜索R
yum search R | grep "^R" 
#查看R版本
yum info R.x86_64
yum install -y R.x86_64 --skip-broken

#利用bioconda安装R
mamba install -y r-base

#=========================
#    2 安装rstudio-server   #
#=========================
#安装rstudio-server
wget https://download2.rstudio.org/server/centos8/x86_64/rstudio-server-rhel-2021.09.2-382-x86_64.rpm
yum install rstudio-server-rhel-2021.09.2-382-x86_64.rpm

#检查是否安装成功
rstudio-server verify-installation
rstudio-server status

#修改rstudio-server配置
# vim /etc/rstudio/rserver.conf
# Server Configuration File
#rsession-which-r=/usr/local/bin/R
rsession-which-r=/ifs1/Software/miniconda3/bin/R

#=========================
#     3 修改防火墙        #
#=========================
systemctl start rstudio-server.service
systemctl enable rstudio-server.service
systemctl status rstudio-server
firewall-cmd --permanent --add-port=8787/tcp
firewall-cmd --permanent --add-port=8787/udp
firewall-cmd --reload

#=========================
#      4  配置R包         #
#=========================
mamba search deseq2
mamba install -y bioconductor-deseq2

########################################################
#                     六十九 安装网络应用                #
########################################################
#=========================
#       1  jupyterlab     #
#=========================
# 1打开端口
#打开8888窗口
firewall-cmd --permanent --add-port=8888/tcp
firewall-cmd --permanent --add-port=8888/udp
firewall-cmd --reload
#2 安装 jupyterlab
#juypterlab
mamba install -c conda-forge -y jupyterlab
mamba install -c conda-forge -y voila
#安装中文插件
mamba install -y jupyterlab-language-pack-zh-CN
#设置密码
jupyter notebook password
#生成配置文件
jupyter notebook --generate-config
#修改配置文件
vim /root/.jupyter/jupyter_notebook_config.py
c.ServerApp.port = 8888  
c.ServerApp.allow_remote_access = True
c.ServerApp.ip='*'  
c.ServerApp.open_browser = False
#3 启动 jupyterlab
nohup jupyter lab --allow-root > jupyter.log 2>&1 &

#=========================
#          2 vscode      #
#=========================
#安装podman
yum install -y podman
#下载
podman pull  docker.io/jmcdice/vscode-server:latest
#开通8080端口
firewall-cmd --add-port=8080/tcp --zone=public --permanent
firewall-cmd --reload
#运行
podman run -d -p 8080:8080 --restart=always --name=vscodeserver jmcdice/vscode-server

#=========================
#       3 在线blast   #
#=========================
#1 安装网络环境
yum install -y httpd httpd-devel
systemctl start  httpd
systemctl enable  httpd
yum install --skip-broken -y php*

#2 开启80端口
firewall-cmd --permanent --zone=public --add-service=http
firewall-cmd --permanent --zone=public --add-service=https
firewall-cmd --reload

#3 安装viroBlast
# 下载 viroblast
git clone  https://github.com/MullinsLab/ViroblastStandalone.git

#4 将整个目录拷贝到网络浏览目录下
mv viroblast /var/www/html/

#3 构建 blast 数据库
cd /var/www/html/viroblast/db

#构建nt库
cd nucleotide/
ln -s /ifs1/Database/nt/nt ./ 
makeblastdb -in nt -out nt -dbtype nucl -parse_seqids

#构建nr库
cd protein/
ln -s /ifs1/Database/nr/nr ./ 
makeblastdb -in nr -out nr -dbtype prot -parse_seqids

#修改配置文件
vim viroblast.ini
# 配置 blast+ 路径
blast+:/ifs1/Software/miniconda3/bin/
#配置数据库
blastn: nt => Nucleotide NT database
blastp: nr => Protein NR database
blastx: nr => Protein NR database
tblastn: nt => Nucleotide NT database
tblastx: nt => Nucleotide NT database

#5 通过浏览器访问IP地址测试


########################################################
#                七十 ubuntu系统配置                        #
########################################################
#创建root账户
sudo passwd root
#输入当前用户密码：
#输入root密码：
#再次输入root密码：

#以下操作使用root账户完成
su -

#apt安装生物软件
apt install -y bwa
apt install -y samtools 
apt install -y bcftools
apt install -y blast2
apt install -y bedtools
apt install -y seqtk
apt install -y minimap2
apt install -y bowtie2
apt install -y phylip
apt install -y clustalx
apt install -y canu
apt install -y kraken2
apt install -y hisat2
apt install -y stringtie
apt install -y jellyfish
apt install -y circos
apt install -y nanopolish
apt install -y nanook
apt install -y centrifuge
apt install -y rna-star
apt install -y freebayes
apt install -y cnvkit
apt install -y spades
apt install -y mothur
apt install -y muscle
apt install -y mafft
apt install -y iqtree
apt install -y sniffles
apt install -y last-align
apt install -y augustus
apt install -y bamtools 
apt install -y bedops
apt install -y delly
