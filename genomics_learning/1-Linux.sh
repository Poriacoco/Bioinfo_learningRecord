#############################
#        基础命令       #
#############################
#1 修改密码
passwd
#输入两次密码，不少于八位，大小写数字组合
#2 退出登录
exit
#3 显示当前目录
pwd
#4 显示日期
date
#5 显示日历
cal

#6 清空屏幕
clear
#7 显示历史记录
history
# 其他基础命令
id
whoami
last
clear
uptime
ifconfig

#############################
#           四快捷键       #
#############################
鼠标	#鼠标左键选中，右键粘贴
Tab	    #Tab键补齐
方向键	#上下方向键查找历史记录
Ctrl+A	#光标移动到行首
Ctrl+E	#光标移动到行尾
Ctrl+C	#终止当前操作
Ctrl+D	#关闭进程
Ctrl+L	#清空屏幕
Ctrl+R	#搜索历史记录
Ctrl+W	#剪切光标所在处之前的一个词
Ctrl+Z	#暂停进程

#############################
#       特殊符号             #
#############################
~   #家目录
#   #Linux文件中注释行，表示不起作用。
$   #文件行结尾标识符，变量标识符。
^   #文件行首标识符
&   #任务放到后台
*   #通配符，代表一个字符或者很多个字符。
\   #用来转义，\t表示制表符，\n表示换行符。
<   #数据流的流入方向，表示输入，将数据传入给左侧软件。
|   #管道，改变数据流的方向，将数据传入给另外的软件。
>   #数据流的流出方向，表示输出，将屏幕输出的内容写入一个文件。
2>  #数据流的流出的第二个方向，表示错误输出，报错信息会写入到这个文件中。
>>  #表示追加写入
/   #根目录，目录分隔符
-   #短选项标识符-h
--  #长选项标识符--help

#############################
#       目录结构    #
#############################
#8 切换目录
cd  	#切换到家目录
cd /	#切换到根目录
cd ~	#切换到用户个人目录
cd ~;ls	#切换到用户个人目录并显示目录下文件
cd -	#切换到上次使用目录
cd .	#使用相对路径
cd ..	#使用相对路径，回到上层

#9 显示文件内容
ls      #显示文件和目录名字
ll      #按列显示文件和目录名字

#绝对目录与相对目录
/ifs1/MetaDatabase/genome/human.fa 
../../MetaDatabase/genome/human.fa

#检查文件是否存在
#绝对目录 
ll /ifs1/MetaDatabase/ genome/human.fa 
#相对目录 
ll ../../MetaDatabase/genome/human.fa

#############################
#       文件（夹）操作             #
#############################
#10 创建目录（文件夹）
mkdir 0723
#11 拷贝文件或者目录
cp /ifs1/TestDatas/coronavirus/amplicon/nanopore/SRR14859525.sra ./
cp -r /ifs1/TestDatas/coronavirus/amplicon/nanopore/ ./

#链接文件
ln -s /ifs1/TestDatas/nanopore7/ ./

#12 剪切粘贴文件或者目录
mv SRR14859525.sra nanopore.sra
mv nanopore.sra 0723

#13 重命名文件或者目录
mv 0723 0724

#14 删除文件或者目录
rm  SRR14859525.sra
rm -f SRR14859525.sra
rm -rf 0724 


#15 查看文件
less
cat
more
head
tail

#16 数据流方向
<
>
|
>>
1>
2>

#17 打包压缩
gzip ref.fna
gunzip ref.fna.gz
cp /ifs1/Software/seqkit_linux_amd64.tar.gz ./
tar -zxvf seqkit_linux_amd64.tar.gz
tar -zxvf 0723.tar.gz seqkit_linux_amd64 ref.fna

#18 其他文件操作

#查看文件小
du -sh seqkit_linux_amd64.tar.gz
#查看当前磁盘大小
df ./

#############################
#       编写脚本            #
#############################
#18 文件编辑
vim scripts.sh 
#i a u 切换为插入模式
#ESC切换为命令模式

#如何退出 vim？ 
#首先按 esc 键切换到命令模式 然后按“shift+：”冒号表示可以输入命令了 然后按 
q！#不保存退出 
wq #保持退出或者 x 保存退 
w+文件名 #另存为一个文件

#19 运行脚本
sh prodigal.sh


#############################
#       权限管理             #
#############################
#421模式修改
chmod 644 a1.index.sh

#ugo模式修改  
chmod ug+x a1.index.sh

#给一个软件赋予可执行权限 
chmod u+x bwa 
#文件只对个人可见
chmod -R 700 *

#############################
#       进程管理        #
#############################
#查看进程
top  #press "q" to exit 查看系统运行状态   
top -b  # press "Ctrl +C" to exit   
top -u $USER  #查看自己进程状态  
htop #htop查看系统运行状态  

#显示进程信息，包括无终端的
ps -aux 
#显示所有进程信息，连同命令行
ps -ef 
#根据 CPU 使用来升序排序
ps -aux --sort -pcpu | less  
#根据内存使用 来升序排序
ps -aux --sort -pmem | less 
#消耗CPU和内存前十名用户
ps -aux --sort -pcpu,+pmem | head -n 10 

# kill杀死进程
kill -9 "process number"  #杀死进程   

#休眠
sleep 60       #休眠60秒   

#后台程序 forehead   
fg 

#前台程序 background
bg  
#查看后台进程   
jobs

#不挂起运行程序，关闭登录窗口后程序继续运行 
nohup  sh scripts.sh &

#忘记使用nohup之后，将后台任务转换为nohup
disown 

#screen不间断进程
#1  新建会话，命名为wget
screen -S  wget

# 2 运行命令
wget -c ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/swissprot.gz

# 3 按ctrl+a，然后字母d，退出会话，任务仍在运行
$ screen -S  wget
[detached from 283349.wget]

# 4 screen -ls查看任务

#5 重新进入wget终端，任务正在运行
$ screen -r wget

# 6 关闭会话任务,如果在会话中使用exit，就会在退出会话，也关闭了该会话，或者按ctrl+a，k
# screen ls查看会话
$ screen -r wget


#############################
#       环境设置             #
#############################

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

#4 更新设置
#首先备份一下
cp ~/.bashrc ~/.bashrc.bak
#打开vim修改
vim ~/.bashrc
#保存退出，刷新设置
source ~/.bashrc

#5 恢复初始配置
#恢复备份
cp ~/.bashrc.bak ~/.bashrc

#也可以恢复到系统初始化
cp /etc/skel/.bashrc ~/.bashrc
#选择覆盖原始文件
