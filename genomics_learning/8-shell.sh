
########################################################
#                五十五 管理数据流                      #
########################################################
mkdir 55.pipe
#重定向
zcat ../data/demo.fastq.gz >demo.fastq

#追加写入
ps -aux >>ps.txt

#写入黑洞文件
ll ../data/demo.sam 1>/dev/null 

#使用错误输出
ll ../data/demo.sam 1>/dev/null 2>stderr.txt

#使用管道
ps -aux | grep root >root.ps

#=========================
#          xargs         #
#=========================
#1 将内容拆成多列显示
cat ../data/xargs.txt | xargs
cat ../data/xargs.txt | xargs -n 3

#2分割符
cat /etc/passwd | xargs -d ":"

#将一行内容分成三列
echo {a..z} | xargs -n 3

#4 将目录下全部fa找出来，拷贝到当前目录下
find ../corona/ -name "*.fa" | xargs -t -I{} cp {} ./

#搜索名字为sleep的进程，然后kill掉
sleep 20 &
ps -u $USER | awk '/sleep/ {print $1}' | xargs echo kill

########################################################
#                   五十六 正则表达式                    #
########################################################





########################################################
#                   五十七 文件搜索                     #
########################################################
#=========================
#           1 locate     #
#=========================
#查找文件名
locate SRX5299464

#使用正则表达式
locate -r "samtools$"

#=========================
#          2 whereis     #
#=========================
#whereis 搜索路径
whereis -l

#查找程序文件
whereis less
whereis -b less

#=========================
#          3 which       #
#=========================
#查找软件
which -a bwa

#显示全部内容
which -a bwa

#=========================
#          4 find        #
#=========================
#1 搜索目录下以点fna结尾的文件；
find  ~/ -name *.fna

#2 搜索系统中最近5分钟内编辑过的文件；
find / -amin 5

#3 查找大于100M的文件
find ./ -size 100M

#4 按照文件类型搜索；
find  ./ -type d ;文件类型 #c :档案,d: 目录,b: 区块装置档案 ，p: 具名贮列,f: 一般档案,l: 符号连结,s: socket

#5 查找完进行筛选
find /ifs1/VipData/  -name "*.sh" | grep "[0-9]-.*"

#6 正则表达式筛选
find ~/ -name "[A-Z]*"        #查以大写字母开头的文件
find  /etc -name  "host*"     #查以host开头的文件
find  ~/ -name "[a-z][a-z][0-9][0-9].fa"  #查以两个小写字母和两个数字开头的fa文件

#7 限制目录层次
find  ~ -maxdepth 4 -name *.fna   #查时深度最多为4层
find  ~ -mindepth 3 -name *.fna   #查时深度最多为3层

#8 搜索文件，直接处理
find . -type f -exec ls -l {} \; #一定要加分号

#9 查找并拷贝
mkdir fna
find  ~ -maxdepth 3 -name  "*.fna" | xargs -I {}  cp {} ./fna

#10 查找并删除
find ./fna -name *.fna 
find ./fna -name *.fna -ok rm '{}' \;
find ./fna -name *.fna -exec rm '{}' \;

#搜索当前目录下所有一点fna结尾的文件，然后删除掉。
find ./fna  -name *.fna | xargs rm


########################################################
#                 五十八 文本处理                       #
########################################################
#=========================
#          1 grep        #
#=========================
#1 统计fasta条数
grep ">" ../data/soapdenovo.fa | wc

#2 去除#开头的行
grep -v "^#" ../data/demo.gff | head

#3 根据关键字搜索
grep  "C2381"  ../data/soapdenovo.fa

#4 关键字上下内容
grep -A 1 "C2381"  ../data/soapdenovo.fa
seqkit seq -w 0 ../data/soapdenovo.fa  | grep -A 1 "C2381"
cat ../data/demo.gff | grep "lnc_RNA" 

#5 使用正则表达式
seqkit seq -w 0 ../data/demo.fasta | grep "A\{7,10\}"

#6 筛选关键字，并输出行号
grep -n "EGFR" ../data/demo.bed 

#7 删除掉空行
grep -v "^$" ../data/test.bed 

#8 设定锚定符
locate artic | grep "artic"
locate artic | grep "artic$"

#9 计算匹配字符行数
grep -c "EGFR" ../data/demo.bed

#10 计算数目并排序
grep -v "^#" ../data/demo.gff | awk '{print $3}' | sort-uniq-count-rank

#11 同时满足多个条件
grep -e "ncRNA" -e "ncRNA_gene" ../data/demo.gff

#12 显示包含关键字的文件
grep -l aspera  /ifs1/VipData/*/*.sh

#=========================
#          2 sed         #
#=========================
#1 输出固定的行
cat -n ../data/demo.fasta | sed -n '1307p'
cat -n ../data/demo.fasta | sed -n '100,200p'

#2 替换操作
grep ">" ../data/demo.fasta | sed -e 's/gi/GI/' | head 
sed -i 's/gi/GI/g' ../data/demo.fasta  
sed -i.bak 's#GI#gi#' ../data/demo.fasta  
grep ">" ../data/demo.fasta | sed -e 's/|/#/2;s/ref/REF/' | head

#3 打印发生替换的行
sed -n 's/gi/GI/p' ../data/demo.fasta 

#4 同时进行多条件替换；
sed -f sed.list ../data/demo.fasta  

#5 使用正则表达式替换
grep ">" ../data/demo.fasta | sed -e 's/ .*//g' | head

#6 行首添加内容
sed -e 's/^/time /g' ../data/test.bed 

#7 行尾追加内容
sed -e '$a \the end of file' ../data/test.bed

#8 行寻址
sed -n '/ref/p' ../data/demo.fasta  
grep ">" ../data/demo.fasta | cat -n | sed -n '100,200 s/gi/GI/gp'
grep ">" ../data/demo.fasta | cat -n | sed -n '100,200！ s/gi/GI/gp'

#9 删除操作
sed -e '/>/d' ../data/demo.fasta #删除包含ref的行；

#10 删除空白行；
sed -e '/^\s*$/d'  ../data/demo.fasta

#11 对应替换
sed -e 'y/ATCG/atcg/' ../data/demo.fasta 

#12 DNA序列反向互补配对，并修改大小写
sed -e '/>/!y/ATCG/atcg/' ../data/demo.fasta  

#13 fastq转换为fasta
zcat ../data/demo.fastq.gz | sed '0~4d' | sed '0~3d;s/^@/>/1' 

########################################################
#                     五十九 表格处理                   #
########################################################
#=========================
#          1 cut         #
#=========================
#1 分割文件并输出
cut -d : -f 1,3 /etc/passwd
cut -d: -f 2- /etc/passwd 

#2 选取每个文件前两个字符
ls -1 /ifs1/VipData/ | cut -c 1-2

#=========================
#          2 sort         #
#=========================
#1 排序
awk '{print $2,$3}' ../data/scores.txt | sort

#2 按第二列数字大小排序
awk '{print $2,$3}' ../data/scores.txt | sort -n -k 2

#3 逆序排序
awk '{print $2,$3}' ../data/scores.txt | sort -n -r -k 2 

#4 计算特异项，类似uniq
awk '{print $2,$3}' ../data/scores.txt | sort -u

#5 按多值排序
cat ../data/scores.txt | sort -t $'\t' -k 2 -k 3 

#6 按照第二列中第三个字母排序
cat ../data/scores.txt | sort -t $'\t' -k 2.3    |head

#7 复杂排序
cat ../data/scores.txt | sort -t $'\t' -k 2,2  -k 3nr,3 


#=========================
#         3 uniq         #
#=========================
#1 计算特异
cat ../data/scores.txt | cut  -f 2 | uniq 

#2 计算频数
cat ../data/scores.txt | cut  -f 2 | uniq -c

#3 找出重复项
cat ../data/scores.txt | cut  -f 2,3 | uniq -D

#4 忽略固定列 
cat ../data/scores.txt | uniq -f 1 -D


#=========================
#          4 awk         #
#=========================
#1：输出一个列表任意行；
awk '{print $2}' ../data/blast6.out  | head
awk '{print $NF}' ../data/blast6.out | head

#2 修改分隔符以及输出分隔符
awk -F ":" '{print $1,$NF}' /etc/passwd
awk -F ":" 'OFS="," {print $1,$NF}' /etc/passwd

#3 过滤blast结果
awk  '{if ($3>=80 && $4>=100) print $0}'  ../data/blast6.out 

#4 统计数目
awk  '{if ($3>=80 && $4>=100) print $2}'  ../data/blast6.out | sort | uniq | wc

#5 输出固定行内容
awk 'NR>=20 && NR<=80' ../data/blast6.out 

#6 格式转换
samtools view ../data/demo.bam  | awk '{print"@" $1"\n"$10"\n""+\n"$11""}' | gzip >demo.fq.gz

#7 fastq转换fasta
zcat demo.fq.gz | awk '{print NR":"$0}' | head
zcat demo.fq.gz |awk '{getline l2;getline l3;getline l4;print $0 "\n" l2}' | head
zcat demo.fq.gz |awk '{getline l2;getline l3;getline l4;sub("@",">",$0);print $0 "\n" l2}' | head

#8 模式匹配
awk '$0~ /wang/ {print $0}' /etc/passwd
last -w | awk '$0 ~ /in/ {print $1}' 

#9 BEGIN与END功能
awk 'BEGIN{print "The Program Begin\n"} $0 ~ /wang/ {print $0} END {print "The Program End\n"}' /etc/passwd

#10 转为bed文件格式
cat ../data/test.bed
awk 'BEGIN{OFS="\t"} {print $1,$2,$3}' ../data/test.bed | cat -A

#11 替换
awk '{sub(/Escherichia coli str. K-12 substr. MG1655, complete genome/, " ", $0);print}' ../data/demo.fasta | grep ">" |head

#12 计算占用CPU情况
ps hux -U $USER | awk -v user=$USER '{ sum += $6} END { printf "%s %.2f\n", user,sum/100;}'

#=========================
#          5 bioawk       #
#=========================
#安装
mamba install -y bioawk

#显示文件格式
bioawk -c help

#输出fasta/q名字部分
bioawk -c fastx ' { print $name } ' ../data/demo.fastq.gz
#计算gc含量
bioawk -c fastx ' { print $name, gc($seq) } ' ../data/demo.fastq.gz
bioawk -c fastx ' { print $name, gc($seq) } ' ../data/demo.fasta

#输出CIGAR为deletions的列
samtools view -f 2 demo.bam | awk '$6 ~ /D/ { print $6 }' | head
samtools view -f 2 demo.bam | bioawk -c sam '$cigar ~ /D/ { print $cigar }' | head

#打印vcf文件中的CHROM与POS列
grep -v "^##" demo.vcf | bioawk -tc hdr '{print $_CHROM,$POS}'

#输出比对上的行
bioawk -Hc sam '!and($flag,4)'
 
#反向互补fasta
bioawk -c fastx '{print ">"$name;print revcomp($seq)}' ../data/demo.fasta

#输出vcf中特定genotypes类型
grep -v "^##" in.vcf | bioawk -tc hdr '{print $foo,$bar}'

#=========================
#          6 csvtk       #
#=========================
#安装软件
mamba install -y csvtk

#下载案例文件
git clone https://github.com/shenwei356/csvtk.git

#1 显示csv文件
cat ../data/testdata/names.csv | csvtk pretty

#2 转为markdown
cat ../data/testdata/names.csv  | csvtk csv2md

#3 用列或列名来选择指定列，可改变列的顺序
cat ../data/testdata/names.csv
cat ../data/testdata/names.csv  | csvtk cut -f 3,1 
cat ../data/testdata/names.csv  | csvtk cut -f last_name,id

#4 用通配符选择多列
cat ../data/testdata/names.csv  | csvtk cut -F -f '*name,id' | csvtk pretty

#5 按指定列搜索，默认精确匹配
cat ../data/testdata/names.csv | csvtk grep -f id -p 1 | csvtk pretty

#6 使用正则表达式匹配
cat ../data/testdata/names.csv  | csvtk grep -f id -p 1 -r | csvtk pretty

#7 分组计算
cat ../data/testdata/digitals2.csv | csvtk summary -i -g f1,f2 -f f4:sum,f5:sum | csvtk pretty

#8 合并文件
cat ../data/testdata/names.csv 
cat ../data/testdata/phones.csv 
csvtk join -f 'username;username' --keep-unmatched names.csv phones.csv
csvtk join -f 'username;username' --keep-unmatched ../data/testdata/names.csv ../data/testdata/phones.csv

#9 绘制直方图
csvtk -t plot hist ../data/testdata/grouped_data.tsv.gz -f 2 -o histogram.png

#10 绘制箱线图
csvtk -t plot box ../data/testdata/grouped_data.tsv.gz -g "Group" -f "GC Content" --width 3 -o boxplot.png


#=========================
#          7 datamash    #
#=========================
#安装软件
mamba install -y datamash

#查找案例文件
locate scores_h.txt
cp /ifs1/Software/miniconda3/share/datamash/examples/* ./

# 1 计算1-10的和与平均值
seq 10 | datamash sum 1 mean 1

#2 将数据进行转置
seq 10 | paste - - | datamash transpose
seq 10 | paste - - -  | datamash transpose

#3 按列计算，使用不同的语法
seq 100 | paste - - - - | datamash sum 1 sum 2 sum 3 sum 4

seq 100 | paste - - - - | datamash sum 1,2,3,4

seq 100 | paste - - - - | datamash sum 1-4

seq 100 | paste - - - - | datamash sum 1-3,4

#4 调整分隔符 
seq 10 | xargs -n 5 |datamash -W sum 2

#5 分组计算频数，根据第二列进行分组。如果计算其他值，只需更换函数就行
cat ../data/scores.txt | datamash  groupby 2 count 2

#6 根据第二列进行分组，计算第三列的最大值和最小值
cat ../data/scores.txt  | datamash --sort groupby 2 min 3 max 3

#6 输出表头
cat ../data/scores.txt | datamash --header-out groupby 2 min  3 max 3

#7 跳过第一行表头
cat ../data/scores_h.txt | datamash  groupby 2 min  3 max 3
cat ../scores_h.txt | datamash  --header-in groupby 2 mean 3

#8 使用列名代替列号
cat ../scores_h.txt  |  datamash --headers groupby Major mean Score 

#9 计算频数
cat ../data/genes_h.txt | datamash -s -g 3 count 3


#10 就算特异的频数
cat ../data/genes_h.txt | datamash -s -g 3 countunique 2


########################################################
#             六十  Linux Shell脚本                    #
########################################################
#=========================
#          1 变量         #
#=========================
#定义变量
a=1
b=2
echo $a $b

#引用变量
a=illumina
echo ${a}_1.fq.gz ${a}_2.fq.gz 

#查看内置变量
echo $PATH
echo $PS1
echo $LANG
echo $HISTSIZE

#修改内置变量
export HISTSIZE=20000

#=========================
#          2 for循环      #
#=========================
#连续数字
echo {1..10}

#用户控制循环次数
for i in {1..10};do echo $i;done;

#连续字符串
echo {a..k}
for i in {a..k};do echo "$i";done;

#进度条
for i in {1..100};do  printf "\rComplete: %d%%" $i;sleep 0.1;done;

#批量生成文件
for i in {A..J};do echo touch sample${i}_1.fq.gz sample${i}_2.fq.gz;done;
for i in {A..J};do touch sample${i}_1.fq.gz sample${i}_2.fq.gz;done;

#批量操作
ls -1 *.fq.gz 
for i in *.fq.gz ;do echo "fastqc -f fastq "${i};done;

#=========================
#        2 while循环      #
#=========================
ls -1 *.fq.gz | xargs -n 2
ls -1 *.fq.gz | xargs -n 2 | while read {i,j};do echo $i,$j;done;

# 生成脚本
ls -1 *.fq.gz | xargs -n 2 | while read {i,j};do echo "spades.py -1 $i -2 $j -o spades";done;
# 生成更加完美的脚本
ls -1 *.fq.gz | xargs -n 2 | while read {i,j};do echo "spades.py -1 $i -2 $j -o ${i%_*}";done;

#生成数据列表
ls *.fq.gz | xargs  -I{} echo "$PWD/{}" | xargs -n 2 | awk -F "/" '{print $7,$0}' | sed -e 's/_1.fq.gz//1'
ls *.fq.gz | xargs  -I{} echo "$PWD/{}" | xargs -n 2 | awk -F "/" '{print $7,$0}' | sed -e 's/_1.fq.gz//1' >reads.list

#生成脚本
cat reads.list | while read {i,j,k};do echo $i,$j,$k;done;


#=========================
#          3 判断        #
#=========================
if [ -f $PWD/reads.list ];then echo "file exist";else echo "no such file or dirctory";fi;

#1 一个简单脚本
#/bin/bash

a=$1
b=$2
if [ $a == $b ]
then
    echo "a = b"
elif [ $a -gt $b ]
then
    echo "a > b"
elif [ $a -lt $b ]
then
    echo "a < b"
else
    echo "no condition"
fi
