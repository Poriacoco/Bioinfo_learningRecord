
########################################################
#                   49 全局比对                         #
########################################################
mkdir 07.alignment;cd 07.alignment/;

#下载基因组序列
#ncbi搜索  Klebsiella pneumoniae
mkdir data;cd data;
#将拼接好的基因组文件拷贝至目录下命名为mgh78578.fasta
efetch -db nuccore -format fasta -id "NC_016845.1 NC_016838.1 NC_016846.1 NC_016839.1 NC_016840.1 NC_016847.1 NC_016841.1" >ref.fna
efetch -db nuccore -format fasta -id "CP042858.1" >qd23.fna
efetch -db nuccore -format fasta -id "CP087123.1" >th12908.fna
#将拼接好的基因组文件拷贝至目录下；
#安装软件
mamba install -y mummer
mamba install -y gnuplot


#2 mummer比对
mkdir 49.mummer;cd 49.mummer;

#nucmer比对
nucmer --mum --maxgap=500 --mincluster=100 --prefix=nucmer ../data/ref.fna ../data/mgh78578.fasta
delta-filter -1 -q -r nucmer.delta > nucmer.filter

#显示比对结果
grep ">" nucmer.delta
show-aligns nucmer.filter NC_016846.1 contig_1_pilon

#显示差别
show-diff nucmer.filter -q
show-diff nucmer.filter -r

#显示突变位点
show-snps -C -H -I -T -r -l nucmer.filter >nucmer.snp

#显示坐标
show-coords nucmer.filter -r >nucmer.coords

#show-tiling 轨迹
cp ../../05.assembly/35.illumina/4.soapdenovo/kmer45/kmer45.scafSeq .
nucmer --mum --maxgap=500 --mincluster=100 --prefix=kmer45 ../data/ref.fna kmer45.scafSeq
delta-filter -1 -q -r kmer45.delta > kmer45.filter
show-tiling kmer45.filter -a
show-tiling  kmer45.filter -l 10000 >kmer45.tiling

#mummerplot绘图
mummerplot -p p1 nucmer.filter --png 
mummerplot -p p2 nucmer.filter --png --medium
mummerplot -p kmer45 kmer45.tiling --png --medium

#promer比对
promer --mum --maxgap=500 --mincluster=100 --prefix=promer ../data/ref.fna ../data/mgh78578.fasta

#dnadiff比对
dnadiff ../data/ref.fna ../data/mgh78578.fasta

########################################################
#                   50 共线性分析                     #
########################################################
#安装软件
#dotplotly
wget https://github.com/tpoorten/dotPlotly/archive/refs/heads/master.zip
#配置R包
install.packages(c("optparse", "ggplot2", "plotly"))

mkdir 50.synteny;cd 50.synteny;

#利用nucmer+dotPlotly绘图
nucmer --maxmatch -l 80 -c 100 ../data/ref.fna ../data/mgh78578.fasta -p nucmer
delta-filter -r nucmer.delta >filter.delta
show-coords -c  filter.delta >filter.delta.coords
./dotPlotly-master/mummerCoordsDotPlotly.R -i filter.delta.coords -o filter.delta.plot  -m 1000 -q 300000 -k 10 -s -t -l -p 12

用Minimap2 + dotPlotly
minimap2 -x asm5 ../data/ref.fna ../data/mgh78578.fasta >minimap2.paf
./dotPlotly-master/pafCoordsDotPlotly.R -i minimap2.paf -o minimap2.plot -m 2000 -q 500000 -k 10 -s -t -l -p 12


#安装mscanx
git clone https://github.com/wyp1125/MCScanX.git
cd MCScanx
make

#准备输入文件
#下载参考序列与gff文件
#https://www.ncbi.nlm.nih.gov/genome/?term=Klebsiella%20pneumoniae
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/240/185/GCF_000240185.1_ASM24018v2/GCF_000240185.1_ASM24018v2_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/240/185/GCF_000240185.1_ASM24018v2/GCF_000240185.1_ASM24018v2_genomic.gff.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/240/185/GCF_000240185.1_ASM24018v2/GCF_000240185.1_ASM24018v2_protein.faa.gz

#处理gff或者bed格式文件
perl -F'\t' -lane 'next unless $F[2] eq "CDS";print join qq{\t},$F[0],$F[-1]=~s/ID=cds-([^;]+).Parent=.*$/$1/r,$F[3],$F[4]' GCF_000240185.1_ASM24018v2_genomic.gff >ref.gff
grep ">" mgh78578.faa | awk '{print "contig\t"$1"\t"$3"\t"$5}' | sed -e 's/>//' >mgh78578.gff

#合并文件
cat GCF_000240185.1_ASM24018v2_protein.faa mgh78578.faa >all.faa
cat cat ref.gff mgh78578.gff >test/all.gff


#blast比对
makeblastdb -in all.faa -dbtype prot -out all -parse_seqids
blastp -query all.faa -db all -out all.blast -evalue 1e-10 -num_threads 12 -outfmt 6 -num_alignments 5

#运行结果
mkdir all;
mv all.blast all.gff all
~/biosoft/MCScanX/MCScanX all/all

#结果
all.collinearity
all.html/
all.tandem

#过滤duplicate_gene
~/biosoft/MCScanX/duplicate_gene_classifier all/all

#下游分析，将downstream_analyses程序拷贝至结果目录
#案例地址：https://github.com/wyp1125/MCScanx#the-following-is-the-list-of-executable-programs
cp ~/biosoft/MCScanX/downstream_analyses/* .

#1 Detect_syntenic_tandem_arrays
detect_collinear_tandem_arrays -g os_sb.gff -b os_sb.blast -c os_sb.collinearity -o exam1.txt

#2 Dissect_multiple_alignment
dissect_multiple_alignment -g os_sb.gff -c os_sb.collinearity -o exam2.txt

#3 dot_plotter
java dot_plotter -g os_sb.gff -s os_sb.collinearity -c dot.ctl -o exam3.png

#4 dual_synteny_plotter
java dual_synteny_plotter -g os_sb.gff -s os_sb.collinearity -c dot.ctl -o exam4.png

#5 Circle_plotter
java circle_plotter -g os_sb.gff -s os_sb.collinearity -c circle.ctl -o exam5.png

#6 Bar_plotter
java bar_plotter -g os_sb.gff -s os_sb.collinearity -c bar.ctl -o exam6.png

#7 add_kaks_to_synteny.pl
perl add_kaks_to_synteny.pl -i os_sb.collinearity -d cds_file -o exam7

#8 group_collinear_genes.pl
perl group_collinear_genes.pl -i os_sb.collinearity -o exam8.cluster


########################################################
#                   51 基因家族分析                     #
########################################################
mkdir 51.genefamily；cd 51.genefamily;
#序列下载：下载水稻全部氨基酸序列以及GFF文件
https://phytozome.jgi.doe.gov/pz/portal.html 
http://plants.ensembl.org/index.html  

#下载水稻Dynamin_N家族hmm文件
http://pfam.xfam.org/family/PF00350#curationBlock

#hmmersearch
hmmsearch -o hmmer.out --domE 1E-5 -E 1E-5 Dynamin_N.hmm Oryza_sativa.IRGSP-1.0.pep.all.fa

#复制比对上的基因ID，保存到ids.txt文件中

#根据比对上的ID提取序列
samtools faidx Oryza_sativa.IRGSP-1.0.pep.all.fa
cat ids.txt  |xargs -n 1 samtools faidx Oryza_sativa.IRGSP-1.0.pep.all.fa >Oryza_family.fa


########################################################
#                   52 测序数据比对                     #
########################################################
#安装软件
mamba install -y bwa
mamba install -y subread
mamba install -y hisat2
mamba install -y minimap2
mamba install -y bwa-mem2
mamba install -y bowtie2

mkdir 52.bwa
#1 bwa比对
#建立索引
ln -s ../data/mgh78578.fasta .

#bwa比对
bwa mem mgh78578.fasta ../../05.assembly/data/illumina.sra_1.fastq.gz ../../05.assembly/data/illumina.sra_2.fastq.gz >mgh78578.sam

#bwa-mem2比对
bwa-mem2 index ../data/mgh78578.fasta
bwa-mem2 mem mgh78578.fasta ../../05.assembly/data/illumina.sra_1.fastq.gz ../../05.assembly/data/illumina.sra_2.fastq.gz >mgh78578.sam

#bwa+samtools一条命令
bwa-mem2 mem -t 12 mgh78578.fasta ../../../05.assembly/data/illumina.sra_1.fastq.gz ../../../05.assembly/data/illumina.sra_2.fastq.gz | samtools sort -O bam - >mgh78578.sorted.bam

#拟南芥比对
ln -s /ifs1/VipData/05.assembly/ninanjie/tair10.fna
bwa-mem2 index tair10.fna
bwa-mem2 mem -t 12 tair10.fna /ifs1/VipData/05.assembly/ninanjie/illumina/il_1.fq.gz /ifs1/VipData/05.assembly/ninanjie/illumina/il_2.fq.gz >tair10.sam 2>bwa.log

#bowtie2与hisat2比较
ln -s /ifs1/VipData/07.aligment/data/chrX.fa .
bowtie2-build chrX.fa chrX
bowtie2 -x chrX -p 12 -1 /ifs1/VipData/07.aligment/data/ERR188044_chrX_1.fastq.gz -2 /ifs1/VipData/07.aligment/data/ERR188044_chrX_2.fastq.gz >bowtie2.sam 2> bowtie2.log 

hisat2-build chrX.fa chrX
hisat2 -p 12 --dta -x  chrX -1 /ifs1/VipData/07.aligment/data/ERR188044_chrX_1.fastq.gz -2 /ifs1/VipData/07.aligment/data/ERR188044_chrX_2.fastq.gz -S hisat2.sam 2> hisat2.log 

#2 minimap2比对
#minimap2建立索引
minimap2 mgh78578.fasta -d mgh78578.min
#minimap2比对
minimap2 -ax map-ont mgh78578.fasta ../../05.assembly/data/nanopore.sra.fastq.gz -t 12 -o mgh78578_ont.sam

#pacbio比对
minimap2 -ax map-pb mgh78578.fasta ../../05.assembly/data/pacbio.sra.fastq.gz -t 12 -o mgh78578.pb.sam

#短序列比对
minimap2 -sr mgh78578.fasta ../../05.assembly/data/illumina.sra_1.fastq.gz ../../05.assembly/data/illumina.sra_2.fastq.gz -t 12 -o mgh78578.ill.sam

#reads直接比对，找overlap
minimap2 -x ava-ont ../../05.assembly/data/nanopore.sra.fastq.gz ../../05.assembly/data/nanopore.sra.fastq.gz > ovlp.paf 

#序列比对
minimap2 -x asm5 ref.fna mgh78578.fasta -o mgh78578_ref.paf

#paf文件可视化
https://tom-poorten.shinyapps.io/dotplotly_shiny/


########################################################
#                   53 samtools使用                    #
########################################################
#安装软件
mamba install -y samtools=1.14
mamba install -y bamtools
#查看文件路径
which samtools
#查看选项
samtools

#2 案例
mkdir 52.samtools;cd 52.samtools;
cp ../52.bwa/mgh78578.sam all.sam

#1 sam文件验证
samtools quickcheck *.sam  && echo 'all ok' || echo 'fail!'

#2 sam和bam格式转换
samtools view -O bam -o all.bam all.sam  
samtools view all.bam -o all.sam

#3 bam排序
samtools sort -@ 4 -m 12G -O bam -o all.sorted.bam all.sam  

#4 排序后建立索引
samtools index all.sorted.bam  

#5 比对结果统计
samtools stats all.sorted.bam

#6 按照flag值进行统计
samtools flagstat all.sorted.bam

#7 按不同染色体统计
samtools idxstats all.sorted.bam 

#8 统计目标区域覆盖度
samtools bedcov test.bed all.sorted.bam

#9 统计每个位点测序深度depth
samtools depth all.sorted.bam

#10计算覆盖比率
samtools coverage all.sorted.bam

#11 按照染色体拆分bam
bamtools split -in all.sorted.bam -reference

#12 合并bam
samtools cat -o cat.bam all.sorted.REF_contig_1_pilon.bam all.sorted.REF_contig_2_pilon.bam all.sorted.REF_contig_3_pilon.bam all.sorted.REF_unmapped.bam
samtools merge merge.bam all.sorted.REF_contig_1_pilon.bam all.sorted.REF_contig_2_pilon.bam all.sorted.REF_contig_3_pilon.bam all.sorted.REF_unmapped.bam
#13 统计bam并绘图  
samtools stats all.sorted.bam >all.stats  
plot-bamstats -p test all.stats 

#14 重新打乱
samtools collate all.sorted.bam -o all.shuffle.bam

#15 过滤数据
#将没有比对上的reads筛选出来
samtools view -f 4 all.sorted.bam

#16 保留比对上的reads输出出来
samtools veiw -F 4 all.sorted.bam

#17 输出fastq格式 
samtools fastq all.sorted.bam -1 A.1.fq.gz -2 A.2.fq.gz
#输出fasta格式
samtools fasta all.sorted.bam -1 A.1.fa.gz -2 A.2.fa.gz

#18 文本方式查看结果tview  
samtools tview all.sorted.bam  
samtools tview all.sorted.bam ref.fna  

#19 bam转换为bed格式
bedtools bamtobed -i all.sorted.bam

#20 mpileup
samtools mpileup --reference ../data/mgh78578.fasta  all.sorted.bam



