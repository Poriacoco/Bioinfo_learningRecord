#====================================
#        七十三、开始使用R         #
#====================================
1+1
1+2
1+3
# 设置工作目录
#windows "C:/Users/xxx/Desktop"
getwd()    
setwd("~/")
dir()

x <- plot(x = women$height,y = women$weight,col="red",pch=c(16,17,18),ann = F)

a <- 1
b <- 2
c <- a+b
c
print(c)
1.02 ** 365
0.98 ** 365
ls()
rm(a)
x <- runif(10000000,min = 1,max = 100)
length(x)
sum(x)
mean(x)
var(x)
sd(x)
median(x)
range(x)
abs(x)
sqrt(x)

#====================================
#       七十四、R包管理             #
#====================================

# R 语言包管理
.libPaths()  
#将4.0 的目录复制到.libPaths()得到的结果中  
# R 语言包管理
#选择国内镜像站点
chooseCRANmirror()
chooseBioCmirror()
setRepositories()

#查看R包安装目录
.libPaths()  

#安装vcd包
install.packages("vcd")
install.packages("plotrix")
install.packages("vioplot")
install.packages("venn")
install.packages("RColorBrewer")
install.packages("ggplot2")
install.packages("ggExtra")
install.packages("gcookbook")
install.packages(c("carData","gridExtra"))
install.packages("ggsci")
install.packages("ggpubr")
install.packages("ggthemes")
install.packages(c("Rwordseg","wordcloud2"))
install.packages("pheatmap")
install.packages("qqman") 
install.packages("maps")
install.packages("factoextra")
install.packages("circlize")
install.packages(c("FactoMineR", "ggfortify"))


#查看已安装的包
installed.packages()

#加载R包
library(vcd)
require(vcd)

#升级软件包
update.packages()

#删除扩展包
remove.packages("vcd")

#查看当前环境哪些包加载
find.package()
path.package()

#取消加载包
detach("package:vcd")

#其余R包管理函数
.Library
.Library.site
search()
(.packages())#列出当前包
path.package() #列出加载包的路径
(.packages(all.available=TRUE))#列出有效包

#查看帮助文档
  help()  
?ggplot2  
??heatmap  
help( package="ggplot2" )  
example("heatmap")  
browseVignettes()  
ls(package:base)


# bioconductor包安装

chooseBioCmirror()  
if (!requireNamespace("BiocManager"))
  install.packages("BiocManager")
BiocManager::install()
library()

#R包迁移
#保存A设备上的R包名字列表；
oldip <- installed.packages()[,1]
save(oldip,file = "installedPacckages.Rdata")
#在B设备上进行安装；
load("installedPacckages.Rdata")
newip <- installed.packages()[,1]
for (i in setdiff(oldip,newip)) {
  install.packages(i)
}


#====================================
#        七十五、文件操作           #
#====================================

#3 读写文件
dir()  
#设置工作目录
dta <- read.table("Rdata/CountMatrix.csv")  

dta <- read.table(file = "Rdata/CountMatrix.csv",header = T,sep = ",")  
dta <- read.csv(file = "Rdata/CountMatrix.csv",header = T,row.names = 1)  
dta <- read.csv(file = "Rdata/CountMatrix.csv",header = T,sep = ",",row.names = 1,na.strings = "NA",stringsAsFactors = F)  
class(dta)  
head(dta)  
tail(dta)  
rownames(dta)  
colnames(dta) 
View(dta)

#4 写入文件
write.csv(x = dta,file = "matrix.csv",sep = ",",append = F) 

#5 读写excel文件
#安装openxlsx包
#install.packages("openxlsx")
library(openxlsx)
x <- read.xlsx("Rdata/vlookup.xlsx",sheet = 1)
write.xlsx(x = dta,sheetName = "dta",file = "dta.xlsx",append = F)


#6 读写R文件
?iris
head(iris)
getwd()
saveRDS(iris,file="iris.RDS")
rdsdata <- readRDS("iris.RDS")
#Write RData file
save(iris,iris3,file = "iris.Rdata")
save.image()
load(file = "iris.Rdata")


#====================================
#       七十六、数据结构              #
#====================================
# R内置数据集
# 显示所有内置数据集
data()
#加载扩展包数据集
data(package = "MASS") 
data(package = "ggplot2")

#R内置向量
rivers
euro
plot(rivers)  
hist(rivers)  
names(euro)  
length(euro)  
a <- c(1,2,3)  
#字符型向量
b <- c("red","green","blue")  
#生成连续型向量
c <- seq(1,100,2)  
d <- rep(c(1,2),5)  
class(a)  
typeof(a)  


# 向量化操作  
a  
a*2  
paste(b,a,sep = "-")  

# 向量化操作
a
a*2
paste(b,a,sep = "-")

#向量计算
x <- runif(10000000,min = 1,max = 100)  
length(x)  
sum(x)  
mean(x)  
var(x)  
sd(x)  
median(x)  
range(x)  
abs(x)  
sqrt(x)  
a <- 1:5
b <- 10:1
a+1
a*2
a+b

# 实践：绘制人染色体长度分布图
hg19_len <- read.csv(file = "Rdata/homo_length.csv",header = T,row.names = 1)
hg19_24 <- hg19_len[1:24,]   
barplot(height = hg19_24)  
barplot(height = hg19_24,names.arg = names(hg19_24))  
library(RColorBrewer)  
cols <- brewer.pal(n = 6,name = "Set1")  
barplot(height = hg19_24,names.arg = names(hg19_24),col = cols) 

#5矩阵（Matrix）是一个按照长方阵列排列的复数或实数集合。向量是一维的，而矩阵是二维的，需要有行和列。
state.x77
m <- matrix(1:20,nrow = 4,ncol = 5)
rowSums(state.x77)
colMeans(state.x77)
cor(state.x77)

#====================================
#        2、热图                    #
#====================================
#1 heatmap()
#2 gplots 包 heatmap.2()
#3 pheatmap包 pheatmap
#4 ComplexHeatmap包

example("heatmap")
m <- read.csv("Rdata/heatmap.csv",header = T,row.names = 1)
class(m)
x <- as.matrix(m)
heatmap(x)
heatmap(t(x))
heatmap(x,col=c("green","red"))
yanse <- colorRampPalette(c("red","black","green"))(nrow(x))
heatmap(x,col=yanse)
heatmap(x,col=yanse,RowSideColors = yanse)
heatmap(x,col=yanse,ColSideColors = colorRampPalette(c("red","black","green"))(ncol(x)))
heatmap(x,col=yanse,Rowv = NA)
heatmap(x,col=yanse,Rowv = NA,Colv = NA)
heatmap(state.x77)
heatmap(state.x77,scale = "col")

#利用pheatmap绘制热图
#install.packages("pheatmap")
library(pheatmap)
x <- read.csv("Rdata/heatmap.csv",header = T,row.names = 1)
class(x)
pheatmap(x)
annotation_col <- data.frame(CellType=factor(rep(c("N1","T1","N2","T2"),each=5)))
rownames(annotation_col) <- colnames(x)
pheatmap(x,annotation_col = annotation_col)
pheatmap(x,annotation_col = annotation_col,display_numbers = T)
pheatmap(x,annotation_col = annotation_col,display_numbers = T,number_format = "%.2f")
pheatmap(x,annotation_col = annotation_col,display_numbers = T,number_format = "%.1f",number_color = "black")
#6 数据框
mtcars
iris
adata <- data.frame(state.name,state.division,state.area)
str(adata)

#6.2 数据框索引
euro[1:3]  
euro[-(1:3)]  
euro["ATS"]  
mtcars[1:4,1:5]  
mtcars[c(2,3,4,5)]  
mtcars$mpg  
mtcars["Mazda RX4",]  
alist <- list(euro,iris,mtcars)  
alist[1]  
alist[[1]]  
x <- runif(100,min = 1,max = 100)  
x[x>50]  
x [x>20 & x <40] 

#6.3 利用R实现vlookup功能
#=VLOOKUP($A2,Sheet1!$A$1:$I$201,COLUMN(Sheet1!A2),FALSE)

#读入200个基因的列表 bigtable  
genes200 <- read.csv("Rdata/200genes.csv",header = T,stringsAsFactors = F)  

#读入121个基因的list small list  
genes121 <-  read.csv("Rdata/121genes.csv",header = T,stringsAsFactors = F)  
head(genes200)  
head(genes121)  
genes200["CLU",]  

#修改gene table行名  
rownames(genes200) <- genes200$gene  
genes200["CLU",]  
genes200[c("CLU","DCN"),]  
#去除重复gene ID  
gene93 <- unique(genes121$gene)  
#利用数据框的访问功能，一次实现Excel Vlookup功能  
dta <- genes200[gene93,]  
dta  
#数据中包含没有检索到的gene ID，返回值为NA，利用na.omit删除包含NA的行  
dta <- na.omit(dta)  
rownames(dta) <- 1:nrow(dta)  
write.csv(dta,file = "Rdata/gene93.csv",row.names = F) 


#7 因子
state.division  
state.region  
table(state.division)  
afactor <- factor(c("blue","red","red","green","red","blue"),levels = c("red","green","blue","yellow"))  
bfactor <- factor(letters[1:5], labels = c("one","two","three","four","five"))

#8 列表
state.center
alist <- list(euro,state.x77,iris)

#9 时间序列
co2  
presidents 

#10 缺失数据
x <- c(1:10,NA)
sum(x)
sum(x,na.rm = T)
is.na(x)
#install.packages("VIM")  
library(VIM)  
data(sleep,package = "VIM")  
na.omit(sleep) 

#====================================
#        七十七、数据处理           #
#====================================
#1 判断数据类型
x <- c(1:100)  
#判断数据类型  
is.vector(x)  
is.character(x)  

#转换数据类型  
as.character(x)  
is.numeric(mtcars$cyl)  
as.factor(mtcars$cyl)  

#不同数据类型绘图结果不同  
plot(as.numeric(mtcars$cyl))  
plot(as.factor(mtcars$cyl)) 

#2 修改数据内容
#rbind,cbind  
state <- cbind(as.data.frame(state.x77),state.region)  
tapply(state$Population, state$state.region, mean)  
a <- head(mtcars,20)  
b <- tail(mtcars,20)  
c <- rbind(a,b)  
length(rownames(c))  
duplicated(c)  
c[duplicated(c),]  
c[!duplicated(c),]  
unique (c)  

#取子集subset
dta <- read.csv("Rdata/gene93.csv",header = T)  
dta  
#过滤掉没有表达的转录本  
dta1 <- subset(dta,rowSums(dta[2:9])>1)


#3 不同数类型转换
#向量和矩阵之间相互转换：给向量加上维度。
x <- c(1:10)
dim(x) <- c(2,5)
#向量和数据框之间相互转换：data.frame，cbind和rbind将向量转换为数据框，取出数据框的每一列为一个向量。
state <- data.frame(state.name,state.abb,state.division,state.area)
#数据框和矩阵之间相互转换：
#as.matrix()将数据框转换为矩阵
iris.mat <- as.matrix(iris[1:4])
#as.data.frame()矩阵转换为数据框。
state <- as.data.frame(state.x77)
#向量和因子之间相互转换；as.factor()函数。
mtcars$cyl <- as.factor(mtcars$cyl)
#修改因子的水平和标签；
mtcars$cyl <- factor(mtcars$cyl,levels = c(4,6,8),labels = c("four","six","eight"))

#4 随机抽样
x <- runif(n=100,min = 1,max = 100)
round(x,digits = 0)
sample(x,10,replace = F)
sample(x,10,replace = T)
set.seed(1234)
sample(x,10,replace = F)

#5 利用R语言斗地主
type <- c("red","spades","cube","plum")
amount <- c("A",2:10,"J","Q","K")
group <- expand.grid(type,amount)
poker <- paste(group$Var2,group$Var1,sep = "-")
poker
poker[c(53,54)] <- c("black Joker)","red Joker")
shuffle <- sample(poker,54,replace = F)
one <- shuffle[c(T,F,F)]
two <- shuffle[c(F,T,F)]
three <- shuffle[c(F,F,T)]

#6 排序
#一维数据排序
sort(rivers)  
order(rivers)  
rivers[order(rivers)]  
sort(rivers,decreasing = T)
rev(sort(rivers))  

#对二维数据进行排序  
#sort(mtcars)  
#order(mtcars)  
order(mtcars$mpg)  
mtcars[order(mtcars$mpg),]  

#同时两个排序条件  
mtcars[order(mtcars$mpg,mtcars$cyl),]  
#反向排序  
mtcars[order(mtcars$mpg,-mtcars$gear),]  

#7 筛选
rivers > 100
rivers[rivers>100]
rivers[rivers >100 & rivers < 500]
rivers[rivers==900]
which.max(rivers)
which.min(rivers)
which(rivers==320)
who <- read.csv("Rdata/WHO.csv",header = T)
who1 <- who[c(1:50),c(1:10)]
who2 <- who[c(1,3,5,8),c(2,14,16,18)]
who3 <- who[which(who$Continent==7),]
who4 <- who[which(who$CountryID>50 & who$CountryID<=100),]
#subset
?subset
who3=subset(who,who$CountryID>50 & who$CountryID<=100)

#8 数据透视表
data("mtcars")
apply(mtcars,2,sum)  
apply(iris[1:4], 2, sd)  
lapply(state.center,mean)  
state.name  
state.division  
tapply(state.name,state.division,length)

