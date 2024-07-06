#====================================
#       八十二、R语言基础绘图        #
#====================================
#安装R包
install.packages("plotrix")
install.packages("vioplot")
install.packages("venn")
install.packages("RColorBrewer")
#切换工作目录
#setwd("C:/Users/xxx/Desktop/Rcourse")

#====================================
#        绘图设备                  #
#====================================
x11()
pdf()
dev.list()
dev.off(3)
dev.list()
dev.off(4)
dev.list()

#====================================
#        1、点图                      #
#====================================
m <- read.table("Rdata/prok_representative.csv",sep = ",",header = T);  
x <- m[,2]
y <- m[,4]  
plot(x,y,pch=16,xlab="Genome Size",ylab="Genes");  
fit <- lm(y~x);  
abline( fit,col="blue",lwd=1.8 );  
rr <- round( summary(fit)$adj.r.squared,2);  
intercept <- round( summary(fit)$coefficients[1],2);  
slope <- round( summary(fit)$coefficients[2],2);  
eq <- bquote( atop( "y = " * .(slope) * " x + " * .(intercept), R^2 == .(rr) ) );  
text(12,6e3,eq); 


#====================================
#        2、直方图                      #
#====================================
#基因长度分布图
x <- read.table("Rdata/H37Rv.gff",sep = "\t",header = F,skip = 7,quote = "")  
x <- x[x$V3=="gene",]  
x <- abs(x$V5-x$V4)+1  
length(x)  
range(x)  
hist(x)  
hist(x,breaks = 80)  
hist(x,breaks = c(0,500,1000,1500,2000,2500,15000))  
hist(x,breaks = 80,freq = F)  
hist(x,breaks = 80,density = T)  
hist(rivers,density = T,breaks = 10)  
?hist  
h <- hist(x,nclass=80,col="pink",xlab="Gene Length (bp)",main="Histogram of Gene Length");  
rug(x);
xfit<-seq(min(x),max(x),length=100);
yfit<-dnorm(xfit,mean=mean(x),sd=sd(x));
yfit <- yfit*diff(h$mids[1:2])*length(x);
lines(xfit, yfit, col="blue", lwd=2);

#====================================
#        3、条形图                      #
#====================================
# 实践：绘制人染色体长度分布图
hg19_len <- read.csv(file = "Rdata/homo_length.csv",header = T,row.names = 1)
hg19_24 <- hg19_len[1:24,]   
barplot(height = hg19_24)  
barplot(height = hg19_24,names.arg = names(hg19_24))  
library(RColorBrewer)  
cols <- brewer.pal(n = 6,name = "Set1")  
barplot(height = hg19_24,names.arg = names(hg19_24),col = cols) 
#绘制分组条形图
x <- read.csv("Rdata/sv_distrubution.csv",header = T,row.names = 1)  
x  
#barplot(x)  
barplot(as.matrix(x))  
barplot(t(as.matrix(x)))  
barplot(t(as.matrix(x)),col = rainbow(4))  
barplot(t(as.matrix(x)),col = rainbow(4),beside = T)  
barplot(t(as.matrix(x)),col = rainbow(4),legend.text = colnames(x))  
barplot(t(as.matrix(x)),col = rainbow(4),legend.text = colnames(x),ylim = c(0,800))  
barplot(t(as.matrix(x)),col = rainbow(4),legend.text = colnames(x),ylim = c(0,800),  
        main = "SV Distribution",xlab="Chromosome Number",ylab="SV Numbers")  

#====================================
#        4、饼图                     #
#====================================
x <- read.csv("Rdata/homo_length.csv",header = T)
x <- x[1:24,]
barplot(height = x$length,names.arg = x$chr)
pie(x$length/sum(x$length))

m <- read.table("Rdata/Species.txt");
m
x <- m[,3]
pie(x);
pie(x,col=rainbow(length(x)))
lbls <- paste(m[,1],m[,2],"\n",m[,3],"%" )
pie(x,col=rainbow(length(x)),labels = lbls)
pie(x,col=rainbow(length(x)),labels = lbls)
pie(x,col=rainbow(length(x)),labels = lbls,radius = 1)
pie(x,col=rainbow(length(x)),labels = lbls,radius = 1,cex=0.8)

#3D饼图
#install.packages("plotrix")
library(plotrix)
pie3D(x,col=rainbow(length(x)),labels = lbls)
pie3D(x,col=rainbow(length(x)),labels = lbls,cex=0.8)
pieplot <- pie3D(x,col=rainbow(length(x)),radius = 1,explode = 0.1)
pie3D.labels(pieplot,labels = lbls,labelcex = 0.8,height = 0.1,labelrad = 1.75)

#扇形图
fan.plot(x,col=rainbow(length(x)),labels = lbls,cex=0.8,radius = 1)
#====================================
#        5、箱线图                     #
#====================================

boxplot(mpg ~cyl,data = mtcars)
boxplot(len ~ dose, data = ToothGrowth)
boxplot(len ~ dose:supp, data = ToothGrowth)
boxplot(len ~ dose:supp, data = ToothGrowth,col = c("orange", "red"))
boxplot(len ~ dose:supp, data = ToothGrowth,col = c("orange", "red"),notch=T)
boxplot(len ~ dose:supp, data = ToothGrowth,col = c("orange", "red"),width=c(1,0.5,1,0.5,1,0.5))
boxplot(len ~ dose:supp, data = ToothGrowth,col = c("orange", "red"),varwidth=T)
boxplot(len ~ dose:supp, data = ToothGrowth,col = c("orange", "red"),boxwex=0.5)
boxplot(len ~ dose:supp, data = ToothGrowth,col = c("orange", "red"),staplewex=0.5)
boxplot(len ~ dose:supp, data = ToothGrowth,col = c("orange", "red"),sep = ":",lex.order = T)


#====================================
#        6、小提琴图                       #
#====================================
#install.packages("vioplot")
library(vioplot)
vioplot(len ~ supp+dose,data = ToothGrowth)
vioplot(len ~ supp+dose,data = ToothGrowth,col=rep(c("cyan","violet"),3))

#====================================
#        7、韦恩图                  #
#====================================
#install.packages("venn")
library(venn)
listA <- read.csv("Rdata/genes_list_A.txt",header=FALSE)
A <- listA$V1
listB <- read.csv("Rdata/genes_list_B.txt",header=FALSE)
B <- listB$V1
listC <- read.csv("Rdata/genes_list_C.txt",header=FALSE)
C <- listC$V1
listD <- read.csv("Rdata/genes_list_D.txt",header=FALSE)
D <- listD$V1
listE <- read.csv("Rdata/genes_list_E.txt",header=FALSE)
E <- listE$V1
alist <- list(A,B,C,D,E)
alist
venn(alist)
venn(alist[1:3])
venn(alist,snames = "A,B,C,D,E", counts = NULL, ilabels = FALSE, ellipse = FALSE,
     zcolor = "bw", opacity = 0.3, size = 15, cexil = 0.6, cexsn = 0.85,
     borders = TRUE)

venn(alist,col = "red",zcolor = "blue")
venn(alist,col = c("red","blue"),zcolor = c("blue","green"))
venn(alist[1:4],col = c("red","blue"),zcolor = c("blue","green"),ellipse = T)
venn(alist[1:4],zcolor = rainbow(5),ellipse = T,ilabels =T )

#====================================
#        布局                 #
#====================================

#mfrow或者mfcol
opar <- par(no.readonly = TRUE)
par(mfrow=c(2,2))
plot(pressure,col="red",main="Pic 1")
barplot(table(mtcars$cyl),col = c("red","cyan","orange"),main = "Pic 2")
hist(rivers,breaks = 30,col = "pink",main = "Pic 3")
pie(c(1,3,4,2),labels = c("A","B","C","D"),main = "Pic 4")

#layout布局
layout(matrix(c(1,2)),heights = c(2,1))
layout.show(2)
#绘图
barplot(table(mtcars$cyl),col = c("red","cyan","orange"),main = "Pic 2")
plot(pressure,col="red",type="l",main="Pic 1")
layout(matrix(c(0,2,0,0,1,3),2,3,byrow=T),widths = c(0.5,3,1),heights =  c(1,3,0.5),TRUE);
layout.show(3)
plot(pressure,col="red",main="Pic 1")
barplot(table(mtcars$cyl),col = c("red","cyan","orange"),main = "Pic 2")
hist(rivers,breaks = 30,col = "pink",main = "Pic 3")
#恢复
par(opar)


#====================================
#       八十三、ggplot2        #
#====================================

#本案例需要安装扩展包
install.packages("ggplot2")
install.packages("ggExtra")
install.packages("gcookbook")
install.packages(c("carData","gridExtra"))
install.packages("ggsci")
install.packages("ggpubr")
install.packages("ggthemes")

#====================================
#        ggplot2基础                #
#====================================
#加载包
library(ggplot2)
plot(mtcars$wt,mtcars$mpg)
#1数据
ggplot(data=mtcars)
#2映射
ggplot(data=mtcars, mapping = aes(x=wt, y=mpg)) 
#3几何图形
ggplot(data=mtcars, mapping = aes(x=wt, y=mpg)) + geom_point() 
ggplot(data=mtcars, mapping = aes(x=wt, y=mpg,size=mpg)) + geom_point() 
ggplot(data=mtcars, mapping = aes(x=wt, y=mpg,size=mpg,shape=cyl)) + geom_point() 
ggplot(data=mtcars, mapping = aes(x=wt, y=mpg,size=mpg,shape=cyl,color=cyl)) + geom_point() 

#4标尺（Scale）
ggplot(data=mtcars, mapping = aes(x=wt, y=mpg,color=mpg)) + geom_point()+
  scale_color_gradient(low = "orange",high = "red")
#5统计变换（Statistics）
ggplot(data=mtcars, mapping = aes(x=wt, y=mpg)) + geom_point()+
  stat_smooth( method = 'loess' ,formula =  'y ~ x')

#6坐标（Coordinate） 
ggplot(data=mtcars, mapping = aes(x=wt, y=mpg)) + geom_point()+coord_flip()

#7图层（Layer）
ggplot(data=mtcars, mapping = aes(x=cyl, y=mpg)) + geom_point()+geom_boxplot()

#8分面
ggplot(data=mtcars, mapping = aes(x=wt, y=mpg)) + geom_point()+facet_grid(. ~ cyl)

#9主题（Theme）
ggplot(data=mtcars, mapping = aes(x=wt, y=mpg)) + geom_point()+theme_bw()

#保存绘图
p <- ggplot(data=mtcars, mapping = aes(x=wt, y=mpg)) + geom_point()+theme_bw()
ggsave(filename = "mtcars.pdf",plot = p)

#====================================
#        4.1 点图                      #
#====================================
x <- read.table("Rdata/prok_representative.csv",sep = ",",header = T);  
head(x)
ggplot(data = x,aes(x=Size,y=Genes))+geom_point()
ggplot(data = x,aes(x=Size,y=Genes))+geom_point(size=1,color="blue")
fit <- lm(data = x,Genes~ Size)
fit
ggplot(data = x,aes(x=Size,y=Genes))+geom_point(size=1,color="blue")+geom_abline(intercept = 286.6,slope = 843.7,col="red",lwd=1)
p <- ggplot(data = x,aes(x=Size,y=Genes))+geom_point(size=1,color="blue")+geom_abline(intercept = 286.6,slope = 843.7,col="red",lwd=1)
p+annotate(geom = "text",x=4,y=10000,label="y=286x+842.7\nR2=0.9676")
p+annotate(geom = "text",x=4,y=10000,label="y=286x+842.7\nR2=0.9676")+labs(title="Genome Size vs Gene Number",x="Genome Size",y="Genes")

#====================================
#       4.2 直方图                      #
#====================================
#基因长度分布图
x <- read.table("Rdata/H37Rv.gff",sep = "\t",header = F,skip = 7,quote = "")  
x <- x[x$V3=="gene",]  
x <- abs(x$V5-x$V4+1)  
length(x)  
range(x)  
ggplot(data = NULL,aes(x=x))
ggplot(data = NULL,aes(x=x))+geom_histogram(bins = 80)
ggplot(data = NULL,aes(x=x))+geom_histogram(bins = 80)+geom_rug()


#====================================
#        4.3 条形图                      #
#====================================
# 实践：绘制人染色体长度分布图
hg19_len <- read.csv(file = "Rdata/homo_length.csv",header = T)
x <- hg19_len[1:24,]   

ggplot(data = x,aes(x=chr,y=length,fill=chr))+geom_bar(stat = "identity")
p <- ggplot(data = x,aes(x=chr,y=length,fill=chr))+geom_bar(stat = "identity")
p+scale_x_discrete(limits=x$chr)
p+scale_x_discrete(limits=x$chr)+coord_flip()
p+scale_x_discrete(limits=x$chr)+coord_flip()+guides(fill=FALSE)

#绘制分组条形图
library(tidyverse)
x <- read.csv("Rdata/sv_distrubution.csv",header = T)  
x  
svs <- x %>% gather(key = Variation,value =Number,-X)
p <- ggplot(data = svs,aes(x=X,y=Number,fill=Variation))+geom_bar(stat = "identity") 
p
p+labs(title ="SV Distribution",x="Chromosome Number",y="SV Numbers") 

#====================================
#        4.4 饼图                     #
#====================================
m <- read.table("Rdata/Species.txt");
y <- paste(m[,1],m[,2])
x <- data.frame(name=y,values=m$V3/sum(m$V3))
p <- ggplot(data = x,aes(x = "",y=values,fill=name))+geom_bar(stat = "identity",width = 1)+guides(fill=FALSE)
p
p+coord_polar(theta = 'y')+labs(x = '', y = '', title = '')

#====================================
#        4.5 箱线图                   #
#====================================
head(ToothGrowth)
ToothGrowth$dose <- as.factor(ToothGrowth$dose)
#按提供药物种类分组
ggplot(data = ToothGrowth,aes(x=supp,y=len,fill=supp))+geom_boxplot()
#按剂量分组
ggplot(data = ToothGrowth,aes(x=dose,y=len,group=dose,fill=dose))+geom_boxplot()
#两组
ggplot(data = ToothGrowth,aes(x=dose,y=len,group=supp:dose,fill=supp:dose))+geom_boxplot()

#box图加抖动点
ggplot(data = ToothGrowth,aes(x=supp,y=len,fill=supp))+geom_boxplot()+geom_jitter(aes(x=supp,y=len),width = 0.1)
ggplot(data = ToothGrowth,aes(x=dose,y=len,group=dose,fill=dose))+geom_boxplot()+geom_jitter(aes(x=dose,y=len),width = 0.1)

#====================================
#       4.6 小提琴图                #
#====================================
#小提琴图
ggplot(data = ToothGrowth,aes(x=supp,y=len,fill=supp))+geom_violin()
#小提琴图+箱线图
ggplot(data = ToothGrowth,aes(x=supp,y=len,fill=supp))+geom_violin()+geom_boxplot(width=0.1,fill="white")
#两组小提琴图
ggplot(data = ToothGrowth,aes(x=dose,y=len,group=supp:dose,fill=supp:dose))+geom_violin()


#====================================
#        5 修改主题                 #
#====================================
#5.1 修改默认主题
p <- ggplot(mtcars, aes(x=cyl, y=mpg, fill=cyl)) +geom_boxplot()
p
p+theme_bw()
p+theme_classic()
p+theme_void()
p+theme_light()
p+theme_linedraw()

#5.2 自定义主题
p+theme(panel.background = element_blank())
p+theme(panel.background = element_rect(colour = "red"))
p+theme(panel.grid.major = element_line(colour = "blue"))
p+theme(panel.grid.minor = element_line(colour = "green"))


#5.3 去除图例三种方法
p <- ggplot(PlantGrowth,aes(x=group,y=weight,fill=group))+geom_boxplot()
p
#方法一：使用guides()
p+guides(fill=FALSE)
#方法二
p+theme(legend.position = "none")
#方法三：
p+scale_fill_discrete(guide=FALSE)


#====================================
#        6 修改颜色                 #
#====================================
#离散型
mtcars$cyl <- as.factor(mtcars$cyl)
p <- ggplot(mtcars, aes(x=cyl, y=mpg, fill=cyl)) +geom_boxplot()
p+scale_fill_brewer(palette = "Set2")
p+scale_fill_manual(values = c("red","green","blue"))

#连续型
p <- ggplot(mtcars, aes(x=wt, y=mpg,color=mpg)) +geom_point()
#两种渐变色
p+scale_color_gradient(low = "yellow",high = "red")
#三种渐变色
p+scale_color_gradient2(low = "yellow",mid = "orange",high = "red")

#====================================
#        科学文献配色              #
#====================================
#install.packages("ggsci")
library(ggsci)
mtcars$cyl <- as.factor(mtcars$cyl)
p <- ggplot(mtcars, aes(x=cyl, y=mpg, fill=cyl)) +geom_boxplot()
p
#查看帮助文档
help(package="ggsci")
p+scale_fill_aaas() 
p+scale_fill_npg() 
p+scale_fill_nejm() 
p+scale_fill_jama()
p+scale_fill_lancet()

#====================================
#        常见期刊主题               #
#====================================
#install.packages("ggthemes")
library(ggthemes)
mtcars$cyl <- as.factor(mtcars$cyl)
p <- ggplot(mtcars, aes(x=cyl, y=mpg, fill=cyl)) +geom_boxplot()
p
p+theme_economist()
p+theme_wsj()
p+theme_excel()
p+theme_stata()

#====================================
#        7 坐标系                   #
#====================================
x <- read.csv("Rdata/sv_distrubution.csv",header = T)  
x  
svs <- x %>% gather(key = Variation,value =Number,-X)
p <- ggplot(data = svs,aes(x=X,y=Number,fill=Variation))+geom_bar(stat = "identity")
p
#交换xy轴
p+coord_flip()
#修改极坐标：玫瑰风向图
p+coord_polar()
p+coord_polar()+guides(fill=FALSE)

#====================================
#        8 分面                      #
#====================================
library(gcookbook)
p <- ggplot(mpg,aes(x=displ,y=hwy,color=drv))+geom_point(size=3)
#横向分面
p+facet_grid(drv ~ .)

#纵向分面
p+facet_grid(.~cyl)
#两个条件分面
p+facet_grid(drv ~ cyl)

#facet_wrap()分面
p+facet_wrap( ~ class)
p+facet_wrap(~ class,nrow = 2)

#使用不同坐标
ggplot(mpg,aes(x=displ,y=hwy))+geom_point()+facet_grid(drv ~ cyl,scales = "free_y")

#====================================
#        9 布局与组合               #
#====================================
#9.1 利用ggExtra布局
library(ggplot2)
library(ggExtra)
head(mtcars)

# classic plot :
p <- ggplot(mtcars, aes(x=wt, y=mpg, color=cyl, size=cyl)) +geom_point() + theme(legend.position="none")

# Set relative size of marginal plots (main plot 10x bigger than marginals)
p1 <- ggMarginal(p, type="histogram", size=10)
p1
# Custom marginal plots:
p2 <- ggMarginal(p, type="histogram", fill = "slateblue", xparams = list(  bins=10))
p2
# Show only marginal plot for x axis
p3 <- ggMarginal(p, margins = 'x', color="purple", size=4)
p3

#9.2 利用gridExtra组合图形
# libraries
library(ggplot2)
library(gridExtra)

# Make 3 simple graphics:
g1 <- ggplot(mtcars, aes(x=qsec)) + geom_density(fill="slateblue")
g2 <- ggplot(mtcars, aes(x=drat, y=qsec, color=cyl)) + geom_point(size=5) + theme(legend.position="none")
g3 <- ggplot(mtcars, aes(x=factor(cyl), y=qsec, fill=cyl)) + geom_boxplot() + theme(legend.position="none")
g4 <- ggplot(mtcars , aes(x=factor(cyl), fill=factor(cyl))) +  geom_bar()

# Plots
grid.arrange(g2, arrangeGrob(g3, g4, ncol=2), nrow = 2)
grid.arrange(g1, g2, g3, nrow = 3)
grid.arrange(g2, arrangeGrob(g3, g4, ncol=2), nrow = 1)
grid.arrange(g2, arrangeGrob(g3, g4, nrow=2), nrow = 1)


#====================================
#        10 更多内容               #
#====================================
#10.1 ggpubr绘图
#install.packages("ggpubr")
library(ggpubr)
help(package="ggpubr")
?ggviolin
ggviolin(ToothGrowth, x = "dose", y = "len",add = "jitter", shape = "dose")
ggviolin(ToothGrowth, "dose", "len", fill = "dose",palette = c("#00AFBB", "#E7B800", "#FC4E07"),
         add = "boxplot", add.params = list(fill = "white"))

#10.2 ploty交互式绘图
#install.packages('plotly')
library(plotly)
#查看版本
packageVersion('plotly')
p <- ggplot(mtcars, aes(x=cyl, y=mpg, fill=cyl)) +geom_boxplot()
ggplotly(p)

#====================================
#       八十四、科学文献绘图        #
#====================================

#安装R包
install.packages(c("Rwordseg","wordcloud2"))
install.packages("pheatmap")
install.packages("qqman") 
install.packages("maps")
install.packages("factoextra")
install.packages("circlize")
install.packages(c("FactoMineR", "ggfortify"))
install.packages("Rcircos")

#切换工作目录
#setwd("C:/Users/xxx/Desktop/Rcourse")


#====================================
#        1 星云图                   #
#====================================
#install.packages(c("Rwordseg","wordcloud2"))

library(Rwordseg)
library(wordcloud2)

#读入文件
x <- readLines("zfgz.txt",encoding = 'UTF-8',)
head(x)
#开始分词,可以使用system.time()函数计时
y <- segmentCN(strwords = x,analyzer = "hmm",returnType = "vector")
#查看数据是否乱码
y[1:3]
#system.time(y <- segmentCN(strwords = x,analyzer = "hmm",returnType = "vector"))
#拆分列表为向量
y <- unlist(y)
#过滤数字
y <- y[!grepl('[0-9]',y)]
#过滤空白以及单个词
y <- y[nchar(y)>=2]
#统计频数
table(y)
#排序获得前50个关键字
top50 <- sort(table(y),decreasing = TRUE)[1:50]
top50
#绘图
wordcloud2(top50)
#修改形状和配色
wordcloud2(top50,shape = "star",color = rep_len(c("red","darkred"),length(top50)))

#====================================
#        2 相关性图                 #
#====================================
install.packages("corrplot")
library(corrplot)
help(package="corrplot")
corrplot(as.matrix(mtcars),is.corr = F)
M <- cor(mtcars)
corrplot(M)
# method = c("circle", "square", "ellipse", "number", "shade", "color", "pie"),
corrplot(M,method = "ellipse")
corrplot(M,method = "shade")
corrplot(M,method = "color")
corrplot(M,method = "shade")
corrplot(M,method = "circle")
corrplot(M,method = "pie",type = "upper")
corrplot(M,method = "circle",type = "lower")
corrplot(M, order = "AOE", type = "upper", tl.pos = "d")

corrplot(M, p.mat = res1$p, insig = "p-value")

#====================================
#        3、曼哈顿图                #
#====================================
#install.packages("qqman")  
library(qqman)  
library(RColorBrewer)  
str(gwasResults)  
head(gwasResults)  
manhattan(gwasResults)  
manhattan(gwasResults, main = "Manhattan Plot", ylim = c(0, 6), cex = 0.6,cex.axis = 0.9, col = c("blue4", "orange3"), suggestiveline =  
            F, genomewideline = F,chrlabs = c(1:20, "P", "Q"))  
unique(gwasResults$CHR)  
number <- length(unique(gwasResults$CHR))  
yanse <- brewer.pal(n = 4,name = "Set1")  
manhattan(gwasResults, main = "Manhattan Plot", ylim = c(0, 6), cex = 0.6,cex.axis = 0.9, col = yanse, suggestiveline =    F, genomewideline = F,chrlabs = c(1:20, "P", "Q"))  

manhattan(subset(gwasResults,CHR==3))  
#高亮显示部分SNP结果  
snpsOfInterest  
manhattan(gwasResults, highlight = snpsOfInterest)  

#注释SNP结果  
manhattan(gwasResults, annotatePval = 0.001)  
manhattan(gwasResults, annotatePval = 0.001, annotateTop = FALSE)  
#更多内容可以查看manhattan与qqman的帮助文档  
help("manhattan")  
vignette("qqman") 

#====================================
#        4、地图                    #
#====================================
#install.packages("maps")
library(maps)
library(mapdata)
library(ggplot2)
world_map <- map_data("world")
#查看全部区域
world_map$region
unique(world_map$region)
sort(unique(world_map$region))
#获取某一地区地图数据
states_map <- map_data("state")
#绘制地图，使用geom_polygon()或者geom_path()
ggplot(states_map,aes(x=long,y=lat,group=group))+geom_polygon(fill="white",color="black")+coord_map("mercator")

ggplot(states_map,aes(x=long,y=lat,group=group))+geom_path()+coord_map("mercator")


#绘制中国地图
china <- world_map[world_map$region=="China",]
ggplot(china,aes(x=long,y=lat,group=group))+geom_polygon(fill="white",color="black")
head(china)
china <- subset(x = world_map,subset = region==c("China","Taiwan"))
china
ggplot(china,aes(x=long,y=lat,group=group))+geom_polygon(fill="white",color="black")

#mapdata包中的worldHires提供高分辨率地图数据
install.packages("mapdata")
library(mapdata)
map('worldHires', col=1:10)
map('worldHires', 'China')
map_data(map = "china")
china <- map_data(map = "china")
ggplot(china,aes(x=long,y=lat,group=group,fill=region))+geom_polygon(color="black")
m_polygon(color="black")
ggplot(china,aes(x=long,y=lat,group=group))+geom_polygon()
ggplot(china,aes(x=long,y=lat,group=group))+geom_polygon()+coord_map("mercator")

#数据映射到地图
crimes <- data.frame(state=tolower(rownames(USArrests)),USArrests)
crimes
states_map <- map_data(map="state")
crime_map <- merge(states_map,crimes,by.x = "region",by.y = "state")
library(dplyr)
dplyr::arrange(crime_map,group,order)
crime_map <-  dplyr::arrange(crime_map,group,order)
ggplot(crime_map,aes(x=long,y=lat,group=group,fill=Assault))+geom_polygon(color="black")+coord_map("mercator")


#====================================
#        5、树形图                      #
#====================================
#install.packages("factoextra")
library(factoextra)
dd <- dist(mtcars,method = "euclidean")
dd
hc <- hclust(dd,method = "ward.D2")
plot(hc)
fviz_dend(hc)
fviz_dend(hc,k=4)
fviz_dend(hc,k=4,cex = 0.8,k_colors = rainbow(4))
fviz_dend(hc,k=4,cex = 0.8,k_colors = rainbow(4),ccolor_labels_by_k = FALSE,rect_border = rainbow(4))
fviz_dend(hc,k=4,cex = 0.8,ccolor_labels_by_k = FALSE,rect_border = rainbow(4),rect = TRUE,rect_fill = TRUE)
fviz_dend(hc,k=4,horiz = TRUE,type = c("circular"))
fviz_dend(hc,k=4,horiz = TRUE,type = c("rectangle"))
fviz_dend(hc,k=4,horiz = TRUE,type = c("phylogenic"))


#====================================
#        6、基因组圈图                     #
#====================================
install.packages("circlize")
library(circlize)
library(RColorBrewer)
help(package="circlize")
set.seed(999)
mat<-matrix(sample(18, 18), 3, 6)
rownames(mat) <- paste0("S", 1:3)
colnames(mat) <- paste0("E", 1:6)
df<- data.frame(from = rep(rownames(mat), times = ncol(mat)),
                to = rep(colnames(mat), each = nrow(mat)),
                value = as.vector(mat),
                stringsAsFactors = FALSE)

chordDiagram(df,grid.col = brewer.pal(9,"Set1")[1:9],link.border="grey")
circos.clear()


chordDiagram(mat,grid.col = brewer.pal(9,"Set1")[1:9],link.border="grey")
circos.clear()




