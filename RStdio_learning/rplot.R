#====================================
#       ��ʮ����R���Ի�����ͼ        #
#====================================
#��װR��
install.packages("plotrix")
install.packages("vioplot")
install.packages("venn")
install.packages("RColorBrewer")
#�л�����Ŀ¼
#setwd("C:/Users/xxx/Desktop/Rcourse")

#====================================
#        ��ͼ�豸                  #
#====================================
x11()
pdf()
dev.list()
dev.off(3)
dev.list()
dev.off(4)
dev.list()

#====================================
#        1����ͼ                      #
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
#        2��ֱ��ͼ                      #
#====================================
#���򳤶ȷֲ�ͼ
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
#        3������ͼ                      #
#====================================
# ʵ����������Ⱦɫ�峤�ȷֲ�ͼ
hg19_len <- read.csv(file = "Rdata/homo_length.csv",header = T,row.names = 1)
hg19_24 <- hg19_len[1:24,]   
barplot(height = hg19_24)  
barplot(height = hg19_24,names.arg = names(hg19_24))  
library(RColorBrewer)  
cols <- brewer.pal(n = 6,name = "Set1")  
barplot(height = hg19_24,names.arg = names(hg19_24),col = cols) 
#���Ʒ�������ͼ
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
#        4����ͼ                     #
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

#3D��ͼ
#install.packages("plotrix")
library(plotrix)
pie3D(x,col=rainbow(length(x)),labels = lbls)
pie3D(x,col=rainbow(length(x)),labels = lbls,cex=0.8)
pieplot <- pie3D(x,col=rainbow(length(x)),radius = 1,explode = 0.1)
pie3D.labels(pieplot,labels = lbls,labelcex = 0.8,height = 0.1,labelrad = 1.75)

#����ͼ
fan.plot(x,col=rainbow(length(x)),labels = lbls,cex=0.8,radius = 1)
#====================================
#        5������ͼ                     #
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
#        6��С����ͼ                       #
#====================================
#install.packages("vioplot")
library(vioplot)
vioplot(len ~ supp+dose,data = ToothGrowth)
vioplot(len ~ supp+dose,data = ToothGrowth,col=rep(c("cyan","violet"),3))

#====================================
#        7��Τ��ͼ                  #
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
#        ����                 #
#====================================

#mfrow����mfcol
opar <- par(no.readonly = TRUE)
par(mfrow=c(2,2))
plot(pressure,col="red",main="Pic 1")
barplot(table(mtcars$cyl),col = c("red","cyan","orange"),main = "Pic 2")
hist(rivers,breaks = 30,col = "pink",main = "Pic 3")
pie(c(1,3,4,2),labels = c("A","B","C","D"),main = "Pic 4")

#layout����
layout(matrix(c(1,2)),heights = c(2,1))
layout.show(2)
#��ͼ
barplot(table(mtcars$cyl),col = c("red","cyan","orange"),main = "Pic 2")
plot(pressure,col="red",type="l",main="Pic 1")
layout(matrix(c(0,2,0,0,1,3),2,3,byrow=T),widths = c(0.5,3,1),heights =  c(1,3,0.5),TRUE);
layout.show(3)
plot(pressure,col="red",main="Pic 1")
barplot(table(mtcars$cyl),col = c("red","cyan","orange"),main = "Pic 2")
hist(rivers,breaks = 30,col = "pink",main = "Pic 3")
#�ָ�
par(opar)


#====================================
#       ��ʮ����ggplot2        #
#====================================

#��������Ҫ��װ��չ��
install.packages("ggplot2")
install.packages("ggExtra")
install.packages("gcookbook")
install.packages(c("carData","gridExtra"))
install.packages("ggsci")
install.packages("ggpubr")
install.packages("ggthemes")

#====================================
#        ggplot2����                #
#====================================
#���ذ�
library(ggplot2)
plot(mtcars$wt,mtcars$mpg)
#1����
ggplot(data=mtcars)
#2ӳ��
ggplot(data=mtcars, mapping = aes(x=wt, y=mpg)) 
#3����ͼ��
ggplot(data=mtcars, mapping = aes(x=wt, y=mpg)) + geom_point() 
ggplot(data=mtcars, mapping = aes(x=wt, y=mpg,size=mpg)) + geom_point() 
ggplot(data=mtcars, mapping = aes(x=wt, y=mpg,size=mpg,shape=cyl)) + geom_point() 
ggplot(data=mtcars, mapping = aes(x=wt, y=mpg,size=mpg,shape=cyl,color=cyl)) + geom_point() 

#4��ߣ�Scale��
ggplot(data=mtcars, mapping = aes(x=wt, y=mpg,color=mpg)) + geom_point()+
  scale_color_gradient(low = "orange",high = "red")
#5ͳ�Ʊ任��Statistics��
ggplot(data=mtcars, mapping = aes(x=wt, y=mpg)) + geom_point()+
  stat_smooth( method = 'loess' ,formula =  'y ~ x')

#6���꣨Coordinate�� 
ggplot(data=mtcars, mapping = aes(x=wt, y=mpg)) + geom_point()+coord_flip()

#7ͼ�㣨Layer��
ggplot(data=mtcars, mapping = aes(x=cyl, y=mpg)) + geom_point()+geom_boxplot()

#8����
ggplot(data=mtcars, mapping = aes(x=wt, y=mpg)) + geom_point()+facet_grid(. ~ cyl)

#9���⣨Theme��
ggplot(data=mtcars, mapping = aes(x=wt, y=mpg)) + geom_point()+theme_bw()

#�����ͼ
p <- ggplot(data=mtcars, mapping = aes(x=wt, y=mpg)) + geom_point()+theme_bw()
ggsave(filename = "mtcars.pdf",plot = p)

#====================================
#        4.1 ��ͼ                      #
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
#       4.2 ֱ��ͼ                      #
#====================================
#���򳤶ȷֲ�ͼ
x <- read.table("Rdata/H37Rv.gff",sep = "\t",header = F,skip = 7,quote = "")  
x <- x[x$V3=="gene",]  
x <- abs(x$V5-x$V4+1)  
length(x)  
range(x)  
ggplot(data = NULL,aes(x=x))
ggplot(data = NULL,aes(x=x))+geom_histogram(bins = 80)
ggplot(data = NULL,aes(x=x))+geom_histogram(bins = 80)+geom_rug()


#====================================
#        4.3 ����ͼ                      #
#====================================
# ʵ����������Ⱦɫ�峤�ȷֲ�ͼ
hg19_len <- read.csv(file = "Rdata/homo_length.csv",header = T)
x <- hg19_len[1:24,]   

ggplot(data = x,aes(x=chr,y=length,fill=chr))+geom_bar(stat = "identity")
p <- ggplot(data = x,aes(x=chr,y=length,fill=chr))+geom_bar(stat = "identity")
p+scale_x_discrete(limits=x$chr)
p+scale_x_discrete(limits=x$chr)+coord_flip()
p+scale_x_discrete(limits=x$chr)+coord_flip()+guides(fill=FALSE)

#���Ʒ�������ͼ
library(tidyverse)
x <- read.csv("Rdata/sv_distrubution.csv",header = T)  
x  
svs <- x %>% gather(key = Variation,value =Number,-X)
p <- ggplot(data = svs,aes(x=X,y=Number,fill=Variation))+geom_bar(stat = "identity") 
p
p+labs(title ="SV Distribution",x="Chromosome Number",y="SV Numbers") 

#====================================
#        4.4 ��ͼ                     #
#====================================
m <- read.table("Rdata/Species.txt");
y <- paste(m[,1],m[,2])
x <- data.frame(name=y,values=m$V3/sum(m$V3))
p <- ggplot(data = x,aes(x = "",y=values,fill=name))+geom_bar(stat = "identity",width = 1)+guides(fill=FALSE)
p
p+coord_polar(theta = 'y')+labs(x = '', y = '', title = '')

#====================================
#        4.5 ����ͼ                   #
#====================================
head(ToothGrowth)
ToothGrowth$dose <- as.factor(ToothGrowth$dose)
#���ṩҩ���������
ggplot(data = ToothGrowth,aes(x=supp,y=len,fill=supp))+geom_boxplot()
#����������
ggplot(data = ToothGrowth,aes(x=dose,y=len,group=dose,fill=dose))+geom_boxplot()
#����
ggplot(data = ToothGrowth,aes(x=dose,y=len,group=supp:dose,fill=supp:dose))+geom_boxplot()

#boxͼ�Ӷ�����
ggplot(data = ToothGrowth,aes(x=supp,y=len,fill=supp))+geom_boxplot()+geom_jitter(aes(x=supp,y=len),width = 0.1)
ggplot(data = ToothGrowth,aes(x=dose,y=len,group=dose,fill=dose))+geom_boxplot()+geom_jitter(aes(x=dose,y=len),width = 0.1)

#====================================
#       4.6 С����ͼ                #
#====================================
#С����ͼ
ggplot(data = ToothGrowth,aes(x=supp,y=len,fill=supp))+geom_violin()
#С����ͼ+����ͼ
ggplot(data = ToothGrowth,aes(x=supp,y=len,fill=supp))+geom_violin()+geom_boxplot(width=0.1,fill="white")
#����С����ͼ
ggplot(data = ToothGrowth,aes(x=dose,y=len,group=supp:dose,fill=supp:dose))+geom_violin()


#====================================
#        5 �޸�����                 #
#====================================
#5.1 �޸�Ĭ������
p <- ggplot(mtcars, aes(x=cyl, y=mpg, fill=cyl)) +geom_boxplot()
p
p+theme_bw()
p+theme_classic()
p+theme_void()
p+theme_light()
p+theme_linedraw()

#5.2 �Զ�������
p+theme(panel.background = element_blank())
p+theme(panel.background = element_rect(colour = "red"))
p+theme(panel.grid.major = element_line(colour = "blue"))
p+theme(panel.grid.minor = element_line(colour = "green"))


#5.3 ȥ��ͼ�����ַ���
p <- ggplot(PlantGrowth,aes(x=group,y=weight,fill=group))+geom_boxplot()
p
#����һ��ʹ��guides()
p+guides(fill=FALSE)
#������
p+theme(legend.position = "none")
#��������
p+scale_fill_discrete(guide=FALSE)


#====================================
#        6 �޸���ɫ                 #
#====================================
#��ɢ��
mtcars$cyl <- as.factor(mtcars$cyl)
p <- ggplot(mtcars, aes(x=cyl, y=mpg, fill=cyl)) +geom_boxplot()
p+scale_fill_brewer(palette = "Set2")
p+scale_fill_manual(values = c("red","green","blue"))

#������
p <- ggplot(mtcars, aes(x=wt, y=mpg,color=mpg)) +geom_point()
#���ֽ���ɫ
p+scale_color_gradient(low = "yellow",high = "red")
#���ֽ���ɫ
p+scale_color_gradient2(low = "yellow",mid = "orange",high = "red")

#====================================
#        ��ѧ������ɫ              #
#====================================
#install.packages("ggsci")
library(ggsci)
mtcars$cyl <- as.factor(mtcars$cyl)
p <- ggplot(mtcars, aes(x=cyl, y=mpg, fill=cyl)) +geom_boxplot()
p
#�鿴�����ĵ�
help(package="ggsci")
p+scale_fill_aaas() 
p+scale_fill_npg() 
p+scale_fill_nejm() 
p+scale_fill_jama()
p+scale_fill_lancet()

#====================================
#        �����ڿ�����               #
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
#        7 ����ϵ                   #
#====================================
x <- read.csv("Rdata/sv_distrubution.csv",header = T)  
x  
svs <- x %>% gather(key = Variation,value =Number,-X)
p <- ggplot(data = svs,aes(x=X,y=Number,fill=Variation))+geom_bar(stat = "identity")
p
#����xy��
p+coord_flip()
#�޸ļ����꣺õ�����ͼ
p+coord_polar()
p+coord_polar()+guides(fill=FALSE)

#====================================
#        8 ����                      #
#====================================
library(gcookbook)
p <- ggplot(mpg,aes(x=displ,y=hwy,color=drv))+geom_point(size=3)
#�������
p+facet_grid(drv ~ .)

#�������
p+facet_grid(.~cyl)
#������������
p+facet_grid(drv ~ cyl)

#facet_wrap()����
p+facet_wrap( ~ class)
p+facet_wrap(~ class,nrow = 2)

#ʹ�ò�ͬ����
ggplot(mpg,aes(x=displ,y=hwy))+geom_point()+facet_grid(drv ~ cyl,scales = "free_y")

#====================================
#        9 ���������               #
#====================================
#9.1 ����ggExtra����
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

#9.2 ����gridExtra���ͼ��
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
#        10 ��������               #
#====================================
#10.1 ggpubr��ͼ
#install.packages("ggpubr")
library(ggpubr)
help(package="ggpubr")
?ggviolin
ggviolin(ToothGrowth, x = "dose", y = "len",add = "jitter", shape = "dose")
ggviolin(ToothGrowth, "dose", "len", fill = "dose",palette = c("#00AFBB", "#E7B800", "#FC4E07"),
         add = "boxplot", add.params = list(fill = "white"))

#10.2 ploty����ʽ��ͼ
#install.packages('plotly')
library(plotly)
#�鿴�汾
packageVersion('plotly')
p <- ggplot(mtcars, aes(x=cyl, y=mpg, fill=cyl)) +geom_boxplot()
ggplotly(p)

#====================================
#       ��ʮ�ġ���ѧ���׻�ͼ        #
#====================================

#��װR��
install.packages(c("Rwordseg","wordcloud2"))
install.packages("pheatmap")
install.packages("qqman") 
install.packages("maps")
install.packages("factoextra")
install.packages("circlize")
install.packages(c("FactoMineR", "ggfortify"))
install.packages("Rcircos")

#�л�����Ŀ¼
#setwd("C:/Users/xxx/Desktop/Rcourse")


#====================================
#        1 ����ͼ                   #
#====================================
#install.packages(c("Rwordseg","wordcloud2"))

library(Rwordseg)
library(wordcloud2)

#�����ļ�
x <- readLines("zfgz.txt",encoding = 'UTF-8',)
head(x)
#��ʼ�ִ�,����ʹ��system.time()������ʱ
y <- segmentCN(strwords = x,analyzer = "hmm",returnType = "vector")
#�鿴�����Ƿ�����
y[1:3]
#system.time(y <- segmentCN(strwords = x,analyzer = "hmm",returnType = "vector"))
#����б�Ϊ����
y <- unlist(y)
#��������
y <- y[!grepl('[0-9]',y)]
#���˿հ��Լ�������
y <- y[nchar(y)>=2]
#ͳ��Ƶ��
table(y)
#������ǰ50���ؼ���
top50 <- sort(table(y),decreasing = TRUE)[1:50]
top50
#��ͼ
wordcloud2(top50)
#�޸���״����ɫ
wordcloud2(top50,shape = "star",color = rep_len(c("red","darkred"),length(top50)))

#====================================
#        2 �����ͼ                 #
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
#        3��������ͼ                #
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
#������ʾ����SNP���  
snpsOfInterest  
manhattan(gwasResults, highlight = snpsOfInterest)  

#ע��SNP���  
manhattan(gwasResults, annotatePval = 0.001)  
manhattan(gwasResults, annotatePval = 0.001, annotateTop = FALSE)  
#�������ݿ��Բ鿴manhattan��qqman�İ����ĵ�  
help("manhattan")  
vignette("qqman") 

#====================================
#        4����ͼ                    #
#====================================
#install.packages("maps")
library(maps)
library(mapdata)
library(ggplot2)
world_map <- map_data("world")
#�鿴ȫ������
world_map$region
unique(world_map$region)
sort(unique(world_map$region))
#��ȡĳһ������ͼ����
states_map <- map_data("state")
#���Ƶ�ͼ��ʹ��geom_polygon()����geom_path()
ggplot(states_map,aes(x=long,y=lat,group=group))+geom_polygon(fill="white",color="black")+coord_map("mercator")

ggplot(states_map,aes(x=long,y=lat,group=group))+geom_path()+coord_map("mercator")


#�����й���ͼ
china <- world_map[world_map$region=="China",]
ggplot(china,aes(x=long,y=lat,group=group))+geom_polygon(fill="white",color="black")
head(china)
china <- subset(x = world_map,subset = region==c("China","Taiwan"))
china
ggplot(china,aes(x=long,y=lat,group=group))+geom_polygon(fill="white",color="black")

#mapdata���е�worldHires�ṩ�߷ֱ��ʵ�ͼ����
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

#����ӳ�䵽��ͼ
crimes <- data.frame(state=tolower(rownames(USArrests)),USArrests)
crimes
states_map <- map_data(map="state")
crime_map <- merge(states_map,crimes,by.x = "region",by.y = "state")
library(dplyr)
dplyr::arrange(crime_map,group,order)
crime_map <-  dplyr::arrange(crime_map,group,order)
ggplot(crime_map,aes(x=long,y=lat,group=group,fill=Assault))+geom_polygon(color="black")+coord_map("mercator")


#====================================
#        5������ͼ                      #
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
#        6��������Ȧͼ                     #
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



