library(lattice)
library(plyr)
library(Rmisc)
library(ggplot2)
library(RColorBrewer)

theme <- theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = 'white', colour = 'white'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

setwd("/Users/yifeng/Library/Mobile Documents/com~apple~CloudDocs/project/81worm/ragoo_assembly/snpdata/fe_m/")
f <- read.table("Sman_f.fastx.fastq.gz.depth.50k-win.txt",header = TRUE,sep="",stringsAsFactors = FALSE)
m <- read.table("Sman_m.fastx.fastq.gz.depth.50k-win.txt",header = TRUE,sep="",stringsAsFactors = FALSE)
f <- read.table("Sman_f.fastx.fastq.gz.depth.5k-win.txt",header = TRUE,sep="",stringsAsFactors = FALSE)
m <- read.table("Sman_m.fastx.fastq.gz.depth.5k-win.txt",header = TRUE,sep="",stringsAsFactors = FALSE)
head(f)
head(m)
# 用一个常染色体片段计算测序深度d
f_auto <- f[f$scaf=="SM_V7_1",]
m_auto <- m[m$scaf=="SM_V7_1",]
f_mapsum <- as.numeric(f_auto$sum)
mapsum_f <- sum(f_mapsum)
m_mapsum <- as.numeric(m_auto$sum)
mapsum_m <- sum(m_mapsum)
mapsum_m
mapsum_f
dm <- mapsum_m/(2*87785169) #88881357-1096188=87785169 TMUE_LG2,35622338bp是这个scaffold除去N的长度
df <- mapsum_f/(2*87785169)
dm#22.31994
df#22.69224
f$coverage <- f$sum/(f$count*df)
m$coverage <- m$sum/(m$count*dm)
f$coverage[ f$coverage=="NaN" ]<- 0
m$coverage[ m$coverage=="NaN" ]<- 0
f$ratio <- log2((m$coverage+0.01)/(f$coverage+0.01))
f$pos<-(f$start+f$end)/2
head(f)
pSmancov<-ggplot(f[f$scaf=="SM_V7_ZW",],aes(x=pos,y=ratio))+geom_point(size=0.01,alpha=0.2)+
  geom_hline(aes(yintercept=0),linetype="dashed",size = 0.2,colour="red")+
  #geom_smooth(fill = NA, size=0.2,colour="black")+
  labs(y="Log2(F/M coverage)")+scale_y_continuous(limits = c(-1,1))+scale_x_continuous(limits = c(0,88385488))+
  theme+theme(axis.title.x = element_text(size = 0))+theme(axis.text.x = element_text(size = 0))+theme(axis.title.y = element_text(size = 0))+theme(axis.text.y = element_text(size = 0))
pSmancov

###change point test,
#install.packages("cpm") #初次使用需安装，以后就不需要了
library(cpm)
fz<-f[f$scaf=="SM_V7_ZW",c(2,9)]
names(fz) <- c("x", "y")
par(mfrow = c(2,1))
plot(fz1, type = "l", col = "steelblue", lwd = 2)
shapiro.test(fz1$y)#检验数据是否服从高斯分布，原假设H0是数据服从正态分布，这里p<0.05，拒绝原假设，数据服从正态分布不成立。发现不服从。所以选择一个非高斯分布的方法使用
cpm.res = processStream(fz1$y, cpmType = "Kolmogorov-Smirnov")
# 可视化变点
plot(fz1, type = "l", col = "steelblue", lwd = 2)
abline(v = cpm.res$changePoints, lwd = 3.5, col = "red")
# 变点坐标信息提取
print(cpm.res$changePoints)
########################
library(changepoint)
results <- cpt.mean(fz1$y,method="PELT",penalty = "AIC",Q=3)
cpts(results)
param.est(results)
plot(results,cpt.col="red",xlab="Index",cpt.width=1)

################################################################################
f <- read.table("Bmal_f.fastq.gz.depth.50k-win.txt",header = TRUE,sep="",stringsAsFactors = FALSE)
m <- read.table("Bmal_m.fastq.gz.depth.50k-win.txt",header = TRUE,sep="",stringsAsFactors = FALSE)
head(f)
head(m)
# 用一个常染色体片段计算测序深度d
f_auto <- f[f$scaf=="Bm_v4_Chr1_scaffold_001",]
m_auto <- m[m$scaf=="Bm_v4_Chr1_scaffold_001",]
f_mapsum <- as.numeric(f_auto$sum)
mapsum_f <- sum(f_mapsum)
m_mapsum <- as.numeric(m_auto$sum)
mapsum_m <- sum(m_mapsum)
mapsum_m
mapsum_f
dm <- mapsum_m/(2*14590985) #14701151-110166=14590985 Bm_v4_Chr1_scaffold_001,14590985 bp是这个scaffold除去N的长度
df <- mapsum_f/(2*14590985)
dm#50.76692
df#71.98064
f$coverage <- f$sum/(f$count*df)
m$coverage <- m$sum/(m$count*dm)
f$ratio <- log2((f$coverage+0.01)/(m$coverage+0.01))
f$pos<-(f$start+f$end)/2
head(f)
pBmalcov<-ggplot(f[f$scaf=="Bm_v4_ChrX_scaffold_001",],aes(x=pos,y=ratio))+
  geom_line(size=0.3,alpha=0.5,color="#6a51a3")+
  #geom_point(size=0.01,alpha=0.3)+
  geom_hline(aes(yintercept=0),linetype="dashed",size = 0.2,colour="grey")+geom_hline(aes(yintercept=1),linetype="dashed",size = 0.2,colour="grey")+geom_hline(aes(yintercept=-1),linetype="dashed",size = 0.2,colour="grey")+
  #geom_smooth(fill = NA, size=0.2,colour="black")+
  labs(y="Log2(F/M coverage)")+scale_y_continuous(limits = c(-1,1.4))+scale_x_continuous(limits = c(0,24943668))+
  theme+theme(axis.title.x = element_text(size = 0))+theme(axis.title.y = element_text(size = 0))+theme(axis.text.x = element_text(size = 0))+theme(axis.text.y = element_text(size = 0))
pBmalcov

##########change poit##############
fz1<-f[f$scaf=="Bm_v4_ChrX_scaffold_001",c(2,7)]
names(fz1) <- c("x1", "y")
fz1$x <- 1:nrow(fz1)
fz2<-fz1[c(3,2)]
#########cool!!!!
library(changepoint)
results <- cpt.mean(fz2$y,method="BinSeg",penalty = "AIC",Q=2)  ##single changepoint (AMOC) or multiple changepoints using exact (PELT or SegNeigh) or approximate (BinSeg) methods. default: AIC
cpts(results)
param.est(results)
plot(results,cpt.col="red",xlab="Index",cpt.width=1)

####mask####heatmap
fSNP <- read.table("Bmal.fe.masked.SNP.50k-win.txt",header = TRUE,sep="",stringsAsFactors = FALSE)
mSNP <- read.table("Bmal.m.masked.SNP.50k-win.txt",header = TRUE,sep="",stringsAsFactors = FALSE)
head(fSNP)
head(mSNP)
SNP<-data.frame(scaf = fSNP$scaf,start=fSNP$start,end=fSNP$end,fe=fSNP$count,m=mSNP$count)
###for change point
SNP[SNP$m+SNP$fe<20, "m"] <- 0
SNP[SNP$m+SNP$fe<20, "fe"] <- 0 
###
SNP$ratio<-log2((SNP$m+1)/(SNP$fe+1))
SNP<-SNP[SNP$scaf=="Bm_v4_ChrX_scaffold_001",]
head(SNP)
interval <- (max(SNP$ratio)-min(SNP$ratio))/200
SNP$range <- 0
SNP$range <- round((SNP$ratio-min(SNP$ratio))/interval)+1
min(SNP$range);max(SNP$range)
colfunc<-colorRampPalette(brewer.pal(9,"Reds"))
colour <- colfunc(201)
names(colour) <- c(seq(1,201))
SNP$color <- colour[as.character(SNP$range)]
#guide
plot(rep(1,201),col=colour, pch=15,cex=2)
plot(rep(1,201),col=SNP$color, pch=15,cex=2)
###x轴与length相关？
pSNPrepeat<-ggplot(SNP)+guides(color=FALSE)+ 
  scale_y_continuous(limits = c(0.35,0.45))+scale_x_continuous(limits = c(0,24943668))+
  geom_rect(data=SNP,aes(xmin=start,xmax=end),ymin=0.35,ymax=0.45,fill=SNP$color)###fill=SNP$color,color=SNP$color同时加会变粗；color=SNP$color放大会有黑线
pSNPrepeat<-pSNPrepeat+theme+theme(axis.text.y = element_text(size = 0))+theme(axis.text.x = element_text(size = 0))

##########change poit##############
head(SNP)
SNP[SNP$ratio<0, "ratio"] <- 0 
SNP1<-SNP[SNP$scaf=="Bm_v4_ChrX_scaffold_001",c(2,6)]
names(SNP1) <- c("x1", "y")
SNP1$x <- 1:nrow(SNP1)
SNP2<-SNP1[c(3,2)]
#################
results <- cpt.mean(SNP2$y,method="BinSeg",penalty = "AIC",Q=8)  #PELT or SegNeigh) or approximate (BinSeg) methods
results <- cpt.mean(SNP2$y,method="BinSeg",penalty = "Manual", pen.value = "1 * log(n)",Q=13)
results <- cpt.mean(SNP2$y,method="BinSeg",penalty = "AIC",Q=12)
cpts(results)  #0.5 * log(n),pen.value="0"
param.est(results)
plot(results,cpt.col="red",xlab="Index",cpt.width=1)

########sig PAR S3
head(SNP1) ##two column, start and ratio
SNP1S3<-subset(SNP1,SNP1$start>=22400000 & SNP1$start<23370000)
SNP1PAR<-subset(SNP1,SNP1$start>=23370000)
SNP1S3$group<-"S3"
SNP1PAR$group<-"PAR"
SNP1S3PAR<-rbind(SNP1S3,SNP1PAR)
library(ggsignif)
ggplot(SNP1S3PAR,aes(x=group,y=ratio,fill=group))+labs(title="SNP")+
  geom_boxplot()+geom_signif(comparisons = list(c("S3", "PAR")), test = wilcox.test, step_increase = 0.2)


####################################################################################
setwd("/Users/yifeng/Desktop/81worm/ragoo_assembly/snpdata/fe_m/")
f <- read.table("Ovol.f.merge.depth.50k-win.txt",header = TRUE,sep="",stringsAsFactors = FALSE)
m <- read.table("Ovol.m.merge.depth.50k-win.txt",header = TRUE,sep="",stringsAsFactors = FALSE)
head(f)
head(m)
# 用一个常染色体片段计算测序深度d
f_auto <- f[f$scaf=="OVOC_OM1b",]
m_auto <- m[m$scaf=="OVOC_OM1b",]
f_mapsum <- as.numeric(f_auto$sum)
mapsum_f <- sum(f_mapsum)
m_mapsum <- as.numeric(m_auto$sum)
mapsum_m <- sum(m_mapsum)
mapsum_m
mapsum_f
dm <- mapsum_m/(2*27520521) #28345163-824642=27520521 OVOC_OM1b,27520521 bp是这个scaffold除去N的长度
df <- mapsum_f/(2*27520521)
dm#12.96822
df#8.905521
dm/df#1.4562
f$coverage <- f$sum/((f$end-f$start)*df)
m$coverage <- m$sum/((m$end-m$start)*dm)
f$coverage[ f$coverage=="NaN" ]<- 0
m$coverage[ m$coverage=="NaN" ]<- 0

f$ratio <- log2((f$coverage+0.01)/(m$coverage+0.01))
f$pos<-(f$start+f$end)/2
head(f)
pOvolcov<-ggplot(f[f$scaf=="OVOC_OM2",],aes(x=pos,y=ratio))+
  geom_line(size=0.3,alpha=0.5,color="#6a51a3")+
  #geom_point(size=0.01,alpha=0.3)+
  geom_hline(aes(yintercept=0),linetype="dashed",size = 0.2,colour="grey")+geom_hline(aes(yintercept=1),linetype="dashed",size = 0.2,colour="grey")+geom_hline(aes(yintercept=-1),linetype="dashed",size = 0.2,colour="grey")+
  #geom_smooth(fill = NA, size=0.2,colour="black")+
  labs(y="Log2(F/M coverage)")+scale_y_continuous(limits = c(-1,1.4))+scale_x_continuous(limits = c(0,25485961))+
  theme+theme(axis.title.x = element_text(size = 0))+theme(axis.title.y = element_text(size = 0))+theme(axis.text.x = element_text(size = 0))+theme(axis.text.y = element_text(size = 0))
pOvolcov##par 22200000
pOvolcov+scale_x_continuous(limits = c(22200000,22970000))
pdf("Ovolcov.pdf",width=9,height=1.7)
pOvolcov
dev.off()

##########change poit##############
fz1<-f[f$scaf=="OVOC_OM2",c(2,9)]
names(fz1) <- c("x1", "y")
fz1$x <- 1:nrow(fz1)
fz2<-fz1[c(3,2)]
par(mfrow = c(2,1))
plot(fz2, type = "l", col = "steelblue", lwd = 2)
shapiro.test(fz2$y)#检验数据是否服从高斯分布，原假设H0是数据服从正态分布，这里p<0.05，拒绝原假设，数据服从正态分布不成立。发现不服从。所以选择一个非高斯分布的方法使用
cpm.res = processStream(fz2$y, cpmType = "Kolmogorov-Smirnov")
# 可视化变点
plot(fz2, type = "l", col = "steelblue", lwd = 2)
abline(v = cpm.res$changePoints, lwd = 1, col = "red")
# 变点坐标信息提取
print(cpm.res$changePoints)
########################
#########cool!!!!
results <- cpt.mean(fz2$y,method="PELT")  ##single changepoint (AMOC) or multiple changepoints using exact (PELT or SegNeigh) or approximate (BinSeg) methods.
results <- cpt.mean(fz2$y,method="BinSeg",penalty = "AIC")
cpts(results)
param.est(results)
plot(results,cpt.col="red",xlab="Index",cpt.width=1)+scale_y_continuous(limits = c(-1,1.4))



##############
####heatmap
fSNP <- read.table("Ovol.f.merge.SNP.50k-win.txt",header = TRUE,sep="",stringsAsFactors = FALSE)
mSNP <- read.table("Ovol.m.merge.SNP.50k-win.txt",header = TRUE,sep="",stringsAsFactors = FALSE)
head(fSNP)
head(mSNP)
SNP<-data.frame(scaf = mSNP$scaf,start=mSNP$start,end=mSNP$end,fe=fSNP$count,m=mSNP$count)
###for change point
SNP[SNP$m+SNP$fe<20, "m"] <- 0
SNP[SNP$m+SNP$fe<20, "fe"] <- 0 
###
SNP$ratio<-log2((SNP$m+1)/(SNP$fe+1))
SNP<-SNP[SNP$scaf=="OVOC_OM2",]
head(SNP)
interval <- (max(SNP$ratio)-min(SNP$ratio))/200
SNP$range <- 0
SNP$range <- round((SNP$ratio-min(SNP$ratio))/interval)+1
min(SNP$range);max(SNP$range)
colfunc<-colorRampPalette(brewer.pal(9,"Reds"))
colour <- colfunc(201)
names(colour) <- c(seq(1,201))
SNP$color <- colour[as.character(SNP$range)]
#guide
plot(rep(1,201),col=colour, pch=15,cex=2)
plot(rep(1,201),col=SNP$color, pch=15,cex=2)
###x轴与length相关？
pSNPrepeat<-ggplot(SNP)+guides(color=FALSE)+ 
  scale_y_continuous(limits = c(0.35,0.45))+scale_x_continuous(limits = c(0,25485961))+
  geom_rect(data=SNP,aes(xmin=start,xmax=end),ymin=0.35,ymax=0.45,fill=SNP$color)
pSNPrepeat<-pSNPrepeat+theme+theme(axis.text.y = element_text(size = 0))+theme(axis.text.x = element_text(size = 0))
pSNPrepeat
pSNPrepeat+scale_x_continuous(limits = c(22800000,23490000))
pSNPrepeat+scale_x_continuous(limits = c(22200000,22970000))

pdf("Ovol-maskSNPrepeat-heatmap.pdf",width=9,height=2)
pSNPrepeat  #ratio  0-9	
dev.off()

##########change poit##############
head(SNP)
SNP[SNP$ratio<0, "ratio"] <- 0 
SNP1<-SNP[SNP$scaf=="OVOC_OM2",c(2,6)]
names(SNP1) <- c("x1", "y")
SNP1$x <- 1:nrow(SNP1)
SNP2<-SNP1[c(3,2)]
#################
results <- cpt.mean(SNP2$y,method="PELT")  #PELT or SegNeigh) or approximate (BinSeg) methods
results <- cpt.mean(SNP2$y,method="BinSeg",penalty = "AIC",Q=9)
cpts(results)
param.est(results)
plot(results,cpt.col="red",xlab="Index",cpt.width=1)



###sig S4 par
head(SNP1)
SNP1S4<-subset(SNP1,SNP1$start>=22000000 & SNP1$start<22960000)
SNP1PAR<-subset(SNP1,SNP1$start>=22960000)
SNP1S4$group<-"S4"
SNP1PAR$group<-"PAR"
SNP1S4PAR<-rbind(SNP1S4,SNP1PAR)
library(ggsignif)
ggplot(SNP1S4PAR,aes(x=group,y=ratio,fill=group))+labs(title="SNP")+
  geom_boxplot()+geom_signif(comparisons = list(c("S4", "PAR")), test = wilcox.test, step_increase = 0.2)

