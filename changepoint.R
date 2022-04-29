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
f <- read.table("Bmal_f.fastq.gz.depth.50k-win.txt",header = TRUE,sep="",stringsAsFactors = FALSE)
m <- read.table("Bmal_m.fastq.gz.depth.50k-win.txt",header = TRUE,sep="",stringsAsFactors = FALSE)
head(f)
head(m)
# coverage calculation
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
##
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

