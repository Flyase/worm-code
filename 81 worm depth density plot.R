#install.packages(c("ggolot2","Rmisc","lattice","plyr"))
#"dark grey"-hemo  salmon-female "forest green"-male  tan-mix
#10k window
library(lattice)
library(plyr)
library(Rmisc)
library(ggplot2)
theme <- theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = 'white', colour = 'white'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
setwd("/Users/yifeng/Library/Mobile Documents/com~apple~CloudDocs/project/81worm/depthfile/schi/")
mix <- read.table("Shae.10k-win.txt")
mix <- subset(mix, mix$V5 > 500)
mix_v4 <- as.numeric(mix$V4)
sum_mix_v4 <- sum(mix_v4)
sum_mix_v4
head(mix)
mix_v5 <- as.numeric(mix$V5)
sum_mix_v5 <- sum(mix_v5)
sum_mix_v5
mix$dmix <- mix_v4/(mix_v5*2)
dmix <- mix[c(1,8)]
#install.packages("ggplot2")

#dmix <- subset(dmix,dmix<250)
dmix$foldd <- dmix$dmix/6.5
p2 <- ggplot(dmix,aes(x=foldd))+geom_density(colour="tan")+
  labs(x="depth",title="Schistosoma haematobium")+scale_x_continuous(limits = c(0,2))+
  theme(axis.text.x = element_text(size = 6))+
  theme(axis.text.y = element_text(size = 6))+
  theme(axis.title.y = element_text(size = 0))+
  theme(axis.title.x = element_text(size = 0))+
  theme(title=element_text(size = 4.5,face = "bold.italic",color = "#e6ab02"))
p2

mix <- read.table("Sjap_1.fastq.gz.depth.10k-win.txt",header=FALSE)
mix <- subset(mix, mix$V5 > 500)
mix_v4 <- as.numeric(mix$V4)
mix_v5 <- as.numeric(mix$V5)
mix$dmix <- mix_v4/mix_v5
head(mix)
dmix <- mix[c(1,8)]
#install.packages("ggplot2")
library(ggplot2)
dmix$foldd <- dmix$dmix/29.7
p3 <- ggplot(dmix,aes(x=foldd))+geom_density(colour="salmon")+
  labs(x="depth",title="Schistosoma japonicum")+scale_x_continuous(limits = c(0,2))+  
  geom_vline(aes(xintercept = 0.632), linetype = "dashed", size = 0.4)+
  geom_vline(aes(xintercept = 0.33), linetype = "dashed", size = 0.4)+
  theme(axis.text.x = element_text(size = 6))+
  theme(axis.text.y = element_text(size = 6))+
  theme(axis.title.y = element_text(size = 0))+
  theme(axis.title.x = element_text(size = 0))+
  theme(title=element_text(size = 4.5,face = "bold.italic",color = "#e6ab02"))
p3


mix <- read.table("Sman.depth.bed.10k-win.txt")
mix <- subset(mix, mix$V5 > 500)
mix_v4 <- as.numeric(mix$V4)
sum_mix_v4 <- sum(mix_v4)
sum_mix_v4
head(mix)
mix_v5 <- as.numeric(mix$V5)
sum_mix_v5 <- sum(mix_v5)
sum_mix_v5
mix$dmix <- mix_v4/(mix_v5*2)
#install.packages("ggplot2")
library(ggplot2)
dmix <- mix[c(1,8)]

#dmix <- subset(dmix,dmix<600)
dmix$foldd <- dmix$dmix/15
p4 <- ggplot(dmix,aes(x=foldd))+geom_density(colour="salmon")+
  labs(x="depth",title="Schistosoma mansoni")+scale_x_continuous(limits = c(0,2))+  
  geom_vline(aes(xintercept = 0.2), linetype = "dashed", size = 0.4)+
  geom_vline(aes(xintercept = 0.77), linetype = "dashed", size = 0.4)+
  theme(axis.text.x = element_text(size = 6))+
  theme(axis.text.y = element_text(size = 6))+
  theme(axis.title.y = element_text(size = 0))+
  theme(axis.title.x = element_text(size = 0))+
  theme(title=element_text(size = 4.5,face = "bold.italic",color = "#e6ab02"))
p4
pdf("Smanfe.pdf",width=1.4,height=0.9)
p4
dev.off()

mix <- read.table("Smargrebowiei_1.fastq.gz.depth.10k-win.txt")
mix <- subset(mix, mix$V5 > 500)
mix_v4 <- as.numeric(mix$V4)
sum_mix_v4 <- sum(mix_v4)
sum_mix_v4
head(mix)
mix_v5 <- as.numeric(mix$V5)
sum_mix_v5 <- sum(mix_v5)
sum_mix_v5
mix$dmix <- mix_v4/(mix_v5*2)
dmix <- mix[c(1,8)]
library(ggplot2)

##dmix <- subset(dmix,dmix<130)
dmix$foldd <- dmix$dmix/25.6
dmix$sex <- "male"
p5 <- ggplot(dmix,aes(x=foldd))+geom_density(colour="forest green")+
  labs(x="depth",title="Schistosoma margrebowiei")+scale_x_continuous(limits = c(0,2))+  
  theme(axis.text.x = element_text(size = 6))+
  theme(axis.text.y = element_text(size = 6))+
  theme(axis.title.y = element_text(size = 0))+
  theme(axis.title.x = element_text(size = 0))+
  theme(title=element_text(size = 4.5,face = "bold.italic",color = "#e6ab02"))
  #legend("topleft",col='green',bty="n",horiz=TRUE)
  #scale_fill_manual(#name='??????', #?????????????????????scale_fill_discrete)
                    #labels=c('1-SL','2-SW','3-PL','4-PW'), #????????????
                    #values=c('darkgreen','red','darkblue','black'))+ #??????
                    #labels=dmix$sex, #????????????
                    #values='green')
p5


mix <- read.table("Smattheei_1.fastq.gz.depth.10k-win.txt")
mix <- subset(mix, mix$V5 > 500)
mix_v4 <- as.numeric(mix$V4)
sum_mix_v4 <- sum(mix_v4)
sum_mix_v4
head(mix)
mix_v5 <- as.numeric(mix$V5)
sum_mix_v5 <- sum(mix_v5)
sum_mix_v5
mix$dmix <- mix_v4/(mix_v5*2)
dmix <- mix[c(1,8)]
library(ggplot2)

dmix$foldd <- dmix$dmix/18.7
p6 <- ggplot(dmix,aes(x=foldd))+geom_density(colour="forest green")+
  labs(x="depth",title="Schistosoma mattheei")+scale_x_continuous(limits = c(0,2))+  
  theme(axis.text.x = element_text(size = 6))+
  theme(axis.text.y = element_text(size = 6))+
  theme(axis.title.y = element_text(size = 0))+
  theme(axis.title.x = element_text(size = 0))+
  theme(title=element_text(size = 4.5,face = "bold.italic",color = "#e6ab02"))
p6
pdf("Smatm.pdf",width=1.4,height=0.9)
p6
dev.off()

mix <- read.table("Srod_1.fastq.gz.depth.10k-win.txt")
mix <- subset(mix, mix$V5 > 500)
mix_v4 <- as.numeric(mix$V4)
sum_mix_v4 <- sum(mix_v4)
sum_mix_v4
head(mix)
mix_v5 <- as.numeric(mix$V5)
sum_mix_v5 <- sum(mix_v5)
sum_mix_v5
mix$dmix <- mix_v4/(mix_v5*2)
dmix <- mix[c(1,8)]
library(ggplot2)

dmix$foldd <- dmix$dmix/4.4
psrod <- ggplot(dmix,aes(x=foldd))+geom_density(colour="forest green")+
  labs(x="depth",title="Schistosoma rodhaini")+scale_x_continuous(limits = c(0,2))+  
  theme(axis.text.x = element_text(size = 6))+
  theme(axis.text.y = element_text(size = 6))+
  theme(axis.title.y = element_text(size = 0))+
  theme(axis.title.x = element_text(size = 0))+
  theme(title=element_text(size = 4.5,face = "bold.italic",color = "#e6ab02"))
psrod

mix <- read.table("Tregenti_1.fastq.gz.depth.10k-win.txt")
mix <- subset(mix, mix$V5 > 500)
mix_v4 <- as.numeric(mix$V4)
sum_mix_v4 <- sum(mix_v4)
sum_mix_v4
head(mix)
mix_v5 <- as.numeric(mix$V5)
sum_mix_v5 <- sum(mix_v5)
sum_mix_v5
mix$dmix <- mix_v4/(mix_v5*2)
dmix <- mix[c(1,2,3,8)]
library(ggplot2)

dmix$foldd <- dmix$dmix/7.7
p7 <- ggplot(dmix,aes(x=foldd))+geom_density(colour="tan")+
  labs(x="depth",title="Trichobilharzia regenti")+scale_x_continuous(limits = c(0,2))+
  theme(axis.text.x = element_text(size = 6))+
  theme(axis.text.y = element_text(size = 6))+
  theme(axis.title.y = element_text(size = 0))+
  theme(axis.title.x = element_text(size = 0))+
  theme(title=element_text(size = 4.5,face = "bold.italic",color = "#e6ab02"))
p7

head(dmix)
ptregdepth1<-ggplot(dmix[dmix$V1=="TRE_scaffold0000001",],aes(x=V3,y=foldd))+geom_point(size=0.2,alpha=0.4,color="tan")+
  geom_smooth(fill = NA, size=0.5,color="tan") +labs(x="coords",y="depth",title="TRE_scaffold0000001")+theme+scale_y_continuous(limits = c(0,1.3))+theme(axis.title.x = element_text(size = 0))
ptregdepth2<-ggplot(dmix[dmix$V1=="TRE_scaffold0000002",],aes(x=V3,y=foldd))+geom_point(size=0.2,alpha=0.4,color="tan")+
  geom_smooth(fill = NA, size=0.5,color="tan") +labs(x="coords",y="depth",title="TRE_scaffold0000002")+theme+scale_y_continuous(limits = c(0,1.3))+theme(axis.title.x = element_text(size = 0))
ptregdepth3<-ggplot(dmix[dmix$V1=="TRE_scaffold0000003",],aes(x=V3,y=foldd))+geom_point(size=0.2,alpha=0.4,color="tan")+
  geom_smooth(fill = NA, size=0.5,color="tan") +labs(x="coords",y="depth",title="TRE_scaffold0000003")+theme+scale_y_continuous(limits = c(0,1.3))+theme(axis.title.x = element_text(size = 0))
ptregdepth4<-ggplot(dmix[dmix$V1=="TRE_scaffold0000004",],aes(x=V3,y=foldd))+geom_point(size=0.2,alpha=0.4,color="tan")+
  geom_smooth(fill = NA, size=0.5,color="tan") +labs(x="coords",y="depth",title="TRE_scaffold0000004")+theme+scale_y_continuous(limits = c(0,1.3))+theme(axis.title.x = element_text(size = 0))
ptregdepth5<-ggplot(dmix[dmix$V1=="TRE_scaffold0000005",],aes(x=V3,y=foldd))+geom_point(size=0.2,alpha=0.4,color="tan")+
  geom_smooth(fill = NA, size=0.5,color="tan") +labs(x="coords",y="depth",title="TRE_scaffold0000005")+theme+scale_y_continuous(limits = c(0,1.3))+theme(axis.title.x = element_text(size = 0))
ptregdepth6<-ggplot(dmix[dmix$V1=="TRE_scaffold0000488",],aes(x=V3,y=foldd))+geom_point(size=0.2,alpha=0.4,color="tan")+
  geom_smooth(fill = NA, size=0.5,color="tan") +labs(x="coords",y="depth",title="TRE_scaffold0000488")+theme+scale_y_continuous(limits = c(0,1.3))+theme(axis.title.x = element_text(size = 0))
ptregdepth7<-ggplot(dmix[dmix$V1=="TRE_scaffold0001157",],aes(x=V3,y=foldd))+geom_point(size=0.2,alpha=0.4,color="tan")+
  geom_smooth(fill = NA, size=0.5,color="tan") +labs(x="coords",y="depth",title="TRE_scaffold0001157")+theme+scale_y_continuous(limits = c(0,1.3))+theme(axis.title.x = element_text(size = 0))
ptregdepth8<-ggplot(dmix[dmix$V1=="TRE_scaffold0001354",],aes(x=V3,y=foldd))+geom_point(size=0.2,alpha=0.4,color="tan")+
  geom_smooth(fill = NA, size=0.5,color="tan") +labs(x="coords",y="depth",title="TRE_scaffold0001354")+theme+scale_y_continuous(limits = c(0,1.3))+theme(axis.title.x = element_text(size = 0))
ptregdepth9<-ggplot(dmix[dmix$V1=="TRE_scaffold0002501",],aes(x=V3,y=foldd))+geom_point(size=0.2,alpha=0.4,color="tan")+
  geom_smooth(fill = NA, size=0.5,color="tan") +labs(x="coords",y="depth",title="TRE_scaffold0002501")+theme+scale_y_continuous(limits = c(0,1.3))+theme(axis.title.x = element_text(size = 0))
ptregdepth10<-ggplot(dmix[dmix$V1=="TRE_scaffold0002715",],aes(x=V3,y=foldd))+geom_point(size=0.2,alpha=0.4,color="tan")+
  geom_smooth(fill = NA, size=0.5,color="tan") +labs(x="coords",y="depth",title="TRE_scaffold0002715")+theme+scale_y_continuous(limits = c(0,1.3))+theme(axis.title.x = element_text(size = 0))
ptregdepth11<-ggplot(dmix[dmix$V1=="TRE_scaffold0003810",],aes(x=V3,y=foldd))+geom_point(size=0.2,alpha=0.4,color="tan")+
  geom_smooth(fill = NA, size=0.5,color="tan") +labs(x="coords",y="depth",title="TRE_scaffold0003810")+theme+scale_y_continuous(limits = c(0,1.3))+theme(axis.title.x = element_text(size = 0))
ptregdepth12<-ggplot(dmix[dmix$V1=="TRE_scaffold0005590",],aes(x=V3,y=foldd))+geom_point(size=0.2,alpha=0.4,color="tan")+
  geom_smooth(fill = NA, size=0.5,color="tan") +labs(x="coords",y="depth",title="TRE_scaffold0005590")+theme+scale_y_continuous(limits = c(0,1.3))+theme(axis.title.x = element_text(size = 0))
ptregdepth13<-ggplot(dmix[dmix$V1=="TRE_scaffold0010507",],aes(x=V3,y=foldd))+geom_point(size=0.2,alpha=0.4,color="tan")+
  geom_smooth(fill = NA, size=0.5,color="tan") +labs(x="coords",y="depth",title="TRE_scaffold0010507")+theme+scale_y_continuous(limits = c(0,1.3))+theme(axis.title.x = element_text(size = 0))
ptregdepth14<-ggplot(dmix[dmix$V1=="TRE_scaffold0011037",],aes(x=V3,y=foldd))+geom_point(size=0.2,alpha=0.4,color="tan")+
  geom_smooth(fill = NA, size=0.5,color="tan") +labs(x="coords",y="depth",title="TRE_scaffold0011037")+theme+scale_y_continuous(limits = c(0,1.3))+theme(axis.title.x = element_text(size = 0))
ptregdepth15<-ggplot(dmix[dmix$V1=="TRE_scaffold0011886",],aes(x=V3,y=foldd))+geom_point(size=0.2,alpha=0.4,color="tan")+
  geom_smooth(fill = NA, size=0.5,color="tan") +labs(x="coords",y="depth",title="TRE_scaffold0011886")+theme+scale_y_continuous(limits = c(0,1.3))+theme(axis.title.x = element_text(size = 0))
ptregdepth16<-ggplot(dmix[dmix$V1=="TRE_scaffold0014813",],aes(x=V3,y=foldd))+geom_point(size=0.2,alpha=0.4,color="tan")+
  geom_smooth(fill = NA, size=0.5,color="tan") +labs(x="coords",y="depth",title="TRE_scaffold0014813")+theme+scale_y_continuous(limits = c(0,1.3))+theme(axis.title.x = element_text(size = 0))
ptregdepth17<-ggplot(dmix[dmix$V1=="TRE_scaffold0015391",],aes(x=V3,y=foldd))+geom_point(size=0.2,alpha=0.4,color="tan")+
  geom_smooth(fill = NA, size=0.5,color="tan") +labs(x="coords",y="depth",title="TRE_scaffold0015391")+theme+scale_y_continuous(limits = c(0,1.3))+theme(axis.title.x = element_text(size = 0))
ptregdepth18<-ggplot(dmix[dmix$V1=="TRE_scaffold0016254",],aes(x=V3,y=foldd))+geom_point(size=0.2,alpha=0.4,color="tan")+
  geom_smooth(fill = NA, size=0.5,color="tan") +labs(x="coords",y="depth",title="TRE_scaffold0016254")+theme+scale_y_continuous(limits = c(0,1.3))+theme(axis.title.x = element_text(size = 0))
ptregdepth19<-ggplot(dmix[dmix$V1=="TRE_scaffold0017642",],aes(x=V3,y=foldd))+geom_point(size=0.2,alpha=0.4,color="tan")+
  geom_smooth(fill = NA, size=0.5,color="tan") +labs(x="coords",y="depth",title="TRE_scaffold0017642")+theme+scale_y_continuous(limits = c(0,1.3))+theme(axis.title.x = element_text(size = 0))
ptregdepth20<-ggplot(dmix[dmix$V1=="TRE_scaffold0018945",],aes(x=V3,y=foldd))+geom_point(size=0.2,alpha=0.4,color="tan")+
  geom_smooth(fill = NA, size=0.5,color="tan") +labs(x="coords",y="depth",title="TRE_scaffold0018945")+theme+scale_y_continuous(limits = c(0,1.3))+theme(axis.title.x = element_text(size = 0))
ptregdepth21<-ggplot(dmix[dmix$V1=="TRE_scaffold0019044",],aes(x=V3,y=foldd))+geom_point(size=0.2,alpha=0.4,color="tan")+
  geom_smooth(fill = NA, size=0.5,color="tan") +labs(x="coords",y="depth",title="TRE_scaffold0019044")+theme+scale_y_continuous(limits = c(0,1.3))+theme(axis.title.x = element_text(size = 0))
ptregdepth22<-ggplot(dmix[dmix$V1=="TRE_scaffold0019511",],aes(x=V3,y=foldd))+geom_point(size=0.2,alpha=0.4,color="tan")+
  geom_smooth(fill = NA, size=0.5,color="tan") +labs(x="coords",y="depth",title="TRE_scaffold0019511")+theme+scale_y_continuous(limits = c(0,1.3))+theme(axis.title.x = element_text(size = 0))
ptregdepth23<-ggplot(dmix[dmix$V1=="TRE_scaffold0022670",],aes(x=V3,y=foldd))+geom_point(size=0.2,alpha=0.4,color="tan")+
  geom_smooth(fill = NA, size=0.5,color="tan") +labs(x="coords",y="depth",title="TRE_scaffold0022670")+theme+scale_y_continuous(limits = c(0,1.3))+theme(axis.title.x = element_text(size = 0))
ptregdepth24<-ggplot(dmix[dmix$V1=="TRE_contig0000078",],aes(x=V3,y=foldd))+geom_point(size=0.2,alpha=0.4,color="tan")+
  geom_smooth(fill = NA, size=0.5,color="tan") +labs(x="coords",y="depth",title="TRE_contig0000078")+theme+scale_y_continuous(limits = c(0,1.3))+theme(axis.title.x = element_text(size = 0))
ptregdepth25<-ggplot(dmix[dmix$V1=="TRE_scaffold0053261",],aes(x=V3,y=foldd))+geom_point(size=0.2,alpha=0.4,color="tan")+
  geom_smooth(fill = NA, size=0.5,color="tan") +labs(x="coords",y="depth",title="TRE_scaffold0053261")+theme+scale_y_continuous(limits = c(0,1.3))+theme(axis.title.x = element_text(size = 0))
ptregdepth1
ptmurdepth<-multiplot(ptregdepth1,ptregdepth2,ptregdepth3,ptregdepth4,ptregdepth5,ptregdepth6,ptregdepth7,ptregdepth8,ptregdepth9,ptregdepth10,ptregdepth11,ptregdepth12,ptregdepth13,ptregdepth14,ptregdepth15,ptregdepth16,ptregdepth17,ptregdepth18,ptregdepth19,ptregdepth20,ptregdepth21,ptregdepth22,ptregdepth23,ptregdepth24,ptregdepth25,cols=5)

#pschi<-multiplot(p2,p3,p4,p5,p6,pa,p7,p8,p9,p23,cols=2)

mix <- read.table("Csin_1.fastq.gz.depth.10k-win.txt")
mix <- subset(mix, mix$V5 > 500)
mix_v4 <- as.numeric(mix$V4)
sum_mix_v4 <- sum(mix_v4)
sum_mix_v4
head(mix)
mix_v5 <- as.numeric(mix$V5)
sum_mix_v5 <- sum(mix_v5)
sum_mix_v5
mix$dmix <- mix_v4/(mix_v5*2)
#install.packages("ggplot2")
library(ggplot2)
dmix <- mix[c(1,8)]
dmix$foldd <- dmix$dmix/5.4
p8 <- ggplot(dmix,aes(x=foldd))+geom_density(colour="dark grey")+
  labs(x="depth",title="Clonorchis sinensis")+scale_x_continuous(limits = c(0,2))+
  theme(axis.text.x = element_text(size = 6))+
  theme(axis.text.y = element_text(size = 6))+
  theme(axis.title.y = element_text(size = 0))+
  theme(axis.title.x = element_text(size = 0))+
  theme(title=element_text(size = 4.5,face = "bold.italic",color = "#a6761d"))
p8

mix <- read.table("Ecaproni_1.fastq.gz.depth.10k-win.txt")
mix <- subset(mix, mix$V5 > 500)
mix_v4 <- as.numeric(mix$V4)
sum_mix_v4 <- sum(mix_v4)
sum_mix_v4
head(mix)
mix_v5 <- as.numeric(mix$V5)
sum_mix_v5 <- sum(mix_v5)
sum_mix_v5
mix$dmix <- mix_v4/(mix_v5*2)
dmix <- mix[c(1,8)]
library(ggplot2)

dmix$foldd <- dmix$dmix/6.2
p9 <- ggplot(dmix,aes(x=foldd))+geom_density(colour="dark grey")+
  labs(x="depth",title="Echinostoma caproni")+scale_x_continuous(limits = c(0,2))+  
  #geom_vline(aes(xintercept = 0.87), linetype = "dashed", size = 0.4)+
  theme(axis.text.x = element_text(size = 6))+
  theme(axis.text.y = element_text(size = 6))+
  theme(axis.title.y = element_text(size = 0))+
  theme(axis.title.x = element_text(size = 0))+
  theme(title=element_text(size = 4.5,face = "bold.italic",color = "#a6761d"))
p9

mix <- read.table("Fhepatica_1.fastq.gz.depth.10k-win.txt")
mix <- subset(mix, mix$V5 > 500)
mix_v4 <- as.numeric(mix$V4)
sum_mix_v4 <- sum(mix_v4)
sum_mix_v4
head(mix)
mix_v5 <- as.numeric(mix$V5)
sum_mix_v5 <- sum(mix_v5)
sum_mix_v5
mix$dmix <- mix_v4/(mix_v5*2)
#install.packages("ggplot2")
library(ggplot2)
dmix <- mix[c(1,8)]

dmix$foldd <- dmix$dmix/14.5
p10 <- ggplot(dmix,aes(x=foldd))+geom_density(colour="dark grey")+
  labs(x="depth",title="Fasciola hepatica")+scale_x_continuous(limits = c(0,2))+
  theme(axis.text.x = element_text(size = 6))+
  theme(axis.text.y = element_text(size = 6))+
  theme(axis.title.y = element_text(size = 0))+
  theme(axis.title.x = element_text(size = 0))+
  theme(title=element_text(size = 4.5,face = "bold.italic",color = "#a6761d"))
p10

#ptrematoda<-multiplot(p2,p3,p4,p5,p6,psrod,p7,p8,p9,p10,cols=2)

mix <- read.table("Pxenopodis_1.fastq.gz.depth.10k-win.txt")
mix <- subset(mix, mix$V5 > 500)
mix_v4 <- as.numeric(mix$V4)
sum_mix_v4 <- sum(mix_v4)
sum_mix_v4
head(mix)
mix_v5 <- as.numeric(mix$V5)
sum_mix_v5 <- sum(mix_v5)
sum_mix_v5
mix$dmix <- mix_v4/(mix_v5*2)
#install.packages("ggplot2")
library(ggplot2)
dmix <- mix[c(1,8)]

dmix$foldd <- dmix$dmix/18.5
*p23 <- ggplot(dmix,aes(x=foldd))+geom_density(colour="dark grey")+
  labs(x="depth",title="Protopolystoma xenopodis")+scale_x_continuous(limits = c(0,2))+
  theme(axis.text.x = element_text(size = 6))+
  theme(axis.text.y = element_text(size = 6))+
  theme(axis.title.y = element_text(size = 0))+
  theme(axis.title.x = element_text(size = 0))+
  theme(title=element_text(size = 4.5,face = "bold.italic",color = "#e7298a"))
*p23


mix <- read.table("Smed_1.fastq.gz.depth.10k-win.txt")
mix <- subset(mix, mix$V5 > 500)
mix_v4 <- as.numeric(mix$V4)
sum_mix_v4 <- sum(mix_v4)
sum_mix_v4
head(mix)
mix_v5 <- as.numeric(mix$V5)
sum_mix_v5 <- sum(mix_v5)
sum_mix_v5
mix$dmix <- mix_v4/(mix_v5*2)
#install.packages("ggplot2")
library(ggplot2)
dmix <- mix[c(1,8)]

dmix$foldd <- dmix$dmix/30
p24 <- ggplot(dmix,aes(x=foldd))+geom_density(colour="dark grey")+
  labs(x="depth",title="Schmidtea mediterranea")+scale_x_continuous(limits = c(0,2))+
  #geom_vline(aes(xintercept = 0.125), linetype = "dashed", size = 0.4)+
  theme(axis.text.x = element_text(size = 6))+
  theme(axis.text.y = element_text(size = 6))+
  theme(axis.title.y = element_text(size = 0))+
  theme(axis.title.x = element_text(size = 0))+
  theme(title=element_text(size = 4.5,face = "bold.italic",color = "#e7298a"))
p24

#pplathy<-multiplot(p23,p24,cols=2)

mix <- read.table("Dlatus.depth.bed.txt")
mix <- subset(mix, mix$V5 > 500)
mix_v4 <- as.numeric(mix$V4)
sum_mix_v4 <- sum(mix_v4)
sum_mix_v4
head(mix)
mix_v5 <- as.numeric(mix$V5)
sum_mix_v5 <- sum(mix_v5)
sum_mix_v5
mix$dmix <- mix_v4/(mix_v5*2)
#install.packages("ggplot2")
library(ggplot2)
dmix <- mix[c(1,8)]

dmix$foldd <- dmix$dmix/8.6
p11 <- ggplot(dmix,aes(x=foldd))+geom_density(colour="dark grey")+
  labs(x="depth",title="Dibothriocephalus latus")+scale_x_continuous(limits = c(0,2))+  
  geom_vline(aes(xintercept = 0.7), linetype = "dashed", size = 0.4)+
  theme(axis.text.x = element_text(size = 6))+
  theme(axis.text.y = element_text(size = 6))+
  theme(axis.title.y = element_text(size = 0))+
  theme(axis.title.x = element_text(size = 0))+
  theme(title=element_text(size = 4.5,face = "bold.italic"))
p11

mix <- read.table("Egranulosus.depth.bed.txt")
mix <- subset(mix, mix$V5 > 500)
mix_v4 <- as.numeric(mix$V4)
sum_mix_v4 <- sum(mix_v4)
sum_mix_v4
head(mix)
mix_v5 <- as.numeric(mix$V5)
sum_mix_v5 <- sum(mix_v5)
sum_mix_v5
mix$dmix <- mix_v4/(mix_v5*2)
#install.packages("ggplot2")
library(ggplot2)
dmix <- mix[c(1,8)]

dmix$foldd <- dmix$dmix/14.5
p12 <- ggplot(dmix,aes(x=foldd))+geom_density(colour="dark grey")+
  labs(x="depth",title="Echinococcus granulosus")+scale_x_continuous(limits = c(0,2))+
  theme(axis.text.x = element_text(size = 6))+
  theme(axis.text.y = element_text(size = 6))+
  theme(axis.title.y = element_text(size = 0))+
  theme(axis.title.x = element_text(size = 0))+
  theme(title=element_text(size = 4.5,face = "bold.italic"))
p12

setwd("/Users/yifeng/Library/Mobile Documents/com~apple~CloudDocs/project/81worm/depthfile/cestoda/")#"original wbps12 masked depth"
mix <- read.table("Emultilocularis.depth.bed.txt")
mix <- subset(mix, mix$V5 > 500)
mix_v4 <- as.numeric(mix$V4)
sum_mix_v4 <- sum(mix_v4)
sum_mix_v4
head(mix)
mix_v5 <- as.numeric(mix$V5)
sum_mix_v5 <- sum(mix_v5)
sum_mix_v5
mix$dmix <- mix_v4/(mix_v5*2)
#install.packages("ggplot2")
library(ggplot2)
dmix <- mix[c(1,8)]

dmix$foldd <- dmix$dmix/80
p13 <- ggplot(dmix,aes(x=foldd))+geom_density(colour="dark grey")+
  labs(x="depth",title="Echinococcus multilocularis")+scale_x_continuous(limits = c(0,2))+
  theme(axis.text.x = element_text(size = 6))+
  theme(axis.text.y = element_text(size = 6))+
  theme(axis.title.y = element_text(size = 0))+
  theme(axis.title.x = element_text(size = 0))+
  theme(title=element_text(size = 4.5,face = "bold.italic",color="#4d4d4d"))
p13
pdf("Emul.depth.peak.pdf",width=1.4,height=0.9)
p13
dev.off()

mix <- read.table("H.taeni.bed1.depth.bed.txt")
mix <- subset(mix, mix$V5 > 500)
mix_v4 <- as.numeric(mix$V4)
sum_mix_v4 <- sum(mix_v4)
sum_mix_v4
head(mix)
mix_v5 <- as.numeric(mix$V5)
sum_mix_v5 <- sum(mix_v5)
sum_mix_v5
mix$dmix <- mix_v4/(mix_v5*2)
dmix <- mix[c(1,8)]
library(ggplot2)

dmix$foldd <- dmix$dmix/80
p14 <- ggplot(dmix,aes(x=foldd))+geom_density(colour="dark grey")+
  labs(x="depth",title="Hydatigera taeniaeformis")+scale_x_continuous(limits = c(0,2))+
  theme(axis.text.x = element_text(size = 6))+
  theme(axis.text.y = element_text(size = 6))+
  theme(axis.title.y = element_text(size = 0))+
  theme(axis.title.x = element_text(size = 0))+
  theme(title=element_text(size = 4.5,face = "bold.italic"))
p14

mix <- read.table("Hdiminuta.depth.bed.txt")
mix <- subset(mix, mix$V5 > 500)
mix_v4 <- as.numeric(mix$V4)
sum_mix_v4 <- sum(mix_v4)
sum_mix_v4
head(mix)
mix_v5 <- as.numeric(mix$V5)
sum_mix_v5 <- sum(mix_v5)
sum_mix_v5
mix$dmix <- mix_v4/(mix_v5*2)
#install.packages("ggplot2")
library(ggplot2)
dmix <- mix[c(1,8)]

dmix$foldd <- dmix$dmix/64
p15 <- ggplot(dmix,aes(x=foldd))+geom_density(colour="dark grey")+
  labs(x="depth",title="Hymenolepis diminuta")+scale_x_continuous(limits = c(0,2))+
  theme(axis.text.x = element_text(size = 6))+
  theme(axis.text.y = element_text(size = 6))+
  theme(axis.title.y = element_text(size = 0))+
  theme(axis.title.x = element_text(size = 0))+
  theme(title=element_text(size = 4.5,face = "bold.italic"))
p15


mix <- read.table("Hmicrostoma.depth.bed.txt")
mix <- subset(mix, mix$V5 > 500)
mix_v4 <- as.numeric(mix$V4)
sum_mix_v4 <- sum(mix_v4)
sum_mix_v4
head(mix)
mix_v5 <- as.numeric(mix$V5)
sum_mix_v5 <- sum(mix_v5)
sum_mix_v5
mix$dmix <- mix_v4/(mix_v5*2)
#install.packages("ggplot2")
library(ggplot2)
dmix <- mix[c(1,8)]

dmix$foldd <- dmix$dmix/50.5
p16 <- ggplot(dmix,aes(x=foldd))+geom_density(colour="dark grey")+
  labs(x="depth",title="Hymenolepis microstoma")+scale_x_continuous(limits = c(0,2))+
  theme(axis.text.x = element_text(size = 6))+
  theme(axis.text.y = element_text(size = 6))+
  theme(axis.title.y = element_text(size = 0))+
  theme(axis.title.x = element_text(size = 0))+
  theme(title=element_text(size = 4.5,face = "bold.italic"))
p16

mix <- read.table("Hnana.depth.bed.txt")
mix <- subset(mix, mix$V5 > 500)
mix_v4 <- as.numeric(mix$V4)
sum_mix_v4 <- sum(mix_v4)
sum_mix_v4
head(mix)
mix_v5 <- as.numeric(mix$V5)
sum_mix_v5 <- sum(mix_v5)
sum_mix_v5
mix$dmix <- mix_v4/(mix_v5*2)
#install.packages("ggplot2")
library(ggplot2)
dmix <- mix[c(1,8)]

dmix$foldd <- dmix$dmix/75
p17 <- ggplot(dmix,aes(x=foldd))+geom_density(colour="dark grey")+
  labs(x="depth",title="Hymenolepis nana")+scale_x_continuous(limits = c(0,2))+  
  #geom_vline(aes(xintercept = 0.77), linetype = "dashed", size = 0.4)+
  theme(axis.text.x = element_text(size = 6))+
  theme(axis.text.y = element_text(size = 6))+
  theme(axis.title.y = element_text(size = 0))+
  theme(axis.title.x = element_text(size = 0))+
  theme(title=element_text(size = 4.5,face = "bold.italic"))
p17

mix <- read.table("Mcorti.depth.bed.txt")
mix <- subset(mix, mix$V5 > 500)
mix_v4 <- as.numeric(mix$V4)
sum_mix_v4 <- sum(mix_v4)
sum_mix_v4
head(mix)
mix_v5 <- as.numeric(mix$V5)
sum_mix_v5 <- sum(mix_v5)
sum_mix_v5
mix$dmix <- mix_v4/(mix_v5*2)
#install.packages("ggplot2")
library(ggplot2)
dmix <- mix[c(1,8)]

dmix$foldd <- dmix$dmix/102
p18 <- ggplot(dmix,aes(x=foldd))+geom_density(colour="dark grey")+
  labs(x="depth",title="Mesocestoides corti")+scale_x_continuous(limits = c(0,2))+
  theme(axis.text.x = element_text(size = 6))+
  theme(axis.text.y = element_text(size = 6))+
  theme(axis.title.y = element_text(size = 0))+
  theme(axis.title.x = element_text(size = 0))+
  theme(title=element_text(size = 4.5,face = "bold.italic"))
p18


mix <- read.table("Ssolidus.depth.bed.txt")
mix <- subset(mix, mix$V5 > 500)
mix_v4 <- as.numeric(mix$V4)
sum_mix_v4 <- sum(mix_v4)
sum_mix_v4
head(mix)
mix_v5 <- as.numeric(mix$V5)
sum_mix_v5 <- sum(mix_v5)
sum_mix_v5
mix$dmix <- mix_v4/(mix_v5*2)
#install.packages("ggplot2")
library(ggplot2)
dmix <- mix[c(1,8)]

dmix$foldd <- dmix$dmix/19
p19 <- ggplot(dmix,aes(x=foldd))+geom_density(colour="dark grey")+
  labs(x="depth",title="Schistocephalus solidus")+scale_x_continuous(limits = c(0,2))+
  geom_vline(aes(xintercept = 0.85), linetype = "dashed", size = 0.4)+
  theme(axis.text.x = element_text(size = 6))+
  theme(axis.text.y = element_text(size = 6))+
  theme(axis.title.y = element_text(size = 0))+
  theme(axis.title.x = element_text(size = 0))+
  theme(title=element_text(size = 4.5,face = "bold.italic"))
p19


mix <- read.table("Serinaceieuropaei.depth.bed.txt")
mix <- subset(mix, mix$V5 > 500)
mix_v4 <- as.numeric(mix$V4)
sum_mix_v4 <- sum(mix_v4)
sum_mix_v4
head(mix)
mix_v5 <- as.numeric(mix$V5)
sum_mix_v5 <- sum(mix_v5)
sum_mix_v5
mix$dmix <- mix_v4/(mix_v5*2)
#install.packages("ggplot2")
library(ggplot2)
dmix <- mix[c(1,8)]

dmix$foldd <- dmix$dmix/4.6
p20<-ggplot(dmix,aes(x=foldd))+geom_density(colour="dark grey")+
  labs(x="depth",title="Spirometra erinaceieuropaei")+scale_x_continuous(limits = c(0,2))+
  theme(axis.text.x = element_text(size = 6))+
  theme(axis.text.y = element_text(size = 6))+
  theme(axis.title.y = element_text(size = 0))+
  theme(axis.title.x = element_text(size = 0))+
  theme(title=element_text(size = 4.5,face = "bold.italic"))
p20

mix <- read.table("Tasiatica.depth.bed.txt")
mix <- subset(mix, mix$V5 > 500)
mix_v4 <- as.numeric(mix$V4)
sum_mix_v4 <- sum(mix_v4)
sum_mix_v4
head(mix)
mix_v5 <- as.numeric(mix$V5)
sum_mix_v5 <- sum(mix_v5)
sum_mix_v5
mix$dmix <- mix_v4/(mix_v5*2)
#install.packages("ggplot2")
library(ggplot2)
dmix <- mix[c(1,8)]

dmix$foldd <- dmix$dmix/55
p21 <- ggplot(dmix,aes(x=foldd))+geom_density(colour="dark grey")+
  labs(x="depth",title="Taenia asiatica")+scale_x_continuous(limits = c(0,2))+
  #geom_vline(aes(xintercept = 0.84), linetype = "dashed", size = 0.4)+
  theme(axis.text.x = element_text(size = 6))+
  theme(axis.text.y = element_text(size = 6))+
  theme(axis.title.y = element_text(size = 0))+
  theme(axis.title.x = element_text(size = 0))+
  theme(title=element_text(size = 4.5,face = "bold.italic"))
p21

mix <- read.table("Tsolium.depth.bed.txt")
mix <- subset(mix, mix$V5 > 500)
mix_v4 <- as.numeric(mix$V4)
sum_mix_v4 <- sum(mix_v4)
sum_mix_v4
head(mix)
mix_v5 <- as.numeric(mix$V5)
sum_mix_v5 <- sum(mix_v5)
sum_mix_v5
mix$dmix <- mix_v4/(mix_v5*2)
#install.packages("ggplot2")
library(ggplot2)
dmix <- mix[c(1,8)]

dmix$foldd <- dmix$dmix/15.5
p22 <- ggplot(dmix,aes(x=foldd))+geom_density(colour="dark grey")+
  labs(x="depth",title="Taenia solium")+scale_x_continuous(limits = c(0,2))+
  theme(axis.text.x = element_text(size = 6))+
  theme(axis.text.y = element_text(size = 6))+
  theme(axis.title.y = element_text(size = 0))+
  theme(axis.title.x = element_text(size = 0))+
  theme(title=element_text(size = 4.5,face = "bold.italic"))
p22



#pceswith2other<-multiplot(p8,p9,p10,p12,p13,p14,p15,p16,p17,p18,p20,p21,p22,p23,p11,p19,p24,cols=3)

##########
setwd("/Users/yifeng/Library/Mobile Documents/com~apple~CloudDocs/project/81worm/depthfile/I/")
mix <- read.table("Rcul_1.fastq.gz.depth.10k-win.txt")
mix <- subset(mix, mix$V5 > 500)
mix_v4 <- as.numeric(mix$V4)
sum_mix_v4 <- sum(mix_v4)
sum_mix_v4
head(mix)
mix_v5 <- as.numeric(mix$V5)
sum_mix_v5 <- sum(mix_v5)
sum_mix_v5
mix$dmix <- mix_v4/(mix_v5*2)
#install.packages("ggplot2")
library(ggplot2)
dmix <- mix[c(1,8)]

dmix$foldd <- dmix$dmix/5
p25 <- ggplot(dmix,aes(x=foldd))+geom_density(colour="tan")+
  labs(x="depth",title="Romanomermis culicivorax")+scale_x_continuous(limits = c(0,2))+
  geom_vline(aes(xintercept = 0.73), linetype = "dashed", size = 0.4)+
  theme(axis.text.x = element_text(size = 6))+
  theme(axis.text.y = element_text(size = 6))+
  theme(axis.title.y = element_text(size = 0))+
  theme(axis.title.x = element_text(size = 0))+
  theme(title=element_text(size = 4.5,face = "bold.italic",color="#66a61e"))
p25
dmix<-subset(dmix,foldd<0.73 & foldd>0.25)

mix <- read.table("Tmuri.m_1.fastq.gz.depth.10k-win.txt")
mix <- subset(mix, mix$V5 > 500)
mix_v4 <- as.numeric(mix$V4)
sum_mix_v4 <- sum(mix_v4)
sum_mix_v4
head(mix)
mix_v5 <- as.numeric(mix$V5)
sum_mix_v5 <- sum(mix_v5)
sum_mix_v5
mix$dmix <- mix_v4/(mix_v5*2)
#install.packages("ggplot2")
library(ggplot2)
dmix <- mix[c(1,2,3,8)]

##dmix <- subset(dmix,dmix<360)
dmix$foldd <- dmix$dmix/75
p26 <- ggplot(dmix,aes(x=foldd))+geom_density(colour="forest green")+
  labs(x="depth",title="Trichuris muris")+scale_x_continuous(limits = c(0,2))+
  geom_vline(aes(xintercept = 0.66), linetype = "dashed", size = 0.4)+
  #geom_vline(aes(xintercept = 0.27), linetype = "dashed", size = 0.4)+
  theme(axis.text.x = element_text(size = 6))+
  theme(axis.text.y = element_text(size = 6))+
  theme(axis.title.y = element_text(size = 0))+
  theme(axis.title.x = element_text(size = 0))+
  theme(title=element_text(size = 4.5,face = "bold.italic",color="#66a61e"))
p26
pdf("Tmur.depth.peak.pdf",width=2,height=2)
p26
dev.off()
dmix<-subset(dmix,foldd>0.23 & foldd<0.66)
write.csv(dmix,file = "Tmur.1stpeak.window.csv")

dmix<-subset(dmix,foldd<0.23)
write.csv(dmix,file = "Tmur.1stpeakabout0.window.csv")
#
head(dmix)
ptmurdepth1<-ggplot(dmix[dmix$V1=="TMUE_LG1",],aes(x=V3,y=foldd))+geom_point(size=0.2,alpha=0.4,color="forest green")+
  geom_smooth(fill = NA, size=0.5,color="forest green") +labs(x="coords",y="depth",title="Tmur_chr1")+theme+
  scale_y_continuous(limits = c(0,1.3))+theme(axis.title.x = element_text(size = 0))
ptmurdepth1
ptmurdepth2<-ggplot(dmix[dmix$V1=="TMUE_LG2",],aes(x=V3,y=foldd))+geom_point(size=0.2,alpha=0.4,color="forest green")+
  geom_smooth(fill = NA, size=0.5,color="forest green") +labs(x="coords",y="depth",title="Tmur_chr2")+theme+
  scale_y_continuous(limits = c(0,1.3))+theme(axis.title.x = element_text(size = 0))+theme(axis.title.y = element_text(size = 0))
ptmurdepth2
ptmurdepth3<-ggplot(dmix[dmix$V1=="TMUE_LG3",],aes(x=V3,y=foldd))+geom_point(size=0.2,alpha=0.4,color="forest green")+
  geom_smooth(fill = NA, size=0.5,color="forest green") +labs(x="coords",y="depth",title="Tmur_chr3")+theme+
  scale_y_continuous(limits = c(0,1.3))+theme(axis.title.x = element_text(size = 0))+theme(axis.title.y = element_text(size = 0))
ptmurdepth3
ptmurdepth<-multiplot(ptmurdepth1,ptmurdepth1,ptmurdepth2,ptmurdepth1,ptmurdepth3,cols=3)



mix <- read.table("Tsui.m_1.fastq.gz.depth.10k-win.txt")
mix <- subset(mix, mix$V5 > 500)
mix_v4 <- as.numeric(mix$V4)
sum_mix_v4 <- sum(mix_v4)
sum_mix_v4
head(mix)
mix_v5 <- as.numeric(mix$V5)
sum_mix_v5 <- sum(mix_v5)
sum_mix_v5
mix$dmix <- mix_v4/(mix_v5*2)
dmix <- mix[c(1,8)]
#install.packages("ggplot2")
library(ggplot2)

#dmix <- subset(dmix,dmix<480)
dmix$foldd <- dmix$dmix/20
p27 <- ggplot(dmix,aes(x=foldd))+geom_density(colour="forest green")+
  labs(x="depth",title="Trichuris suis")+scale_x_continuous(limits = c(0,2))+
  geom_vline(aes(xintercept = 0.73), linetype = "dashed", size = 0.4)+
  theme(axis.text.x = element_text(size = 6))+
  theme(axis.text.y = element_text(size = 6))+
  theme(axis.title.y = element_text(size = 0))+
  theme(axis.title.x = element_text(size = 0))+
  theme(title=element_text(size = 4.5,face = "bold.italic",color="#66a61e"))
p27

mix <- read.table("Ttri_1.fastq.gz.depth.10k-win.txt")
mix <- subset(mix, mix$V5 > 500)
mix_v4 <- as.numeric(mix$V4)
sum_mix_v4 <- sum(mix_v4)
sum_mix_v4
head(mix)
mix_v5 <- as.numeric(mix$V5)
sum_mix_v5 <- sum(mix_v5)
sum_mix_v5
mix$dmix <- mix_v4/(mix_v5*2)
#install.packages("ggplot2")
library(ggplot2)
dmix <- mix[c(1,8)]

#dmix <- subset(dmix,dmix<600)
dmix$foldd <- dmix$dmix/94
p28<-ggplot(dmix,aes(x=foldd))+geom_density(colour="forest green")+
  labs(x="depth",title="Trichuris trichiura")+scale_x_continuous(limits = c(0,2))+
  geom_vline(aes(xintercept = 0.665), linetype = "dashed", size = 0.4)+
  theme(axis.text.x = element_text(size = 6))+
  theme(axis.text.y = element_text(size = 6))+
  theme(axis.title.y = element_text(size = 0))+
  theme(axis.title.x = element_text(size = 0))+
  theme(title=element_text(size = 4.5,face = "bold.italic",color="#66a61e"))
p28
#pI<-multiplot(p25,p26,p27,p28,cols=2)


setwd("/Users/yifeng/Library/Mobile Documents/com~apple~CloudDocs/project/81worm/depthfile/IV/")
mix <- read.table("Ptri_1.fastq.gz.depth.10k-win.txt")
mix <- subset(mix, mix$V5 > 500)
mix_v4 <- as.numeric(mix$V4)
sum_mix_v4 <- sum(mix_v4)
sum_mix_v4
head(mix)
mix_v5 <- as.numeric(mix$V5)
sum_mix_v5 <- sum(mix_v5)
sum_mix_v5
mix$dmix <- mix_v4/(mix_v5*2)
dmix <- mix[c(1,2,3,8)]
#install.packages("ggplot2")
library(ggplot2)

#dmix <- subset(dmix,dmix<360)
dmix$foldd <- dmix$dmix/46.5
p29 <- ggplot(dmix,aes(x=foldd))+geom_density(colour="forest green")+
  labs(x="depth",title="Parastrongyloides trichosuri")+scale_x_continuous(limits = c(0,2))+
  geom_vline(aes(xintercept = 0.6), linetype = "dashed", size = 0.4)+
  #geom_vline(aes(xintercept = 0.19), linetype = "dashed", size = 0.4)+
  theme(axis.text.x = element_text(size = 6))+
  theme(axis.text.y = element_text(size = 6))+
  theme(axis.title.y = element_text(size = 0))+
  theme(axis.title.x = element_text(size = 0))+
  theme(title=element_text(size = 4.5,face = "bold.italic",color="#7570b3"))
p29
pdf("Ptri.depth.peak.pdf",width=2,height=2)
p29
dev.off()

dmix<-subset(dmix,foldd>0.2 & foldd<0.6)
write.csv(dmix,file = "Ptri.1stpeak.window.csv")
#
pptridepth1<-ggplot(dmix[dmix$V1=="PTRK_scaffold0000001",],aes(x=V3,y=foldd))+geom_point(size=0.2,alpha=0.4,color="forest green")+
  geom_smooth(fill = NA, size=0.5,color="forest green") +labs(x="coords",y="depth",title="Ptri_scaf1")+theme+
  scale_y_continuous(limits = c(0,1.3))+theme(axis.title.x = element_text(size = 0))
pptridepth1
pptridepth2<-ggplot(dmix[dmix$V1=="PTRK_scaffold0000002",],aes(x=V3,y=foldd))+geom_point(size=0.2,alpha=0.4,color="forest green")+
  geom_smooth(fill = NA, size=0.5,color="forest green") +labs(x="coords",y="depth",title="Ptri_scaf2")+theme+
  scale_y_continuous(limits = c(0,1.3))+theme(axis.title.x = element_text(size = 0))+theme(axis.title.y = element_text(size = 0))
pptridepth2
pptridepth3<-ggplot(dmix[dmix$V1=="PTRK_scaffold0000003",],aes(x=V3,y=foldd))+geom_point(size=0.2,alpha=0.4,color="forest green")+
  geom_smooth(fill = NA, size=0.5,color="forest green") +labs(x="coords",y="depth",title="Ptri_scaf3")+theme+
  scale_y_continuous(limits = c(0,1.3))+theme(axis.title.x = element_text(size = 0))+theme(axis.title.y = element_text(size = 0))
pptridepth3
pptridepth4<-ggplot(dmix[dmix$V1=="PTRK_scaffold0000004",],aes(x=V3,y=foldd))+geom_point(size=0.2,alpha=0.4,color="forest green")+
  geom_smooth(fill = NA, size=0.5,color="forest green") +labs(x="coords",y="depth",title="Ptri_scaf4")+theme+
  scale_y_continuous(limits = c(0,1.3))+theme(axis.title.x = element_text(size = 0))+theme(axis.title.y = element_text(size = 0))
pptridepth4
pptridepth5<-ggplot(dmix[dmix$V1=="PTRK_scaffold0000005",],aes(x=V3,y=foldd))+geom_point(size=0.2,alpha=0.4,color="forest green")+
  geom_smooth(fill = NA, size=0.5,color="forest green") +labs(x="coords",y="depth",title="Ptri_scaf5")+theme+
  scale_y_continuous(limits = c(0,1.3))+theme(axis.title.x = element_text(size = 0))+theme(axis.title.y = element_text(size = 0))
pptridepth5
pptridepth<-multiplot(pptridepth1,pptridepth1,pptridepth2,pptridepth1,pptridepth3,pptridepth1,pptridepth4,pptridepth1,pptridepth5,cols=5)





mix <- read.table("Rhab_1.fastq.gz.depth.10k-win.txt")
mix <- subset(mix, mix$V5 > 500)
mix_v4 <- as.numeric(mix$V4)
sum_mix_v4 <- sum(mix_v4)
sum_mix_v4
head(mix)
mix_v5 <- as.numeric(mix$V5)
sum_mix_v5 <- sum(mix_v5)
sum_mix_v5
mix$dmix <- mix_v4/(mix_v5*2)
#install.packages("ggplot2")
library(ggplot2)
dmix <- mix[c(1,8)]

#dmix <- subset(dmix,dmix<1500)
dmix$foldd <- dmix$dmix/298
p30 <- ggplot(dmix,aes(x=foldd))+geom_density(colour="salmon")+
  labs(x="depth",title="Rhabditophanes sp. KR3021")+scale_x_continuous(limits = c(0,2))+
  theme(axis.text.x = element_text(size = 6))+
  theme(axis.text.y = element_text(size = 6))+
  theme(axis.title.y = element_text(size = 0))+
  theme(axis.title.x = element_text(size = 0))+
  theme(title=element_text(size = 4.5,face = "bold.italic",color="#7570b3"))
p30

##Spap new
mix <- read.table("Spap_1.fastq.gz.depth.10k-win.txt")
mix <- subset(mix, mix$V5 > 500)
mix_v4 <- as.numeric(mix$V4)
sum_mix_v4 <- sum(mix_v4)
sum_mix_v4
head(mix)
mix_v5 <- as.numeric(mix$V5)
sum_mix_v5 <- sum(mix_v5)
sum_mix_v5
mix$dmix <- mix_v4/(mix_v5*2)
#install.packages("ggplot2")
library(ggplot2)
dmix <- mix[c(1,8)]

#dmix <- subset(dmix,dmix<200)
dmix$foldd <- dmix$dmix/40.5
p31 <- ggplot(dmix,aes(x=foldd))+geom_density(colour="forest green")+
  labs(x="depth",title="Strongyloides papillosus")+scale_x_continuous(limits = c(0,2))+
  #geom_vline(aes(xintercept = 0.315), linetype = "dashed", size = 0.4)+
  geom_vline(aes(xintercept = 0.67), linetype = "dashed", size = 0.4)+
  theme(axis.text.x = element_text(size = 6))+
  theme(axis.text.y = element_text(size = 6))+
  theme(axis.title.y = element_text(size = 0))+
  theme(axis.title.x = element_text(size = 0))+
  theme(title=element_text(size = 4.5,face = "bold.italic",color="#7570b3"))
p31

###new
mix <- read.table("Srat.m_1.fastq.gz.depth.10k-win.txt")
mix <- subset(mix, mix$V5 > 500)
mix_v4 <- as.numeric(mix$V4)
sum_mix_v4 <- sum(mix_v4)
sum_mix_v4
head(mix)
mix_v5 <- as.numeric(mix$V5)
sum_mix_v5 <- sum(mix_v5)
sum_mix_v5
mix$dmix <- mix_v4/(mix_v5*2)
dmix <- mix[c(1,8)]
#install.packages("ggplot2")
library(ggplot2)

#dmix <- subset(dmix,dmix<120)
dmix$foldd <- dmix$dmix/51
p32 <- ggplot(dmix,aes(x=foldd))+geom_density(colour="forest green")+
  labs(x="depth",title="Strongyloides ratti")+scale_x_continuous(limits = c(0,2))+
  theme(axis.text.x = element_text(size = 6))+
  theme(axis.text.y = element_text(size = 6))+
  theme(axis.title.y = element_text(size = 0))+
  theme(axis.title.x = element_text(size = 0))+
  geom_vline(aes(xintercept = 0.56), linetype = "dashed", size = 0.4)+
  theme(title=element_text(size = 4.5,face = "bold.italic",color="#7570b3"))
p32

mix <- read.table("Sste_1.fastq.gz.depth.10k-win.txt")
mix <- subset(mix, mix$V5 > 1000)
mix_v4 <- as.numeric(mix$V4)
sum_mix_v4 <- sum(mix_v4)
sum_mix_v4
head(mix)
mix_v5 <- as.numeric(mix$V5)
sum_mix_v5 <- sum(mix_v5)
sum_mix_v5
mix$dmix <- mix_v4/(mix_v5*2)
dmix <- mix[c(1,8)]
#install.packages("ggplot2")
library(ggplot2)

#dmix <- subset(dmix,dmix<1300)
dmix$foldd <- dmix$dmix/57
p33 <- ggplot(dmix,aes(x=foldd))+geom_density(colour="forest green")+
  labs(x="depth",title="Strongyloides stercoralis")+scale_x_continuous(limits = c(0,2))+
  theme(axis.text.x = element_text(size = 6))+
  theme(axis.text.y = element_text(size = 6))+
  theme(axis.title.y = element_text(size = 0))+
  theme(axis.title.x = element_text(size = 0))+
  geom_vline(aes(xintercept = 0.67), linetype = "dashed", size = 0.4)+
  theme(title=element_text(size = 4.5,face = "bold.italic",color="#7570b3"))
p33

mix <- read.table("Pred_1.fastq.gz.depth.10k-win.txt")
mix <- subset(mix, mix$V5 > 500)
mix_v4 <- as.numeric(mix$V4)
sum_mix_v4 <- sum(mix_v4)
sum_mix_v4
head(mix)
mix_v5 <- as.numeric(mix$V5)
sum_mix_v5 <- sum(mix_v5)
sum_mix_v5
mix$dmix <- mix_v4/(mix_v5*2)
#install.packages("ggplot2")
library(ggplot2)
dmix <- mix[c(1,8)]

#dmix <- subset(dmix,dmix<110)
dmix$foldd <- dmix$dmix/24
p34 <- ggplot(dmix,aes(x=foldd))+geom_density(colour="tan")+
  labs(x="depth",title="Panagrellus redivivus")+scale_x_continuous(limits = c(0,2))+
  theme(axis.text.x = element_text(size = 6))+
  theme(axis.text.y = element_text(size = 6))+
  theme(axis.title.y = element_text(size = 0))+
  theme(axis.title.x = element_text(size = 0))+
  theme(title=element_text(size = 4.5,face = "bold.italic",color="#7570b3"))
p34

mix <- read.table("Sven_1.fastq.gz.depth.10k-win.txt")
mix <- subset(mix, mix$V5 > 500)
mix_v4 <- as.numeric(mix$V4)
sum_mix_v4 <- sum(mix_v4)
sum_mix_v4
head(mix)
mix_v5 <- as.numeric(mix$V5)
sum_mix_v5 <- sum(mix_v5)
sum_mix_v5
mix$dmix <- mix_v4/(mix_v5*2)
#install.packages("ggplot2")
library(ggplot2)
dmix <- mix[c(1,8)]

#dmix <- subset(dmix,dmix<110)
dmix$foldd <- dmix$dmix/48
psven <- ggplot(dmix,aes(x=foldd))+geom_density(colour="salmon")+
  labs(x="depth",title="Strongyloides venezuelensis")+scale_x_continuous(limits = c(0,2))+
  theme(axis.text.x = element_text(size = 6))+
  theme(axis.text.y = element_text(size = 6))+
  theme(axis.title.y = element_text(size = 0))+
  theme(axis.title.x = element_text(size = 0))+
  theme(title=element_text(size = 4.5,face = "bold.italic",color="#7570b3"))
psven
#pIV<-multiplot(p29,p30,p31,p32,p33,p34,psven,cols=2)

III
setwd("/Users/yifeng/Library/Mobile Documents/com~apple~CloudDocs/project/81worm/depthfile/III/")
mix <- read.table("Bmal_1.fastq.gz.depth.10k-win.txt")
mix <- subset(mix, mix$V5 > 500)
mix_v4 <- as.numeric(mix$V4)
sum_mix_v4 <- sum(mix_v4)
sum_mix_v4
head(mix)
mix_v5 <- as.numeric(mix$V5)
sum_mix_v5 <- sum(mix_v5)
sum_mix_v5
mix$dmix <- mix_v4/(mix_v5*2)
#install.packages("ggplot2")
library(ggplot2)
dmix <- mix[c(1,8)]

#dmix <- subset(dmix,dmix<200)
dmix$foldd <- dmix$dmix/45.5
p35 <- ggplot(dmix,aes(x=foldd))+geom_density(colour="forest green")+
  labs(x="depth",title="Brugia malayi")+scale_x_continuous(limits = c(0,2))+
  geom_vline(aes(xintercept = 0.8), linetype = "dashed", size = 0.4)+
  theme(axis.text.x = element_text(size = 6))+
  theme(axis.text.y = element_text(size = 6))+
  theme(axis.title.y = element_text(size = 0))+
  theme(axis.title.x = element_text(size = 0))+
  theme(title=element_text(size = 4.5,face = "bold.italic",color="#d95f02"))
p35


mix <- read.table("Bpah_1.fastq.gz.depth.10k-win.txt")
mix <- subset(mix, mix$V5 > 500)
mix_v4 <- as.numeric(mix$V4)
sum_mix_v4 <- sum(mix_v4)
sum_mix_v4
head(mix)
mix_v5 <- as.numeric(mix$V5)
sum_mix_v5 <- sum(mix_v5)
sum_mix_v5
mix$dmix <- mix_v4/(mix_v5*2)
dmix <- mix[c(1,2,3,8)]
#install.packages("ggplot2")
library(ggplot2)

#dmix <- subset(dmix,dmix<790)
dmix$foldd <- dmix$dmix/193
p36 <- ggplot(dmix,aes(x=foldd))+geom_density(colour="tan")+
  labs(x="depth",title="Brugia pahangi")+scale_x_continuous(limits = c(0,2))+
  #geom_vline(aes(xintercept = 0.63), linetype = "dashed", size = 0.4)+
  geom_vline(aes(xintercept = 0.86), linetype = "dashed", size = 0.4)+
  theme(axis.text.x = element_text(size = 6))+
  theme(axis.text.y = element_text(size = 6))+
  theme(axis.title.y = element_text(size = 0))+
  theme(axis.title.x = element_text(size = 0))+
  theme(title=element_text(size = 4.5,face = "bold.italic",color="#d95f02"))
p36
pBpahdepth1<-ggplot(dmix[dmix$V1=="BPAG_scaffold0000001",],aes(x=V3,y=foldd))+geom_point(size=0.2,alpha=0.4,color="tan")+
  geom_smooth(fill = NA, size=0.5,color="tan") +labs(x="coords",y="depth",title="Bpah_scaf1")+theme+
  #scale_y_continuous(limits = c(0,1.3))+
  theme(axis.title.x = element_text(size = 0))
pBpahdepth1
pBpahdepth2<-ggplot(dmix[dmix$V1=="BPAG_scaffold0000002",],aes(x=V3,y=foldd))+geom_point(size=0.2,alpha=0.4,color="tan")+
  geom_smooth(fill = NA, size=0.5,color="tan") +labs(x="coords",y="depth",title="Btim_scaf2")+theme+
  #scale_y_continuous(limits = c(0,1.3))+
  theme(axis.title.x = element_text(size = 0))
pBpahdepth2
pBpahdepth3<-ggplot(dmix[dmix$V1=="BPAG_scaffold0000003",],aes(x=V3,y=foldd))+geom_point(size=0.2,alpha=0.4,color="tan")+
  geom_smooth(fill = NA, size=0.5,color="tan") +labs(x="coords",y="depth",title="Btim_scaf3")+theme+
  #scale_y_continuous(limits = c(0,1.3))+
  theme(axis.title.x = element_text(size = 0))
pBpahdepth3
pBpahdepth4<-ggplot(dmix[dmix$V1=="BPAG_scaffold0000004",],aes(x=V3,y=foldd))+geom_point(size=0.2,alpha=0.4,color="tan")+
  geom_smooth(fill = NA, size=0.5,color="tan") +labs(x="coords",y="depth",title="Btim_scaf3")+theme+
  #scale_y_continuous(limits = c(0,1.3))+
  theme(axis.title.x = element_text(size = 0))
pBpahdepth4
pBpahdepth5<-ggplot(dmix[dmix$V1=="BPAG_scaffold0000005",],aes(x=V3,y=foldd))+geom_point(size=0.2,alpha=0.4,color="tan")+
  geom_smooth(fill = NA, size=0.5,color="tan") +labs(x="coords",y="depth",title="Btim_scaf3")+theme+
  #scale_y_continuous(limits = c(0,1.3))+
  theme(axis.title.x = element_text(size = 0))
pBpahdepth5
pBpahdepth6<-ggplot(dmix[dmix$V1=="BPAG_contig0000007",],aes(x=V3,y=foldd))+geom_point(size=0.2,alpha=0.4,color="tan")+
  geom_smooth(fill = NA, size=0.5,color="tan") +labs(x="coords",y="depth",title="Btim_scaf35")+theme+
  #scale_y_continuous(limits = c(0,1.3))+
  theme(axis.title.x = element_text(size = 0))
pBpahdepth6
pBpahdepth6<-ggplot(dmix[dmix$V1=="BPAG_scaffold0000035",],aes(x=V3,y=foldd))+geom_point(size=0.2,alpha=0.4,color="tan")+
  geom_smooth(fill = NA, size=0.5,color="tan") +labs(x="coords",y="depth",title="Btim_scaf35")+theme+
  #scale_y_continuous(limits = c(0,1.3))+
  theme(axis.title.x = element_text(size = 0))
pBpahdepth6



mix <- read.table("Btim_1.fastq.gz.depth.10k-win.txt")
mix <- subset(mix, mix$V5 > 500)
mix_v4 <- as.numeric(mix$V4)
sum_mix_v4 <- sum(mix_v4)
sum_mix_v4
head(mix)
mix_v5 <- as.numeric(mix$V5)
sum_mix_v5 <- sum(mix_v5)
sum_mix_v5
mix$dmix <- mix_v4/(mix_v5*2)
#install.packages("ggplot2")
library(ggplot2)
dmix <- mix[c(1,2,3,8)]

#dmix <- subset(dmix,dmix<70)
dmix$foldd <- dmix$dmix/(16)
#dmix$foldd <- dmix$dmix/(16*1.8)##画scaf depth
p37 <- ggplot(dmix,aes(x=foldd))+geom_density(colour="tan")+
  labs(x="depth",title="Brugia timori")+scale_x_continuous(limits = c(0,2))+
  theme(axis.text.x = element_text(size = 6))+
  theme(axis.text.y = element_text(size = 6))+
  theme(axis.title.y = element_text(size = 0))+
  theme(axis.title.x = element_text(size = 0))+
  theme(title=element_text(size = 4.5,face = "bold.italic",color="#d95f02"))
p37
pptridepth1<-ggplot(dmix[dmix$V1=="BTMF_scaffold0000001",],aes(x=V3,y=foldd))+geom_point(size=0.2,alpha=0.4,color="forest green")+
  geom_smooth(fill = NA, size=0.5,color="forest green") +labs(x="coords",y="depth",title="Btim_scaf1")+theme+
  #scale_y_continuous(limits = c(0,1.3))+
  theme(axis.title.x = element_text(size = 0))
pptridepth1
pptridepth2<-ggplot(dmix[dmix$V1=="BTMF_scaffold0000002",],aes(x=V3,y=foldd))+geom_point(size=0.2,alpha=0.4,color="forest green")+
  geom_smooth(fill = NA, size=0.5,color="forest green") +labs(x="coords",y="depth",title="Btim_scaf1")+theme+
  #scale_y_continuous(limits = c(0,1.3))+
  theme(axis.title.x = element_text(size = 0))
pptridepth2
pptridepth3<-ggplot(dmix[dmix$V1=="BTMF_scaffold0000003",],aes(x=V3,y=foldd))+geom_point(size=0.2,alpha=0.4,color="forest green")+
  geom_smooth(fill = NA, size=0.5,color="forest green") +labs(x="coords",y="depth",title="Btim_scaf1")+theme+
  #scale_y_continuous(limits = c(0,1.3))+
  theme(axis.title.x = element_text(size = 0))
pptridepth3


mix <- read.table("Dimmi_1.fastq.gz.depth.10k-win.txt")
mix <- subset(mix, mix$V5 > 500)
mix_v4 <- as.numeric(mix$V4)
sum_mix_v4 <- sum(mix_v4)
sum_mix_v4
head(mix)
mix_v5 <- as.numeric(mix$V5)
sum_mix_v5 <- sum(mix_v5)
sum_mix_v5
mix$dmix <- mix_v4/(mix_v5*2)
dmix <- mix[c(1,2,3,8)]
#install.packages("ggplot2")
library(ggplot2)

#dmix <- subset(dmix,dmix<460)
dmix$foldd <- dmix$dmix/105
p38 <- ggplot(dmix,aes(x=foldd))+geom_density(colour="tan")+
  labs(x="depth",title="Dirofilaria immitis")+scale_x_continuous(limits = c(0,2))+
  theme(axis.text.x = element_text(size = 6))+
  theme(axis.text.y = element_text(size = 6))+
  theme(axis.title.y = element_text(size = 0))+
  theme(axis.title.x = element_text(size = 0))+
  theme(title=element_text(size = 4.5,face = "bold.italic",color="#d95f02"))
p38

pdimmdepth1<-ggplot(dmix[dmix$V1=="nDi.2.2.scaf00001",],aes(x=V3,y=foldd))+geom_point(size=0.2,alpha=0.4,color="tan")+
  geom_smooth(fill = NA, size=0.5,color="tan") +labs(x="coords",y="depth",title="Dimm_scaf1")+theme+
  scale_y_continuous(limits = c(0,1.3))+theme(axis.title.x = element_text(size = 0))
pdimmdepth1
pdimmdepth1<-ggplot(dmix[dmix$V1=="nDi.2.2.scaf00002",],aes(x=V3,y=foldd))+geom_point(size=0.2,alpha=0.4,color="tan")+
  geom_smooth(fill = NA, size=0.5,color="tan") +labs(x="coords",y="depth",title="Dimm_scaf1")+theme+
  scale_y_continuous(limits = c(0,1.3))+theme(axis.title.x = element_text(size = 0))
pdimmdepth1
pdimmdepth1<-ggplot(dmix[dmix$V1=="nDi.2.2.scaf00003",],aes(x=V3,y=foldd))+geom_point(size=0.2,alpha=0.4,color="tan")+
  geom_smooth(fill = NA, size=0.5,color="tan") +labs(x="coords",y="depth",title="Dimm_scaf1")+theme+
  scale_y_continuous(limits = c(0,1.3))+theme(axis.title.x = element_text(size = 0))
pdimmdepth1
pdimmdepth1<-ggplot(dmix[dmix$V1=="nDi.2.2.scaf00007",],aes(x=V3,y=foldd))+geom_point(size=0.2,alpha=0.4,color="tan")+
  geom_smooth(fill = NA, size=0.5,color="tan") +labs(x="coords",y="depth",title="Dimm_scaf1")+theme+
  scale_y_continuous(limits = c(0,1.3))+theme(axis.title.x = element_text(size = 0))
pdimmdepth1
pdimmdepth1<-ggplot(dmix[dmix$V1=="nDi.2.2.scaf00021",],aes(x=V3,y=foldd))+geom_point(size=0.2,alpha=0.4,color="tan")+
  geom_smooth(fill = NA, size=0.5,color="tan") +labs(x="coords",y="depth",title="Dimm_scaf21")+theme+
  scale_y_continuous(limits = c(0,1.3))+theme(axis.title.x = element_text(size = 0))
pdimmdepth1


mix <- read.table("Dmed.m_1.fastq.gz.depth.10k-win.txt")
mix <- subset(mix, mix$V5 > 500)
mix_v4 <- as.numeric(mix$V4)
sum_mix_v4 <- sum(mix_v4)
sum_mix_v4
head(mix)
mix_v5 <- as.numeric(mix$V5)
sum_mix_v5 <- sum(mix_v5)
sum_mix_v5
mix$dmix <- mix_v4/(mix_v5*2)
dmix <- mix[c(1,2,3,8)]
#install.packages("ggplot2")
library(ggplot2)

#dmix <- subset(dmix,dmix<500)
dmix$foldd <- dmix$dmix/191
p39 <- ggplot(dmix,aes(x=foldd))+geom_density(colour="forestgreen")+
  labs(x="depth",title="Dracunculus medinensis")+scale_x_continuous(limits = c(0,2))+
  geom_vline(aes(xintercept = 0.4), linetype = "dashed", size = 0.4)+
  geom_vline(aes(xintercept = 0.67), linetype = "dashed", size = 0.4)+
  theme(axis.text.x = element_text(size = 6))+
  theme(axis.text.y = element_text(size = 6))+
  theme(axis.title.y = element_text(size = 0))+
  theme(axis.title.x = element_text(size = 0))+
  theme(title=element_text(size = 4.5,face = "bold.italic",color="#d95f02"))
p39
pdf("Dmed.depth.peak.pdf",width=1.4,height=0.9)
p39
dev.off()
#
foldd <- as.numeric(dmix$foldd)
dmix <- subset(dmix,foldd<0.85 & foldd>0.5)
write.csv(dmix,file = "Dmed.window.csv")
#
head(dmix)
pdmeddepth1<-ggplot(dmix[dmix$V1=="DME_scaffold0000001",],aes(x=V3,y=foldd))+geom_point(size=0.2,alpha=0.4,color="salmon")+
  geom_smooth(fill = NA, size=0.5,color="salmon") +labs(x="coords",y="depth",title="Dmed_scaf1")+theme+
  scale_y_continuous(limits = c(0,1.3))+theme(axis.title.x = element_text(size = 0))
pdmeddepth1
pdmeddepth2<-ggplot(dmix[dmix$V1=="DME_scaffold0000002",],aes(x=V3,y=foldd))+geom_point(size=0.2,alpha=0.4,color="salmon")+
  geom_smooth(fill = NA, size=0.5,color="salmon") +labs(x="coords",y="depth",title="Dmed_scaf2")+theme+
  scale_y_continuous(limits = c(0,1.3))+theme(axis.title.x = element_text(size = 0))+theme(axis.title.y = element_text(size = 0))
pdmeddepth2
pdmeddepth3<-ggplot(dmix[dmix$V1=="DME_scaffold0000003",],aes(x=V3,y=foldd))+geom_point(size=0.2,alpha=0.4,color="salmon")+
  geom_smooth(fill = NA, size=0.5,color="salmon") +labs(x="coords",y="depth",title="Dmed_scaf3")+theme+
  scale_y_continuous(limits = c(0,1.3))+theme(axis.title.x = element_text(size = 0))+theme(axis.title.y = element_text(size = 0))
pdmeddepth3
pdmeddepth4<-ggplot(dmix[dmix$V1=="DME_scaffold0000004",],aes(x=V3,y=foldd))+geom_point(size=0.2,alpha=0.4,color="salmon")+
  geom_smooth(fill = NA, size=0.5,color="salmon") +labs(x="coords",y="depth",title="Dmed_scaf4")+theme+
  scale_y_continuous(limits = c(0,1.3))+theme(axis.title.x = element_text(size = 0))+theme(axis.title.y = element_text(size = 0))
pdmeddepth4
pdmeddepth5<-ggplot(dmix[dmix$V1=="DME_scaffold0000005",],aes(x=V3,y=foldd))+geom_point(size=0.2,alpha=0.4,color="salmon")+
  geom_smooth(fill = NA, size=0.5,color="salmon") +labs(x="coords",y="depth",title="Dmed_scaf5")+theme+
  scale_y_continuous(limits = c(0,1.3))+theme(axis.title.x = element_text(size = 0))+theme(axis.title.y = element_text(size = 0))
pdmeddepth5
pdmeddepth<-multiplot(pdmeddepth1,pdmeddepth1,pdmeddepth2,pdmeddepth1,pdmeddepth3,pdmeddepth1,pdmeddepth4,pdmeddepth1,pdmeddepth5,cols=5)

mix <- read.table("Gpul_1.fastq.gz.depth.10k-win.txt")
mix <- subset(mix, mix$V5 > 500)
mix_v4 <- as.numeric(mix$V4)
sum_mix_v4 <- sum(mix_v4)
sum_mix_v4
head(mix)
mix_v5 <- as.numeric(mix$V5)
sum_mix_v5 <- sum(mix_v5)
sum_mix_v5
mix$dmix <- mix_v4/(mix_v5*2)
dmix <- mix[c(1,2,3,8)]
#install.packages("ggplot2")
library(ggplot2)

#dmix <- subset(dmix,dmix<120)
dmix$foldd <- dmix$dmix/32
p40 <- ggplot(dmix,aes(x=foldd))+geom_density(colour="tan")+
  labs(x="depth",title="Gongylonema pulchrum")+scale_x_continuous(limits = c(0,2))+
  geom_vline(aes(xintercept = 0.73), linetype = "dashed", size = 0.4)+
  theme(axis.text.x = element_text(size = 6))+
  theme(axis.text.y = element_text(size = 6))+
  theme(axis.title.y = element_text(size = 0))+
  theme(axis.title.x = element_text(size = 0))+
  theme(title=element_text(size = 4.5,face = "bold.italic",color="#d95f02"))
p40
pdepth1<-ggplot(dmix[dmix$V1=="GPUH_scaffold0000001",],aes(x=V3,y=foldd))+geom_point(size=0.2,alpha=0.4,color="tan")+
  geom_smooth(fill = NA, size=0.5,color="tan") +labs(x="coords",y="depth",title="Dmed_scaf1")+theme+
  scale_y_continuous(limits = c(0,1.3))+theme(axis.title.x = element_text(size = 0))
pdepth1
pdepth1<-ggplot(dmix[dmix$V1=="GPUH_scaffold0000002",],aes(x=V3,y=foldd))+geom_point(size=0.2,alpha=0.4,color="tan")+
  geom_smooth(fill = NA, size=0.5,color="tan") +labs(x="coords",y="depth",title="Dmed_scaf1")+theme+
  scale_y_continuous(limits = c(0,1.3))+theme(axis.title.x = element_text(size = 0))
pdepth1
pdepth1<-ggplot(dmix[dmix$V1=="GPUH_scaffold0000003",],aes(x=V3,y=foldd))+geom_point(size=0.2,alpha=0.4,color="tan")+
  geom_smooth(fill = NA, size=0.5,color="tan") +labs(x="coords",y="depth",title="Dmed_scaf1")+theme+
  scale_y_continuous(limits = c(0,1.3))+theme(axis.title.x = element_text(size = 0))
pdepth1
pdepth1<-ggplot(dmix[dmix$V1=="GPUH_scaffold0000004",],aes(x=V3,y=foldd))+geom_point(size=0.2,alpha=0.4,color="tan")+
  geom_smooth(fill = NA, size=0.5,color="tan") +labs(x="coords",y="depth",title="Dmed_scaf1")+theme+
  scale_y_continuous(limits = c(0,1.3))+theme(axis.title.x = element_text(size = 0))
pdepth1
pdepth1<-ggplot(dmix[dmix$V1=="GPUH_scaffold0000011",],aes(x=V3,y=foldd))+geom_point(size=0.2,alpha=0.4,color="tan")+
  geom_smooth(fill = NA, size=0.5,color="tan") +labs(x="coords",y="depth",title="Dmed_scaf1")+theme+
  scale_y_continuous(limits = c(0,1.3))+theme(axis.title.x = element_text(size = 0))
pdepth1

mix <- read.table("Lsig_1.fastq.gz.depth.10k-win.txt")
mix <- subset(mix, mix$V5 > 500)
mix_v4 <- as.numeric(mix$V4)
sum_mix_v4 <- sum(mix_v4)
sum_mix_v4
head(mix)
mix_v5 <- as.numeric(mix$V5)
sum_mix_v5 <- sum(mix_v5)
sum_mix_v5
mix$dmix <- mix_v4/(mix_v5*2)
#install.packages("ggplot2")
library(ggplot2)
dmix <- mix[c(1,8)]

#dmix <- subset(dmix,dmix<250)
dmix$foldd <- dmix$dmix/60
p41 <- ggplot(dmix,aes(x=foldd))+geom_density(colour="salmon")+
  labs(x="depth",title="Litomosoides sigmodontis")+scale_x_continuous(limits = c(0,2))+
  theme(axis.text.x = element_text(size = 6))+
  theme(axis.text.y = element_text(size = 6))+
  theme(axis.title.y = element_text(size = 0))+
  theme(axis.title.x = element_text(size = 0))+
  theme(title=element_text(size = 4.5,face = "bold.italic",color="#d95f02"))
p41


mix <- read.table("Lloa.fastq.gz.depth.10k-win.txt")
mix <- subset(mix, mix$V5 > 500)
mix_v4 <- as.numeric(mix$V4)
sum_mix_v4 <- sum(mix_v4)
sum_mix_v4
head(mix)
mix_v5 <- as.numeric(mix$V5)
sum_mix_v5 <- sum(mix_v5)
sum_mix_v5
mix$dmix <- mix_v4/(mix_v5*2)
#install.packages("ggplot2")
library(ggplot2)
dmix <- mix[c(1,8)]

#dmix <- subset(dmix,dmix<6)
dmix$foldd <- dmix$dmix/1.26
p42 <- ggplot(dmix,aes(x=foldd))+geom_density(colour="tan")+
  labs(x="depth",title="Loa loa")+scale_x_continuous(limits = c(0,2))+
  geom_vline(aes(xintercept = 0.49), linetype = "dashed", size = 0.4)+
  theme(axis.text.x = element_text(size = 6))+
  theme(axis.text.y = element_text(size = 6))+
  theme(axis.title.y = element_text(size = 0))+
  theme(axis.title.x = element_text(size = 0))+
  theme(title=element_text(size = 4.5,face = "bold.italic",color="#d95f02"))
p42


mix <- read.table("Ofle_1.fastq.gz.depth.10k-win.txt")
mix <- subset(mix, mix$V5 > 500)
mix_v4 <- as.numeric(mix$V4)
sum_mix_v4 <- sum(mix_v4)
sum_mix_v4
head(mix)
mix_v5 <- as.numeric(mix$V5)
sum_mix_v5 <- sum(mix_v5)
sum_mix_v5
mix$dmix <- mix_v4/(mix_v5*2)
#install.packages("ggplot2")
library(ggplot2)
dmix <- mix[c(1,2,3,8)]

#dmix <- subset(dmix,dmix<400)
dmix$foldd <- dmix$dmix/84.5
p43 <- ggplot(dmix,aes(x=foldd))+geom_density(colour="tan")+
  labs(x="depth",title="Onchocerca flexuosa")+scale_x_continuous(limits = c(0,2))+
  #geom_vline(aes(xintercept = 0.95), linetype = "dashed", size = 0.4)+
  theme(axis.text.x = element_text(size = 6))+
  theme(axis.text.y = element_text(size = 6))+
  theme(axis.title.y = element_text(size = 0))+
  theme(axis.title.x = element_text(size = 0))+
  theme(title=element_text(size = 4.5,face = "bold.italic",color="#d95f02"))
p43
pofledepth1<-ggplot(dmix[dmix$V1=="OFLC_scaffold0000001",],aes(x=V3,y=foldd))+geom_point(size=0.2,alpha=0.4,color="tan")+
  geom_smooth(fill = NA, size=0.5,color="tan") +labs(x="coords",y="depth",title="ofle_scaf1")+theme+
  scale_y_continuous(limits = c(0,1.3))+theme(axis.title.x = element_text(size = 0))
pofledepth1
pofledepth2<-ggplot(dmix[dmix$V1=="OFLC_scaffold0000002",],aes(x=V3,y=foldd))+geom_point(size=0.2,alpha=0.4,color="tan")+
  geom_smooth(fill = NA, size=0.5,color="tan") +labs(x="coords",y="depth",title="ofle_scaf2")+theme+
  scale_y_continuous(limits = c(0,1.3))+theme(axis.title.x = element_text(size = 0))+theme(axis.title.y = element_text(size = 0))
pofledepth2
pofledepth3<-ggplot(dmix[dmix$V1=="OFLC_scaffold0000003",],aes(x=V3,y=foldd))+geom_point(size=0.2,alpha=0.4,color="tan")+
  geom_smooth(fill = NA, size=0.5,color="tan") +labs(x="coords",y="depth",title="ofle_scaf3")+theme+
  scale_y_continuous(limits = c(0,1.3))+theme(axis.title.x = element_text(size = 0))+theme(axis.title.y = element_text(size = 0))
pofledepth3
pofledepth4<-ggplot(dmix[dmix$V1=="OFLC_scaffold0000004",],aes(x=V3,y=foldd))+geom_point(size=0.2,alpha=0.4,color="tan")+
  geom_smooth(fill = NA, size=0.5,color="tan") +labs(x="coords",y="depth",title="ofle_scaf4")+theme+
  scale_y_continuous(limits = c(0,1.3))+theme(axis.title.x = element_text(size = 0))+theme(axis.title.y = element_text(size = 0))
pofledepth4
pofledepth5<-ggplot(dmix[dmix$V1=="OFLC_scaffold0000005",],aes(x=V3,y=foldd))+geom_point(size=0.2,alpha=0.4,color="tan")+
  geom_smooth(fill = NA, size=0.5,color="tan") +labs(x="coords",y="depth",title="ofle_scaf5")+theme+
  scale_y_continuous(limits = c(0,1.3))+theme(axis.title.x = element_text(size = 0))+theme(axis.title.y = element_text(size = 0))
pofledepth5




mix <- read.table("Ooch_1.fastq.gz.depth.10k-win.txt") ###sure is male data
mix <- subset(mix, mix$V5 > 500)
mix_v4 <- as.numeric(mix$V4)
sum_mix_v4 <- sum(mix_v4)
sum_mix_v4
head(mix)
mix_v5 <- as.numeric(mix$V5)
sum_mix_v5 <- sum(mix_v5)
sum_mix_v5
mix$dmix <- mix_v4/(mix_v5*2)
#install.packages("ggplot2")
library(ggplot2)
dmix <- mix[c(1,2,3,8)]

#dmix <- subset(dmix,dmix<70)
dmix$foldd <- dmix$dmix/13
p44 <- ggplot(dmix,aes(x=foldd))+geom_density(colour="forestgreen")+
  labs(x="depth",title="Onchocerca ochengi")+scale_x_continuous(limits = c(0,2))+
  theme(axis.text.x = element_text(size = 6))+
  theme(axis.text.y = element_text(size = 6))+
  theme(axis.title.y = element_text(size = 0))+
  theme(axis.title.x = element_text(size = 0))+
  theme(title=element_text(size = 4.5,face = "bold.italic",color="#d95f02"))
p44

####Ovol new
mix <- read.table("Ovol_1.fastq.gz.depth.10k-win.txt")
mix <- subset(mix, mix$V5 > 500)
mix_v4 <- as.numeric(mix$V4)
sum_mix_v4 <- sum(mix_v4)
sum_mix_v4
head(mix)
mix_v5 <- as.numeric(mix$V5)
sum_mix_v5 <- sum(mix_v5)
sum_mix_v5
mix$dmix <- mix_v4/(mix_v5*2)
#install.packages("ggplot2")
library(ggplot2)
dmix <- mix[c(1,8)]

#dmix <- subset(dmix,dmix<180)
dmix$foldd <- dmix$dmix/3.1
p45 <- ggplot(dmix,aes(x=foldd))+geom_density(colour="forest green")+
  labs(x="depth",title="Onchocerca volvulus")+scale_x_continuous(limits = c(0,2))+
  geom_vline(aes(xintercept = 0.71), linetype = "dashed", size = 0.4)+
  theme(axis.text.x = element_text(size = 6))+
  theme(axis.text.y = element_text(size = 6))+
  theme(axis.title.y = element_text(size = 0))+
  theme(axis.title.x = element_text(size = 0))+
  theme(title=element_text(size = 4.5,face = "bold.italic",color="#d95f02"))
p45
foldd <- as.numeric(dmix$foldd)
dmix <- subset(dmix,foldd<0.83)
write.csv(dmix,file = "Ovol.window.csv")

mix <- read.table("Tcal_1.fastq.gz.depth.10k-win.txt")
mix <- subset(mix, mix$V5 > 500)
mix_v4 <- as.numeric(mix$V4)
sum_mix_v4 <- sum(mix_v4)
sum_mix_v4
head(mix)
mix_v5 <- as.numeric(mix$V5)
sum_mix_v5 <- sum(mix_v5)
sum_mix_v5
mix$dmix <- mix_v4/(mix_v5*2)
#install.packages("ggplot2")
library(ggplot2)
dmix <- mix[c(1,2,3,8)]

#dmix <- subset(dmix,dmix<550)
dmix$foldd <- dmix$dmix/115
p46 <- ggplot(dmix,aes(x=foldd))+geom_density(colour="forest green")+
  labs(x="depth",title="Thelazia callipedia")+scale_x_continuous(limits = c(0,2))+
  geom_vline(aes(xintercept = 0.85), linetype = "dashed", size = 0.4)+
  theme(axis.text.x = element_text(size = 6))+
  theme(axis.text.y = element_text(size = 6))+
  theme(axis.title.y = element_text(size = 0))+
  theme(axis.title.x = element_text(size = 0))+
  theme(title=element_text(size = 4.5,face = "bold.italic",color="#d95f02"))
p46


mix <- read.table("Wban.fastq.gz.depth.10k-win.txt")
mix <- subset(mix, mix$V5 > 500)
mix_v4 <- as.numeric(mix$V4)
sum_mix_v4 <- sum(mix_v4)
sum_mix_v4
head(mix)
mix_v5 <- as.numeric(mix$V5)
sum_mix_v5 <- sum(mix_v5)
sum_mix_v5
mix$dmix <- mix_v4/(mix_v5*2)
#install.packages("ggplot2")
library(ggplot2)
dmix <- mix[c(1,8)]

#dmix <- subset(dmix,dmix<150)
dmix$foldd <- dmix$dmix/65
p47 <- ggplot(dmix,aes(x=foldd))+geom_density(colour="tan")+
  labs(x="depth",title="Wuchereria bancrofti")+scale_x_continuous(limits = c(0,2))+
  geom_vline(aes(xintercept = 0.6), linetype = "dashed", size = 0.4)+
  theme(axis.text.x = element_text(size = 6))+
  theme(axis.text.y = element_text(size = 6))+
  theme(axis.title.y = element_text(size = 0))+
  theme(axis.title.x = element_text(size = 0))+
  theme(title=element_text(size = 4.5,face = "bold.italic",color="#d95f02"))
p47


mix <- read.table("Asim_1.fastq.gz.depth.10k-win.txt")
mix <- subset(mix, mix$V5 > 500)
mix_v4 <- as.numeric(mix$V4)
sum_mix_v4 <- sum(mix_v4)
sum_mix_v4
head(mix)
mix_v5 <- as.numeric(mix$V5)
sum_mix_v5 <- sum(mix_v5)
sum_mix_v5
mix$dmix <- mix_v4/(mix_v5*2)
#install.packages("ggplot2")
library(ggplot2)
dmix <- mix[c(1,2,3,8)]

#dmix <- subset(dmix,dmix<50)
dmix$foldd <- dmix$dmix/19
p48 <- ggplot(dmix,aes(x=foldd))+geom_density(colour="tan")+
  labs(x="depth",title="Anisakis simplex")+scale_x_continuous(limits = c(0,2))+
  geom_vline(aes(xintercept = 0.54), linetype = "dashed", size = 0.4)+
  theme(axis.text.x = element_text(size = 6))+
  theme(axis.text.y = element_text(size = 6))+
  theme(axis.title.y = element_text(size = 0))+
  theme(axis.title.x = element_text(size = 0))+
  theme(title=element_text(size = 4.5,face = "bold.italic",color="#d95f02"))
p48
pdepth1<-ggplot(dmix[dmix$V1=="ASIM_scaffold0000001",],aes(x=V3,y=foldd))+geom_point(size=0.2,alpha=0.4,color="tan")+
  geom_smooth(fill = NA, size=0.5,color="tan") +labs(x="coords",y="depth",title="asim_scaf1")+theme+
  scale_y_continuous(limits = c(0,1.3))+theme(axis.title.x = element_text(size = 0))
pdepth1
pdepth1<-ggplot(dmix[dmix$V1=="ASIM_scaffold0000002",],aes(x=V3,y=foldd))+geom_point(size=0.2,alpha=0.4,color="tan")+
  geom_smooth(fill = NA, size=0.5,color="tan") +labs(x="coords",y="depth",title="asim_scaf1")+theme+
  scale_y_continuous(limits = c(0,1.3))+theme(axis.title.x = element_text(size = 0))
pdepth1
pdepth1<-ggplot(dmix[dmix$V1=="ASIM_scaffold0000003",],aes(x=V3,y=foldd))+geom_point(size=0.2,alpha=0.4,color="tan")+
  geom_smooth(fill = NA, size=0.5,color="tan") +labs(x="coords",y="depth",title="asim_scaf1")+theme+
  scale_y_continuous(limits = c(0,1.3))+theme(axis.title.x = element_text(size = 0))
pdepth1
pdepth1<-ggplot(dmix[dmix$V1=="ASIM_scaffold0000004",],aes(x=V3,y=foldd))+geom_point(size=0.2,alpha=0.4,color="tan")+
  geom_smooth(fill = NA, size=0.5,color="tan") +labs(x="coords",y="depth",title="asim_scaf1")+theme+
  scale_y_continuous(limits = c(0,1.3))+theme(axis.title.x = element_text(size = 0))
pdepth1
pdepth1<-ggplot(dmix[dmix$V1=="ASIM_scaffold0000110",],aes(x=V3,y=foldd))+geom_point(size=0.2,alpha=0.4,color="tan")+
  geom_smooth(fill = NA, size=0.5,color="tan") +labs(x="coords",y="depth",title="asim_scaf1")+theme+
  scale_y_continuous(limits = c(0,1.3))+theme(axis.title.x = element_text(size = 0))
pdepth1


mix <- read.table("Alum.fastq.gz.depth.10k-win.txt")
mix <- subset(mix, mix$V5 > 500)
mix_v4 <- as.numeric(mix$V4)
sum_mix_v4 <- sum(mix_v4)
sum_mix_v4
head(mix)
mix_v5 <- as.numeric(mix$V5)
sum_mix_v5 <- sum(mix_v5)
sum_mix_v5
mix$dmix <- mix_v4/(mix_v5*2)
#install.packages("ggplot2")
library(ggplot2)
dmix <- mix[c(1,8)]
#dmix <- subset(dmix,dmix<63)
dmix$foldd <- dmix$dmix/26
p49 <- ggplot(dmix,aes(x=foldd))+geom_density(colour="forest green")+
  labs(x="depth",title="Ascaris lumbricoides")+scale_x_continuous(limits = c(0,2))+
  geom_vline(aes(xintercept = 0.4), linetype = "dashed", size = 0.4)+
  geom_vline(aes(xintercept = 0.63), linetype = "dashed", size = 0.4)+
  theme(axis.text.x = element_text(size = 6))+
  theme(axis.text.y = element_text(size = 6))+
  theme(axis.title.y = element_text(size = 0))+
  theme(axis.title.x = element_text(size = 0))+
  theme(title=element_text(size = 4.5,face = "bold.italic",color="#d95f02"))
p49


mix <- read.table("Asuum.m.fastq.gz.depth.10k-win.txt")
mix <- subset(mix, mix$V5 > 500)
mix_v4 <- as.numeric(mix$V4)
sum_mix_v4 <- sum(mix_v4)
sum_mix_v4
head(mix)
mix_v5 <- as.numeric(mix$V5)
sum_mix_v5 <- sum(mix_v5)
sum_mix_v5
mix$dmix <- mix_v4/(mix_v5*2)
#install.packages("ggplot2")
library(ggplot2)
dmix <- mix[c(1,2,3,8)]

#dmix <- subset(dmix,dmix<200)
dmix$foldd <- dmix$dmix/20
pasuum <- ggplot(dmix,aes(x=foldd))+geom_density(colour="forest green")+
  labs(x="depth",title="Ascaris suum")+scale_x_continuous(limits = c(0,2))+
  geom_vline(aes(xintercept = 0.68), linetype = "dashed", size = 0.4)+
  theme(axis.text.x = element_text(size = 6))+
  theme(axis.text.y = element_text(size = 6))+
  theme(axis.title.y = element_text(size = 0))+
  theme(axis.title.x = element_text(size = 0))+
  theme(title=element_text(size = 4.5,face = "bold.italic",color="#d95f02"))
pasuum
pasuumdepth1<-ggplot(dmix[dmix$V1=="AgB01",],aes(x=V3,y=foldd))+geom_point(size=0.2,alpha=0.4,color="forest green")+
  geom_smooth(fill = NA, size=0.5,color="forest green") +labs(x="coords",y="depth",title="asuum_scaf1")+theme+
  scale_y_continuous(limits = c(0,1.3))+theme(axis.title.x = element_text(size = 0))
pasuumdepth1
pasuumdepth2<-ggplot(dmix[dmix$V1=="AgB02",],aes(x=V3,y=foldd))+geom_point(size=0.2,alpha=0.4,color="forest green")+
  geom_smooth(fill = NA, size=0.5,color="forest green") +labs(x="coords",y="depth",title="asuum_scaf2")+theme+
  scale_y_continuous(limits = c(0,1.3))+theme(axis.title.x = element_text(size = 0))+theme(axis.title.y = element_text(size = 0))
pasuumdepth2
pasuumdepth3<-ggplot(dmix[dmix$V1=="AgB03",],aes(x=V3,y=foldd))+geom_point(size=0.2,alpha=0.4,color="forest green")+
  geom_smooth(fill = NA, size=0.5,color="forest green") +labs(x="coords",y="depth",title="asuum_scaf3")+theme+
  scale_y_continuous(limits = c(0,1.3))+theme(axis.title.x = element_text(size = 0))+theme(axis.title.y = element_text(size = 0))
pasuumdepth3
pasuumdepth4<-ggplot(dmix[dmix$V1=="AgB04",],aes(x=V3,y=foldd))+geom_point(size=0.2,alpha=0.4,color="forest green")+
  geom_smooth(fill = NA, size=0.5,color="forest green") +labs(x="coords",y="depth",title="asuum_scaf4")+theme+
  scale_y_continuous(limits = c(0,1.3))+theme(axis.title.x = element_text(size = 0))+theme(axis.title.y = element_text(size = 0))
pasuumdepth4

pdepth1 <-ggplot(dmix[dmix$V1=="AgB08X",],aes(x=V3,y=foldd))+geom_point(size=0.2,alpha=0.4,color="forest green")+geom_smooth(fill = NA, size=0.5,color="forest green") +labs(x="coords",y="depth",title="AgB08X")+theme+scale_y_continuous(limits = c(0,1.3))+theme(axis.title.x = element_text(size = 0))
pdepth2 <-ggplot(dmix[dmix$V1=="AgB09X",],aes(x=V3,y=foldd))+geom_point(size=0.2,alpha=0.4,color="forest green")+geom_smooth(fill = NA, size=0.5,color="forest green") +labs(x="coords",y="depth",title="AgB09X")+theme+scale_y_continuous(limits = c(0,1.3))+theme(axis.title.x = element_text(size = 0))
pdepth3 <-ggplot(dmix[dmix$V1=="AgB13X",],aes(x=V3,y=foldd))+geom_point(size=0.2,alpha=0.4,color="forest green")+geom_smooth(fill = NA, size=0.5,color="forest green") +labs(x="coords",y="depth",title="AgB13X")+theme+scale_y_continuous(limits = c(0,1.3))+theme(axis.title.x = element_text(size = 0))
pdepth4 <-ggplot(dmix[dmix$V1=="AgB15X",],aes(x=V3,y=foldd))+geom_point(size=0.2,alpha=0.4,color="forest green")+geom_smooth(fill = NA, size=0.5,color="forest green") +labs(x="coords",y="depth",title="AgB15X")+theme+scale_y_continuous(limits = c(0,1.3))+theme(axis.title.x = element_text(size = 0))
pdepth5 <-ggplot(dmix[dmix$V1=="AgB16X",],aes(x=V3,y=foldd))+geom_point(size=0.2,alpha=0.4,color="forest green")+geom_smooth(fill = NA, size=0.5,color="forest green") +labs(x="coords",y="depth",title="AgB16X")+theme+scale_y_continuous(limits = c(0,1.3))+theme(axis.title.x = element_text(size = 0))
pdepth6 <-ggplot(dmix[dmix$V1=="AgB22X",],aes(x=V3,y=foldd))+geom_point(size=0.2,alpha=0.4,color="forest green")+geom_smooth(fill = NA, size=0.5,color="forest green") +labs(x="coords",y="depth",title="AgB22X")+theme+scale_y_continuous(limits = c(0,1.3))+theme(axis.title.x = element_text(size = 0))
pdepth7 <-ggplot(dmix[dmix$V1=="AgB27X",],aes(x=V3,y=foldd))+geom_point(size=0.2,alpha=0.4,color="forest green")+geom_smooth(fill = NA, size=0.5,color="forest green") +labs(x="coords",y="depth",title="AgB27X")+theme+scale_y_continuous(limits = c(0,1.3))+theme(axis.title.x = element_text(size = 0))
pdepth8 <-ggplot(dmix[dmix$V1=="AgB28X",],aes(x=V3,y=foldd))+geom_point(size=0.2,alpha=0.4,color="forest green")+geom_smooth(fill = NA, size=0.5,color="forest green") +labs(x="coords",y="depth",title="AgB28X")+theme+scale_y_continuous(limits = c(0,1.3))+theme(axis.title.x = element_text(size = 0))
pdepth9 <-ggplot(dmix[dmix$V1=="AgB31X",],aes(x=V3,y=foldd))+geom_point(size=0.2,alpha=0.4,color="forest green")+geom_smooth(fill = NA, size=0.5,color="forest green") +labs(x="coords",y="depth",title="AgB31X")+theme+scale_y_continuous(limits = c(0,1.3))+theme(axis.title.x = element_text(size = 0))
pdepth10 <-ggplot(dmix[dmix$V1=="AgB35X",],aes(x=V3,y=foldd))+geom_point(size=0.2,alpha=0.4,color="forest green")+geom_smooth(fill = NA, size=0.5,color="forest green") +labs(x="coords",y="depth",title="AgB35X")+theme+scale_y_continuous(limits = c(0,1.3))+theme(axis.title.x = element_text(size = 0))
pdepth11 <-ggplot(dmix[dmix$V1=="AgB36X",],aes(x=V3,y=foldd))+geom_point(size=0.2,alpha=0.4,color="forest green")+geom_smooth(fill = NA, size=0.5,color="forest green") +labs(x="coords",y="depth",title="AgB36X")+theme+scale_y_continuous(limits = c(0,1.3))+theme(axis.title.x = element_text(size = 0))
pdepth12 <-ggplot(dmix[dmix$V1=="AgR009X",],aes(x=V3,y=foldd))+geom_point(size=0.2,alpha=0.4,color="forest green")+geom_smooth(fill = NA, size=0.5,color="forest green") +labs(x="coords",y="depth",title="AgR009X")+theme+scale_y_continuous(limits = c(0,1.3))+theme(axis.title.x = element_text(size = 0))
pdepth13 <-ggplot(dmix[dmix$V1=="AgR013X",],aes(x=V3,y=foldd))+geom_point(size=0.2,alpha=0.4,color="forest green")+geom_smooth(fill = NA, size=0.5,color="forest green") +labs(x="coords",y="depth",title="AgR013X")+theme+scale_y_continuous(limits = c(0,1.3))+theme(axis.title.x = element_text(size = 0))
pdepth14 <-ggplot(dmix[dmix$V1=="AgR014X",],aes(x=V3,y=foldd))+geom_point(size=0.2,alpha=0.4,color="forest green")+geom_smooth(fill = NA, size=0.5,color="forest green") +labs(x="coords",y="depth",title="AgR014X")+theme+scale_y_continuous(limits = c(0,1.3))+theme(axis.title.x = element_text(size = 0))
pdepth15 <-ggplot(dmix[dmix$V1=="AgR018X",],aes(x=V3,y=foldd))+geom_point(size=0.2,alpha=0.4,color="forest green")+geom_smooth(fill = NA, size=0.5,color="forest green") +labs(x="coords",y="depth",title="AgR018X")+theme+scale_y_continuous(limits = c(0,1.3))+theme(axis.title.x = element_text(size = 0))
pdepth16 <-ggplot(dmix[dmix$V1=="AgR019X",],aes(x=V3,y=foldd))+geom_point(size=0.2,alpha=0.4,color="forest green")+geom_smooth(fill = NA, size=0.5,color="forest green") +labs(x="coords",y="depth",title="AgR019X")+theme+scale_y_continuous(limits = c(0,1.3))+theme(axis.title.x = element_text(size = 0))
pdepth17 <-ggplot(dmix[dmix$V1=="AgR021X",],aes(x=V3,y=foldd))+geom_point(size=0.2,alpha=0.4,color="forest green")+geom_smooth(fill = NA, size=0.5,color="forest green") +labs(x="coords",y="depth",title="AgR021X")+theme+scale_y_continuous(limits = c(0,1.3))+theme(axis.title.x = element_text(size = 0))
pdepth1
pdepth2
pdepth3
pdepth4
pdepth5
pdepth6
pdepth7
pdepth8
pdepth9
pdepth10
pdepth11
pdepth12
pdepth13
pdepth14
pdepth15
pdepth16
pdepth17
pasuumdepth<-multiplot(pdepth1,pdepth2,pdepth3,pdepth4,pdepth5,pdepth6,pdepth7,pdepth8,pdepth9,pdepth10,pdepth11,pdepth12,pdepth13,pdepth14,pdepth15,pdepth16,cols=4)

mix <- read.table("Avit_1.fastq.gz.depth.10k-win.txt")
mix <- subset(mix, mix$V5 > 500)
mix_v4 <- as.numeric(mix$V4)
sum_mix_v4 <- sum(mix_v4)
sum_mix_v4
head(mix)
mix_v5 <- as.numeric(mix$V5)
sum_mix_v5 <- sum(mix_v5)
sum_mix_v5
mix$dmix <- mix_v4/(mix_v5*2)
#install.packages("ggplot2")
library(ggplot2)
dmix <- mix[c(1,8)]

#dmix <- subset(dmix,dmix<200)
dmix$foldd <- dmix$dmix/75
p50 <- ggplot(dmix,aes(x=foldd))+geom_density(colour="salmon")+
  labs(x="depth",title="Acanthocheilonema viteae")+scale_x_continuous(limits = c(0,2))+
  #geom_vline(aes(xintercept = 0.75), linetype = "dashed", size = 0.4)+
  theme(axis.text.x = element_text(size = 6))+
  theme(axis.text.y = element_text(size = 6))+
  theme(axis.title.y = element_text(size = 0))+
  theme(axis.title.x = element_text(size = 0))+
  theme(title=element_text(size = 4.5,face = "bold.italic",color="#d95f02"))
p50



mix <- read.table("Tcan_1.fastq.gz.depth.10k-win.txt")
mix <- subset(mix, mix$V5 > 1500)
mix_v4 <- as.numeric(mix$V4)
sum_mix_v4 <- sum(mix_v4)
sum_mix_v4
head(mix)
mix_v5 <- as.numeric(mix$V5)
sum_mix_v5 <- sum(mix_v5)
sum_mix_v5
mix$dmix <- mix_v4/(mix_v5*2)
#install.packages("ggplot2")
library(ggplot2)
dmix <- mix[c(1,8)]

#dmix <- subset(dmix,dmix<40)
dmix$foldd <- dmix$dmix/36
p51 <- ggplot(dmix,aes(x=foldd))+geom_density(colour="forest green")+
  labs(x="depth",title="Toxocara canis")+scale_x_continuous(limits = c(0,2))+
  geom_vline(aes(xintercept = 0.44), linetype = "dashed", size = 0.4)+
  geom_vline(aes(xintercept = 0.64), linetype = "dashed", size = 0.4)+
  theme(axis.text.x = element_text(size = 6))+
  theme(axis.text.y = element_text(size = 6))+
  theme(axis.title.y = element_text(size = 0))+
  theme(axis.title.x = element_text(size = 0))+
  theme(title=element_text(size = 4.5,face = "bold.italic",color="#d95f02"))
p51


mix <- read.table("Ever.fastq.gz.depth.10k-win.txt")
mix <- subset(mix, mix$V5 > 500)
mix_v4 <- as.numeric(mix$V4)
sum_mix_v4 <- sum(mix_v4)
sum_mix_v4
head(mix)
mix_v5 <- as.numeric(mix$V5)
sum_mix_v5 <- sum(mix_v5)
sum_mix_v5
mix$dmix <- mix_v4/(mix_v5*2)
#install.packages("ggplot2")
library(ggplot2)
dmix <- mix[c(1,2,3,8)]

#dmix <- subset(dmix,dmix<60)
dmix$foldd <- dmix$dmix/13.2
p52 <- ggplot(dmix,aes(x=foldd))+geom_density(colour="tan")+
  labs(x="depth",title="Enterobius vermicularis")+scale_x_continuous(limits = c(0,2))+
  theme(axis.text.x = element_text(size = 6))+
  theme(axis.text.y = element_text(size = 6))+
  theme(axis.title.y = element_text(size = 0))+
  theme(axis.title.x = element_text(size = 0))+
  theme(title=element_text(size = 4.5,face = "bold.italic",color="#d95f02"))
p52


mix <- read.table("Smur_1.fastq.gz.depth.10k-win.txt")
mix <- subset(mix, mix$V5 > 500)
mix_v4 <- as.numeric(mix$V4)
sum_mix_v4 <- sum(mix_v4)
sum_mix_v4
head(mix)
mix_v5 <- as.numeric(mix$V5)
sum_mix_v5 <- sum(mix_v5)
sum_mix_v5
mix$dmix <- mix_v4/(mix_v5*2)
#install.packages("ggplot2")
library(ggplot2)
dmix <- mix[c(1,8)]

#dmix <- subset(dmix,dmix<150)
dmix$foldd <- dmix$dmix/40
p53 <- ggplot(dmix,aes(x=foldd))+geom_density(colour="tan")+
  labs(x="depth",title="Syphacia muris")+scale_x_continuous(limits = c(0,2))+
  theme(axis.text.x = element_text(size = 6))+
  theme(axis.text.y = element_text(size = 6))+
  theme(axis.title.y = element_text(size = 0))+
  theme(axis.title.x = element_text(size = 0))+
  theme(title=element_text(size = 4.5,face = "bold.italic",color="#d95f02"))
p53

#pIII<-multiplot(p35,p36,p37,p38,p39,p40,p41,p42,p43,p44,p45,p46,p47,p48,p49,pasuum,p50,p51,p52,p53,cols=4)

#pIII<-multiplot(p35,p46,p50,p51,p36,p42,p48,p47,p41,p44,p49,p37,p53,p38,p43,cols=4)

V
setwd("/Users/yifeng/Library/Mobile Documents/com~apple~CloudDocs/project/81worm/depthfile/V/")
mix <- read.table("Cgol_1.fastq.gz.depth.10k-win.txt")
mix <- subset(mix, mix$V5 > 500)
mix_v4 <- as.numeric(mix$V4)
sum_mix_v4 <- sum(mix_v4)
sum_mix_v4
head(mix)
mix_v5 <- as.numeric(mix$V5)
sum_mix_v5 <- sum(mix_v5)
sum_mix_v5
mix$dmix <- mix_v4/(mix_v5*2)
#install.packages("ggplot2")
library(ggplot2)
dmix <- mix[c(1,8)]

#dmix <- subset(dmix,dmix<600)
dmix$foldd <- dmix$dmix/14.5
p53_1 <- ggplot(dmix,aes(x=foldd))+geom_density(colour="tan")+
  labs(x="depth",title="Cylicostephanus goldi")+scale_x_continuous(limits = c(0,2))+
  theme(axis.text.x = element_text(size = 6))+
  theme(axis.text.y = element_text(size = 6))+
  theme(axis.title.y = element_text(size = 0))+
  theme(axis.title.x = element_text(size = 0))+
  theme(title=element_text(size = 4.5,face = "bold.italic",color="#1b9e77"))
p53_1

mix <- read.table("Name_1.fastq.gz.depth.10k-win.txt")
mix <- subset(mix, mix$V5 > 500)
mix_v4 <- as.numeric(mix$V4)
sum_mix_v4 <- sum(mix_v4)
sum_mix_v4
head(mix)
mix_v5 <- as.numeric(mix$V5)
sum_mix_v5 <- sum(mix_v5)
sum_mix_v5
mix$dmix <- mix_v4/(mix_v5*2)
#install.packages("ggplot2")
library(ggplot2)
dmix <- mix[c(1,8)]

#dmix <- subset(dmix,dmix<28)
dmix$foldd <- dmix$dmix/4.8
p54 <- ggplot(dmix,aes(x=foldd))+geom_density(colour="tan")+
  labs(x="depth",title="Necator americanus")+scale_x_continuous(limits = c(0,2))+
  theme(axis.text.x = element_text(size = 6))+
  theme(axis.text.y = element_text(size = 6))+
  theme(axis.title.y = element_text(size = 0))+
  theme(axis.title.x = element_text(size = 0))+
  theme(title=element_text(size = 4.5,face = "bold.italic",color="#1b9e77"))
p54

mix <- read.table("Nbra_1.fastq.gz.depth.10k-win.txt")
mix <- subset(mix, mix$V5 > 500)
mix_v4 <- as.numeric(mix$V4)
sum_mix_v4 <- sum(mix_v4)
sum_mix_v4
head(mix)
mix_v5 <- as.numeric(mix$V5)
sum_mix_v5 <- sum(mix_v5)
sum_mix_v5
mix$dmix <- mix_v4/(mix_v5*2)
#install.packages("ggplot2")
library(ggplot2)
dmix <- mix[c(1,8)]

#dmix <- subset(dmix,dmix<4.8)
dmix$foldd <- dmix$dmix/34
p55 <- ggplot(dmix,aes(x=foldd))+geom_density(colour="tan")+
  labs(x="depth",title="Nippostrongylus brasiliensis")+scale_x_continuous(limits = c(0,2))+
  geom_vline(aes(xintercept = 0.7), linetype = "dashed", size = 0.4)+
  theme(axis.text.x = element_text(size = 6))+
  theme(axis.text.y = element_text(size = 6))+
  theme(axis.title.y = element_text(size = 0))+
  theme(axis.title.x = element_text(size = 0))+
  theme(title=element_text(size = 4.5,face = "bold.italic",color="#1b9e77"))
p55

mix <- read.table("Acos_1.fastq.gz.depth.10k-win.txt")
mix <- subset(mix, mix$V5 > 500)
mix_v4 <- as.numeric(mix$V4)
sum_mix_v4 <- sum(mix_v4)
sum_mix_v4
head(mix)
mix_v5 <- as.numeric(mix$V5)
sum_mix_v5 <- sum(mix_v5)
sum_mix_v5
mix$dmix <- mix_v4/(mix_v5*2)
#install.packages("ggplot2")
library(ggplot2)
dmix <- mix[c(1,8)]

#dmix <- subset(dmix,dmix<160)
dmix$foldd <- dmix$dmix/35
p56 <- ggplot(dmix,aes(x=foldd))+geom_density(colour="salmon")+
  labs(x="depth",title="Angiostrongylus costaricensis")+scale_x_continuous(limits = c(0,2))+
  theme(axis.text.x = element_text(size = 6))+
  theme(axis.text.y = element_text(size = 6))+
  theme(axis.title.y = element_text(size = 0))+
  theme(axis.title.x = element_text(size = 0))+
  theme(title=element_text(size = 4.5,face = "bold.italic",color="#1b9e77"))
p56

mix <- read.table("Dviv_1.fastq.gz.depth.10k-win.txt")
mix <- subset(mix, mix$V5 > 500)
mix_v4 <- as.numeric(mix$V4)
sum_mix_v4 <- sum(mix_v4)
sum_mix_v4
head(mix)
mix_v5 <- as.numeric(mix$V5)
sum_mix_v5 <- sum(mix_v5)
sum_mix_v5
mix$dmix <- mix_v4/(mix_v5*2)
#install.packages("ggplot2")
library(ggplot2)
dmix <- mix[c(1,8)]

#dmix <- subset(dmix,dmix<18)
dmix$foldd <- dmix$dmix/5
p57 <- ggplot(dmix,aes(x=foldd))+geom_density(colour="forest green")+
  labs(x="depth",title="Dictyocaulus viviparus")+scale_x_continuous(limits = c(0,2))+
  geom_vline(aes(xintercept = 0.67), linetype = "dashed", size = 0.4)+
  theme(axis.text.x = element_text(size = 6))+
  theme(axis.text.y = element_text(size = 6))+
  theme(axis.title.y = element_text(size = 0))+
  theme(axis.title.x = element_text(size = 0))+
  theme(title=element_text(size = 4.5,face = "bold.italic",color="#1b9e77"))
p57
dmix1 <- subset(dmix,foldd<0.67 & foldd>0.37)


mix <- read.table("Hcon_1.fastq.gz.depth.10k-win.txt")
mix <- subset(mix, mix$V5 > 500)
mix_v4 <- as.numeric(mix$V4)
sum_mix_v4 <- sum(mix_v4)
sum_mix_v4
head(mix)
mix_v5 <- as.numeric(mix$V5)
sum_mix_v5 <- sum(mix_v5)
sum_mix_v5
mix$dmix <- mix_v4/(mix_v5*2)
#install.packages("ggplot2")
library(ggplot2)
dmix <- mix[c(1,8)]

#dmix <- subset(dmix,dmix<130)
dmix$foldd <- dmix$dmix/26.5
p58 <- ggplot(dmix,aes(x=foldd))+geom_density(colour="forest green")+
  labs(x="depth",title="Haemonchus contortus")+scale_x_continuous(limits = c(0,2))+
  geom_vline(aes(xintercept = 0.61), linetype = "dashed", size = 0.4)+
  theme(axis.text.x = element_text(size = 6))+
  theme(axis.text.y = element_text(size = 6))+
  theme(axis.title.y = element_text(size = 0))+
  theme(axis.title.x = element_text(size = 0))+
  theme(title=element_text(size = 4.5,face = "bold.italic",color="#1b9e77"))
p58

mix <- read.table("Acan_1.fastq.gz.depth.10k-win.txt")
mix <- subset(mix, mix$V5 > 500)
mix_v4 <- as.numeric(mix$V4)
sum_mix_v4 <- sum(mix_v4)
sum_mix_v4
head(mix)
mix_v5 <- as.numeric(mix$V5)
sum_mix_v5 <- sum(mix_v5)
sum_mix_v5
mix$dmix <- mix_v4/(mix_v5*2)
#install.packages("ggplot2")
library(ggplot2)
dmix <- mix[c(1,8)]

#dmix <- subset(dmix,dmix<100)
dmix$foldd <- dmix$dmix/23.8
p59 <- ggplot(dmix,aes(x=foldd))+geom_density(colour="forest green")+
  labs(x="depth",title="Angiostrongylus cantonensis")+scale_x_continuous(limits = c(0,2))+
  geom_vline(aes(xintercept = 0.65), linetype = "dashed", size = 0.4)+  ##0.35-0.65
  theme(axis.text.x = element_text(size = 6))+
  theme(axis.text.y = element_text(size = 6))+
  theme(axis.title.y = element_text(size = 0))+
  theme(axis.title.x = element_text(size = 0))+
  theme(title=element_text(size = 4.5,face = "bold.italic",color="#1b9e77"))
p59

mix <- read.table("Cele_1.fastq.gz.depth.10k-win.txt")
mix <- subset(mix, mix$V5 > 500)
mix_v4 <- as.numeric(mix$V4)
sum_mix_v4 <- sum(mix_v4)
sum_mix_v4
head(mix)
mix_v5 <- as.numeric(mix$V5)
sum_mix_v5 <- sum(mix_v5)
sum_mix_v5
mix$dmix <- mix_v4/(mix_v5*2)
#install.packages("ggplot2")
library(ggplot2)
dmix <- mix[c(1,8)]

#dmix <- subset(dmix,dmix<700)
dmix$foldd <- dmix$dmix/16
p60 <- ggplot(dmix,aes(x=foldd))+geom_density(colour="forest green")+
  labs(x="depth",title="Caenorhabditis elegans")+scale_x_continuous(limits = c(0,2))+
  theme(axis.text.x = element_text(size = 6))+
  theme(axis.text.y = element_text(size = 6))+
  theme(axis.title.y = element_text(size = 0))+
  theme(axis.title.x = element_text(size = 0))+
  theme(title=element_text(size = 4.5,face = "bold.italic",color="#1b9e77"))
p60
dmix1 <- subset(dmix,foldd>0.2 & foldd<0.5)


mix <- read.table("Svul_1.fastq.gz.depth.10k-win.txt")
mix <- subset(mix, mix$V5 > 500)
mix_v4 <- as.numeric(mix$V4)
sum_mix_v4 <- sum(mix_v4)
sum_mix_v4
head(mix)
mix_v5 <- as.numeric(mix$V5)
sum_mix_v5 <- sum(mix_v5)
sum_mix_v5
mix$dmix <- mix_v4/(mix_v5*2)
#install.packages("ggplot2")
library(ggplot2)
dmix <- mix[c(1,8)]

#dmix <- subset(dmix,dmix<700)
dmix$foldd <- dmix$dmix/27
psvul <- ggplot(dmix,aes(x=foldd))+geom_density(colour="forest green")+
  labs(x="depth",title="Strongylus vulgaris")+scale_x_continuous(limits = c(0,2))+
  geom_vline(aes(xintercept = 0.73), linetype = "dashed", size = 0.4)+
  theme(axis.text.x = element_text(size = 6))+
  theme(axis.text.y = element_text(size = 6))+
  theme(axis.title.y = element_text(size = 0))+
  theme(axis.title.x = element_text(size = 0))+
  theme(title=element_text(size = 4.5,face = "bold.italic",color="#1b9e77"))
psvul

mix <- read.table("Ppac_1.fastq.gz.depth.10k-win.txt")
mix <- subset(mix, mix$V5 > 500)
mix_v4 <- as.numeric(mix$V4)
sum_mix_v4 <- sum(mix_v4)
sum_mix_v4
head(mix)
mix_v5 <- as.numeric(mix$V5)
sum_mix_v5 <- sum(mix_v5)
sum_mix_v5
mix$dmix <- mix_v4/(mix_v5*2)
#install.packages("ggplot2")
library(ggplot2)
dmix <- mix[c(1,8)]

#dmix <- subset(dmix,dmix<700)
dmix$foldd <- dmix$dmix/17.9
p61 <- ggplot(dmix,aes(x=foldd))+geom_density(colour="forest green")+
  labs(x="depth",title="Pristionchus pacificus")+scale_x_continuous(limits = c(0,2))+
  geom_vline(aes(xintercept = 0.74), linetype = "dashed", size = 0.4)+
  theme(axis.text.x = element_text(size = 6))+
  theme(axis.text.y = element_text(size = 6))+
  theme(axis.title.y = element_text(size = 0))+
  theme(axis.title.x = element_text(size = 0))+
  theme(title=element_text(size = 4.5,face = "bold.italic", color="#1b9e77"))
p61

#pV <- multiplot(p53_1,p54,p55,p56, p57, p58,p59,psvul, p60,p61,cols=2)



      
pall<-multiplot(p2,p3,p4,p5,p6,psrod,p7,p8,p9,p10,p13,p24,p25,p26,p27,p28,p35,p36,p37,p38,p39,p40,p41,p42,p43,p44,p45,p46,p47,p48,p49,pasuum,p50,p51,p52,p53,p29,p30,p31,p32,p33,p34,psven,p53_1,p54,p55,p56, p57, p58,p59,psvul, p60,p61,cols=6)        

p100 <- multiplot(p1, p2, p3, p4 ,p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20, p21, p22 ,p23, p24, p25, p26, p27, p28, p29, p30, p31, p32, p33, p34, p35, p36, 
          p37, p38, p39, p40, p41, p42, p43, p44, p45, p46, p47, p48, p49, p50, p51, p52, p53, p54, p55, p56, p57, p58, p59, p60, cols=6)
pr <- multiplot(p31,p39,p45,p33_1,p33,p6,p30,p29_1,p32_1,cols=5)
