######排序！！！！
rm(list = ls())
library(ggplot2)
library(ggprism)

theme <- theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = 'white', colour = 'white'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
setwd("/Users/yifeng/Library/Mobile Documents/com~apple~CloudDocs/project/81worm/RNA/csin-sman/synteny/")
Smancsin<-read.table("Csin_Sman.ZW.CM030370_3.1.txt", header = FALSE)
head(Smancsin)

p<-ggplot(Smancsin)+
  geom_segment(data=(Smancsin[Smancsin$V1=="CM030370.1",]),aes(x=(V3+V2)/2,y=20.5,xend=(V5+V6)/2,yend=19),color="#EF3B2C",alpha=0.6,size=0.1)+
  geom_segment(data=(Smancsin[Smancsin$V1=="CM030373.1",]),aes(x=(V3+V2)/2+169711085,y=20.5,xend=(V5+V6)/2,yend=19),color="#3288BD",alpha=0.6,size=0.1)+
  xlab("")+ylab("")+
  theme_classic()
p
#8da0cb

p +
  geom_rect(xmin=0,xmax=168711085,ymin=20.6,ymax=20.7,fill="#EF3B2C")+#CM030370
  geom_rect(xmin=169711085,xmax=169711085+32400773,ymin=20.6,ymax=20.7,fill="#3288BD")+   ##CM030373
  geom_rect(xmin=0,xmax=10800000,ymin=18.8,ymax=18.9,fill="#fb6a4a")+
  geom_rect(xmin=10800000,xmax=23402947,ymin=18.8,ymax=18.9,fill="#ef3b2c")+
  geom_rect(xmin=23402947,xmax=44200000,ymin=18.8,ymax=18.9,fill="#cb181d")+ 
  geom_rect(xmin=44200000,xmax=88385488,ymin=18.8,ymax=18.9,fill="#fb6a4a")+  ##Sman
  annotate(geom="text", x=-3000, y=21.3,label = "CM030370")+annotate(geom="text", x=80086874, y=21.3,label = "CM030373")+
  annotate(geom="text", x=-3000, y=18.3,label = "SmanZW")+
  geom_vline(xintercept = 27767564)+ ##(fox-1)
  geom_vline(xintercept = 35898301)+ ##(mag-1)
  geom_vline(xintercept = 29797040)+ ##(U2AF2)
  theme_classic()

