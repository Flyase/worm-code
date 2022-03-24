library(ggplot2)
library(ggsignif)
library(ggprism)
library(reshape2)
library(ggpubr)
library(dplyr)

setwd("/Users/yifeng/Library/Mobile Documents/com~apple~CloudDocs/project/81worm/depthfile/III/ortholog/normalized cov/")
###
folder <- "/Users/yifeng/Library/Mobile Documents/com~apple~CloudDocs/project/81worm/depthfile/III/ortholog/normalized cov/"      
file_list <- list.files(path=folder, pattern="*.nigon.nor.depth")                              
# read in each .txt file in file_list and rbind them into a data frame called data 
combine <- 
  do.call("rbind", 
          lapply(list.files(pattern="*.nigon.nor.depth"),
                 function(x) 
                   read.table(paste(folder, x, sep=''), 
                              header = FALSE, 
                              stringsAsFactors = FALSE)))

combine<-combine[c(1,11,9,2,10)]
combine<-data.frame(combine)
colnames(combine)<-c("gene","norcov","group","scaf","spe")
combine$group= factor(combine$group, levels=c('NAA','NB','NC','ND','NE','NX','NN'))
##对spe顺序排序
combine$spe= factor(combine$spe, levels=c('Asim','Asuum','Alum','Tcan','Dmed','Tcal','Dimm','Ofle','Ooch','Ovol','Wban','Bpah','Btim','Bmal','Ever','Smur','Gpul'))
pdepth <- ggplot(combine,aes(x=group,y=log2(norcov),colour=group))+geom_boxplot(size=0.8,alpha=0.8,outlier.shape= NA)+
  #labs(y="Log2(M/F depth)")+
  #geom_hline(aes(yintercept=5))+
  labs(title="Depth")+
  scale_color_manual(values=c("#9ECAE1","#74C476","#FEE5D9","#BCBDDC","#C7E9C0","#FCBBA1","#6BAED6"))+theme(axis.title.x = element_text(size = 0))+
  theme_classic()+theme(legend.position = "none")+
  #facet_grid(spe ~ .)
  facet_wrap( ~ spe, ncol=5)
pdepth
##sig
my_comparisons <- list(c("NAA","NC"),c("NB","NC"),c("ND","NC"),c("NE","NC"),c("NX","NC"),c("NN","NC"))
pdepth+
  geom_signif(comparisons = my_comparisons, y_position=8,test = wilcox.test,tip_length = 0, size = 0.4,map_signif_level=TRUE)


###rm DNA elimination
#combine<-subset(combine,norcov>0.2)
combine<-subset(combine,norcov>0.3)

###根据不同的coverage分布，分成四行，避免box被压缩的太短
##NC
my_comparisons <- list(c("NAA","NC"),c("NB","NC"),c("ND","NC"),c("NE","NC"),c("NX","NC"),c("NN","NC"))
pdepth1<-pdepth+
  geom_signif(comparisons = my_comparisons, y_position=1.6,test = wilcox.test,tip_length = 0, size = 0.4,map_signif_level=TRUE)+coord_cartesian(ylim = c(-1.3,1.8))

#######中位数，效果不明显
depth_median<-aggregate(combine$norcov,by = list(combine$spe,combine$group),FUN = median)
depth_median<-depth_median[grepl("NC", depth_median$Group.2),]
depth_median<-as.data.frame(depth_median)
colnames(depth_median)<-c("spe","group","norcov")
depth_median$group= factor(depth_median$group, levels=c('NAA','NB','NC','ND','NE','NX','NN'))
depth_median$spe= factor(depth_median$spe, levels=c('Asim','Asuum','Alum','Tcan','Dmed','Tcal','Dimm','Ofle','Ooch','Ovol','Wban','Bpah','Btim','Bmal','Ever','Smur','Gpul'))
pdepth1+geom_hline(aes(yintercept = log2(norcov)), depth_median,linetype="dotted",alpha=0.5,size=0.5)+
  #geom_hline(aes(yintercept = (log2(norcov))/2), depth_median,linetype="dotted",alpha=0.5)+
  theme(legend.position = "none")+coord_cartesian(ylim = c(-1.3,1.8))
#############













##NN
my_comparisons <- list(c("NAA","NN"),c("NB","NN"),c("NC","NN"),c("ND","NN"),c("NE","NN"),c("NX","NN"))
pdepth+
  geom_signif(comparisons = my_comparisons, y_position=1,test = wilcox.test,tip_length = 0, size = 0.4,map_signif_level=TRUE)+coord_cartesian(ylim = c(-1.5,1.5))
##NE
my_comparisons <- list(c("NAA","NE"),c("NB","NE"),c("NC","NE"),c("ND","NE"),c("NX","NE"),c("NN","NE"))
pdepth+
  geom_signif(comparisons = my_comparisons, y_position=1,test = wilcox.test,tip_length = 0, size = 0.4,map_signif_level=TRUE)+coord_cartesian(ylim = c(-1.5,1.5))
##NX
my_comparisons <- list(c("NAA","NX"),c("NB","NX"),c("NC","NX"),c("ND","NX"),c("NE","NX"),c("NN","NX"))
pdepth+
  geom_signif(comparisons = my_comparisons, y_position=1,test = wilcox.test,tip_length = 0, size = 0.4,map_signif_level=TRUE)+coord_cartesian(ylim = c(-1.5,1.5))
##NA
my_comparisons <- list(c("NB","NAA"),c("NC","NAA"),c("ND","NAA"),c("NE","NAA"),c("NX","NAA"),c("NN","NAA"))
pdepth+
  geom_signif(comparisons = my_comparisons, y_position=1,test = wilcox.test,tip_length = 0, size = 0.4,map_signif_level=TRUE)+coord_cartesian(ylim = c(-1.5,1.5))
##NA smur
my_comparisons <- list(c("NB","NAA"),c("NC","NAA"),c("ND","NAA"),c("NE","NAA"),c("NX","NAA"),c("NN","NAA"))
pdepth+
  geom_signif(comparisons = my_comparisons, y_position=0.2,test = wilcox.test,tip_length = 0, size = 0.4,map_signif_level=TRUE)+coord_cartesian(ylim = c(-0.3,0.3))
##NB
my_comparisons <- list(c("NAA","NB"),c("NC","NB"),c("ND","NB"),c("NE","NB"),c("NX","NB"),c("NN","NB"))
pdepth+
  geom_signif(comparisons = my_comparisons, y_position=1,test = wilcox.test,tip_length = 0, size = 0.4,map_signif_level=TRUE)+coord_cartesian(ylim = c(-1.5,1.5))
#######中位数，效果不明显
depth_median<-aggregate(combine$cov,by = list(combine$spe,combine$group),FUN = median)
depth_median<-depth_median[grepl("NC", depth_median$Group.2),]
depth_median<-as.data.frame(depth_median)
colnames(depth_median)<-c("spe","group","cov")
depth_median$group= factor(depth_median$group, levels=c('NAA','NB','NC','ND','NE','NX','NN'))
depth_median$spe= factor(depth_median$spe, levels=c('Alum','Dimm','Ooch','Asuum','Asim','Gpul','Btim','Tcan','Tcal','Ever','Bpah','Dmed','Ofle','Wban','Smur'))
pdepth1+geom_hline(aes(yintercept = log2(cov+1)), depth_median,linetype="dotted",alpha=0.5)+
  geom_hline(aes(yintercept = log2((cov+1)/2)), depth_median,linetype="dotted",alpha=0.5)+theme(legend.position = "none")
#############

pdepth <- ggplot(combine,aes(x=group,y=log2(cov+1),colour=group))+geom_boxplot(size=0.8,alpha=0.8,outlier.shape= NA)+
  scale_y_continuous(limits = c(4,8.5))+
  #labs(y="Log2(M/F depth)")+
  #geom_hline(aes(yintercept=5))+
  stat_summary(fun.y = "mean",geom = "point", shape = 8,size = 0.01)+labs(title="Depth")+
  scale_color_manual(values=c("#FCBBA1","#BCBDDC","#C7E9C0","#6BAED6","#969696"))+theme(axis.title.x = element_text(size = 0))+
  theme_classic()+theme(legend.position = "none")+
  #facet_grid(spe ~ .)
  facet_wrap( ~ spe, ncol=4)
pdepth
##sig
my_comparisons <- list(c("NX","Others"),c("ND","Others"),c("NE","Others"),c("NN","Others"))
pdepth+
  geom_signif(comparisons = my_comparisons, test = wilcox.test,y_position=c(8.2,7.8,7.4,7),tip_length = 0, size = 0.4,map_signif_level=TRUE)


pdepth <- ggplot(combine,aes(x=group,y=log2(cov+1),colour=group))+geom_boxplot(size=0.8,alpha=0.8,outlier.shape= NA)+
  scale_y_continuous(limits = c(7,10.5))+
  #labs(y="Log2(M/F depth)")+
  #geom_hline(aes(yintercept=5))+
  stat_summary(fun.y = "mean",geom = "point", shape = 8,size = 0.01)+labs(title="Depth")+
  scale_color_manual(values=c("#FCBBA1","#BCBDDC","#C7E9C0","#6BAED6","#969696"))+theme(axis.title.x = element_text(size = 0))+
  theme_classic()+theme(legend.position = "none")+
  #facet_grid(spe ~ .)
  facet_wrap( ~ spe, ncol=4)
pdepth
##sig
my_comparisons <- list(c("NX","Others"),c("ND","Others"),c("NE","Others"),c("NN","Others"))
pdepth+
  geom_signif(comparisons = my_comparisons, test = wilcox.test,y_position=c(10.2,9.8,9.4,9),tip_length = 0, size = 0.4,map_signif_level=TRUE)



pdepth <- ggplot(combine,aes(x=group,y=log2(cov+1),colour=group))+geom_boxplot(size=0.8,alpha=0.8,outlier.shape= NA)+
  scale_y_continuous(limits = c(6.2,8.5))+
  #labs(y="Log2(M/F depth)")+
  #geom_hline(aes(yintercept=5))+
  stat_summary(fun.y = "mean",geom = "point", shape = 8,size = 0.01)+labs(title="Depth")+
  scale_color_manual(values=c("#FCBBA1","#BCBDDC","#C7E9C0","#6BAED6","#969696"))+theme(axis.title.x = element_text(size = 0))+
  theme_classic()+theme(legend.position = "none")+
  #facet_grid(spe ~ .)
  facet_wrap( ~ spe, ncol=4)
pdepth
##sig
my_comparisons <- list(c("NX","Others"),c("ND","Others"),c("NE","Others"),c("NN","Others"))
pdepth+
  geom_signif(comparisons = my_comparisons, test = wilcox.test,y_position=c(8.3,8.1,7.9,7.7),tip_length = 0, size = 0.4,map_signif_level=TRUE)


pdepth <- ggplot(combine,aes(x=group,y=log2(cov+1),colour=group))+geom_boxplot(size=0.8,alpha=0.8,outlier.shape= NA)+
  scale_y_continuous(limits = c(6.2,6.82))+
  #labs(y="Log2(M/F depth)")+
  #geom_hline(aes(yintercept=5))+
  stat_summary(fun.y = "mean",geom = "point", shape = 8,size = 0.01)+labs(title="Depth")+
  scale_color_manual(values=c("#FCBBA1","#BCBDDC","#C7E9C0","#6BAED6","#969696"))+theme(axis.title.x = element_text(size = 0))+
  theme_classic()+theme(legend.position = "none")+
  #facet_grid(spe ~ .)
  facet_wrap( ~ spe, ncol=4)
pdepth
##sig
my_comparisons <- list(c("NX","Others"),c("ND","Others"),c("NE","Others"),c("NN","Others"))
pdepth+
  geom_signif(comparisons = my_comparisons, test = wilcox.test,y_position=c(6.8,6.73,6.66,6.6),tip_length = 0, size = 0.4,map_signif_level=TRUE)



###计数
library(dplyr)
head(combine)
str(combine)
combine$nigoncount<-1
nigon_count<-aggregate(combine$nigoncount,by = list(combine$group,combine$spe),FUN = sum)
