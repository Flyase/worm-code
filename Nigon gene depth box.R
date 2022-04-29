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
##sort spe
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

##NC
my_comparisons <- list(c("NAA","NC"),c("NB","NC"),c("ND","NC"),c("NE","NC"),c("NX","NC"),c("NN","NC"))
pdepth1<-pdepth+
  geom_signif(comparisons = my_comparisons, y_position=1.6,test = wilcox.test,tip_length = 0, size = 0.4,map_signif_level=TRUE)+coord_cartesian(ylim = c(-1.3,1.8))

#######median value
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
