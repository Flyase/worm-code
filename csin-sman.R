library(ggplot2)
library(lattice)
library(plyr)
library(Rmisc)
library(ggsignif)
library(pheatmap)
library(patchwork)
library(reshape2)


theme <- theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = 'white', colour = 'white'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
setwd("/Users/yifeng/Library/Mobile Documents/com~apple~CloudDocs/project/81worm/RNA/csin-sman/")
Smanso<-read.table("so.sort.bam.featurecount.txt", header = TRUE)
Smanst<-read.table("st.sort.bam.featurecount.txt", header = TRUE)
Smanbo<-read.table("bo.sort.bam.featurecount.txt", header = TRUE)
Smanbt<-read.table("bt.sort.bam.featurecount.txt", header = TRUE)
Smansf<-read.table("sf.sort.bam.featurecount.txt", header = TRUE)
Smansm<-read.table("sm.sort.bam.featurecount.txt", header = TRUE)
Smanbf<-read.table("bf.sort.bam.featurecount.txt", header = TRUE)
Smanbm<-read.table("bm.sort.bam.featurecount.txt", header = TRUE)
Csinmucle<-read.table("ERR604979_1.fastq.gz.featurecount.txt", header = FALSE)
Csinmg<-read.table("ERR604978_1.fastq.gz.featurecount.txt", header = FALSE)
Csinfg<-read.table("ERR604981_1.fastq.gz.featurecount.txt", header = FALSE)
#Cele<-read.table("Cele.fem.count.length.txt", header = FALSE)


Sman<-data.frame(row.names=Smanso$Geneid, so=Smanso$so.sort.bam, st=Smanst$st.sort.bam,bo=Smanbo$bo.sort.bam, bt=Smanbt$bt.sort.bam,sf=Smansf$sf.sort.bam, sm=Smansm$sm.sort.bam,bf=Smanbf$bf.sort.bam, bm=Smanbm$bm.sort.bam)
tail(Sman)
head(Sman)
nrow(Sman)
genelength1<-read.table("Sman.genelength.transcript.txt",header = FALSE)
genelength1<-data.frame(row.names=genelength1$V1,L=genelength1$V2)
head(genelength1)
nrow(genelength1)
nrow(Sman)
Smanlength<-merge(genelength1,Sman,by=c("row.names"),all.x = TRUE)
Smanlength<-subset(Smanlength,Smanlength$L>0)
nrow(Smanlength)
head(Smanlength)
####################
rpkm <- function(counts, lengths) {
  rate <- counts / lengths 
  rate / sum(counts) * 1e9
}

genes <- data.frame(Gene=Smanlength[,1],
                    Length=Smanlength[,2]
)

counts <- data.frame(Smanlength[,3:10]
)

rpkms <- apply(counts, 2, function(x) rpkm(x, genes$Length))
RPKM1 <- cbind(genes, rpkms)

tpm <- function(counts, lengths) {
  rate <- counts / lengths
  rate / sum(rate) * 1e6
}
tpms <- apply(counts, 2, function(x) tpm(x, genes$Length))
TPM1 <- cbind(genes, tpms)
#write.table(TPM, file="CD_musclegene_TPM.txt",sep = '\t',row.names = F,quote = F)

#################################################################################
Csin<-data.frame(row.names=Csinmucle$V1, mucle=Csinmucle$V2,mg=Csinmg$V2, fg=Csinfg$V2)
genelength1<-read.table("Csin.longest.trans.genelength",header = FALSE)
genelength1<-data.frame(row.names=genelength1$V1,L=genelength1$V2)
head(genelength1)
nrow(genelength1)
nrow(Csin)
Csinlength<-merge(genelength1,Csin,by=c("row.names"),all.x = TRUE)
Csinlength<-subset(Csinlength,Csinlength$L>0)
nrow(Csinlength)
head(Csinlength)

genes <- data.frame(Gene=Csinlength[,1],
                    Length=Csinlength[,2]
)

counts <- data.frame(Csinlength[,3:5]
)

rpkms <- apply(counts, 2, function(x) rpkm(x, genes$Length))
RPKM2 <- cbind(genes, rpkms)

tpms <- apply(counts, 2, function(x) tpm(x, genes$Length))
TPM2 <- cbind(genes, tpms)



#######subtract Csin candidate gonad-high exp gene
head(RPKM2)
rownames(RPKM2)<-RPKM2$Gene
Csin_f_high<-subset(RPKM2,2*(RPKM2$mucle)<(RPKM2$fg) & RPKM2$fg>2)
Csin_m_high<-subset(RPKM2,2*(RPKM2$mucle)<(RPKM2$mg) & RPKM2$mg>2)
Csingonadhigh1<-rbind(Csin_f_high,Csin_m_high)
Csingonadhigh <- Csingonadhigh1[!duplicated(Csingonadhigh1$Gene),]  

#######subtract Csin candidate gonad-low exp gene
Csin_f_low<-subset(RPKM2,(RPKM2$mucle)>2*(RPKM2$fg) & RPKM2$mucle>2 & RPKM2$fg<1)
Csin_m_low<-subset(RPKM2,(RPKM2$mucle)>2*(RPKM2$mg) & RPKM2$mucle>2 & RPKM2$mg<1)
Csingonadlow1<-rbind(Csin_f_low,Csin_m_low)
Csingonadlow <- Csingonadlow1[!duplicated(Csingonadlow1$Gene),]  #
Csingonadlow = Csingonadlow[!rownames(Csingonadlow) %in% rownames(Csingonadhigh),] 

#######subtract Csin candidate gonad-soma similiar exp gene
Csingonad<-rbind(Csingonadhigh,Csingonadlow)
Csingonad_dup <- Csingonad[duplicated(Csingonad$Gene),]  ##check dup gene, must be 0
Csingonad_soma = RPKM2[!rownames(RPKM2) %in% rownames(Csingonad),] 

#######merge with ortholog#########
A<-read.table("Csin__v__Sman.11.csv", header = TRUE)
head(A)
SmanA<-data.frame(Gene=A[,2])
CsinA<-data.frame(Gene=A[,1])
rownames(SmanA)<-SmanA$Gene
rownames(CsinA)<-CsinA$Gene
rownames(TPM1)<-TPM1$Gene
rownames(TPM2)<-TPM2$Gene
SmanAtpm<-merge(SmanA,TPM1,by=c("row.names"),sort=FALSE)
CsinAtpm<-merge(CsinA,TPM2,by=c("row.names"),sort=FALSE)
combine<-cbind(SmanAtpm,CsinAtpm)
head(combine)

combine1<-combine[c(1,5,6,7,8,9,10,11,12,13,18,19,17)]
#names(combine)<-c("smangene","smanlength","Smanso","Smanst","Smanbo","Smanbt","Smanbf","Smanbm","Csingene","Csinlength","Csinon","Csinmucle")
names(combine1)<-c("smangene","Smanso","Smanst","Smanbo","Smanbt","Smansf","Smansm","Smanbf","Smanbm","Csingene","Csinmg","Csinfg","Csinmucle")
rownames(combine1)<-combine1$Csingene


#################################################！！！！！3 parts##################################################
##csin gonad high
combine_csingonad_high = combine1[rownames(combine1) %in% rownames(Csingonadhigh),] 
head(combine_csingonad_high)

##csin gonad low
combine_csingonad_low = combine1[rownames(combine1) %in% rownames(Csingonadlow),] 
head(combine_csingonad_low)

##rest, similiar with muscle
combine_csingonad_soma = combine1[rownames(combine1) %in% rownames(Csingonad_soma),] 
head(combine_csingonad_soma)

#######################density####################

############not just Csingonadhigh gene， density plot all genes,not just ortholog##########first step, check, selective step####
head(TPM1)
head(TPM2)
bo<-TPM1[c(1,5)]
bt<-TPM1[c(1,6)]
mg<-TPM2[c(1,4)]
fg<-TPM2[c(1,5)]
bo$group<-"bo"
bt$group<-"bt"
mg$group<-"mg"
fg$group<-"fg"
names(bo)<-c("gene","tpm","group")
names(bt)<-c("gene","tpm","group")
names(mg)<-c("gene","tpm","group")
names(fg)<-c("gene","tpm","group")

combine11<-rbind(bo,bt,mg,fg)
head(combine11)
combine11$group=factor(combine11$group, levels=c('bo','bt','fg','mg'))
pexp1<-ggplot(combine11,aes(x=log2(tpm+1),color=group))+ geom_density(alpha=0.2,size=1)+ scale_color_manual(values=c("#EF3B2C","#3288BD","#636363","#000000"))+
  theme_classic()
pexp1   ####sman mf both up


############Csingonadhigh gene，ortholog############
head(combine1)
bo<-combine1[c(4)]
bt<-combine1[c(5)]
mg<-combine1[c(11)]
fg<-combine1[c(12)]

##
head(combine_csingonad_high)
bo<-combine_csingonad_high[c(4)]
bt<-combine_csingonad_high[c(5)]
mg<-combine_csingonad_high[c(11)]
fg<-combine_csingonad_high[c(12)]
##

##
head(combine_csingonad_low)
bo<-combine_csingonad_low[c(4)]
bt<-combine_csingonad_low[c(5)]
mg<-combine_csingonad_low[c(11)]
fg<-combine_csingonad_low[c(12)]
##

##
head(combine_csingonad_soma)
bo<-combine_csingonad_soma[c(4)]
bt<-combine_csingonad_soma[c(5)]
mg<-combine_csingonad_soma[c(11)]
fg<-combine_csingonad_soma[c(12)]
##

bo$group<-"bo"
bt$group<-"bt"
mg$group<-"mg"
fg$group<-"fg"
names(bo)<-c("tpm","group")
names(bt)<-c("tpm","group")
names(mg)<-c("tpm","group")
names(fg)<-c("tpm","group")

combine111<-rbind(bo,bt,mg,fg)
head(combine111)
combine111$group=factor(combine111$group, levels=c('bo','bt','fg','mg'))
pexp1<-ggplot(combine111,aes(x=log2(tpm+1),color=group))+ geom_density(alpha=0.2,size=1)+scale_color_manual(values=c("#EF3B2C","#3288BD","#636363","#000000"))+
  theme_classic()
pexp1   ####sman f down, m both slightly down & up

########soma ###ortholog###
head(combine1)
bf<-combine1[c(8)]
bm<-combine1[c(9)]
muscle<-combine1[c(13)]
bf$group<-"bf"
bm$group<-"bm"
muscle$group<-"muscle"
names(bf)<-c("tpm","group")
names(bm)<-c("tpm","group")
names(muscle)<-c("tpm","group")
combine111soma<-rbind(bf,bm,muscle)
head(combine111soma)
combine111soma$group=factor(combine111soma$group, levels=c('bf','bm','muscle'))
pexp1<-ggplot(combine111soma,aes(x=log2(tpm+1),color=group))+ geom_density(alpha=0.2,size=1)+scale_color_manual(values=c("#EF3B2C","#3288BD","#000000"))+
  theme_classic()
pexp1 
####soma #####gonad high#####
head(combine_csingonad_high)
bf<-combine_csingonad_high[c(8)]
bm<-combine_csingonad_high[c(9)]
muscle<-combine_csingonad_high[c(13)]
bf$group<-"bf"
bm$group<-"bm"
muscle$group<-"muscle"
names(bf)<-c("tpm","group")
names(bm)<-c("tpm","group")
names(muscle)<-c("tpm","group")
combine111soma<-rbind(bf,bm,muscle)
head(combine111soma)
combine111soma$group=factor(combine111soma$group, levels=c('bf','bm','muscle'))
pexp1<-ggplot(combine111soma,aes(x=log2(tpm+1),color=group))+ geom_density(alpha=0.2,size=1)+scale_color_manual(values=c("#EF3B2C","#3288BD","#000000"))+
  theme_classic()
pexp1   ####sman f down, m both slightly down & up

#######gonad ratio #######
combine1$sosf<-log2((combine1$Smanso+0.01)/(combine1$Smansf+0.01))
combine1$stsm<-log2((combine1$Smanst+0.01)/(combine1$Smansm+0.01))
combine1$bobf<-log2((combine1$Smanbo+0.01)/(combine1$Smanbf+0.01))
combine1$btbm<-log2((combine1$Smanbt+0.01)/(combine1$Smanbm+0.01))
combine1$fgmu<-log2((combine1$Csinfg+0.01)/(combine1$Csinmucle+0.01))
combine1$mgmu<-log2((combine1$Csinmg+0.01)/(combine1$Csinmucle+0.01))
combine12<-combine1[c(1,16,17,18,19)]

#
combine12 = combine12[rownames(combine12) %in% rownames(Csingonadhigh),] 
#

#
combine12 = combine12[rownames(combine12) %in% rownames(Csingonadlow),] 
#

#
combine12 = combine12[rownames(combine12) %in% rownames(Csingonad_soma),] 
#

combine12 <- melt(combine12, id.vars=c("smangene"))
combine12$variable=factor(combine12$variable, levels=c('bobf','btbm','fgmu','mgmu'))
pexp2<-ggplot(combine12,aes(x=value,color=variable))+ geom_density(alpha=0.2,size=1)+ scale_color_manual(values=c("#EF3B2C","#3288BD","#636363","#000000"))+
  theme_classic()+coord_cartesian(xlim = c(-8,8))
pexp2


############################################pheatmap############################################

##########################******Sman gene#######################
head(combine1)
row.names(combine1)<-combine1$smangene
S0<-read.table("Sman.strata.list.txt", header = FALSE,row.names = 1)
head(S0)
combines0<-merge(combine1,S0,by= "row.names",sort=FALSE)
head(combines0)
#combines0<-combines0[c(2,3,4,5,6,7,8,10,11,12,13,14)]
combines01<-combines0[c(2,5,6,12,13,15,11)]##计算ratio后15应改为21
row.names(combines01)<-combines01$Csingene 

########calcu fold change percentage (mas,demas,fem,defem)#######
head(combines01)
combines01_demas<-subset(combines01,Csinmg>2*Smanbt)
combines01_mas<-subset(combines01,Csinmg*2<Smanbt)
combines01_defem<-subset(combines01,Csinfg>2*Smanbo)
combines01_fem<-subset(combines01,Csinfg*2<Smanbo)
combines01_demas$g<-"c5" #9ECAE1
combines01_defem$g<-"c6" ##FCBBA1
combines01_mas$g<-"c7" ##084594
combines01_fem$g<-"c8" #A65628
combines01_formassig<-rbind(combines01_demas,combines01_mas,combines01_defem,combines01_fem)
write.table(combines01_formassig,file="smangene.csin.expchange.formassig.txt",sep = '\t',row.names = F,quote = F) ##calcu mas sig, combined data
########################

#
combines0_gonadhigh = combines01[rownames(combines01) %in% rownames(Csingonadhigh),] 
#write.table(combines0_gonadhigh,file="smangene.csingonadhigh.strata.txt",sep = '\t',row.names = F,quote = F) ##not add mas.demas.fem.defem filtration 
#

#
combines0_gonadlow = combines01[rownames(combines01) %in% rownames(Csingonadlow),] 
#

#
combines0_gonad_soma = combines01[rownames(combines01) %in% rownames(Csingonad_soma),] 
#

######add cele######no#####
head(TPM3)
row.names(TPM3)<-TPM3$Gene
celesman<-read.table("Cele__v__Sman.11.csv", header = TRUE,row.names = 2)
celesman1<-merge(TPM3,celesman,by= "row.names",sort=FALSE)
row.names(celesman1)<-celesman1$Sman
celesmantmul<-merge(combines0,celesman1,by= "row.names",,all.x = TRUE,sort=FALSE)
head(celesmantmul)
S1<-data.frame(row.names=celesmantmul$smangene,Smanbt=log2(celesmantmul$Smanbt+1),Smanbo=log2(celesmantmul$Smanbo+1),Tmulimma=log2(celesmantmul$Tmulimma+1),Tmulpro=log2(celesmantmul$Tmulpro+1),Smanbm=log2(celesmantmul$Smanbm+1),Smanbf=log2(celesmantmul$Smanbf+1),Celeherma=log2(celesmantmul$fe+1),Celem=log2(celesmantmul$m+1))
S1[is.na(S1)] <- 0
##################no####
annotation_row <- data.frame(row.names=combines01$smangene,  group=factor(combines01$V2, levels=c('A','PAR','S1','S0')))
ann_colors = list(group =c(A = "#4daf4a", PAR = "#377eb8", S1 = "#ff7f00",S0="#e41a1c")
)

S1<-data.frame(row.names=combines01$smangene,Csinfg=log2(combines01$Csinfg+1),Csinmg=log2(combines01$Csinmg+1),Smanbo=log2(combines01$Smanbo+1),Smanbt=log2(combines01$Smanbt+1))
#S1<-S1[order(S1[,2],decreasing=TRUE),]
pheatmap(S1,annotation_row = annotation_row,show_rownames = FALSE,annotation_colors = ann_colors,
         cluster_cols=FALSE,border=FALSE,
         color = colorRampPalette(colors = c("darkblue","#ffffbf","red"))(100))


###############more significant pattern 2-fold,final used; no 3 parts###############
head(combines01)
combines0_1<-subset(combines01,Csinmg>2*Smanbt & Csinfg>2*Smanbo)
combines0_1<-combines0_1[order(-combines0_1[,4]),]
combines0_1$g<-"c1" #green
combines0_2<-subset(combines01,Csinmg>2*Smanbt & Csinfg*2<Smanbo) 
combines0_2<-combines0_2[order(-combines0_2[,2]),]
combines0_2$g<-"c2" ##blue
combines0_3<-subset(combines01,Csinmg*2<Smanbt & Csinfg>2*Smanbo) 
combines0_3<-combines0_3[order(-combines0_3[,3]),]
combines0_3$g<-"c3" ##orange
combines0_4<-subset(combines01,Csinmg*2<Smanbt & Csinfg*2<Smanbo) 
combines0_4<-combines0_4[order(-combines0_4[,3]),]
combines0_4$g<-"c4" #red
combines0_com<-rbind(combines0_4,combines0_3,combines0_2,combines0_1)
head(combines0_com)

#####
combines0_com = combines0_4[rownames(combines0_4) %in% rownames(Csingonadlow),] 
#

#write.table(combines0_1$smangene, file="smangene.demas_defemi.txt",sep = '\t',row.names = F,quote = F) ##demas and defemi
#write.table(combines0_3$smangene, file="smangene.mas_defemi.txt",sep = '\t',row.names = F,quote = F) ##mas and defemi
#write.table(combines0_4$smangene, file="smangene.mas_femi.txt",sep = '\t',row.names = F,quote = F) ##mas and femi
write.table(combines0_com,file="smangene.csin.expchange.txt",sep = '\t',row.names = F,quote = F) ##mas and femi

S11<-data.frame(row.names=combines0_com$smangene,Csinfg=log2(combines0_com$Csinfg+1),Csinmg=log2(combines0_com$Csinmg+1),Smanbo=log2(combines0_com$Smanbo+1),Smanbt=log2(combines0_com$Smanbt+1))
head(S11)
p<-pheatmap(S11,
            cluster_cols=FALSE,cluster_rows=FALSE,border=FALSE,#show_rownames = F,
            color = colorRampPalette(colors = c("#3288bd","#ffffbf","#EF3B2C"))(100),
            #legend = FALSE,
            fontsize_row=3)  ##not good


annotation_row <- data.frame(row.names=combines0_com$smangene,  group=factor(combines0_com$V2, levels=c('A','PAR','S1','S0')),g=factor(combines0_com$g, levels=c('c1','c2','c3','c4')))
ann_colors = list(group =c(A = "#4daf4a", PAR = "#377eb8", S1 = "#ff7f00",S0="#e41a1c"),
                  g =c(c1 = "#4daf4a", c2 = "#377eb8", c3 = "#ff7f00",c4="#e41a1c")
)

p<-pheatmap(S11,annotation_row = annotation_row,annotation_colors = ann_colors,
            cluster_cols=FALSE,border=FALSE,#show_rownames = F,
            color = colorRampPalette(colors = c("#3288bd","#ffffbf","#EF3B2C"))(100),
            #legend = FALSE,
            fontsize_row=3,
            annotation_legend=F)



#########3 parts##########
###
head(combines0_gonadhigh)
combines0_gonadhigh_1<-subset(combines0_gonadhigh,Csinmg>2*Smanbt & Csinfg>2*Smanbo)
combines0_gonadhigh_1<-combines0_gonadhigh_1[order(-combines0_gonadhigh_1[,4]),]
combines0_gonadhigh_1$g<-"c1" #green
combines0_gonadhigh_2<-subset(combines0_gonadhigh,Csinmg>2*Smanbt & Csinfg*2<Smanbo) 
combines0_gonadhigh_2<-combines0_gonadhigh_2[order(-combines0_gonadhigh_2[,2]),]
combines0_gonadhigh_2$g<-"c2" ##orange
combines0_gonadhigh_3<-subset(combines0_gonadhigh,Csinmg*2<Smanbt & Csinfg>2*Smanbo) 
combines0_gonadhigh_3<-combines0_gonadhigh_3[order(-combines0_gonadhigh_3[,3]),]
combines0_gonadhigh_3$g<-"c3" ##blue
combines0_gonadhigh_4<-subset(combines0_gonadhigh,Csinmg*2<Smanbt & Csinfg*2<Smanbo) 
combines0_gonadhigh_4<-combines0_gonadhigh_4[order(-combines0_gonadhigh_4[,3]),]
combines0_gonadhigh_4$g<-"c4" #red
combines0_gonadhigh_com<-rbind(combines0_gonadhigh_4,combines0_gonadhigh_3,combines0_gonadhigh_2,combines0_gonadhigh_1)
head(combines0_gonadhigh_com)
#write.table(combines0_gonadhigh_com,file="smangene.csingonadhigh.expchange.txt",sep = '\t',row.names = F,quote = F)
##############
annotation_row <- data.frame(row.names=combines0_gonadhigh_com1$smangene,g=factor(combines0_gonadhigh_com1$g, levels=c('c1','c2','c3','c4','c5','c6','c7','c8')),group=factor(combines0_gonadhigh_com1$V2, levels=c('A','PAR','S1','S0')))
ann_colors = list(g =c(c1 = "#4daf4a", c2 = "#ff7f00", c3 = "#3288BD",c4="#e41a1c",c5 = "#9ECAE1", c6 = "#FCBBA1", c7 = "#084594",c8="#A65628"),
                  group =c(A = "#4daf4a", PAR = "#3288BD", S1 = "#ff7f00",S0="#e41a1c")
                  
)

######
S11_1<-data.frame(row.names=combines0_gonadhigh_1$smangene,Csinfg=log2(combines0_gonadhigh_1$Csinfg+1),Csinmg=log2(combines0_gonadhigh_1$Csinmg+1),Smanbo=log2(combines0_gonadhigh_1$Smanbo+1),Smanbt=log2(combines0_gonadhigh_1$Smanbt+1))
p1<-pheatmap(S11_1,annotation_row = annotation_row,annotation_colors = ann_colors,
            cluster_cols=FALSE,border=FALSE,#show_rownames = F,
            color = colorRampPalette(colors = c("#3288bd","#ffffbf","#EF3B2C"))(100),
            #legend = FALSE,
            fontsize_row=3,
            annotation_legend=F)
gn1=rownames(S11_1)[p1$tree_row[["order"]]]
S11_2<-data.frame(row.names=combines0_gonadhigh_2$smangene,Csinfg=log2(combines0_gonadhigh_2$Csinfg+1),Csinmg=log2(combines0_gonadhigh_2$Csinmg+1),Smanbo=log2(combines0_gonadhigh_2$Smanbo+1),Smanbt=log2(combines0_gonadhigh_2$Smanbt+1))
p2<-pheatmap(S11_2,annotation_row = annotation_row,annotation_colors = ann_colors,
            cluster_cols=FALSE,border=FALSE,#show_rownames = F,
            color = colorRampPalette(colors = c("#3288bd","#ffffbf","#EF3B2C"))(100),
            #legend = FALSE,
            fontsize_row=3,
            annotation_legend=F)
gn2=rownames(S11_2)[p2$tree_row[["order"]]]
S11_3<-data.frame(row.names=combines0_gonadhigh_3$smangene,Csinfg=log2(combines0_gonadhigh_3$Csinfg+1),Csinmg=log2(combines0_gonadhigh_3$Csinmg+1),Smanbo=log2(combines0_gonadhigh_3$Smanbo+1),Smanbt=log2(combines0_gonadhigh_3$Smanbt+1))
p3<-pheatmap(S11_3,annotation_row = annotation_row,annotation_colors = ann_colors,
            cluster_cols=FALSE,border=FALSE,#show_rownames = F,
            color = colorRampPalette(colors = c("#3288bd","#ffffbf","#EF3B2C"))(100),
            #legend = FALSE,
            fontsize_row=3,
            annotation_legend=F)
gn3=rownames(S11_3)[p3$tree_row[["order"]]]
S11_4<-data.frame(row.names=combines0_gonadhigh_4$smangene,Csinfg=log2(combines0_gonadhigh_4$Csinfg+1),Csinmg=log2(combines0_gonadhigh_4$Csinmg+1),Smanbo=log2(combines0_gonadhigh_4$Smanbo+1),Smanbt=log2(combines0_gonadhigh_4$Smanbt+1))
p4<-pheatmap(S11_4,annotation_row = annotation_row,annotation_colors = ann_colors,
            cluster_cols=FALSE,border=FALSE,#show_rownames = F,
            color = colorRampPalette(colors = c("#3288bd","#ffffbf","#EF3B2C"))(100),
            #legend = FALSE,
            fontsize_row=3,
            annotation_legend=F)
gn4=rownames(S11_4)[p4$tree_row[["order"]]]
gn<-paste(gn1,gn2,gn3,gn4)

S11<-data.frame(row.names=combines0_gonadhigh_com$smangene,Csinfg=log2(combines0_gonadhigh_com$Csinfg+1),Csinmg=log2(combines0_gonadhigh_com$Csinmg+1),Smanbo=log2(combines0_gonadhigh_com$Smanbo+1),Smanbt=log2(combines0_gonadhigh_com$Smanbt+1))
S11=rbind(S11_4[gn4,],S11_3[gn3,],S11_2[gn2,],S11_1[gn1,])
p1<-pheatmap(S11,annotation_row = annotation_row,annotation_colors = ann_colors,
             cluster_cols=FALSE,cluster_rows=FALSE,border=FALSE,#show_rownames = F,
             color = colorRampPalette(colors = c("#3288bd","#ffffbf","#EF3B2C"))(100),
             #legend = FALSE,
             fontsize_row=3,
             annotation_legend=F)

########more annotation####
####venn####
combines0_gonadhigh_demas<-subset(combines0_gonadhigh,Csinmg>2*Smanbt)
combines0_gonadhigh_mas<-subset(combines0_gonadhigh,Csinmg*2<Smanbt)
combines0_gonadhigh_defem<-subset(combines0_gonadhigh,Csinfg>2*Smanbo)
combines0_gonadhigh_fem<-subset(combines0_gonadhigh,Csinfg*2<Smanbo)
combines0_gonadhigh_demas$g<-"c5" #9ECAE1
combines0_gonadhigh_defem$g<-"c6" ##FCBBA1
combines0_gonadhigh_mas$g<-"c7" ##084594
combines0_gonadhigh_fem$g<-"c8" #A65628
combines0_gonadhigh_formassig<-rbind(combines0_gonadhigh_demas,combines0_gonadhigh_mas,combines0_gonadhigh_defem,combines0_gonadhigh_fem)
#write.table(combines0_gonadhigh_formassig,file="smangene.csingonadhigh.expchange.formassig.txt",sep = '\t',row.names = F,quote = F) ##calcu gonad high mas sig



library(VennDiagram)       
p<-venn.diagram(list("demasculinization"= combines0_gonadhigh_demas$smangene,
                     "masculinization" = combines0_gonadhigh_mas$smangene,
                     "defeminization" = combines0_gonadhigh_defem$smangene,
                     "feminization" = combines0_gonadhigh_fem$smangene),
                fill=c('#9ECAE1','#3288BD','#FCBBA1','#EF3B2C'), 
                resolution = 500, 
                alpha=c(0.8,0.8,0.8,0.8),
                cat.fontface =3, 
                cat.pos = c(-5, 5, -110, 110),
                filename = NULL)
p
pdf(file="combines0_gonadhigh_VennDiagram.pdf")
grid.draw(p)
dev.off()
##########
combines0_gonadhigh_demas = combines0_gonadhigh_demas[!rownames(combines0_gonadhigh_demas) %in% rownames(combines0_gonadhigh_1),]  
combines0_gonadhigh_demas = combines0_gonadhigh_demas[!rownames(combines0_gonadhigh_demas) %in% rownames(combines0_gonadhigh_2),]  
combines0_gonadhigh_defem = combines0_gonadhigh_defem[!rownames(combines0_gonadhigh_defem) %in% rownames(combines0_gonadhigh_1),]  
combines0_gonadhigh_defem = combines0_gonadhigh_defem[!rownames(combines0_gonadhigh_defem) %in% rownames(combines0_gonadhigh_3),] 
combines0_gonadhigh_mas = combines0_gonadhigh_mas[!rownames(combines0_gonadhigh_mas) %in% rownames(combines0_gonadhigh_4),]  
combines0_gonadhigh_mas = combines0_gonadhigh_mas[!rownames(combines0_gonadhigh_mas) %in% rownames(combines0_gonadhigh_3),]  
combines0_gonadhigh_fem = combines0_gonadhigh_fem[!rownames(combines0_gonadhigh_fem) %in% rownames(combines0_gonadhigh_4),]  
combines0_gonadhigh_fem = combines0_gonadhigh_fem[!rownames(combines0_gonadhigh_fem) %in% rownames(combines0_gonadhigh_2),]
combines0_gonadhigh_demas$g<-"c5" #9ECAE1
combines0_gonadhigh_defem$g<-"c6" ##FCBBA1
combines0_gonadhigh_mas$g<-"c7" ##084594
combines0_gonadhigh_fem$g<-"c8" #A65628
combines0_gonadhigh_com1<-rbind(combines0_gonadhigh_4,combines0_gonadhigh_3,combines0_gonadhigh_2,combines0_gonadhigh_1,combines0_gonadhigh_demas,combines0_gonadhigh_defem,combines0_gonadhigh_mas,combines0_gonadhigh_fem)
head(combines0_gonadhigh_com1)

S11_5<-data.frame(row.names=combines0_gonadhigh_demas$smangene,Csinfg=log2(combines0_gonadhigh_demas$Csinfg+1),Csinmg=log2(combines0_gonadhigh_demas$Csinmg+1),Smanbo=log2(combines0_gonadhigh_demas$Smanbo+1),Smanbt=log2(combines0_gonadhigh_demas$Smanbt+1))
p5<-pheatmap(S11_5,annotation_row = annotation_row,annotation_colors = ann_colors,
             cluster_cols=FALSE,border=FALSE,#show_rownames = F,
             color = colorRampPalette(colors = c("#3288bd","#ffffbf","#EF3B2C"))(100),
             #legend = FALSE,
             fontsize_row=3,
             annotation_legend=F)
gn5=rownames(S11_5)[p5$tree_row[["order"]]]
S11_6<-data.frame(row.names=combines0_gonadhigh_defem$smangene,Csinfg=log2(combines0_gonadhigh_defem$Csinfg+1),Csinmg=log2(combines0_gonadhigh_defem$Csinmg+1),Smanbo=log2(combines0_gonadhigh_defem$Smanbo+1),Smanbt=log2(combines0_gonadhigh_defem$Smanbt+1))
p6<-pheatmap(S11_6,annotation_row = annotation_row,annotation_colors = ann_colors,
             cluster_cols=FALSE,border=FALSE,#show_rownames = F,
             color = colorRampPalette(colors = c("#3288bd","#ffffbf","#EF3B2C"))(100),
             #legend = FALSE,
             fontsize_row=3,
             annotation_legend=F)
gn6=rownames(S11_6)[p6$tree_row[["order"]]]
S11_7<-data.frame(row.names=combines0_gonadhigh_mas$smangene,Csinfg=log2(combines0_gonadhigh_mas$Csinfg+1),Csinmg=log2(combines0_gonadhigh_mas$Csinmg+1),Smanbo=log2(combines0_gonadhigh_mas$Smanbo+1),Smanbt=log2(combines0_gonadhigh_mas$Smanbt+1))
p7<-pheatmap(S11_7,annotation_row = annotation_row,annotation_colors = ann_colors,
             cluster_cols=FALSE,border=FALSE,#show_rownames = F,
             color = colorRampPalette(colors = c("#3288bd","#ffffbf","#EF3B2C"))(100),
             #legend = FALSE,
             fontsize_row=3,
             annotation_legend=F)
gn7=rownames(S11_7)[p7$tree_row[["order"]]]
S11_8<-data.frame(row.names=combines0_gonadhigh_fem$smangene,Csinfg=log2(combines0_gonadhigh_fem$Csinfg+1),Csinmg=log2(combines0_gonadhigh_fem$Csinmg+1),Smanbo=log2(combines0_gonadhigh_fem$Smanbo+1),Smanbt=log2(combines0_gonadhigh_fem$Smanbt+1))
p8<-pheatmap(S11_8,annotation_row = annotation_row,annotation_colors = ann_colors,
             cluster_cols=FALSE,border=FALSE,#show_rownames = F,
             color = colorRampPalette(colors = c("#3288bd","#ffffbf","#EF3B2C"))(100),
             #legend = FALSE,
             fontsize_row=3,
             annotation_legend=F)
gn8=rownames(S11_8)[p8$tree_row[["order"]]]

annotation_row <- data.frame(row.names=combines0_gonadhigh_com1$smangene,group=factor(combines0_gonadhigh_com1$V2, levels=c('A','PAR','S1','S0')),g=factor(combines0_gonadhigh_com1$g, levels=c('c1','c2','c3','c4','c5','c6','c7','c8')))
ann_colors = list(group =c(A = "#4daf4a", PAR = "#3288BD", S1 = "#ff7f00",S0="#e41a1c"),
                  g =c(c1 = "#c1c1c1", c2 = "#ff7f00", c3 = "#3288BD",c4="#e41a1c",c5 = "#9ECAE1", c6 = "#FCBBA1", c7 = "#084594",c8="#A65628")
)                  
S11<-data.frame(row.names=combines0_gonadhigh_com1$smangene,Csinfg=log2(combines0_gonadhigh_com1$Csinfg+1),Csinmg=log2(combines0_gonadhigh_com1$Csinmg+1),Smanbo=log2(combines0_gonadhigh_com1$Smanbo+1),Smanbt=log2(combines0_gonadhigh_com1$Smanbt+1))
S11=rbind(S11_4[gn4,],S11_7[gn7,],S11_8[gn8,],S11_3[gn3,],S11_2[gn2,],S11_5[gn5,],S11_6[gn6,],S11_1[gn1,])
S11=rbind(S11_7[gn7,],S11_4[gn4,],S11_3[gn3,],S11_5[gn5,],S11_8[gn8,],S11_2[gn2,],S11_6[gn6,],S11_1[gn1,])
p1<-pheatmap(S11,annotation_row = annotation_row,annotation_colors = ann_colors,
             cluster_cols=FALSE,cluster_rows=FALSE,border=FALSE,#show_rownames = F,
             color = colorRampPalette(colors = c("#3288bd","#ffffbf","#EF3B2C"))(100),
             #legend = FALSE,
             fontsize_row=3,
             annotation_legend=F)
#write.table(S11,file="phatmap_8.11_order",sep = '\t',row.names = T,quote = F)

######add RNAi annatation genes######
setwd("/Users/yifeng/Library/Mobile Documents/com~apple~CloudDocs/project/81worm/RNA/csin-sman/density/01phenotype/")
head(combines0_gonadhigh_com1)
row.names(combines0_gonadhigh_com1)<-combines0_gonadhigh_com1$smangene
celesman<-read.table("phatmap_8.11_order.rnai.reproduction.smancele.selected.txt", header = FALSE,row.names = 3)
head(celesman)
combines0_gonadhigh_com1cele<-merge(combines0_gonadhigh_com1,celesman,by= "row.names",all.x = TRUE,sort=FALSE)
head(combines0_gonadhigh_com1cele)

annotation_row <- data.frame(row.names=combines0_gonadhigh_com1cele$smangene,group=factor(combines0_gonadhigh_com1cele$V2.x, levels=c('A','PAR','S1','S0')),g=factor(combines0_gonadhigh_com1cele$g, levels=c('c1','c2','c3','c4','c5','c6','c7','c8')),cele=combines0_gonadhigh_com1cele$V2.y)
ann_colors = list(group =c(A = "#4daf4a", PAR = "#3288BD", S1 = "#ff7f00",S0="#e41a1c"),
                  g =c(c1 = "#c1c1c1", c2 = "#ff7f00", c3 = "#3288BD",c4="#e41a1c",c5 = "#9ECAE1", c6 = "#FCBBA1", c7 = "#084594",c8="#A65628"),
                  cele=c(Fem = "#e41a1c", B = "#000000", M = "#3288BD")
)      

p1<-pheatmap(S11,annotation_row = annotation_row,annotation_colors = ann_colors,
             cluster_cols=FALSE,cluster_rows=FALSE,border=FALSE,#show_rownames = F,
             color = colorRampPalette(colors = c("#3288bd","#ffffbf","#EF3B2C"))(100),
             #legend = FALSE,
             fontsize_row=3,
             annotation_legend=F)




###
head(combines0_gonadlow)
combines0_gonadlow_1<-subset(combines0_gonadlow,Csinmg>2*Smanbt & Csinfg>2*Smanbo)
combines0_gonadlow_1<-combines0_gonadlow_1[order(-combines0_gonadlow_1[,4]),]
combines0_gonadlow_1$g<-"c1" #green
combines0_gonadlow_2<-subset(combines0_gonadlow,Csinmg>2*Smanbt & Csinfg*2<Smanbo) 
combines0_gonadlow_2<-combines0_gonadlow_2[order(-combines0_gonadlow_2[,2]),]
combines0_gonadlow_2$g<-"c2" ##orange
combines0_gonadlow_3<-subset(combines0_gonadlow,Csinmg*2<Smanbt & Csinfg>2*Smanbo) 
combines0_gonadlow_3<-combines0_gonadlow_3[order(-combines0_gonadlow_3[,3]),]
combines0_gonadlow_3$g<-"c3" ##blue
combines0_gonadlow_4<-subset(combines0_gonadlow,Csinmg*2<Smanbt & Csinfg*2<Smanbo) 
combines0_gonadlow_4<-combines0_gonadlow_4[order(-combines0_gonadlow_4[,3]),]
combines0_gonadlow_4$g<-"c4" #red
combines0_gonadlow_com<-rbind(combines0_gonadlow_4,combines0_gonadlow_3,combines0_gonadlow_2,combines0_gonadlow_1)
head(combines0_gonadlow_com)
write.table(combines0_gonadlow_com,file="smangene.csingonadlow.expchange.txt",sep = '\t',row.names = F,quote = F)

annotation_row <- data.frame(row.names=combines0_gonadlow_com$smangene,g=factor(combines0_gonadlow_com$g, levels=c('c1','c2','c3','c4')),group=factor(combines0_gonadlow_com$V2, levels=c('A','PAR','S1','S0')))
ann_colors = list(g =c(c1 = "#4daf4a", c2 = "#ff7f00", c3 = "#3288BD",c4="#e41a1c"),
                  group =c(A = "#4daf4a", PAR = "#3288BD", S1 = "#ff7f00",S0="#e41a1c")
)


S11_1<-data.frame(row.names=combines0_gonadlow_1$smangene,Csinfg=log2(combines0_gonadlow_1$Csinfg+1),Csinmg=log2(combines0_gonadlow_1$Csinmg+1),Smanbo=log2(combines0_gonadlow_1$Smanbo+1),Smanbt=log2(combines0_gonadlow_1$Smanbt+1))
p1<-pheatmap(S11_1,annotation_row = annotation_row,annotation_colors = ann_colors,
             cluster_cols=FALSE,border=FALSE,#show_rownames = F,
             color = colorRampPalette(colors = c("#3288bd","#ffffbf","#EF3B2C"))(100),
             #legend = FALSE,
             fontsize_row=3,
             annotation_legend=F)
gn1=rownames(S11_1)[p1$tree_row[["order"]]]
S11_2<-data.frame(row.names=combines0_gonadlow_2$smangene,Csinfg=log2(combines0_gonadlow_2$Csinfg+1),Csinmg=log2(combines0_gonadlow_2$Csinmg+1),Smanbo=log2(combines0_gonadlow_2$Smanbo+1),Smanbt=log2(combines0_gonadlow_2$Smanbt+1))
p2<-pheatmap(S11_2,annotation_row = annotation_row,annotation_colors = ann_colors,
             cluster_cols=FALSE,border=FALSE,#show_rownames = F,
             color = colorRampPalette(colors = c("#3288bd","#ffffbf","#EF3B2C"))(100),
             #legend = FALSE,
             fontsize_row=3,
             annotation_legend=F)
gn2=rownames(S11_2)[p2$tree_row[["order"]]]
S11_3<-data.frame(row.names=combines0_gonadlow_3$smangene,Csinfg=log2(combines0_gonadlow_3$Csinfg+1),Csinmg=log2(combines0_gonadlow_3$Csinmg+1),Smanbo=log2(combines0_gonadlow_3$Smanbo+1),Smanbt=log2(combines0_gonadlow_3$Smanbt+1))
p3<-pheatmap(S11_3,annotation_row = annotation_row,annotation_colors = ann_colors,
             cluster_cols=FALSE,border=FALSE,#show_rownames = F,
             color = colorRampPalette(colors = c("#3288bd","#ffffbf","#EF3B2C"))(100),
             #legend = FALSE,
             fontsize_row=3,
             annotation_legend=F)
gn3=rownames(S11_3)[p3$tree_row[["order"]]]
S11_4<-data.frame(row.names=combines0_gonadlow_4$smangene,Csinfg=log2(combines0_gonadlow_4$Csinfg+1),Csinmg=log2(combines0_gonadlow_4$Csinmg+1),Smanbo=log2(combines0_gonadlow_4$Smanbo+1),Smanbt=log2(combines0_gonadlow_4$Smanbt+1))
p4<-pheatmap(S11_4,annotation_row = annotation_row,annotation_colors = ann_colors,
             cluster_cols=FALSE,border=FALSE,#show_rownames = F,
             color = colorRampPalette(colors = c("#3288bd","#ffffbf","#EF3B2C"))(100),
             #legend = FALSE,
             fontsize_row=3,
             annotation_legend=F)
gn4=rownames(S11_4)[p4$tree_row[["order"]]]
gn<-paste(gn1,gn2,gn3,gn4)

S11<-data.frame(row.names=combines0_gonadlow_com$smangene,Csinfg=log2(combines0_gonadlow_com$Csinfg+1),Csinmg=log2(combines0_gonadlow_com$Csinmg+1),Smanbo=log2(combines0_gonadlow_com$Smanbo+1),Smanbt=log2(combines0_gonadlow_com$Smanbt+1))
S11=rbind(S11_4[gn4,],S11_3[gn3,],S11_2[gn2,],S11_1[gn1,])
p1<-pheatmap(S11,annotation_row = annotation_row,annotation_colors = ann_colors,
             cluster_cols=FALSE,cluster_rows=FALSE,border=FALSE,#show_rownames = F,
             color = colorRampPalette(colors = c("#3288bd","#ffffbf","#EF3B2C"))(100),
             #legend = FALSE,
             fontsize_row=3,
             annotation_legend=F)



###
head(combines0_gonad_soma)
combines0_gonad_soma_1<-subset(combines0_gonad_soma,Csinmg>2*Smanbt & Csinfg>2*Smanbo)
combines0_gonad_soma_1<-combines0_gonad_soma_1[order(-combines0_gonad_soma_1[,4]),]
combines0_gonad_soma_1$g<-"c1" #green
combines0_gonad_soma_2<-subset(combines0_gonad_soma,Csinmg>2*Smanbt & Csinfg*2<Smanbo) 
combines0_gonad_soma_2<-combines0_gonad_soma_2[order(-combines0_gonad_soma_2[,2]),]
combines0_gonad_soma_2$g<-"c2" ##orange
combines0_gonad_soma_3<-subset(combines0_gonad_soma,Csinmg*2<Smanbt & Csinfg>2*Smanbo) 
combines0_gonad_soma_3<-combines0_gonad_soma_3[order(-combines0_gonad_soma_3[,3]),]
combines0_gonad_soma_3$g<-"c3" ##blue
combines0_gonad_soma_4<-subset(combines0_gonad_soma,Csinmg*2<Smanbt & Csinfg*2<Smanbo) 
combines0_gonad_soma_4<-combines0_gonad_soma_4[order(-combines0_gonad_soma_4[,3]),]
combines0_gonad_soma_4$g<-"c4" #red
combines0_gonad_soma_com<-rbind(combines0_gonad_soma_4,combines0_gonad_soma_3,combines0_gonad_soma_2,combines0_gonad_soma_1)
head(combines0_gonad_soma_com)
write.table(combines0_gonad_soma_com,file="smangene.csingonadsoma.expchange.txt",sep = '\t',row.names = F,quote = F)

annotation_row <- data.frame(row.names=combines0_gonad_soma_com$smangene,g=factor(combines0_gonad_soma_com$g, levels=c('c1','c2','c3','c4')),group=factor(combines0_gonad_soma_com$V2, levels=c('A','PAR','S1','S0')))
ann_colors = list(g =c(c1 = "#4daf4a", c2 = "#ff7f00", c3 = "#3288BD",c4="#e41a1c"),
                  group =c(A = "#4daf4a", PAR = "#3288BD", S1 = "#ff7f00",S0="#e41a1c")
)


S11_1<-data.frame(row.names=combines0_gonad_soma_1$smangene,Csinfg=log2(combines0_gonad_soma_1$Csinfg+1),Csinmg=log2(combines0_gonad_soma_1$Csinmg+1),Smanbo=log2(combines0_gonad_soma_1$Smanbo+1),Smanbt=log2(combines0_gonad_soma_1$Smanbt+1))
p1<-pheatmap(S11_1,annotation_row = annotation_row,annotation_colors = ann_colors,
             cluster_cols=FALSE,border=FALSE,#show_rownames = F,
             color = colorRampPalette(colors = c("#3288bd","#ffffbf","#EF3B2C"))(100),
             #legend = FALSE,
             fontsize_row=3,
             annotation_legend=F)
gn1=rownames(S11_1)[p1$tree_row[["order"]]]
S11_2<-data.frame(row.names=combines0_gonad_soma_2$smangene,Csinfg=log2(combines0_gonad_soma_2$Csinfg+1),Csinmg=log2(combines0_gonad_soma_2$Csinmg+1),Smanbo=log2(combines0_gonad_soma_2$Smanbo+1),Smanbt=log2(combines0_gonad_soma_2$Smanbt+1))
p2<-pheatmap(S11_2,annotation_row = annotation_row,annotation_colors = ann_colors,
             cluster_cols=FALSE,border=FALSE,#show_rownames = F,
             color = colorRampPalette(colors = c("#3288bd","#ffffbf","#EF3B2C"))(100),
             #legend = FALSE,
             fontsize_row=3,
             annotation_legend=F)
gn2=rownames(S11_2)[p2$tree_row[["order"]]]
S11_3<-data.frame(row.names=combines0_gonad_soma_3$smangene,Csinfg=log2(combines0_gonad_soma_3$Csinfg+1),Csinmg=log2(combines0_gonad_soma_3$Csinmg+1),Smanbo=log2(combines0_gonad_soma_3$Smanbo+1),Smanbt=log2(combines0_gonad_soma_3$Smanbt+1))
p3<-pheatmap(S11_3,annotation_row = annotation_row,annotation_colors = ann_colors,
             cluster_cols=FALSE,border=FALSE,#show_rownames = F,
             color = colorRampPalette(colors = c("#3288bd","#ffffbf","#EF3B2C"))(100),
             #legend = FALSE,
             fontsize_row=3,
             annotation_legend=F)
gn3=rownames(S11_3)[p3$tree_row[["order"]]]
S11_4<-data.frame(row.names=combines0_gonad_soma_4$smangene,Csinfg=log2(combines0_gonad_soma_4$Csinfg+1),Csinmg=log2(combines0_gonad_soma_4$Csinmg+1),Smanbo=log2(combines0_gonad_soma_4$Smanbo+1),Smanbt=log2(combines0_gonad_soma_4$Smanbt+1))
p4<-pheatmap(S11_4,annotation_row = annotation_row,annotation_colors = ann_colors,
             cluster_cols=FALSE,border=FALSE,#show_rownames = F,
             color = colorRampPalette(colors = c("#3288bd","#ffffbf","#EF3B2C"))(100),
             #legend = FALSE,
             fontsize_row=3,
             annotation_legend=F)
gn4=rownames(S11_4)[p4$tree_row[["order"]]]

S11<-data.frame(row.names=combines0_gonad_soma_com$smangene,Csinfg=log2(combines0_gonad_soma_com$Csinfg+1),Csinmg=log2(combines0_gonad_soma_com$Csinmg+1),Smanbo=log2(combines0_gonad_soma_com$Smanbo+1),Smanbt=log2(combines0_gonad_soma_com$Smanbt+1))
S11=rbind(S11_4[gn4,],S11_3[gn3,],S11_2[gn2,],S11_1[gn1,])
p1<-pheatmap(S11,annotation_row = annotation_row,annotation_colors = ann_colors,
             cluster_cols=FALSE,cluster_rows=FALSE,border=FALSE,#show_rownames = F,
             color = colorRampPalette(colors = c("#3288bd","#ffffbf","#EF3B2C"))(100),
             #legend = FALSE,
             fontsize_row=3,
             annotation_legend=F)



# selective, add cele sd gene --------------------------------------------------------------

##########no############
####with cele sd genes, not involved in tmur high gonad exp genes
csd<-read.table("cele.sexde.gene", header = FALSE)
head(csd)
rownames(csd)<-csd$V1
combine_csd = combine1[rownames(combine1) %in% rownames(csd),] 
rownames(CsinAtpm)<-CsinAtpm$Row.names
CsinAtpm_csd = CsinAtpm[rownames(CsinAtpm) %in% rownames(csd),] 
CsinAtpm_csd
#"Tm1G004760" "Tm1G003109" "Tm1G003937" "Tm1G005504" "Tm1G000698" "Tm1G001327" "Tm1G001612" "Tm1G000937" "Tm4G009474" "Tm2G007350"
########

