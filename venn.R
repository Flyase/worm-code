#######!!!upset only count for uniq overlap, that overlap genes between 2 of 4 species may lower that overlap genes among the 4 species, if genes emerged in other groups, then they will minus it 
setwd("/Users/yifeng/Library/Mobile Documents/com~apple~CloudDocs/project/81worm/ancestralX/1nigon/")
#install.packages("nVennR")
folder <- "/Users/yifeng/Library/Mobile Documents/com~apple~CloudDocs/project/81worm/ancestralX/1nigon/"      
file_list <- list.files(path=folder, pattern="*.gene")                              
# read in each .txt file in file_list and rbind them into a data frame called data 
combine <- 
  do.call("rbind", 
          lapply(list.files(pattern="*.gene"),
                 function(x) 
                   read.table(paste(folder, x, sep=''), 
                              header = FALSE, 
                              stringsAsFactors = FALSE)))
head(combine)

########venn
library(VennDiagram)
p1 <- venn.diagram(list("Asim_NA"= combine[combine$V3=="Asim" & combine$V2=="NAA",]$V1,
                        "Tcan_NA"= combine[combine$V3=="Tcan" & combine$V2=="NAA",]$V1,
                        "Alum_NA"= combine[combine$V3=="Alum" & combine$V2=="NAA",]$V1,
                        "Asuum_NA"= combine[combine$V3=="Asuum" & combine$V2=="NAA",]$V1
                        ),
                   fill=c("#093B5B","#D54332","#9CBB4C","#C7E9C0"),
                   resolution = 500, 
                   alpha=c(0.36,0.6,0.6,0.6),
                   cat.fontface =3, 
                   cat.cex = 2,# tag size
                   cex = 2,# num size
                   fontface = "bold", # num
                   filename = NULL)
pdf(file="NA.nvenn.pdf")
grid.draw(p1)
dev.off()


p2 <- venn.diagram(list("Asim_NB"= combine[combine$V3=="Asim" & combine$V2=="NB",]$V1,
                        "Tcan_NB"= combine[combine$V3=="Tcan" & combine$V2=="NB",]$V1,
                        "Alum_NB"= combine[combine$V3=="Alum" & combine$V2=="NB",]$V1,
                        "Asuum_NB"= combine[combine$V3=="Asuum" & combine$V2=="NB",]$V1
),
fill=c("#093B5B","#D54332","#9CBB4C","#C7E9C0"),
resolution = 500, 
alpha=c(0.36,0.6,0.6,0.6),
cat.fontface =3, 
cat.cex = 2,# tag size
cex = 2,# num size
fontface = "bold", # num
filename = NULL)
pdf(file="NB.nvenn.pdf")
grid.draw(p2)
dev.off()

p3 <- venn.diagram(list("Asim_ND"= combine[combine$V3=="Asim" & combine$V2=="ND",]$V1,
                        "Tcan_ND"= combine[combine$V3=="Tcan" & combine$V2=="ND",]$V1,
                        "Alum_ND"= combine[combine$V3=="Alum" & combine$V2=="ND",]$V1,
                        "Asuum_ND"= combine[combine$V3=="Asuum" & combine$V2=="ND",]$V1
),
fill=c("#093B5B","#D54332","#9CBB4C","#C7E9C0"),
resolution = 500, 
alpha=c(0.36,0.6,0.6,0.6),
cat.fontface =3, 
cat.cex = 2,# tag size
cex = 2,# num size
fontface = "bold", # num
filename = NULL)
pdf(file="ND.Alumclade.nvenn.pdf")
grid.draw(p3)
dev.off()


p4 <- venn.diagram(list("Alum_ND"= combine[combine$V3=="Alum" & combine$V2=="ND",]$V1,
                        "Dimm_ND"= combine[combine$V3=="Dimm" & combine$V2=="ND",]$V1,
                        "Ovol_ND"= combine[combine$V3=="Ovol" & combine$V2=="ND",]$V1,
                        "Bpah_ND"= combine[combine$V3=="Bpah" & combine$V2=="ND",]$V1,
                        "Bmal_ND"= combine[combine$V3=="Bmal" & combine$V2=="ND",]$V1
),
fill=c("#093B5B","#D54332","#9CBB4C","#C7E9C0","#FDB863"),
resolution = 500, 
alpha=c(0.36,0.6,0.6,0.6,0.6),
cat.fontface =3, 
cat.cex = 2,# tag size
cex = 2,# num size
#fontface = "bold", # num
filename = NULL)
pdf(file="ND.cladeiii.nvenn.pdf")
grid.draw(p4)
dev.off()


p5 <- venn.diagram(list("Ovol_NE"= combine[combine$V3=="Ovol" & combine$V2=="NE",]$V1,
                        "Ooch_NE"= combine[combine$V3=="Ooch" & combine$V2=="NE",]$V1,
                        "Ofle_NE"= combine[combine$V3=="Ofle" & combine$V2=="NE",]$V1
),
fill=c("#093B5B","#D54332","#9CBB4C"),
resolution = 500, 
alpha=c(0.36,0.6,0.6),
cat.fontface =3, 
cat.cex = 2,# tag size
cex = 2,# num size
fontface = "bold", # num
filename = NULL)
pdf(file="NE.nvenn.pdf")
grid.draw(p5)
dev.off()

#NN
p5 <- venn.diagram(list("Wban_NN"= combine[combine$V3=="Wban" & combine$V2=="NN",]$V1,
                        "Bpah_NN"= combine[combine$V3=="Bpah" & combine$V2=="NN",]$V1,
                        "Bmal_NN"= combine[combine$V3=="Bmal" & combine$V2=="NN",]$V1
),
fill=c("#093B5B","#D54332","#9CBB4C"),
resolution = 500, 
alpha=c(0.36,0.6,0.6),
cat.fontface =3, 
cat.cex = 2,# tag size
cex = 2,# num size
fontface = "bold", # num
filename = NULL)
pdf(file="NN.cladeiii.nvenn.pdf")
grid.draw(p5)
dev.off()



#######I
p5 <- venn.diagram(list("Rcul_NA"= combine[combine$V3=="Rcul" & combine$V2=="NAA",]$V1,
                        "Tmur_NA"= combine[combine$V3=="Tmur" & combine$V2=="NAA",]$V1
),
fill=c("#093B5B","#9CBB4C"),
resolution = 500, 
alpha=c(0.36,0.6),
cat.fontface =3, 
cat.cex = 2,# tag size
cex = 2,# num size
fontface = "bold", # num
filename = NULL)
pdf(file="NA.I.nvenn.pdf")
grid.draw(p5)
dev.off()

#######IV
p5 <- venn.diagram(list("Pred_ND"= combine[combine$V3=="Pred" & combine$V2=="ND",]$V1,
                        "Srat_ND"= combine[combine$V3=="Srat" & combine$V2=="ND",]$V1
),
fill=c("#093B5B","#9CBB4C"),
resolution = 500, 
alpha=c(0.36,0.6),
cat.fontface =3, 
cat.cex = 2,# tag size
cex = 2,# num size
fontface = "bold", # num
filename = NULL)
pdf(file="ND.IV.nvenn.pdf")
grid.draw(p5)
dev.off()

######V
p4 <- venn.diagram(list("Dviv_NN"= combine[combine$V3=="Dviv" & combine$V2=="NN",]$V1,
                        "Name_NN"= combine[combine$V3=="Name" & combine$V2=="NN",]$V1,
                        "Svul_NN"= combine[combine$V3=="Svul" & combine$V2=="NN",]$V1,
                        "Cele_NN"= combine[combine$V3=="Cele" & combine$V2=="NN",]$V1
),
fill=c("#093B5B","#D54332","#9CBB4C","#C7E9C0"),
resolution = 500, 
alpha=c(0.36,0.6,0.6,0.6),
cat.fontface =3, 
cat.cex = 2,# tag size
cex = 2,# num size
#fontface = "bold", # num
filename = NULL)
pdf(file="NN.cladev.nvenn.pdf")
grid.draw(p4)
dev.off()








library(nVennR)
p <- plotVenn(list("Adult_muscle"=samples[[1]]$id,
                   "1-cell"= samples[[2]]$id,
                   "32-cell" = samples[[3]]$id,
                   "64-cell" = samples[[4]]$id,
                   "Gastrula" = samples[[5]]$id,
                   "Larvae" = samples[[6]]$id),
              setColors = c('#e41a1c', '#377eb8', '#4daf4a','#984ea3','#ff7f00', '#ffff33'),
              labelRegions=F, opacity=0.9, fontScale=1, 
              outFile = "/Users/yifeng/Library/Mobile Documents/com~apple~CloudDocs/project/81worm/ancestralX/1nigon/NB.svg")




p3 <- venn.diagram(list("Asim_ND"= combine[combine$V3=="Asim" & combine$V2=="ND",]$V1,
                        "Tcan_ND"= combine[combine$V3=="Tcan" & combine$V2=="ND",]$V1,
                        "Alum_ND"= combine[combine$V3=="Alum" & combine$V2=="ND",]$V1,
                        "Asuum_ND"= combine[combine$V3=="Asuum" & combine$V2=="ND",]$V1
),
fill=c("#093B5B","#D54332","#9CBB4C","#C7E9C0"),
resolution = 500, 
alpha=c(0.36,0.6,0.6,0.6),
cat.fontface =3, 
cat.cex = 2,# tag size
cex = 2,# num size
fontface = "bold", # num
filename = NULL)

pdf(file="ND.pdf")
grid.draw(p3)
dev.off()


BmalXnx <- read.table("BmalX.NX.gene", header = FALSE,sep="",stringsAsFactors = FALSE)
HconXnx <- read.table("HconX.NX.gene",header = FALSE,sep="",stringsAsFactors = FALSE)
PpacXnx <- read.table("PpacX.NX.gene",header = FALSE,sep="",stringsAsFactors = FALSE)
OvolXnx <- read.table("OvolX.NX.gene",header = FALSE,sep="",stringsAsFactors = FALSE)
SratXnx <- read.table("SratX.NX.gene",header = FALSE,sep="",stringsAsFactors = FALSE)
p1 <- venn.diagram(list("BmalX_CeleNX"= BmalXnx$V1,
                        "OvolX_CeleNX"= OvolXnx$V1,
                        "SratX_CeleNX" = SratXnx$V1,
                        "HconX_CeleNX" = HconXnx$V1,
                        "PpacX_CeleNX" = PpacXnx$V1),
                   #fill=c("#b3e2cd","#fdcdac","#cbd5e8","#f4cae4","#e6f5c9"),
                   fill=c('#1F4E79','#FB8072','#BCBDDC','#C7E9C0','#FCBBA1'),
                   resolution = 500, 
                   alpha=c(0.6,0.6,0.6,0.6,0.6),
                   cat.fontface =3, 
                   cat.pos = c(-5, 5, -110, 110, -10 ),
                   cat.cex = 2,# tag size
                   cex = 2,# num size
                   fontface = "bold", # num
                   filename = NULL)

pdf(file="allNX.pdf")
grid.draw(p1)
dev.off()

BmalXnx <- read.table("BmalX.NX.gene", header = FALSE,sep="",stringsAsFactors = FALSE)
HconXnx <- read.table("HconX.NX.gene",header = FALSE,sep="",stringsAsFactors = FALSE)
PpacXnx <- read.table("PpacX.NX.gene",header = FALSE,sep="",stringsAsFactors = FALSE)
OvolXnx <- read.table("OvolX.NX.gene",header = FALSE,sep="",stringsAsFactors = FALSE)
SratXnx <- read.table("SratX.NX.gene",header = FALSE,sep="",stringsAsFactors = FALSE)
p1 <- venn.diagram(list("BmalX_CeleNX"= BmalXnx$V1,
                        "OvolX_CeleNX"= OvolXnx$V1,
                        "SratX_CeleNX" = SratXnx$V1,
                        "HconX_CeleNX" = HconXnx$V1),
                   #fill=c("#b3e2cd","#fdcdac","#cbd5e8","#f4cae4","#e6f5c9"),
                   fill=c('#1F4E79','#FB8072','#BCBDDC','#C7E9C0'),
                   resolution = 500, 
                   alpha=c(0.6,0.6,0.6,0.6),
                   cat.fontface =3, 
                   #cat.pos = c(-5, 5, -110, 110, -10 ),
                   cat.cex = 2,# tag size
                   cex = 2,# num size
                   fontface = "bold", # num
                   filename = NULL)

pdf(file="allNX_noPpac.pdf")
grid.draw(p1)
dev.off()
