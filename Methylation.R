###########################
# Author: Vi Dam
# Version: 2019 01 15
###########################
#In Unix,
#conda install -c r r-data.table
#In R:
#install.packages("data.table")
#require(data.table)

library(data.table)
library("ggplot2")
library(reshape)
library(ggplot2)
library(plyr)
library(gtools)
require("ggrepel")

#-----------------Download data from HPV list USING BASH
gdc-client download -m manifest.txt
round(memory.limit()/2^20, 2)

#-------------------merge TCGA Methylation tables
#library(data.table)
path <- "/Users/ViDam/Dropbox/HNSC data/"
setwd(path)
#-----------MSET table: beta value of all patients-----------------------------
#install.packages("reshape")
#library(reshape)

path <- "/home/vdam/methylation/TCGAFiles"
filesList <- list.files(path)
filesList

for (file in filesList){
  #extract substring TCGAid from fileName
  TCGAid <- sub(".*lvl-3.","",file)
  TCGAid <- sub("-01A-.*","",TCGAid)

  # Extract Composite Element REF and Beta_value columns
  # Replace Beta_value column name to TCGAid
  if (!exists("TCGAdataset")){
    TCGAdataset <- fread((file.path(path,file)),select = c("Composite Element REF","Beta_value"))
    colnames(TCGAdataset)[colnames(TCGAdataset)=="Beta_value"] <- TCGAid
  }
  else {
    temp_dataset <-fread((file.path(path,file)),select = c("Composite Element REF","Beta_value"))
    TCGAdataset<- merge(TCGAdataset, temp_dataset,by="Composite Element REF")
    colnames(TCGAdataset)[colnames(TCGAdataset)=="Beta_value"] <- TCGAid
    rm(temp_dataset)
  }
}

colnames(TCGAdataset) <- substr(colnames(TCGAdataset), 1,12)
TCGAdataset <- TCGAdataset[,-1]

TCGA_HSNC <- read.table("/home/vdam/methylation/TCGAdataset.csv", header=TRUE, sep="\t", check.names=F,stringsAsFactors=T)

HSNC <- cbind(TCGA_HSNC, TCGAdataset)
write.table(HSNC, file = "TCGAdataset.csv", sep = "\t", row.names = F, quote=FALSE)

dft <- data.frame(colnames(TCGAdataset))
dft$Condition <- "Normal"
colnames(dft) <- c("TCGA","Condition")

HPV <- read.table("/home/vdam/methylation/HNSC_annotations_3_groups.csv", header=TRUE, sep=",", stringsAsFactors=TRUE)
HPV <- rbind(HPV,dft)
colnames(HPV) <- c("TCGA","Condition")
write.table(HPV, file = "HNSC_annotations_4_groups.tsv", sep = "\t", row.names = F, quote=FALSE)


HSNC1 <- melt(HSNC1, id=("Composite Element REF"))
colnames(HSNC1) <- c("Composite Element REF","TCGA","Beta_value")

HSNC2 <- merge(HSNC1, HPV, by.x="TCGA",by.y="TCGA")

#---------GGPLOT-----------------------------
#install.packages("tidyverse")
#library(ggplot2)

#Beta_value Density Plot
plot <- ggplot(HSNC2, aes(x=Beta_value, color=Condition)) + geom_density() +
  theme_bw(base_size = 20) +
  ggtitle("Beta density")
plot <- plot +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(size=20),
        legend.position="bottom") +
  scale_color_hue(name="Condition",guide=guide_legend(nrow=3))
ggsave(paste("BetaDensity.png",sep=""),plot)

#mean DNA methylation
plot<- ggplot(HSNC2, aes(x=Condition,y=Beta_value, color=Condition)) + geom_boxplot()
ggsave(paste("meanDNAMetByCondition.png",sep=""),plot)

rm(plot)
#----------WILCOX P_VALUE---------------------
library(plyr)
#calculate mean for each Composite Element REF grouped by Condition
colnames(HSNC2) <- c("TCGA","cg","Beta_value","Condition")

diffMSet <- ddply(HSNC2, .(cg), summarize,
                  NegBMean=mean(Beta_value[Condition=="Negative_black_or_african_american"]),
                  NegWMean=mean(Beta_value[Condition=="Negative_white"]),
                  PosWMean=mean(Beta_value[Condition=="Positive_white"]),
                  Normal=mean(Beta_value[Condition=="Normal"]))
diffMSet <- subset(diffMSet,subset = rowSums(is.na(diffMSet))==0)  #remove NA

diff_nb_nw <- diffMSet[,c(1,2,3)]
diff_nb_pw <- diffMSet[,c(1,2,4)]
diff_nw_pw <- diffMSet[,c(1,3,4)]
diff_nw_normal <- diffMSet[,c(1,3,5)]
diff_nb_normal <- diffMSet[,c(1,2,5)]
diff_pw_normal <- diffMSet[,c(1,4,5)]

diff_nb_nw[,"Wilcox_pvalue"] <- NA
diff_nb_pw[,"Wilcox_pvalue"] <- NA
diff_nw_pw[,"Wilcox_pvalue"] <- NA
diff_nw_normal[,"Wilcox_pvalue"] <- NA
diff_nb_normal[,"Wilcox_pvalue"] <- NA
diff_pw_normal[,"Wilcox_pvalue"] <- NA

diff_nb_nw[,"p_adjusted"] <- NA
diff_nb_pw[,"p_adjusted"] <- NA
diff_nw_pw[,"p_adjusted"] <- NA
diff_nw_normal[,"p_adjusted"] <- NA
diff_nb_normal[,"p_adjusted"] <- NA
diff_pw_normal[,"p_adjusted"] <- NA

#create 3 individual data set based on Condition
nb <- HSNC2[HSNC2$Condition=="Negative_black_or_african_american",]
nb <- subset(nb,subset = rowSums(is.na(nb))==0)  #remove NA

nw <- HSNC2[HSNC2$Condition=="Negative_white",]
nw <- subset(nw,subset = rowSums(is.na(nw))==0)  #remove NA

pw <- HSNC2[HSNC2$Condition=="Positive_white",]
pw <- subset(pw,subset = rowSums(is.na(pw))==0)  #remove NA

Normal <- HSNC2[HSNC2$Condition=="Normal",]
Normal <- subset(Normal,subset = rowSums(is.na(pw))==0)  #remove NA

nb$Beta_value <- as.numeric(nb$Beta_value)
nw$Beta_value <- as.numeric(nw$Beta_value)
pw$Beta_value <- as.numeric(pw$Beta_value)
Normal$Beta_value <- as.numeric(Normal$Beta_value)

save.image("methylation.RData")

#-----------------------parallel
#library(foreach)
#library(doParallel)
#cores=detectCores()
#cl <- makeCluster(cores[1]-2)
#registerDoParallel(cl)

#gene <- list()
#gene <- foreach (i = iter(diffMSet$cg[1:150000],by='row'), .combine=rbind) %dopar% {
#cgInfo$Gene[cgInfo$cg == i]
#}

#stopCluster(cl)
#registerDoSEQ()

#diffMSet$Gene <- gene

##--------------parallel---------------
##install.packages("foreach")
##install.packages("doParallel")
#library(foreach)
#library(doParallel)
##setup parallel backend to use many processors
#cores=detectCores()
#cl <- makeCluster(cores[1]-2) #not to overload your computer
#registerDoParallel(cl)

#foreach(i=seq_along(unique(diff_nb_nw$cg))) %dopar% {
#  cg <- unique(diff_nb_nw$cg)[[i]]
#  wilcoxp <- as.numeric(wilcox.test(as.numeric(nb$Beta_value[nb$cg==cg]),
#                                    as.numeric(nw$Beta_value[nw$cg==cg]))$p.value)
#  diff_nb_nw$Wilcox_pvalue[diff_nb_nw$cg==cg] <- wilcoxp
#}
##stop cluster
#stopCluster(cl)

#-----------------WILCOX P AND P ADJUSTED
#4-diff_nb_normal
Sys.time()
list <- unique(diff_nb_normal$cg)
for (i in list){
  print(i)
  b1 <- nb$Beta_value[nb$cg==i]
  b2 <- Normal$Beta_value[Normal$cg==i]
  wilcoxp <- as.numeric(wilcox.test(b1,b2)$p.value)
  diff_nb_normal$Wilcox_pvalue[diff_nb_normal$cg==i] <- wilcoxp
}

bh <- p.adjust(diff_nb_normal$Wilcox_pvalue, method="BH")
diff_nb_normal$p_adjusted <- bh
write.table(diff_nb_normal, file = "diff_nb_normal.csv", sep = "\t",row.names = F, quote=FALSE)
Sys.time()

#5-diff_nw_normal
Sys.time()
list <- unique(diff_nw_normal$cg)
for (i in list){
  print(i)
  b3 <- nw$Beta_value[nw$cg==i]
  b4 <- Normal$Beta_value[Normal$cg==i]
  diff_nw_normal$Wilcox_pvalue[diff_nw_normal$cg==i] <- as.numeric(wilcox.test(b3,b4)$p.value)
}
bh <- p.adjust(diff_nw_normal$Wilcox_pvalue, method="BH")
diff_nw_normal$p_adjusted <- bh
write.table(diff_nw_normal, file = "diff_nw_normal.csv", sep = "\t",row.names = F, quote=FALSE)
Sys.time()

#6-diff_pw_normal
Sys.time()
list <- unique(diff_pw_normal$cg)
for (i in list){
  print(i)
  b5 <- pw$Beta_value[pw$cg==i]
  b6 <- Normal$Beta_value[Normal$cg==i]
  diff_pw_normal$Wilcox_pvalue[diff_pw_normal$cg==i] <- as.numeric(wilcox.test(b5,b6)$p.value)
}
bh <- p.adjust(diff_pw_normal$Wilcox_pvalue, method="BH")
diff_pw_normal$p_adjusted <- bh
write.table(diff_pw_normal, file = "diff_pw_normal.csv", sep = "\t",row.names = F, quote=FALSE)
Sys.time()

#1-diff_nb_nw
load("methylation.RData")
library(plyr)
Sys.time()
for (i in unique(diff_nb_nw$cg)){
  wilcoxp <- as.numeric(wilcox.test(as.numeric(nb$Beta_value[nb$cg==i]),
                                    as.numeric(nw$Beta_value[nw$cg==i]))$p.value)
  diff_nb_nw$Wilcox_pvalue[diff_nb_nw$cg==i] <- wilcoxp
  write(paste(i, wilcoxp, sep="\t"), file="diff_nb_nw.out", append=TRUE)
}
bh <- p.adjust(diff_nb_nw$Wilcox_pvalue, method="BH")
diff_nb_nw$p_adjusted <- bh
write.table(diff_nb_nw, file = "diff_nb_nw.csv", sep = "\t",row.names = F, quote=FALSE)
Sys.time()

#2-diff_nb_pw
Sys.time()
for (i in unique(diff_nb_pw$cg)){
  wilcoxp <- as.numeric(wilcox.test(as.numeric(nb$Beta_value[nb$cg==i]),
                                    as.numeric(pw$Beta_value[pw$cg==i]))$p.value)
  diff_nb_pw$Wilcox_pvalue[diff_nb_pw$cg==i] <- wilcoxp
  write(paste(i, wilcoxp, sep="\t"),file="diff_nb_pw.out", append=TRUE)
}
bh <- p.adjust(diff_nb_pw$Wilcox_pvalue, method="BH")
diff_nb_pw$p_adjusted <- bh
write.table(diff_nb_pw, file = "diff_nb_pw.csv", sep = "\t",row.names = F, quote=FALSE)
Sys.time()

#3-diff_nw_pw
Sys.time()
for (i in unique(diff_nw_pw$cg)){
  wilcoxp <- as.numeric(wilcox.test(as.numeric(nw$Beta_value[nw$cg==i]),
                                    as.numeric(pw$Beta_value[pw$cg==i]))$p.value)
  diff_nw_pw$Wilcox_pvalue[diff_nw_pw$cg==i] <- wilcoxp
  write(paste(i, wilcoxp, sep="\t"),file="diff_nw_pw.out", append=TRUE)
}
bh <- p.adjust(diff_nw_pw$Wilcox_pvalue, method="BH")
diff_nw_pw$p_adjusted <- bh
write.table(diff_nw_pw, file = "diff_nw_pw.csv", sep = "\t",row.names = F, quote=FALSE)
Sys.time()
save.image("methylation.RData")

#---------------------------GET THE UNIQUE GENE list for each cg
#cgInfo is one of the text file
cgInfo <- read.table("cgInfo", header=TRUE, sep="\t",stringsAsFactors=TRUE, check.names=F)
cgInfo <- cgInfo[c("Composite Element REF","Gene_Symbol")]

#unique_gene_per_cg <- "NA"
#for (i in cgInfo$Gene_Symbol){
#  genes <- strsplit(i,";")
#  unique_gene_per_cg <- rbind(unique_gene_per_cg,unique(unlist(genes)))
#}
#write.table(unique_gene_per_cg, file="unique_gene_per_cg.csv",sep = "\t",row.names = F, quote=FALSE)
#geneList <- read.table("unique_gene_per_cg",header=TRUE, sep="\t")

#Do this if unique gene per cg has 1 extra length at the end of list
#length(unique_gene_per_cg)
#unique_gene_per_cg <- unique_gene_per_cg[-485578]

cgInfo$Gene <- "NA"
cgInfo$Gene <- unique_gene_per_cg

ind <- match(diffMSet$cg, cgInfo[,"Composite Element REF"])
diffMSet$Gene <- cgInfo[ind,"Gene"]
ind <- match(diff_nb_normal$cg, cgInfo[,"Composite Element REF"])
diff_nb_normal$Gene <- cgInfo[ind,"Gene"]
ind <- match(diff_nw_normal$cg, cgInfo[,"Composite Element REF"])
diff_nw_normal$Gene <- cgInfo[ind,"Gene"]
ind <- match(diff_pw_normal$cg, cgInfo[,"Composite Element REF"])
diff_pw_normal$Gene <- cgInfo[ind,"Gene"]
save.image("methylation.RData")

gene <- list()
for (i in diffMSet$cg){
  print(i)
  gene <- rbind(gene, cgInfo$Gene[cgInfo$Composite.Element.REF == i])
}

colnames(cgInfo) <- c("cg","Gene_Symbol","Gene")
#-----------------------------------VOCANO PLOT- FOLD CHANGE
#install.packages("gtools")
#install.packages("plotly")
library(gtools)

#-------------------nb_normal------
diff_nb_normal$log2FoldChange <- NA
fc <- foldchange(diff_nb_normal[,2],diff_nb_normal[,3])
fc <- foldchange2logratio(fc, base=2)
diff_nb_normal$log2FoldChange <- fc

#GROUP by pvalue and lof2FoldChange value for COLOR
#4 groups: SignificantFC, Significant, FC, NA
diff_nb_normal$SignificantStatus <- NA
diff_nb_normal[which(diff_nb_normal['Wilcox_pvalue'] < 0.05 &
                       abs(diff_nb_normal['log2FoldChange']) > 1),"SignificantStatus"] <- "SignificantFC"
diff_nb_normal[which(diff_nb_normal['Wilcox_pvalue'] < 0.05 &
                       abs(diff_nb_normal['log2FoldChange']) < 1),"SignificantStatus"] <- "Significant"
diff_nb_normal[which(diff_nb_normal['Wilcox_pvalue'] > 0.05 &
                       abs(diff_nb_normal['log2FoldChange']) > 1),"SignificantStatus"] <- "FC"

#significantFC used for plot's label
sigFC <- diff_nb_normal[with(diff_nb_normal, order(log2FoldChange, Wilcox_pvalue)),][1:30,]
sigFC <- rbind(sigFC, diff_nb_normal[with(diff_nb_normal, order(-log2FoldChange, Wilcox_pvalue)),][1:30,])
sigFC <- sigFC[ which(sigFC$SignificantStatus=="SignificantFC"),]

sigFC$Gene <- "NA"
sigFCGeneList<- list()
ind <- match(sigFC$cg, cgInfo$cg)
sigFC$Gene <- cgInfo[ind,"Gene"]

#Plot
#install.packages("ggrepel")
require("ggrepel")

#plot <- ggplot(diff_nb_normal, aes(x=log2FoldChange,y=-log10(Wilcox_pvalue), color=SignificantStatus)) + geom_point()
plot <- ggplot(diff_nb_normal, aes(x=log2FoldChange,y=-log10(Wilcox_pvalue))) +
  theme_bw(base_size = 20) +
  geom_point(aes(color=SignificantStatus), legend = F) +
  ggtitle("Differences between NB and normal") +
  geom_text_repel(
    data = sigFC,
    aes(label = unlist(Gene)),
    size = 4,
    box.padding = unit(0.3, "lines"),
    point.padding = unit(0.3, "lines"),
    show.legend = FALSE
  )
plot <- plot +  scale_color_hue(name="",
                                labels = c("Fold Change (FC >1)", "Significant (p < 0.05)",
                                           "SignificantFC (p< 0.05 & FC >1)","NA"),
                                guide=guide_legend(nrow=2))
plot <- plot +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(size=20),
        legend.position="bottom")
ggsave("volcano_diff_nb_normal.png", plot)

#write top 100 cgs list
top_peaks <- diff_nb_normal[with(diff_nb_normal, order(log2FoldChange, Wilcox_pvalue)),][1:150000,]
top_peaks <- rbind(top_peaks, diff_nb_normal[with(diff_nb_normal, order(-log2FoldChange, Wilcox_pvalue)),][1:150000,])
top_peaks <- top_peaks[ which(top_peaks$SignificantStatus=="SignificantFC"),]

#match(x,y) where is x in y, x = 1,2 y = 3,2,1 the answer is 3,2
index <- match(top_peaks$cg, cgInfo$cg)
top_peaks <- cgInfo[index,]
#combine top peaks with column Gene_Symbol of cgInfo
top_peaks<- merge(top_peaks, diff_nb_normal,by="cg")
#apply(top_peaks,2,as.character)
top_peaks = as.matrix(top_peaks)

write.table(top_peaks,file="diff_nb_normal.csv",sep = "\t",row.names = F, quote=FALSE)


#--------------pw_normal--
diff_pw_normal$log2FoldChange <- NA
fc <- foldchange(diff_pw_normal[,2],diff_pw_normal[,3])
fc <- foldchange2logratio(fc, base=2)
diff_pw_normal$log2FoldChange <- fc

#GROUP by pvalue and lof2FoldChange value for COLOR
#4 groups: SignificantFC, Significant, FC, NA
diff_pw_normal$SignificantStatus <- NA
diff_pw_normal[which(diff_pw_normal['Wilcox_pvalue'] < 0.05 &
                       abs(diff_pw_normal['log2FoldChange']) > 1),"SignificantStatus"] <- "SignificantFC"
diff_pw_normal[which(diff_pw_normal['Wilcox_pvalue'] < 0.05 &
                       abs(diff_pw_normal['log2FoldChange']) < 1),"SignificantStatus"] <- "Significant"
diff_pw_normal[which(diff_pw_normal['Wilcox_pvalue'] > 0.05 &
                       abs(diff_pw_normal['log2FoldChange']) > 1),"SignificantStatus"] <- "FC"

#significantFC used for plot's label
sigFC <- diff_pw_normal[with(diff_pw_normal, order(log2FoldChange, Wilcox_pvalue)),][1:30,]
sigFC <- rbind(sigFC, diff_pw_normal[with(diff_pw_normal, order(-log2FoldChange, Wilcox_pvalue)),][1:30,])
sigFC <- sigFC[ which(sigFC$SignificantStatus=="SignificantFC"),]

sigFC$Gene <- "NA"
sigFCGeneList<- list()
ind <- match(sigFC$cg, cgInfo$cg)
sigFC$Gene <- cgInfo[ind,"Gene"]

#Plot
#install.packages("ggrepel")
require("ggrepel")

#plot <- ggplot(diff_pw_normal, aes(x=log2FoldChange,y=-log10(Wilcox_pvalue), color=SignificantStatus)) + geom_point()
plot <- ggplot(diff_pw_normal, aes(x=log2FoldChange,y=-log10(Wilcox_pvalue))) +
  theme_bw(base_size = 20) +
  geom_point(aes(color=SignificantStatus), legend = F) +
  ggtitle("Differences between PW and normal") +
  geom_text_repel(
    data = sigFC,
    aes(label = unlist(Gene)),
    size = 4,
    box.padding = unit(0.3, "lines"),
    point.padding = unit(0.3, "lines"),
    show.legend = FALSE
  )
plot <- plot +  scale_color_hue(name="",
                                labels = c("Fold Change (FC >1)", "Significant (p < 0.05)",
                                           "SignificantFC (p< 0.05 & FC >1)","NA"),
                                guide=guide_legend(nrow=2))
plot <- plot +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(size=20),
        legend.position="bottom")
ggsave("volcano_diff_pw_normal.png", plot)



#write top 100 cgs list
top_peaks <- diff_pw_normal[with(diff_pw_normal, order(log2FoldChange, Wilcox_pvalue)),][1:150000,]
top_peaks <- rbind(top_peaks, diff_pw_normal[with(diff_pw_normal, order(-log2FoldChange, Wilcox_pvalue)),][1:150000,])
top_peaks <- top_peaks[ which(top_peaks$SignificantStatus=="SignificantFC"),]

#match(x,y) where is x in y, x = 1,2 y = 3,2,1 the answer is 3,2
index <- match(top_peaks$cg, cgInfo$cg)
top_peaks <- cgInfo[index,]
#combine top peaks with column Gene_Symbol of cgInfo
top_peaks<- merge(top_peaks, diff_pw_normal,by="cg")
#apply(top_peaks,2,as.character)
top_peaks = as.matrix(top_peaks)

write.table(top_peaks,file="diff_pw_normal.csv",sep = "\t",row.names = F, quote=FALSE)

#------------NW_normal------
diff_nw_normal$log2FoldChange <- NA
fc <- foldchange(diff_nw_normal[,2],diff_nw_normal[,3])
fc <- foldchange2logratio(fc, base=2)
diff_nw_normal$log2FoldChange <- fc

#GROUP by pvalue and lof2FoldChange value for COLOR
#4 groups: SignificantFC, Significant, FC, NA
diff_nw_normal$SignificantStatus <- NA
diff_nw_normal[which(diff_nw_normal['Wilcox_pvalue'] < 0.05 &
                       abs(diff_nw_normal['log2FoldChange']) > 1),"SignificantStatus"] <- "SignificantFC"
diff_nw_normal[which(diff_nw_normal['Wilcox_pvalue'] < 0.05 &
                       abs(diff_nw_normal['log2FoldChange']) < 1),"SignificantStatus"] <- "Significant"
diff_nw_normal[which(diff_nw_normal['Wilcox_pvalue'] > 0.05 &
                       abs(diff_nw_normal['log2FoldChange']) > 1),"SignificantStatus"] <- "FC"

#significantFC used for plot's label
sigFC <- diff_nw_normal[with(diff_nw_normal, order(log2FoldChange, Wilcox_pvalue)),][1:30,]
sigFC <- rbind(sigFC, diff_nw_normal[with(diff_nw_normal, order(-log2FoldChange, Wilcox_pvalue)),][1:30,])
sigFC <- sigFC[ which(sigFC$SignificantStatus=="SignificantFC"),]

sigFC$Gene <- "NA"
sigFCGeneList<- list()
ind <- match(sigFC$cg, cgInfo$cg)
sigFC$Gene <- cgInfo[ind,"Gene"]

#Plot
#install.packages("ggrepel")
require("ggrepel")

#plot <- ggplot(diff_nw_normal, aes(x=log2FoldChange,y=-log10(Wilcox_pvalue), color=SignificantStatus)) + geom_point()
plot <- ggplot(diff_nw_normal, aes(x=log2FoldChange,y=-log10(Wilcox_pvalue))) +
  theme_bw(base_size = 20) +
  geom_point(aes(color=SignificantStatus), legend = F) +
  ggtitle("Differences between NW and normal") +
  geom_text_repel(
    data = sigFC,
    aes(label = unlist(Gene)),
    size = 4,
    box.padding = unit(0.3, "lines"),
    point.padding = unit(0.3, "lines"),
    show.legend = FALSE
  )
plot <- plot +  scale_color_hue(name="",
                                labels = c("Fold Change (FC >1)", "Significant (p < 0.05)",
                                           "SignificantFC (p< 0.05 & FC >1)","NA"),
                                guide=guide_legend(nrow=2))
plot <- plot +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(size=20),
        legend.position="bottom")
ggsave("volcano_diff_nw_normal.png", plot)

#write top 100 cgs list
top_peaks <- diff_nw_normal[with(diff_nw_normal, order(log2FoldChange, Wilcox_pvalue)),][1:150000,]
top_peaks <- rbind(top_peaks, diff_nw_normal[with(diff_nw_normal, order(-log2FoldChange, Wilcox_pvalue)),][1:150000,])
top_peaks <- top_peaks[ which(top_peaks$SignificantStatus=="SignificantFC"),]

#match(x,y) where is x in y, x = 1,2 y = 3,2,1 the answer is 3,2
index <- match(top_peaks$cg, cgInfo$cg)
top_peaks <- cgInfo[index,]
#combine top peaks with column Gene_Symbol of cgInfo
top_peaks<- merge(top_peaks, diff_nw_normal,by="cg")
#apply(top_peaks,2,as.character)
top_peaks = as.matrix(top_peaks)

write.table(top_peaks,file="diff_nw_normal.csv",sep = "\t",row.names = F, quote=FALSE)


# ------- Mutational burden of all mutation using maf data- visualized by violin plot per MB
# W vs AA
maf <- read.csv("TCGA.HNSC.mutect.1aa33f25-3893-4f37-a6a4-361c9785d07e.DR-10.0.somatic.maf",
                skip = 5 , header = TRUE, sep='\t', stringsAsFactors=TRUE, check.names=F)

Barcode <- maf$Tumor_Sample_Barcode
Barcode <- substring(Barcode, 1, 12)
maf[,"Barcode"] <- NA
maf$Barcode <- Barcode

maf.tb <- maf[,c("Hugo_Symbol","Barcode")]
colnames(maf.tb) <- c("Gene","Barcode")
maf.tb <- table(maf.tb)
maf.tb <- as.data.frame.matrix(maf.tb)

for(i in 2:ncol(maf.tb)) {
  maf.tb[,i] <- as.numeric(maf.tb[,i])
}

maf.tb <- apply(maf.tb, 2, sum)
maf.tb <- as.data.frame(maf.tb)
maf.tb$TCGA <- rownames(maf.tb)
colnames(maf.tb)[colnames(maf.tb)=="maf.tb"] <- "Count"

ind <- match(HPV$TCGA, maf.tb$TCGA)
ind <- ind[!is.na(ind)]
maf.tb <- maf.tb[ind,] #101 patients

ind <- match(maf.tb$TCGA, HPV$TCGA)
ind <- ind[!is.na(ind)]
maf.tb$Condition <- HPV$Condition[ind]

maf.tb$Condition <- as.character(maf.tb$Condition)
maf.tb[which(maf.tb$Condition=="Negative_black_or_african_american"),"Condition"] <- "African_american"
maf.tb$Condition <- as.factor(maf.tb$Condition)

library(ggplot2)

tmp <- maf.tb
#tmp <- maf.tb[maf.tb$Condition!="African_american",]

# tmp <- maf.tb
# tmp$Condition <- as.character(tmp$Condition)
# tmp[which(tmp$Condition!="African_american"),"Condition"] <- "White"
# tmp$Condition <- as.factor(tmp$Condition)

t <- tmp$Count/3200
tmp$Counts_per_MB <- t

plot <- ggplot(tmp, aes(x=Condition, y=Counts_per_MB, color=Condition)) +
  geom_violin(trim=FALSE) +
  theme_bw(base_size = 20) +
  ggtitle("Mutational burden")
plot <- plot + stat_summary(fun.data=mean_sdl, mult=1,
                 geom="pointrange", color="black")
plot <- plot +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(size=20),
        legend.position="bottom") +
  scale_color_hue(name="Condition",guide=guide_legend(nrow=3))
plot <- plot + stat_compare_means()

ggsave(paste("Mutational_burden_AA_vs_W.svg",sep=""),plot)

# ---SURVIVAL ANALYSIS------  AA vs W  ----- Pos vs Neg W  -------  All
library("survival")
library("survminer")

clin <- read.csv("clinical.tsv",header=TRUE, sep="\t", check.names=FALSE, stringsAsFactors=F)
clin2 <- clin[,c("submitter_id","vital_status","days_to_death","days_to_last_follow_up")] #623

# for all alive patients, replace all "--" of days_to_death by value of days_to_last_follow_up
for(i in 1:nrow(clin2)){
  ind <- grep("--", clin2[, "days_to_death"])
  clin2[ind, 3] <- clin2[ind, "days_to_last_follow_up"]
}
clin2 <- clin2[clin2$vital_status!="--",]

clin2 <- subset(clin2, select=-c(days_to_last_follow_up))
colnames(clin2)[3] <- "time"

#match HPV cases and clin cases
ind <- which(clin2$submitter_id %in% HPV$TCGA)
clin2 <- clin2[ind,]

ind <- match(clin2$submitter_id, HPV$TCGA)
clin2$Condition <- HPV$Condition[ind]

clin2$Condition <- as.character(clin2$Condition)
clin2[which(clin2$Condition=="Negative_black_or_african_american"),"Condition"] <- "African_american"
clin2$Condition <- as.factor(clin2$Condition)

clin2$time <- as.numeric(clin2$time)
# create survObj, if alive, add "+" to time as new value; if dead, leave as is
clin2$survObj <- with(clin2, Surv(time, vital_status == "dead"))

tmp <- clin2

# tmp <- clin2
# tmp$Condition <- as.character(tmp$Condition)
# tmp[which(tmp$Condition!="African_american"),"Condition"] <- "White"
# tmp$Condition <- as.factor(tmp$Condition)

#tmp <- clin2[clin2$Condition!="African_american",]

fit <- survfit(tmp$survObj ~ tmp$Condition , data = tmp)
svg(file="survival_pos_vs_neg_W.svg")
ggsurvplot(fit, data= tmp, font.tickslab = 10, pval= T, risk.table = T)
dev.off()
