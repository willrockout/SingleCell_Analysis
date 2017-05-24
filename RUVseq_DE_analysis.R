source("http://bioconductor.org/biocLite.R")
biocLite("RUVSeq")
biocLite("EDASeq")
library(RUVSeq)
library(EDASeq)
library(RColorBrewer)
library(edgeR)
library(lattice)
library(gplots)
library(calibrate)

####Getting HTSeq Gene/CDS Count Data and Combining into 1 table
Gene_count_dir <- list.files("./SC_gene_Counts/Gene_counts",full=TRUE)
Coding_region_dir <- list.files("./SC_gene_Counts/Coding_region_counts",full=TRUE)

FACS_gene_counts <- Gene_count_dir[grep("FACS",Gene_count_dir)]
FACS_coding_counts <- Coding_region_dir[grep("FACS",Gene_count_dir)]

C1_gene_counts <- Gene_count_dir[grep("C1",Gene_count_dir)]
C1_coding_counts <- Coding_region_dir[grep("C1",Gene_count_dir)]

FACS.Stress_gene_counts <- FACS_gene_counts[grep("Stressed",FACS_gene_counts)]
FACS.Control_gene_counts <- FACS_gene_counts[grep("Control",FACS_gene_counts)]
FACS.Stress_coding_counts <- FACS_coding_counts[grep("Stressed",FACS_coding_counts)]
FACS.Control_coding_counts <- FACS_coding_counts[grep("Control",FACS_coding_counts)]

C1.Stress_gene_counts <- C1_gene_counts[grep("Stressed",C1_gene_counts)]
C1.Control_gene_counts <- C1_gene_counts[grep("Control",C1_gene_counts)]
C1.Stress_coding_counts <- C1_coding_counts[grep("Stressed",C1_coding_counts)]
C1.Control_coding_counts <- C1_coding_counts[grep("Control",C1_coding_counts)]

Genecount_FACSS <- do.call(cbind,lapply(FACS.Stress_gene_counts, function(x) read.table(x,
                              col.names=basename(x))))
Genecount_FACSC <- do.call(cbind,lapply(FACS.Control_gene_counts,function(x) read.table(x,
                              col.names=basename(x))))
CDScount_FACSS <- do.call(cbind,lapply(FACS.Stress_coding_counts,function(x) read.table(x,
                             col.names=basename(x))))
CDScount_FACSC <- do.call(cbind,lapply(FACS.Control_coding_counts,function(x)
    read.table(x,col.names=basename(x))))

Genecount_C1S <- do.call(cbind,lapply(C1.Stress_gene_counts,function(x) read.table(x,
                              col.names=basename(x))))
Genecount_C1C <- do.call(cbind,lapply(C1.Control_gene_counts,function(x) read.table(x,
                              col.names=basename(x))))
CDScount_C1S <- do.call(cbind,lapply(C1.Stress_coding_counts,function(x) read.table(x,
                             col.names=basename(x))))
CDScount_C1C <- do.call(cbind,lapply(C1.Control_coding_counts,function(x) read.table(x,
                             col.names=basename(x))))


Total_Facs_Gene <- cbind(Genecount_FACSC,Genecount_FACSS)
Total_Facs_CDS <- cbind(CDScount_FACSS,CDScount_FACSC)

Total_C1_Gene <- cbind(Genecount_C1C,Genecount_C1S)
Total_C1_CDS <- cbind(CDScount_C1S,CDScount_C1C)

#####Selecting For read cutoff
Read_Count_Cut <- Total_Facs_Gene[,]>=5
Tr <- data.frame(Values=colSums(Read_Count_Cut==TRUE,na.rm=TRUE))
Fal <- data.frame(Values=colSums(Read_Count_Cut==FALSE,na.rm=TRUE))
Combined <- cbind(Tr,Fal)
colnames(Combined) <- c("True","False")
Combined$Total <- rowSums(Combined)
Combined_Percent <- Combined$True/Combined$Total
df <- data.frame(Combined_Percent)
rownames(df) <- rownames(Combined)
df <- df * 100

C1_Read_Count_Cut <- Total_C1_Gene[,]>=5
C1Tr <- data.frame(Values=colSums(C1_Read_Count_Cut==TRUE,na.rm=TRUE))
C1Fal <- data.frame(Values=colSums(C1_Read_Count_Cut==FALSE,na.rm=TRUE))
CombinedC1 <- cbind(C1Tr,C1Fal)
colnames(CombinedC1) <- c("True","False")
CombinedC1$Total <- rowSums(CombinedC1)
Combined_PercentC1 <- CombinedC1$True/CombinedC1$Total
C1df <- data.frame(Combined_PercentC1)
rownames(C1df) <- rownames(CombinedC1)
C1df <- C1df * 100


#####Plotting Percentage of read counts
layout(1:2, heights=c(6, 6))
par(mar=c(5,5,3,2))
grange<-range(0,100)
barplot(as.matrix(t(C1df)),ylim=grange,col="cyan",ylab="%",
        main="Percentage of C1 Read Counts/Gene 5 or Greater",
        las=3,cex.names=.5,xlab="",xaxt='n',yaxt='n',
        cex.lab=1.5,cex.main=1.5,yaxt='n')
axis(1,c(51,157),lab=c("C1 Control","C1 Stressed"),cex.axis=1.5)
axis(2,las=1)
abline(v=104,col="black")
par(mar=c(5,5,3,2))
grange<-range(0,100)
barplot(as.matrix(t(df)),ylim=grange,col="cyan",ylab="%",
        main="Percentage of FACS Read Counts/Gene 5 or Greater",
        las=3,cex.names=.5,xlab="",xaxt='n',yaxt='n',
        cex.lab=1.5,cex.main=1.5,yaxt='n')
axis(1,c(45,145),lab=c("FACS Control","FACS Stressed"),cex.axis=1.5)
axis(2,las=1)
abline(v=89,col="black")

######Boxplot of Distributions
layout(1:2, heights=c(6, 6))
par(mar=c(5,5,3,2))
boxplot(log2(Total_Facs_Gene[1:6606,]+1),col="cyan",
        pch=20,xaxt='n',main="FACS Read Count/Gene Distribution",
        ylab="log2(count+1)",cex.lab=1.5)
axis(1,c(48,144),lab=c("FACS Stressed","FACS Control"),cex.axis=1.5)
abline(v=96.5,col="black")
par(mar=c(5,5,3,2))
boxplot(log2(Total_Facs_CDS[1:6606,]+1),col="cyan",
        pch=20,xaxt='n',main="FACS Read Count/CDS Distribution",
        ylab="log2(count+1)",cex.lab=1.5)
axis(1,c(48,144),lab=c("FACS Stressed","FACS Control"),ylab="log2(count+1)",cex.axis=1.5)
abline(v=96.5,col="black")


layout(1:2, heights=c(6, 6))
par(mar=c(5,5,3,2))
boxplot(log2(Total_C1_Gene[1:6606,]+1),col="cyan",
        pch=20,xaxt='n',main="C1 Read Count/Gene Distribution",
        ylab="log2(count+1)",cex.lab=1.5)
axis(1,c(48,144),lab=c("C1 Stressed","C1 Control"),cex.axis=1.5)
abline(v=96.5,col="black")
par(mar=c(5,5,3,2))
boxplot(log2(Total_C1_CDS[1:6606,]+1),col="cyan",
        pch=20,xaxt='n',main="C1 Read Count/CDS Distribution",
        ylab="log2(count+1)",cex.lab=1.5)
axis(1,c(48,144),lab=c("C1 Stressed","C1 Control"),cex.axis=1.5)
abline(v=96.5,col="black")

######Spike in genes counts
SpikeFilesC1 <- list.files("./Spiked_in/Bowtie2_Counts/C1",full=TRUE)
SpikeFilesFACS <- list.files("./Spiked_in/Bowtie2_Counts/FACS",full=TRUE)
ControlsC1 <- lapply(SpikeFilesC1,function(x) read.table(x,row.names=2))
ControlsFACS <- lapply(SpikeFilesFACS,function(x) read.table(x,row.names=2))

C1rowNames <- unique(unlist(lapply(ControlsC1, rownames)))
C1spikeCounts <- data.frame(row.names=C1rowNames)
C1splits <- strsplit(basename(SpikeFilesC1), "\\.")
C1colNames <- unlist(lapply(C1splits, function(x) x[1]))
names(SpikeFilesC1) <- C1colNames
names(ControlsC1) <- C1colNames
for(name in C1colNames){
    dC1 <- ControlsC1[[name]]
    C1spikeCounts[,name] <- 0
    C1spikeCounts[match(rownames(dC1), rownames(C1spikeCounts)), name] <- dC1[,1] 
}

FACSrowNames <- unique(unlist(lapply(ControlsFACS, rownames)))
FACSspikeCounts <- data.frame(row.names=FACSrowNames)
FACSsplits <- strsplit(basename(SpikeFilesFACS), "\\.")
FACScolNames <- unlist(lapply(FACSsplits, function(x) x[1]))
names(SpikeFilesFACS) <- FACScolNames
names(ControlsFACS) <- FACScolNames
for(name in FACScolNames){
    dFACS <- ControlsFACS[[name]]
    FACSspikeCounts[,name] <- 0
    FACSspikeCounts[match(rownames(dFACS), rownames(FACSspikeCounts)), name] <- dFACS[,1] 
}



######Filtered and Exploratory Data Analysis

C1filter <- apply(Total_C1_Gene, 1, function(x) length(x[x>5])>=2)
C1filtered <- Total_C1_Gene[C1filter,]
C1genes <- rownames(C1filtered)[-grep("^__",rownames(C1filtered))]
C1counts <- C1filtered[C1genes,]
C1x <- as.factor(rep(c("Ctl", "Trt"),c(83,89)))
C1colSplit <- strsplit(colnames(C1counts), "_")
colnames(C1counts) <- unlist(lapply(C1colSplit, function(x) x[1]))
c1 <- C1spikeCounts[,match(colnames(C1counts), colnames(C1spikeCounts))]
C1forRUV <- rbind(C1counts, c1)
cIdx <- rownames(C1forRUV) %in% rownames(c1)
C1set <- newSeqExpressionSet(as.matrix(C1forRUV),
phenoData = data.frame(C1x, row.names=colnames(C1forRUV)))
colors <- brewer.pal(3, "Set2")
#sum(FACSforRUV[,155]==0,na.rm=TRUE)
C1pdno <- pData(C1set)

C1designno <- model.matrix(~C1x , data=C1pdno)
C1yno <- DGEList(counts=counts(C1set), group=C1x)
C1yno <- calcNormFactors(C1yno, method="none")
#C1yno <- calcNormFactors(C1yno, method="RLE")
C1yno <- estimateGLMCommonDisp(C1yno, C1designno)
C1yno <- estimateGLMTagwiseDisp(C1yno, C1designno)
C1fitno <- glmFit(C1yno, C1designno)
C1lrtno <- glmLRT(C1fitno,coef=2)

C1deno <- decideTestsDGE(C1lrtno,lfc=2)
topTags(C1lrt1)
summary(de <- decideTestsDGE(C1lrt1))
C1detags1 <- rownames(C1y1)[as.logical(de)]
plotSmear(C1lrt1, de.tags=C1detags1)
abline(h=c(-1, 1), col="blue")


FACSfilter <- apply(Total_Facs_Gene, 1, function(x) length(x[x>5])>=2)
FACSfiltered <- Total_Facs_Gene[FACSfilter,]
FACSgenes <- rownames(FACSfiltered)[-grep("^__",rownames(FACSfiltered))]
FACScounts <- FACSfiltered[FACSgenes,]
FACSx <- as.factor(rep(c("Ctl", "Trt"),c(75,92)))
FACScolSplit <- strsplit(colnames(FACScounts), "_")
colnames(FACScounts) <- unlist(lapply(FACScolSplit, function(x) x[1]))
facs <- FACSspikeCounts[,match(colnames(FACScounts), colnames(FACSspikeCounts))]
FACSforRUV <- rbind(FACScounts, facs)
cIdxF <- rownames(FACSforRUV) %in% rownames(facs)
FACSset <- newSeqExpressionSet(as.matrix(FACSforRUV),
phenoData = data.frame(FACSx, row.names=colnames(FACSforRUV)))
colors <- brewer.pal(3, "Set2")

FACSpdno <- pData(FACSset)

FACSdesignno <- model.matrix(~FACSx , data=FACSpdno)
FACSyno <- DGEList(counts=counts(FACSset), group=FACSx)
#FACSyno <- calcNormFactors(FACSyno, method="none")
FACSyno <- calcNormFactors(FACSyno, method="TMM")
FACSyno <- estimateGLMCommonDisp(FACSyno, FACSdesignno)
FACSyno <- estimateGLMTagwiseDisp(FACSyno, FACSdesignno)
FACSfitno <- glmFit(FACSyno, FACSdesignno)
FACSlrtno <- glmLRT(FACSfitno,coef=2)

FACSdeno <- decideTestsDGE(FACSlrtno,lfc=2)
topTags(FACSlrt1)
summary(de <- decideTestsDGE(FACSlrt1))
FACSdetags1 <- rownames(FACSy1)[as.logical(de)]
plotSmear(FACSlrt1, de.tags=FACSdetags1)
abline(h=c(-1, 1), col="blue")


layout(1:2, heights=c(6, 6))
par(mar=c(5,5,3,2))
plotRLE(C1set, outline=FALSE, ylim=c(-4, 4), col=colors[C1x],main="C1 Raw Data",xaxt='n')
plotRLE(FACSset, outline=FALSE, ylim=c(-4, 4), col=colors[FACSx],
        main="FACS Raw Data",xaxt='n')
axis(1,c(41.5,127.5),lab=c("Control","Stressed"))
axis(1,c(37.5,121.5),lab=c("Control","Stressed"))
plotPCA(C1set, col=colors[C1x], cex=1.2,main="C1 Raw Data")
plotPCA(FACSset, col=colors[FACSx], cex=1.2,main="FACS Raw Data")

#####RUVg

C1set1 <- RUVg(C1set,rownames(c1),k=1)
FACSset1 <- RUVg(FACSset,rownames(facs),k=1)

plotRLE(C1set1, outline=FALSE, ylim=c(-4, 4),
        col=colors[C1x],main="RUVg Normalized C1 Data",xaxt='n')
plotRLE(FACSset1, outline=FALSE, ylim=c(-4, 4),
        col=colors[FACSx],main="RUVg Normalized FACS Data",xaxt='n')
axis(1,c(37.5,121.5),lab=c("Control","Stressed"))
plotPCA(C1set1, col=colors[C1x], cex=1.2,main="C1 Normalized Data")
plotPCA(FACSset1, col=colors[FACSx], cex=1.2,main="FACS Normalized Data")

C1pd1 <- pData(C1set1)
FACSpd1 <- pData(FACSset1)

C1design1 <- model.matrix(~C1x + W_1, data=C1pd1)
C1y1 <- DGEList(counts=counts(C1set), group=C1x)
C1y1 <- calcNormFactors(C1y1, method="none")
#C1y1 <- calcNormFactors(C1y1, method="RLE")
C1y1 <- estimateGLMCommonDisp(C1y1, C1design1)
C1y1 <- estimateGLMTagwiseDisp(C1y1, C1design1)
C1fit1 <- glmFit(C1y1, C1design1)
C1lrt1 <- glmLRT(C1fit1,coef=2)

C1de1 <- decideTestsDGE(C1lrt1,lfc=2)
topTags(C1lrt1)
summary(de <- decideTestsDGE(C1lrt1))
C1detags1 <- rownames(C1y1)[as.logical(de)]
plotSmear(C1lrt1, de.tags=C1detags1)
abline(h=c(-1, 1), col="blue")

FACSdesign1 <- model.matrix(~FACSx + W_1, data=FACSpd1)
FACSy1 <- DGEList(counts=counts(FACSset), group=FACSx)
FACSy1 <- calcNormFactors(FACSy1, method="none")
#FACSy1 <- calcNormFactors(FACSy1, method="TMM")
FACSy1 <- estimateGLMCommonDisp(FACSy1, FACSdesign1)
FACSy1 <- estimateGLMTagwiseDisp(FACSy1, FACSdesign1)
FACSfit1 <- glmFit(FACSy1, FACSdesign1)
FACSlrt1 <- glmLRT(FACSfit1,coef=2)

FACSde1 <- decideTestsDGE(FACSlrt1,lfc=2)
topTags(FACSlrt1)
summary(de <- decideTestsDGE(FACSlrt1))
FACSdetags1 <- rownames(FACSy1)[as.logical(de)]
plotSmear(FACSlrt1, de.tags=FACSdetags1)
abline(h=c(-1, 1), col="blue")


####RUVs

##Creating differences matrix multplie methods none work


FACSdifferences = matrix(data= seq(1,184,by=1), byrow=TRUE, nrow=2)
FACSdifferences[1,76:92] <- -1
FACSdifferences[2,] <- seq(76,167,by=1)

C1differences = matrix(data= seq(1,178,by=1), byrow=TRUE, nrow=2)
C1differences[1,84:89] <- -1
C1differences[2,] <- seq(84,172,by=1)

C1set2 <- RUVs(C1set,C1genes,k=1,C1differences)
FACSset2 <- RUVs(FACSset,FACSgenes,k=1,FACSdifferences)

plotRLE(C1set2, outline=FALSE, ylim=c(-4, 4),
        col=colors[C1x],xaxt='n',main="C1 Normalized Data")
plotRLE(FACSset2, outline=FALSE, ylim=c(-4, 4),
        col=colors[FACSx],xaxt='n',main="FACS Normalized Data")

plotPCA(C1set2, col=colors[C1x], cex=1.2,main="C1 Normalized Data")
plotPCA(FACSset2, col=colors[FACSx], cex=1.2,main="FACS Normalized Data")

C1pd2 <- pData(C1set2)
FACSpd2 <- pData(FACSset2)

C1design2 <- model.matrix(~C1x + W_1, data=C1pd2)
C1y2 <- DGEList(counts=counts(C1set), group=C1x)
C1y2 <- calcNormFactors(C1y2, method="none")
#C1y2 <- calcNormFactors(C1y2, method="RLE")
C1y2 <- estimateGLMCommonDisp(C1y2, C1design2)
C1y2 <- estimateGLMTagwiseDisp(C1y2, C1design2)
C1fit2 <- glmFit(C1y2, C1design2)
C1lrt2 <- glmLRT(C1fit2,coef=2)

C1de2 <- decideTestsDGE(C1lrt2,lfc=2)
topTags(C1lrt2)
summary(de <- decideTestsDGE(C1lrt2))
C1detags2 <- rownames(C1y2)[as.logical(de)]
plotSmear(C1lrt2, de.tags=C1detags2)
abline(h=c(-1, 1), col="blue")

FACSdesign2 <- model.matrix(~FACSx + W_1, data=FACSpd2)
FACSy2 <- DGEList(counts=counts(FACSset), group=FACSx)
FACSy2 <- calcNormFactors(FACSy2, method="none")
#FACSy2 <- calcNormFactors(FACSy2, method="TMM")
FACSy2 <- estimateGLMCommonDisp(FACSy2, FACSdesign2)
FACSy2 <- estimateGLMTagwiseDisp(FACSy2, FACSdesign2)
FACSfit2 <- glmFit(FACSy2, FACSdesign2)
FACSlrt2 <- glmLRT(FACSfit2,coef=2)

FACSde2 <- decideTestsDGE(FACSlrt2,lfc=2)
topTags(FACSlrt2)
summary(de <- decideTestsDGE(FACSlrt2))
FACSdetags2 <- rownames(FACSy2)[as.logical(de)]
plotSmear(FACSlrt2, de.tags=FACSdetags2)
abline(h=c(-1, 1), col="blue")


###RUVr

C1design <- model.matrix(~C1x,data=pData(C1set))
C1y <- DGEList(counts=counts(C1set), group=C1x)
C1y <- calcNormFactors(C1y, method="none")
C1y <- estimateGLMCommonDisp(C1y, C1design,verbose=FALSE)
C1y <- estimateGLMTagwiseDisp(C1y, C1design)
C1fit <- glmFit(C1y, C1design)
C1res <- residuals(C1fit, type="deviance")
C1set3 <- RUVr(C1set, C1genes, k=1, C1res)

plotRLE(C1set4, outline=FALSE, ylim=c(-4, 4),
        col=colors[C1x],xaxt='n',main="C1 Normalized Data")
plotPCA(C1set4, col=colors[C1x], cex=1.2,main="C1 Normalized Data")


FACSdesign <- model.matrix(~FACSx,data=pData(FACSset))
FACSy <- DGEList(counts=counts(FACSset), group=FACSx)
FACSy <- calcNormFactors(FACSy, method="none")
FACSy <- estimateGLMCommonDisp(FACSy, FACSdesign,verbose=FALSE)
FACSy <- estimateGLMTagwiseDisp(FACSy, FACSdesign)
FACSfit <- glmFit(FACSy, FACSdesign)
FACSres <- residuals(FACSfit, type="deviance")
FACSset3 <- RUVr(FACSset, FACSgenes, k=1, FACSres)

plotRLE(FACSset4, outline=FALSE, ylim=c(-4, 4),
        col=colors[FACSx],xaxt='n',main="FACS Normalized Data")
axis(1,c(37.5,121.5),lab=c("Control","Stressed"))
plotPCA(FACSset4, col=colors[FACSx], cex=1.2,main="FACS Normalized Data")

C1pd3 <- pData(C1set3)
FACSpd3 <- pData(FACSset3)

C1design3 <- model.matrix(~C1x + W_1, data=C1pd3)
C1y3 <- DGEList(counts=counts(C1set), group=C1x)
C1y3 <- calcNormFactors(C1y3, method="none")
#C1y3 <- calcNormFactors(C1y3, method="RLE")
C1y3 <- estimateGLMCommonDisp(C1y3, C1design3)
C1y3 <- estimateGLMTagwiseDisp(C1y3, C1design3)
C1fit3 <- glmFit(C1y3, C1design3)
C1lrt3 <- glmLRT(C1fit3,coef=3)

C1de3 <- decideTestsDGE(C1lrt3,lfc=2)
topTags(C1lrt3)
summary(de <- decideTestsDGE(C1lrt3))
C1detags3 <- rownames(C1y3)[as.logical(de)]
plotSmear(C1lrt3, de.tags=C1detags3)
abline(h=c(-1, 1), col="blue")

FACSdesign3 <- model.matrix(~FACSx + W_1, data=FACSpd3)
FACSy3 <- DGEList(counts=counts(FACSset), group=FACSx)
FACSy3 <- calcNormFactors(FACSy3, method="none")
#FACSy3 <- calcNormFactors(FACSy3, method="TMM")
FACSy3 <- estimateGLMCommonDisp(FACSy3, FACSdesign3)
FACSy3 <- estimateGLMTagwiseDisp(FACSy3, FACSdesign3)
FACSfit3 <- glmFit(FACSy3, FACSdesign3)
FACSlrt3 <- glmLRT(FACSfit3,coef=3)

FACSde3 <- decideTestsDGE(FACSlrt3,lfc=2)
topTags(FACSlrt3)
summary(de <- decideTestsDGE(FACSlrt3))
FACSdetags3 <- rownames(FACSy3)[as.logical(de)]
plotSmear(FACSlrt3, de.tags=FACSdetags3)
abline(h=c(-1, 1), col="blue")


#########Venn Diagram of Differential Expression
DE_C1Matrix <- cbind(C1deno,C1de1,C1de2,C1de3)
colnames(DE_C1Matrix) <- c("NoNorm","RUVg","RUVs","RUVr")
rownames(DE_C1Matrix) <- rownames(C1lrt3$table)
DE_C1Matrix <- abs(DE_C1Matrix)
DE_C1Matrix <- DE_C1Matrix[rowSums(DE_C1Matrix)>=1,]

venn(DE_C1Matrix)

DE_FACSMatrix <- cbind(FACSdeno,FACSde1,FACSde2,FACSde3)
colnames(DE_FACSMatrix) <- c("TMM","RUVg","RUVs","RUVr")
rownames(DE_FACSMatrix) <- rownames(FACSlrt3$table)
DE_FACSMatrix <- abs(DE_FACSMatrix)
DE_FACSMatrix <- DE_FACSMatrix[rowSums(DE_FACSMatrix)>=1,]
colSums(DE_FACSMatrix)
venn(DE_FACSMatrix)

#####Volcano Plot

layout(1:2, heights=c(6, 6))
?layout

Vol_Vals <- C1lrt1$table
Vol_Vals$PValue[Vol_Vals$PValue<1e-20] <- 1e-20
High <- subset(Vol_Vals, PValue<.0005 & abs(logFC)>5)

dim(subset(Vol_Vals, PValue<.05))
dim(subset(Vol_Vals, abs(logFC)>2))
dim(subset(Vol_Vals, PValue<.05 & abs(logFC)>2))

with(Vol_Vals,plot(logFC,-log10(PValue),pch=20,main="RUVg"))
with(subset(Vol_Vals, PValue<.05 ), points(logFC, -log10(PValue), pch=20, col="red"))
with(subset(Vol_Vals, abs(logFC)>2), points(logFC, -log10(PValue), pch=20, col="orange"))
with(subset(Vol_Vals, PValue<.05 & abs(logFC)>2),
     points(logFC, -log10(PValue), pch=20, col="green"))
legend("topleft",legend=c("PValue >= 0.05","logFC >= 2","Both"),
       col=c("red","orange","green"),pch=20,bty="n")

with(subset(Vol_Vals, PValue<.0005 & abs(logFC)>5), textxy(logFC, -log10(PValue), labs=rownames(High), cex=.8))

Vol_Facs <- FACSlrt1$table
Vol_Facs$PValue[Vol_Facs$PValue<1e-10] <- 1e-10
High <- subset(Vol_Facs, PValue<.0005 & abs(logFC)>5)

dim(subset(Vol_Facs, PValue<.05))
dim(subset(Vol_Facs, abs(logFC)>2))
dim(subset(Vol_Facs, PValue<.05 & abs(logFC)>2))

with(Vol_Facs,plot(logFC,-log10(PValue),pch=20,main="RUVg"))
with(subset(Vol_Facs, PValue<.05 ), points(logFC, -log10(PValue), pch=20, col="red"))
with(subset(Vol_Facs, abs(logFC)>2), points(logFC, -log10(PValue), pch=20, col="orange"))
with(subset(Vol_Facs, PValue<.05 & abs(logFC)>2),
     points(logFC, -log10(PValue), pch=20, col="green"))
legend("topleft",legend=c("PValue >= 0.05","logFC >= 2","Both"),
       col=c("red","orange","green"),pch=20,bty="n")


Vol2_Vals <- C1lrt2$table
Vol2_Vals$PValue[Vol2_Vals$PValue<1e-20] <- 1e-20
High2 <- subset(Vol2_Vals, PValue<.0005 & abs(logFC)>5)

dim(subset(Vol2_Vals, PValue<.05))
dim(subset(Vol2_Vals, abs(logFC)>2))
dim(subset(Vol2_Vals, PValue<.05 & abs(logFC)>2))

with(Vol2_Vals,plot(logFC,-log10(PValue),pch=20,main="RUVs"))
with(subset(Vol2_Vals, PValue<.05 ), points(logFC, -log10(PValue), pch=20, col="red"))
with(subset(Vol2_Vals, abs(logFC)>2), points(logFC, -log10(PValue), pch=20, col="orange"))
with(subset(Vol2_Vals, PValue<.05 & abs(logFC)>2),
     points(logFC, -log10(PValue), pch=20, col="green"))
legend("topleft",legend=c("PValue >= 0.05","logFC >= 2","Both"),
       col=c("red","orange","green"),pch=20,bty="n")

with(subset(Vol2_Vals, PValue<.0005 & abs(logFC)>5), textxy(logFC, -log10(PValue), labs=rownames(High2), cex=.8))

Vol_Facs2 <- FACSlrt2$table
Vol_Facs2$PValue[Vol_Facs2$PValue<1e-10] <- 1e-10
Highf2 <- subset(Vol_Facs2, PValue<.0005 & abs(logFC)>5)

dim(subset(Vol_Facs2, PValue<.05))
dim(subset(Vol_Facs2, abs(logFC)>2))
dim(subset(Vol_Facs2, PValue<.05 & abs(logFC)>2))

with(Vol_Facs2,plot(logFC,-log10(PValue),pch=20,main="RUVs"))
with(subset(Vol_Facs2, PValue<.05 ), points(logFC, -log10(PValue), pch=20, col="red"))
with(subset(Vol_Facs2, abs(logFC)>2), points(logFC, -log10(PValue), pch=20, col="orange"))
with(subset(Vol_Facs2, PValue<.05 & abs(logFC)>2),
     points(logFC, -log10(PValue), pch=20, col="green"))
legend("topleft",legend=c("PValue >= 0.05","logFC >= 2","Both"),
       col=c("red","orange","green"),pch=20,bty="n")

Vol3_Vals <- C1lrt3$table
Vol3_Vals$PValue[Vol3_Vals$PValue<1e-10] <- 1e-10
High3 <- subset(Vol3_Vals, PValue<.0005 & abs(logFC)>5)

dim(subset(Vol3_Vals, PValue<.05))
dim(subset(Vol3_Vals, abs(logFC)>2))
dim(subset(Vol3_Vals, PValue<.05 & abs(logFC)>2))

with(Vol3_Vals,plot(logFC,-log10(PValue),pch=20,main="RUVr"))
with(subset(Vol3_Vals, PValue<.05 ), points(logFC, -log10(PValue), pch=20, col="red"))
with(subset(Vol3_Vals, abs(logFC)>2), points(logFC, -log10(PValue), pch=20, col="orange"))
with(subset(Vol3_Vals, PValue<.05 & abs(logFC)>2),
     points(logFC, -log10(PValue), pch=20, col="green"))
legend("topleft",legend=c("PValue >= 0.05","logFC >= 2","Both"),
       col=c("red","orange","green"),pch=20,bty="n")

with(subset(Vol3_Vals, PValue<.0005 & abs(logFC)>5), textxy(logFC, -log10(PValue), labs=rownames(High3), cex=.8))

Vol_Facs3 <- FACSlrt3$table
Vol_Facs3$PValue[Vol_Facs3$PValue<1e-7] <- 1e-7
High <- subset(Vol_Facs3, PValue<.0005 & abs(logFC)>5)

dim(subset(Vol_Facs3, PValue<.05))
dim(subset(Vol_Facs3, abs(logFC)>2))
dim(subset(Vol_Facs3, PValue<.05 & abs(logFC)>2))

with(Vol_Facs3,plot(logFC,-log10(PValue),pch=20,main="RUVr"))
with(subset(Vol_Facs3, PValue<.05 ), points(logFC, -log10(PValue), pch=20, col="red"))
with(subset(Vol_Facs3, abs(logFC)>2), points(logFC, -log10(PValue), pch=20, col="orange"))
with(subset(Vol_Facs3, PValue<.05 & abs(logFC)>2),
     points(logFC, -log10(PValue), pch=20, col="green"))
legend("topleft",legend=c("PValue >= 0.05","logFC >= 2","Both"),
       col=c("red","orange","green"),pch=20,bty="n")

#####Top Differentially Expressed Genes
C1RLE <- rownames(topTags(C1lrtno))
C1RUVg <- rownames(topTags(C1lrt1))
C1RUVs <- rownames(topTags(C1lrt2))
C1RUVr <- rownames(topTags(C1lrt3))
C1TopGenes <- cbind(C1RLE,C1RUVg,C1RUVs,C1RUVr)
rownames(topTags(C1lrt3,n=100))

rownames(subset(Vol_Vals, PValue<.05 & abs(logFC)>2))


FACSRLE <- rownames(topTags(FACSlrtno))
FACSRUVg <- rownames(topTags(FACSlrt1))
FACSRUVs <- rownames(topTags(FACSlrt2))
FACSRUVr <- rownames(topTags(FACSlrt3))
FACSTopGenes <- cbind(FACSRLE,FACSRUVg,FACSRUVs,FACSRUVr)
FACSTopGenes

####Heatmap of topTags

FC_Across_Methods <- cbind(C1lrtno$table$logFC,C1lrt1$table$logFC,
                           C1lrt2$table$logFC,C1lrt3$table$logFC)
rownames(FC_Across_Methods) <- rownames(C1lrt1$table)
colnames(FC_Across_Methods) <- c("RLE","RUVg","RUVs","RUVr")

FC_Across_MethodsFACS <- cbind(FACSlrtno$table$logFC,FACSlrt1$table$logFC,
                           FACSlrt2$table$logFC,FACSlrt3$table$logFC)
rownames(FC_Across_MethodsFACS) <- rownames(FACSlrt1$table)
colnames(FC_Across_MethodsFACS) <- c("RLE","RUVg","RUVs","RUVr")

par(mar=c(5,5,3,2))
heatmap.2(FC_Across_Methods,col=my_palette,breaks=color_ramp,dendrogram="column",labRow="",symkey=FALSE)
color_ramp= c(seq(-10,-3,length=5),seq(-2.9,3,length=10),seq(3.1,10,length=5))
my_palette <- colorRampPalette(c("red","orange","yellow"))(n = 19)
dev.off()
