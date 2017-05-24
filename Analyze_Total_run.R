#### Contamination Plot

library(plyr)
metricsFile <- ("./Total_old_cuff_run/RQCQueryOutput.csv")

##read the information into R
Data <- read.table(metricsFile,header=F,sep=",")

##Consolidate and average fastq contmaination
Names <- subset(Data[,1],!duplicated(Data$V1))
Data<- ddply(Data,"V1",numcolwise(sum))
Libraries <- Data[,1]
Data <- Data[,2:7]/3
Data$V1 <- Libraries
Data <- Data[,c(7,1,2,3,4,5,6)]

##set variables for plotting
colors <- c("black","red","blue","orange","green","yellow")
colnames(Data) <- c("name", "Artifacts","JGI.Contam","rRNA","Chloroplast","Mitochondria","E.coli")
Data <- Data[match(Names,Data[,1]),]

##change the sizes of text and windows for different sizes of input
size <- 1+((96-ifelse(nrow(Data)>96, 96, nrow(Data)))/96)
out.width <- ifelse(nrow(Data)>36, 18, 12)

##write RQC plot to a pdf file

pdf(file="RQCplot.pdf",height=8,width=out.width)
lapply(2:ncol(Data),function(i){
    if(i==2){
        layout(1:2, heights=c(1, 6))
        par(mar=c(0,0,0,0))
        plot(0, 0, type="n", ann=FALSE, axes=FALSE)
        legend("center", fill=colors, colnames(Data)[-1], cex=1.3, ncol=3)
        par(mar=c(5,5,3,2))
        plot(1:nrow(Data), Data[,i], ylim=c(0,100),
             main="RQC Statistics", xlab="", ylab="",
             xaxt='n', yaxt='n', cex.main=2, cex=size/2,
             pch=19, col=colors[i-1], type='o')
        axis(2, at=seq(0,100, by=10),
             labels=paste(seq(0,100, by=10), "%"), las=2, cex.axis=1.3)
        axis(1,c(48,144,240,336),
             labels=c("FACS Control","FACS Stressed","C1 Control","C1 Stressed"),cex.axis=1.5)
        abline(v=c(96.5,192.5,288.5),col="black")
        par(new=TRUE)
        plot(0,0,pch="",ylim=c(0,100), main="", xlab="", ylab="",
             xaxt='n', yaxt='n', abline(h=seq(0,100, by=5), col=rgb(0.5,0.5,1,0.1), lwd=2))
    }
    if(i>2){
        par(new=TRUE)
        plot(1:nrow(Data), Data[,i], ylim=c(0,100),
             main="", xlab="", ylab="",
             xaxt='n', yaxt='n', cex.main=2, cex.lab=1.5, cex=size/2,
             pch=19, col=colors[i-1], type='b')
    }
})   

###Select Directories
C1Dirs <- list.files("./Long_insertSize_run/C1", full=TRUE)
FacsDirs <- list.files("./Long_insertSize_run/FACS", full=TRUE)

###Selecting Stressed Readstat files

##C1
C1Stressed <- C1Dirs[grep("Stressed", C1Dirs)]
StressedC1files <- list.files(C1Stressed,full=TRUE)
C1StressedStatFiles <- StressedC1files[grep("read.stats.txt",StressedC1files)]
Stressed_C1_stats <- do.call(cbind,lapply(C1StressedStatFiles, function(x)
    data.frame(read.table(x))))
Stressed_C1_stats <- Stressed_C1_stats[ , -which(names(Stressed_C1_stats) %in% c("V2"))]
colnames(Stressed_C1_stats) <- basename(dirname(C1StressedStatFiles))
rownames(Stressed_C1_stats) <- c("Raw","Filtered","rRNA_Removed")
Stressed_C1_stats <- Stressed_C1_stats/4

##Facs
FacsStressed <- FacsDirs[grep("Stressed", FacsDirs)]
StressedFacsfiles <- list.files(FacsStressed,full=TRUE)
FacsStressedStatFiles <- StressedFacsfiles[grep("read.stats.txt",StressedFacsfiles)]
Stressed_facs_stats <- do.call(cbind,lapply(FacsStressedStatFiles, function(x)
    data.frame(read.table(x))))
Stressed_facs_stats <- Stressed_facs_stats[ , -which(names(Stressed_facs_stats) %in% c("V2"))]
colnames(Stressed_facs_stats) <- basename(dirname(FacsStressedStatFiles))
rownames(Stressed_facs_stats) <- c("Raw","Filtered","rRNA_Removed")
Stressed_facs_stats <- Stressed_facs_stats/4

###Selecting Control Readstat files

##C1
C1Control <- C1Dirs[grep("Control", C1Dirs)]
ControlC1files <- list.files(C1Control,full=TRUE)
C1ControlStatFiles <- ControlC1files[grep("read.stats.txt",ControlC1files)]
Control_C1_stats <- do.call(cbind,lapply(C1ControlStatFiles, function(x) data.frame(read.table(x))))
Control_C1_stats <- Control_C1_stats[ , -which(names(Control_C1_stats) %in% c("V2"))]
colnames(Control_C1_stats) <- basename(dirname(C1ControlStatFiles))
rownames(Control_C1_stats) <- c("Raw","Filtered","rRNA_Removed")
Control_C1_stats <- Control_C1_stats/4

##Facs
FacsControl <- FacsDirs[grep("Control", FacsDirs)]
ControlFacsfiles <- list.files(FacsControl,full=TRUE)
FacsControlStatFiles <- ControlFacsfiles[grep("read.stats.txt",ControlFacsfiles)]
Control_Facs_stats <- do.call(cbind,lapply(FacsControlStatFiles, function(x)
    data.frame(read.table(x))))
Control_Facs_stats <- Control_Facs_stats[ , -which(names(Control_Facs_stats) %in% c("V2"))]
colnames(Control_Facs_stats) <- basename(dirname(FacsControlStatFiles))
rownames(Control_Facs_stats) <- c("Raw","Filtered","rRNA_Removed")
Control_Facs_stats <- Control_Facs_stats/4

####Plotting Read Stats

##C1 Stressed
layout(1:2, heights=c(1, 6))
par(mar=c(0,0,0,0))
plot(0, 0, type="n", ann=FALSE, axes=FALSE)
legend("center", fill=c("red","blue","green"), rownames(Stressed_C1_stats), cex=1.3, ncol=3)
par(mar=c(10,5,3,2))
plot(1:ncol(Stressed_C1_stats),Stressed_C1_stats[1,],main="C1 Stressed Read Stats",
     xaxt='n',xlab="",col="red",pch=20,type='o',
     ylab="Number of Reads",ylim=c(0,2001000))
axis(1, at=seq(1:ncol(Stressed_C1_stats)), labels=colnames(Stressed_C1_stats),las=3,cex.axis=.85)
par(new=TRUE)
plot(1:ncol(Stressed_C1_stats),Stressed_C1_stats[2,],xaxt='n',yaxt='n',
     xlab="",col="blue",pch=20,type='o',ylab="",ylim=c(0,2001000))
par(new=TRUE)
plot(1:ncol(Stressed_C1_stats),Stressed_C1_stats[3,],xaxt='n',yaxt='n',
     xlab="",col="green",pch=20,type='o',ylab="",ylim=c(0,2001000))
par(new=FALSE)

##C1 Control
layout(1:2, heights=c(1, 6))
par(mar=c(0,0,0,0))
plot(0, 0, type="n", ann=FALSE, axes=FALSE)
legend("center", fill=c("red","blue","green"), rownames(Control_C1_stats), cex=1.3, ncol=3)
par(mar=c(10,5,3,2))
plot(1:ncol(Control_C1_stats),Control_C1_stats[1,],main="C1 Control Read Stats",
     xaxt='n',xlab="",col="red",pch=20,type='o',
     ylab="Number of Reads",ylim=c(0,2400000))
axis(1, at=seq(1:ncol(Control_C1_stats)), labels=colnames(Control_C1_stats),las=3,cex.axis=.85)
par(new=TRUE)
plot(1:ncol(Control_C1_stats),Control_C1_stats[2,],xaxt='n',yaxt='n',
     xlab="",col="blue",pch=20,type='o',ylab="",ylim=c(0,2400000))
par(new=TRUE)
plot(1:ncol(Control_C1_stats),Control_C1_stats[3,],xaxt='n',yaxt='n',
     xlab="",col="green",pch=20,type='o',ylab="",ylim=c(0,2400000))
par(new=FALSE)

##Facs Stressed
layout(1:2, heights=c(1, 6))
par(mar=c(0,0,0,0))
plot(0, 0, type="n", ann=FALSE, axes=FALSE)
legend("center", fill=c("red","blue","green"), rownames(Stressed_facs_stats), cex=1.3, ncol=3)
par(mar=c(10,5,3,2))
plot(1:ncol(Stressed_facs_stats),Stressed_facs_stats[1,],main="facs Stressed Read Stats",
     xaxt='n',xlab="",col="red",pch=20,type='o',
     ylab="Number of Reads",ylim=c(0,6000000))
axis(1, at=seq(1:ncol(Stressed_facs_stats)), labels=colnames(Stressed_facs_stats),las=3,cex.axis=.85)
par(new=TRUE)
plot(1:ncol(Stressed_facs_stats),Stressed_facs_stats[2,],xaxt='n',yaxt='n',
     xlab="",col="blue",pch=20,type='o',ylab="",ylim=c(0,6000000))
par(new=TRUE)
plot(1:ncol(Stressed_facs_stats),Stressed_facs_stats[3,],xaxt='n',yaxt='n',
     xlab="",col="green",pch=20,type='o',ylab="",ylim=c(0,6000000))
par(new=FALSE)
colnames(Stressed_facs_stats[,12])
Stressed_facs_stats
plot(1:ncol(Stressed_facs_stats),Stressed_facs_stats[1,])

##Facs Control
layout(1:2, heights=c(1, 6))
par(mar=c(0,0,0,0))
plot(0, 0, type="n", ann=FALSE, axes=FALSE)
legend("center", fill=c("red","blue","green"), rownames(Control_Facs_stats), cex=1.3, ncol=3)
par(mar=c(10,5,3,2))
plot(1:ncol(Control_Facs_stats),Control_Facs_stats[1,],main="Facs Control Read Stats",
     xaxt='n',xlab="",col="red",pch=20,type='o',
     ylab="Number of Reads",ylim=c(0,2400000))
axis(1, at=seq(1:ncol(Control_Facs_stats)), labels=colnames(Control_Facs_stats),las=3,cex.axis=.85)
par(new=TRUE)
plot(1:ncol(Control_Facs_stats),Control_Facs_stats[2,],xaxt='n',yaxt='n',
     xlab="",col="blue",pch=20,type='o',ylab="",ylim=c(0,2400000))
par(new=TRUE)
plot(1:ncol(Control_Facs_stats),Control_Facs_stats[3,],xaxt='n',yaxt='n',
     xlab="",col="green",pch=20,type='o',ylab="",ylim=c(0,2400000))
par(new=FALSE)

###Combined Plot
Total_data <- cbind(Control_C1_stats,Stressed_C1_stats,
                    Control_Facs_stats,Stressed_facs_stats)
Total_facs <- cbind(Control_Facs_stats,Stressed_facs_stats)
Total_C1 <- cbind(Control_C1_stats,Stressed_C1_stats)

Raw <- Total_data[1,]
RealTotal_lib <- subset(t(Rawq)<1000)

t(RealTotal_lib)

Total_data <- Total_data[,-which(names(Total_data) %in% rownames(RealTotal_lib))]


##############################High FPKM Read Stats
Name_High_Lib <- paste(c(High_Lib_names[1:2],High_Lib_names[3:4]),c(rep("C1_Control",2),rep("C1_Stressed",2)),sep="_")
Total_C1t <- t(Total_C1)
High_fpkm_read <- subset(Total_C1t,rownames(Total_C1t) %in% Name_High_Lib)
High_fpkm_read


pdf("Single_Cell_Read_Stats.pdf",width=25.5,height=11)

layout(1:2, heights=c(1, 6))
par(mar=c(0,0,0,0))
plot(0, 0, type="n", ann=FALSE, axes=FALSE)
legend("center", fill=c("red","blue","green"), rownames(Total_data), cex=1.3, ncol=3)
par(mar=c(5,5,3,2))
plot(1:ncol(Total_data),Total_data[1,],main="Single Cell Read Stats",
     xaxt='n',xlab="",col="red",pch=20,type='o',
     ylab="Number of Reads",ylim=c(0,23000000),cex.lab=1.5)
axis(1, c(48,144,240,336),
     labels=c("C1 Control","C1 Stressed","FACS Control","FACS Stressed"),cex.axis=1.5)
par(new=TRUE)
plot(1:ncol(Total_data),Total_data[2,],xaxt='n',yaxt='n',
     xlab="",col="blue",pch=20,type='o',ylab="",ylim=c(0,23000000))
par(new=TRUE)
plot(1:ncol(Total_data),Total_data[3,],xaxt='n',yaxt='n',
     xlab="",col="green",pch=20,type='o',ylab="",ylim=c(0,23000000))
abline(v=c(96.5,192.5,288.5),col="black",lwd=1.5)
par(new=FALSE)

dev.off()

####Alignment Files

flag <- read.table("./Long_insertSize_run/Compiled_flagstat",header=T,sep=",")
rownames(flag)=flag[,1]
flag=flag[,-c(1,ncol(flag)),drop=F]
flagT <- t(flag)

C1StressedAln <- subset(flagT,
                        paste(rownames(flagT),"C1","Stressed",sep="_") %in% basename(C1Stressed))
C1ControlAln <- subset(flagT,
                        paste(rownames(flagT),"C1","Control",sep="_") %in% basename(C1Control))
FacsStressedAln <- subset(flagT,
                        paste(rownames(flagT),"FACS","Stressed",sep="_") %in% basename(FacsStressed))
FacsControlAln <- subset(flagT,
                        paste(rownames(flagT),"FACS","Control",sep="_") %in% basename(FacsControl))

C1Aln <- rbind(C1StressedAln,C1ControlAln)
FacsAln <- rbind(FacsStressedAln,FacsControlAln)

##################High FPKM ALignment Scores############################
High_fpkm_aln_stats <- subset(C1Aln, rownames(C1Aln) %in% High_Lib_names)
High_fpkm_aln_stats

pdf(file="TotalReadFlagStats.pdf",height=11,width=25.5)

layout(1:3, heights=c(1, 7, 7))
par(mar=c(0,0,0,0))
plot(0, 0, type="n", ann=FALSE, axes=FALSE)
legend("center", col=c("black","green","red","blue"),
       c("%GenMapSingletons","%GenMapPP",
                "%TransMapSingletons","%TransMapPP")
      ,lty=c(1,1,1,1),lwd=c(2,2,2,2), cex=1.7, ncol=2)
par(mar=c(5,5,3,2))
grange<-range(0,100)
colors<-c("black","green","red","blue")
plot(C1Aln[,2],type="l",col=colors[1],ylim=grange,xaxt='n',yaxt='n',ann=FALSE,xlab="",lwd=2)
axis(1,c(48,144),lab=c("C1 Stressed","C1 Control"),cex.axis=2)
axis(2,las=1,cex.axis=1.5)
box();
lines(C1Aln[,4],type="l",col=colors[2],lwd=2)
lines(C1Aln[,10],type="l",col=colors[3],lwd=2)
lines(C1Aln[,12],type="l",col=colors[4],lwd=2)
title(main="Read & Mapping Statistics",col.main="darkblue",font.main=4)
title(ylab="%",cex.lab=1.5)
abline(v=96.5,col="black")

par(mar=c(5,5,3,2))
grange<-range(0,100)
colors<-c("black","green","red","blue")
plot(FacsAln[,2],type="l",col=colors[1],ylim=grange,xaxt='n',yaxt='n',ann=FALSE,xlab="",lwd=2)
axis(1,c(48,144),lab=c("FACS Stressed","FACS Control"),cex.axis=2)
axis(2,las=1,cex.axis=1.5)
box();
lines(FacsAln[,4],type="l",col=colors[2],lwd=2)
lines(FacsAln[,10],type="l",col=colors[3],lwd=2)
lines(FacsAln[,12],type="l",col=colors[4],lwd=2)
title(main="Read & Mapping Statistics",col.main="darkblue",font.main=4)
title(ylab="%",cex.lab=1.5)
abline(v=96.5,col="black")

dev.off()

#### BWA Alignment Stats
info <- read.table("./Total_run/C1/Control/info3.txt")
reads <- read.csv("./Total_run/C1/Control/Compiled_mem_readstat", row.names=1)

libs <- as.character(info[,1])
reads <- reads[,colnames(reads) %in% libs]
read.data <- t(reads)
read.data <- data.frame(read.data)
read.data$Trimmed <- read.data$Raw-read.data$AfterTrimming
data <- t(as.matrix(subset(read.data, select=c(Raw, Trimmed))))
data <- data.frame(data)

alns <- read.csv("./Total_run/C1/Control/Compiled_mem_flagstat", row.names=1)
libs <- as.character(info[,1])
alns <- alns[,colnames(alns) %in% libs]
aln.data <- t(alns)
aln.data <- data.frame(aln.data)
aln.data$total <- read.data[match(rownames(aln.data), rownames(read.data)),]$After.rRNAFiltering

inf <- read.table("./Total_run/C1/Control/info3.txt")
my.files <- as.character(inf[,1])
my.files <- paste("./Total_run/C1/Control",my.files, sep="/")
flags <- paste(my.files,"SORT.mem.flagstat",sep=".")
flagsPP <- paste(my.files,"PP.SORT.mem.flagstat",sep=".")
flags <- unlist(flags)
flagsPP <- unlist(flagsPP)

Allvals <- lapply(flags, function(flag){
    lines <- readLines(flag)
    vals <- as.numeric(unlist(lapply(lines, function(x) strsplit(x, " ")[[1]][1])))
    names(vals) <- c("all", "dup", "map", "pairs", "r1", "r2",
                     "proper", "pairMap", "single", "disChrHigh", "disChrLow")
    vals
})
vals <- do.call(rbind, Allvals)
Genomicvals <- data.frame(vals)
my.data <- Genomicvals

AllvalsPP <- lapply(flagsPP, function(flag){
    lines <- readLines(flag)
    vals <- as.numeric(unlist(lapply(lines, function(x) strsplit(x, " ")[[1]][1])))
    names(vals) <- c("all", "dup", "map", "pairs", "r1", "r2",
                     "proper", "pairMap", "single", "disChrHigh", "disChrLow")
    vals
})
PPvals <- do.call(rbind, AllvalsPP)
GenomicPPvals <- data.frame(PPvals)
my.PPdata <- GenomicPPvals

my.data_Percent <- (my.data[,3]/my.data[,1])*100
my.PPdata_Percent<- (my.PPdata[,3]/my.data[,1])*100
Map_PP_data <- cbind(my.data_Percent,my.PPdata_Percent)


pdf(file="BWAReadFlagStats.pdf",height=11,width=25.5)

layout(1:3, heights=c(1, 6, 6))
par(mar=c(0,0,0,0))
plot(0, 0, type="n", ann=FALSE, axes=FALSE)
legend("center", col=c("red","blue"),
       c("%TransMapSingletons","%TransMapPP")
      ,lty=c(1,1),lwd=c(2,2), cex=1.5, ncol=2)
par(mar=c(5,5,3,2))
grange<-range(0,100)
colors<-c("red","blue")
plot(Map_PP_data[,1],type="l",col=colors[1],ylim=grange,xaxt='n',yaxt='n',ann=FALSE,xlab="",ylab="",lwd=2)
axis(1,48,lab=c("C1 Control"),cex.axis=2)
axis(2,las=1,cex.axis=1.5)
box();
lines(Map_PP_data[,2],type="l",col=colors[2],lty=1,lwd=2)
title(main="BWA MEM Mapping Statistics",col.main="darkblue",font.main=4,cex.main=2)
title(ylab="%",cex.lab=1.5)

grange<-range(0,100)
colors<-c("red","blue")
plot(C1ControlAln[,10],type="l",col=colors[1],ylim=grange,xaxt='n',yaxt='n',ann=FALSE,xlab="",lwd=2)
axis(1,c(48),lab=c("C1 Control"),cex.axis=2)
axis(2,las=1,cex.axis=1.5)
box();
lines(C1ControlAln[,12],type="l",col=colors[2],lwd=2)
title(main="Bowtie2 Mapping Statistics",col.main="darkblue",font.main=4,cex.main=2)
title(ylab="%",cex.lab=1.5)


dev.off()

####Boxplot of FPKM values
#New_cuffs <- read.table("./Total_run/Compiled_cufflinks.1",header=TRUE)
#FACS_cuff <- read.table("./Total_old_cuff_run/Compiled_cufflinks.1",header=TRUE)
#C1_cuff <- read.table("./Total_old_cuff_run/Compiled_cufflinks.2",header=TRUE)
#C1_cuff <- read.table("./Compiled_cufflinks.2",header=TRUE)
#C1_cuff <- read.table("./Compiled_cufflinks.2",header=TRUE)
FACS_cuff <- read.table("./Long_insertSize_run/Compiled_cufflinks.1",header=TRUE)
C1_cuff <- read.table("./Long_insertSize_run/Compiled_cufflinks.2",header=TRUE)
C1_fpkm <- C1_cuff[,grep("fpkm",colnames(C1_cuff))]
Facs_fpkm <- FACS_cuff[,grep("fpkm",colnames(FACS_cuff))]
#test_AUOYY <- read.table("./Total_G_run/C1/AUOYY_C1_Control/cufflink_rerun/genes.fpkm_tracking",header=TRUE)
#test_AUOZA <- read.table("./Total_G_run/C1/AUOZA_C1_Control/cufflink_rerun/genes.fpkm_tracking",header=TRUE)

C1_fpkm[is.na(C1_fpkm)] <- 0
Cufflinks_Positive <- C1_fpkm[,]>=1
Tr <- data.frame(Values=colSums(Cufflinks_Positive==TRUE,na.rm=TRUE))
Fal <- data.frame(Values=colSums(Cufflinks_Positive==FALSE,na.rm=TRUE))
Combined <- cbind(Tr,Fal)
colnames(Combined) <- c("True","False")
Combined$Total <- rowSums(Combined)
Combined_Percent <- Combined$True/Combined$Total
df <- data.frame(Combined_Percent)
rownames(df) <- rownames(Combined)
df <- df * 100
High_FPKM <- subset(df,Combined_Percent > 50)
High_FPKM
High_Lib_names <- rownames(High_FPKM)
High_Lib_names <- gsub("\\..*","",High_Lib_names)
High_Lib_names
###Read Stats
High_fpkm_read
###Align Stats
High_fpkm_aln_stats
###FPKM Stats
High_C1_fpkm_vals <- subset(t(C1_fpkm), rownames(t(C1_fpkm)) %in% High_Lib_names)
head(t(High_C1_fpkm_vals))


Cufflinks_Positive_F <- Facs_fpkm[,]>=1
FTr <- data.frame(Values=colSums(Cufflinks_Positive_F==TRUE,na.rm=TRUE))
FFal <- data.frame(Values=colSums(Cufflinks_Positive_F==FALSE,na.rm=TRUE))
CombinedF <- cbind(FTr,FFal)
colnames(CombinedF) <- c("True","False")
CombinedF$Total <- rowSums(CombinedF)
Combined_PercentF <- CombinedF$True/CombinedF$Total
dfF <- data.frame(Combined_PercentF)
rownames(dfF) <- rownames(CombinedF)
dfF <- dfF * 100

layout(1:2, heights=c(6, 6))
par(mar=c(5,5,3,2))
grange<-range(0,100)
barplot(as.matrix(t(df)),ylim=grange,col="cyan",ylab="%",
        main="Percentage of C1 Sample FPKM's 1 and Above",
        las=3,cex.names=.5,xlab="",xaxt='n',yaxt='n',
        cex.lab=1.5,cex.main=1.5,yaxt='n')
axis(1,c(55,170),lab=c("C1 Control","C1 Stressed"),cex.axis=1.5)
axis(2,las=1)
abline(v=113,col="black")
par(mar=c(5,5,3,2))
grange<-range(0,100)
barplot(as.matrix(t(dfF)),ylim=grange,col="cyan",ylab="%",
        main="Percentage of FACS Sample FPKM's 1 and Above",
        las=3,cex.names=.5,xlab="",xaxt='n',yaxt='n',
        cex.lab=1.5,cex.main=1.5)
axis(1,c(55,170),lab=c("FACS Control","FACS Stressed"),cex.axis=1.5)
axis(2,las=1)
abline(v=113,col="black")

########### Spike in Controls
SpikeInNumbers <- read.table("./Spiked_in/SpikeinNumbers.txt")
SpikeInNumbersMem <- read.table("./Spiked_in/SpikeinNumbersCounts.mem.txt")
SpikeInNumbersBow <- read.table("./Spiked_in/SpikeinNumbersCounts.bow.txt")
Lib_names <- read.table("./Spiked_in/SpikeFiles.txt")
files <- list.files("./Spiked_in/Control_Counts",full=TRUE)
Control <- do.call(cbind,lapply(files, function(x) read.table(x)))
SpikeInNumbers <- cbind(SpikeInNumbersBow,SpikeInNumbersMem)

Alignment_Comp <- cbind(SpikeInNumbersBow,SpikeInNumbersMem)
colnames(Alignment_Comp) <- c("Bow","Mem")
Alignment_Comp$Diff <- Alignment_Comp$Mem - Alignment_Comp$Bow
Alignment_Comp

rownames(SpikeInNumbersMem) <- Lib_names$V1
SpikeInNumbersMem <- t(SpikeInNumbersMem)
SpikeInNumbersMem <- (SpikeInNumbersMem/92) * 100
median(SpikeInNumbersBow[193:384])

FACSSI <- SpikeInNumbersMem[,1:192]
C1SI <- SpikeInNumbersMem[,193:384]

layout(1:2, heights=c(6, 6))
par(mar=c(5,5,3,2))
grange<-range(0,100)
par(mar=c(5,5,3,2))
barplot(as.matrix(t(C1SI)),ylim=grange,col="cyan",ylab="%",
        main="Percentage of RNA Spike-Ins in C1 Samples",
        las=3,cex.names=.5,xlab="",xaxt='n',yaxt='n',
        cex.lab=1.5,cex.main=1.5,yaxt='n')
axis(1,c(55,170),lab=c("C1 Control","C1 Stressed"),cex.axis=1.5)
axis(2,las=1)
abline(v=c(113),col="black")
barplot(as.matrix(t(FACSSI)),ylim=grange,col="cyan",ylab="%",
        main="Percentage of RNA Spike-Ins in FACS Samples",
        las=3,cex.names=.5,xlab="",xaxt='n',yaxt='n',
        cex.lab=1.5,cex.main=1.5,yaxt='n')
axis(1,c(55,170),lab=c("FACS Control","FACS Stressed"),cex.axis=1.5)
axis(2,las=1)
abline(v=c(113),col="black")

###FACS median 10.8% and C1 median 21%
