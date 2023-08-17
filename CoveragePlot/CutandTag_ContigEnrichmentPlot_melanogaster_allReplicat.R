library(karyoploteR)
library(regioneR)
library(GenomicRanges) 
library(rtracklayer) 
library(IRanges) 
library(devtools)
library(stringr)


# Building a Genome
mygenome= read.table("Reference/dmel_scaffold2_plus0310.chrom.sizes",comment.char = "")
mygenome2 = data.frame(mygenome[,1],rep(1,dim(mygenome)[1]),mygenome[,2])
custom.genome <- toGRanges(mygenome2)

# Building a Cytoband
cytobands = read.csv("Reference/dmel_scaffold2_plus0310_Mel_042623.cytoband.mod", sep="\t", header=T)
custom.cytobands <- toGRanges(cytobands)

# Colors Cytoband
color_legend=read.csv("Reference/name-colorv3.csv", sep=";", header=F)
vv=as.character(color_legend$V2)
names(vv)=color_legend$V1

# Peaks
peaks1q0<-read.table("Macs2_Peaks/N25_CID1_bwa_q0_macs2_sort.peaks.narrowPeak",header=F, sep="\t", comment.char = "")
peaks1q0<- toGRanges(peaks1q0)
peaks1q30<-read.table("Macs2_Peaks/N25_CID1_bwa_q30_macs2_sort.peaks.narrowPeak",header=F, sep="\t", comment.char = "")
peaks1q30<- toGRanges(peaks1q30)

peaks2q0<-read.table("Macs2_Peaks/N25_CID2_bwa_q0_macs2_sort.peaks.narrowPeak",header=F, sep="\t", comment.char = "")
peaks2q0<- toGRanges(peaks2q0)
peaks2q30<-read.table("Macs2_Peaks/N25_CID2_bwa_q30_macs2_sort.peaks.narrowPeak",header=F, sep="\t", comment.char = "")
peaks2q30<- toGRanges(peaks2q30)

peaks3q0<-read.table("Macs2_Peaks/N25_CID3_bwa_q0_macs2_sort.peaks.narrowPeak",header=F, sep="\t", comment.char = "")
peaks3q0<- toGRanges(peaks3q0)
peaks3q30<-read.table("Macs2_Peaks/N25_CID3_bwa_q30_macs2_sort.peaks.narrowPeak",header=F, sep="\t", comment.char = "")
peaks3q30<- toGRanges(peaks3q30)



#########  Cen X (Contig79) ##########
chr=c("Contig79")


pdf("Mel_Contig_79_otherReplicate.pdf")
pp <- getDefaultPlotParams(plot.type=1)
pp$leftmargin <- 0.15
pp$topmargin <- 10
pp$bottommargin <- 15
pp$ideogramheight <- 5
pp$data1inmargin <- 10


kp <- plotKaryotype(genome = custom.genome, cytobands = custom.cytobands, ideogram.plotter = NULL,label=NULL,cex=0.5, plot.params = pp, chromosomes=chr)

kpAddBaseNumbers(kp,tick.dist=10000,tick.len=6, cex=0.8, digits=3, units="kb" )
kpAddCytobandsAsLine(kp,color.table=vv,lwd=15)

YY=1600

kp<- kpPlotBigWig(kp, data.panel = 1, r0=0.80, r1=1 , data="Coverage/N25_CID1_bwa_q0_RPM.bw", ymax=YY, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data.panel = 1, r0=0.80, r1=1, data="Coverage/N25_CID1_bwa_q30_RPM.bw", ymax=YY, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=YY, r0=0.80, r1=1 , numticks = 2,cex=0.6)
kpAddLabels(kp, labels = "Cenp-A enrichment (RPM)", r0=0.80, r1=1, srt=90, pos=3, cex=0.4, label.margin = 0.035)
kpPlotRegions(kp, data=peaks1q0, data.panel = 1,  r0=0.76,r1=0.78 ,avoid.overlapping=FALSE,col="grey",border=NA)
kpPlotRegions(kp, data=peaks1q30, data.panel = 1, r0=0.76,r1=0.78,avoid.overlapping=FALSE,col="black",border=NA)

kp<- kpPlotBigWig(kp, data.panel = 1, r0=0.53, r1=0.73, data="Coverage/N25_CID2_bwa_q0_RPM.bw", ymax=YY, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data.panel = 1, r0=0.53, r1=0.73 , data="Coverage/N25_CID2_bwa_q30_RPM.bw", ymax=YY, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=YY, r0=0.53, r1=0.73 , numticks = 2,cex=0.6)
kpAddLabels(kp, labels = "Cenp-A enrichment (RPM)", r0=0.53, r1=0.73, srt=90, pos=3, cex=0.4, label.margin = 0.035)
kpPlotRegions(kp, data=peaks2q0, data.panel = 1, r0=0.49,r1=0.51 ,avoid.overlapping=FALSE,col="grey",border=NA)
kpPlotRegions(kp, data=peaks2q30, data.panel = 1, r0=0.49,r1=0.51,avoid.overlapping=FALSE,col="black",border=NA)

kp<- kpPlotBigWig(kp, data.panel = 1, r0=0.26, r1=0.46 , data="Coverage/N25_CID3_bwa_q0_RPM.bw", ymax=YY, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data.panel = 1, r0=0.26, r1=0.46 , data="Coverage/N25_CID3_bwa_q30_RPM.bw", ymax=YY, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=YY, r0=0.26, r1=0.46 , numticks = 2,cex=0.6)
kpAddLabels(kp, labels = "Cenp-A enrichment (RPM)", r0=0.26, r1=0.46, srt=90, pos=3, cex=0.4, label.margin = 0.035)
kpPlotRegions(kp, data=peaks3q0, data.panel = 1, r0=0.22,r1=0.24,avoid.overlapping=FALSE,col="grey",border=NA)
kpPlotRegions(kp, data=peaks3q30, data.panel = 1, r0=0.22,r1=0.24,avoid.overlapping=FALSE,col="black",border=NA)

kp<- kpPlotBigWig(kp, data.panel = 1,  r0=0, r1=0.19, data="Coverage/N25_IgG_bwa_q0_RPM.bw", ymax=3500, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data.panel = 1,  r0=0, r1=0.19, data="Coverage/N25_IgG_bwa_q30_RPM.bw", ymax=3500, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=3500, r0=0, r1=0.19, numticks = 2,cex=0.6)
kpAddLabels(kp, labels = "Cenp-A enrichment (RPM)", srt=90, pos=3, cex=0.4, r0=0, r1=0.19, label.margin = 0.035)

dev.off()





#########  Cen X (Contig119) ##########
chr=c("Contig119")


pdf("Mel_Contig_119_otherReplicate.pdf")
pp <- getDefaultPlotParams(plot.type=1)
pp$leftmargin <- 0.15
pp$topmargin <- 10
pp$bottommargin <- 15
pp$ideogramheight <- 5
pp$data1inmargin <- 10


kp <- plotKaryotype(genome = custom.genome, cytobands = custom.cytobands, ideogram.plotter = NULL,label=NULL,cex=0.5, plot.params = pp, chromosomes=chr)

kpAddBaseNumbers(kp,tick.dist=10000,tick.len=6, cex=0.8, digits=3, units="kb" )
kpAddCytobandsAsLine(kp,color.table=vv,lwd=15)

YY=3000

kp<- kpPlotBigWig(kp, data.panel = 1, r0=0.80, r1=1 , data="Coverage/N25_CID1_bwa_q0_RPM.bw", ymax=YY, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data.panel = 1, r0=0.80, r1=1, data="Coverage/N25_CID1_bwa_q30_RPM.bw", ymax=YY, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=YY, r0=0.80, r1=1 , numticks = 2,cex=0.6)
kpAddLabels(kp, labels = "Cenp-A enrichment (RPM)", r0=0.80, r1=1, srt=90, pos=3, cex=0.4, label.margin = 0.035)
kpPlotRegions(kp, data=peaks1q0, data.panel = 1,  r0=0.76,r1=0.78 ,avoid.overlapping=FALSE,col="grey",border=NA)
kpPlotRegions(kp, data=peaks1q30, data.panel = 1, r0=0.76,r1=0.78,avoid.overlapping=FALSE,col="black",border=NA)

kp<- kpPlotBigWig(kp, data.panel = 1, r0=0.53, r1=0.73, data="Coverage/N25_CID2_bwa_q0_RPM.bw", ymax=YY, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data.panel = 1, r0=0.53, r1=0.73 , data="Coverage/N25_CID2_bwa_q30_RPM.bw", ymax=YY, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=YY, r0=0.53, r1=0.73 , numticks = 2,cex=0.6)
kpAddLabels(kp, labels = "Cenp-A enrichment (RPM)", r0=0.53, r1=0.73, srt=90, pos=3, cex=0.4, label.margin = 0.035)
kpPlotRegions(kp, data=peaks2q0, data.panel = 1, r0=0.49,r1=0.51 ,avoid.overlapping=FALSE,col="grey",border=NA)
kpPlotRegions(kp, data=peaks2q30, data.panel = 1, r0=0.49,r1=0.51,avoid.overlapping=FALSE,col="black",border=NA)

kp<- kpPlotBigWig(kp, data.panel = 1, r0=0.26, r1=0.46 , data="Coverage/N25_CID3_bwa_q0_RPM.bw", ymax=YY, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data.panel = 1, r0=0.26, r1=0.46 , data="Coverage/N25_CID3_bwa_q30_RPM.bw", ymax=YY, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=YY, r0=0.26, r1=0.46 , numticks = 2,cex=0.6)
kpAddLabels(kp, labels = "Cenp-A enrichment (RPM)", r0=0.26, r1=0.46, srt=90, pos=3, cex=0.4, label.margin = 0.035)
kpPlotRegions(kp, data=peaks3q0, data.panel = 1, r0=0.22,r1=0.24,avoid.overlapping=FALSE,col="grey",border=NA)
kpPlotRegions(kp, data=peaks3q30, data.panel = 1, r0=0.22,r1=0.24,avoid.overlapping=FALSE,col="black",border=NA)

kp<- kpPlotBigWig(kp, data.panel = 1,  r0=0, r1=0.19, data="Coverage/N25_IgG_bwa_q0_RPM.bw", ymax=3500, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data.panel = 1,  r0=0, r1=0.19, data="Coverage/N25_IgG_bwa_q30_RPM.bw", ymax=3500, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=3000, r0=0, r1=0.19, numticks = 2,cex=0.6)
kpAddLabels(kp, labels = "Cenp-A enrichment (RPM)", srt=90, pos=3, cex=0.4, r0=0, r1=0.19, label.margin = 0.035)

dev.off()



#########  Cen 3 (3R_5) ##########
chr=c("3R_5")


pdf("Mel_Contig_3R_5_otherReplicate.pdf")
pp <- getDefaultPlotParams(plot.type=1)
pp$leftmargin <- 0.15
pp$topmargin <- 10
pp$bottommargin <- 15
pp$ideogramheight <- 5
pp$data1inmargin <- 10


kp <- plotKaryotype(genome = custom.genome, cytobands = custom.cytobands, ideogram.plotter = NULL,label=NULL,cex=0.5, plot.params = pp, chromosomes=chr)

kpAddBaseNumbers(kp,tick.dist=10000,tick.len=6, cex=0.8, digits=3, units="kb" )
kpAddCytobandsAsLine(kp,color.table=vv,lwd=15)

YY=1500

kp<- kpPlotBigWig(kp, data.panel = 1, r0=0.80, r1=1 , data="Coverage/N25_CID1_bwa_q0_RPM.bw", ymax=YY, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data.panel = 1, r0=0.80, r1=1, data="Coverage/N25_CID1_bwa_q30_RPM.bw", ymax=YY, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=YY, r0=0.80, r1=1 , numticks = 2,cex=0.6)
kpAddLabels(kp, labels = "Cenp-A enrichment (RPM)", r0=0.80, r1=1, srt=90, pos=3, cex=0.4, label.margin = 0.035)
kpPlotRegions(kp, data=peaks1q0, data.panel = 1,  r0=0.76,r1=0.78 ,avoid.overlapping=FALSE,col="grey",border=NA)
kpPlotRegions(kp, data=peaks1q30, data.panel = 1, r0=0.76,r1=0.78,avoid.overlapping=FALSE,col="black",border=NA)

kp<- kpPlotBigWig(kp, data.panel = 1, r0=0.53, r1=0.73, data="Coverage/N25_CID2_bwa_q0_RPM.bw", ymax=YY, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data.panel = 1, r0=0.53, r1=0.73 , data="Coverage/N25_CID2_bwa_q30_RPM.bw", ymax=YY, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=YY, r0=0.53, r1=0.73 , numticks = 2,cex=0.6)
kpAddLabels(kp, labels = "Cenp-A enrichment (RPM)", r0=0.53, r1=0.73, srt=90, pos=3, cex=0.4, label.margin = 0.035)
kpPlotRegions(kp, data=peaks2q0, data.panel = 1, r0=0.49,r1=0.51 ,avoid.overlapping=FALSE,col="grey",border=NA)
kpPlotRegions(kp, data=peaks2q30, data.panel = 1, r0=0.49,r1=0.51,avoid.overlapping=FALSE,col="black",border=NA)

kp<- kpPlotBigWig(kp, data.panel = 1, r0=0.26, r1=0.46 , data="Coverage/N25_CID3_bwa_q0_RPM.bw", ymax=YY, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data.panel = 1, r0=0.26, r1=0.46 , data="Coverage/N25_CID3_bwa_q30_RPM.bw", ymax=YY, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=YY, r0=0.26, r1=0.46 , numticks = 2,cex=0.6)
kpAddLabels(kp, labels = "Cenp-A enrichment (RPM)", r0=0.26, r1=0.46, srt=90, pos=3, cex=0.4, label.margin = 0.035)
kpPlotRegions(kp, data=peaks3q0, data.panel = 1, r0=0.22,r1=0.24,avoid.overlapping=FALSE,col="grey",border=NA)
kpPlotRegions(kp, data=peaks3q30, data.panel = 1, r0=0.22,r1=0.24,avoid.overlapping=FALSE,col="black",border=NA)

kp<- kpPlotBigWig(kp, data.panel = 1,  r0=0, r1=0.19, data="Coverage/N25_IgG_bwa_q0_RPM.bw", ymax=3500, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data.panel = 1,  r0=0, r1=0.19, data="Coverage/N25_IgG_bwa_q30_RPM.bw", ymax=3500, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=1500, r0=0, r1=0.19, numticks = 2,cex=0.6)
kpAddLabels(kp, labels = "Cenp-A enrichment (RPM)", srt=90, pos=3, cex=0.4, r0=0, r1=0.19, label.margin = 0.035)

dev.off()



#########  Cen 2 (tig00057289) ##########
chr=c("tig00057289")

pdf("Mel_Contig_tig00057289_otherReplicate.pdf")
pp <- getDefaultPlotParams(plot.type=1)
pp$leftmargin <- 0.15
pp$topmargin <- 10
pp$bottommargin <- 15
pp$ideogramheight <- 5
pp$data1inmargin <- 10


kp <- plotKaryotype(genome = custom.genome, cytobands = custom.cytobands, ideogram.plotter = NULL,label=NULL,cex=0.5, plot.params = pp, chromosomes=chr)

kpAddBaseNumbers(kp,tick.dist=5000,tick.len=6, cex=0.8, digits=3, units="kb" )
kpAddCytobandsAsLine(kp,color.table=vv,lwd=15)

YY=1000

kp<- kpPlotBigWig(kp, data.panel = 1, r0=0.80, r1=1 , data="Coverage/N25_CID1_bwa_q0_RPM.bw", ymax=YY, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data.panel = 1, r0=0.80, r1=1, data="Coverage/N25_CID1_bwa_q30_RPM.bw", ymax=YY, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=YY, r0=0.80, r1=1 , numticks = 2,cex=0.6)
kpAddLabels(kp, labels = "Cenp-A enrichment (RPM)", r0=0.80, r1=1, srt=90, pos=3, cex=0.4, label.margin = 0.035)
kpPlotRegions(kp, data=peaks1q0, data.panel = 1,  r0=0.76,r1=0.78 ,avoid.overlapping=FALSE,col="grey",border=NA)
kpPlotRegions(kp, data=peaks1q30, data.panel = 1, r0=0.76,r1=0.78,avoid.overlapping=FALSE,col="black",border=NA)

kp<- kpPlotBigWig(kp, data.panel = 1, r0=0.53, r1=0.73, data="Coverage/N25_CID2_bwa_q0_RPM.bw", ymax=YY, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data.panel = 1, r0=0.53, r1=0.73 , data="Coverage/N25_CID2_bwa_q30_RPM.bw", ymax=YY, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=YY, r0=0.53, r1=0.73 , numticks = 2,cex=0.6)
kpAddLabels(kp, labels = "Cenp-A enrichment (RPM)", r0=0.53, r1=0.73, srt=90, pos=3, cex=0.4, label.margin = 0.035)
kpPlotRegions(kp, data=peaks2q0, data.panel = 1, r0=0.49,r1=0.51 ,avoid.overlapping=FALSE,col="grey",border=NA)
kpPlotRegions(kp, data=peaks2q30, data.panel = 1, r0=0.49,r1=0.51,avoid.overlapping=FALSE,col="black",border=NA)

kp<- kpPlotBigWig(kp, data.panel = 1, r0=0.26, r1=0.46 , data="Coverage/N25_CID3_bwa_q0_RPM.bw", ymax=YY, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data.panel = 1, r0=0.26, r1=0.46 , data="Coverage/N25_CID3_bwa_q30_RPM.bw", ymax=YY, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=YY, r0=0.26, r1=0.46 , numticks = 2,cex=0.6)
kpAddLabels(kp, labels = "Cenp-A enrichment (RPM)", r0=0.26, r1=0.46, srt=90, pos=3, cex=0.4, label.margin = 0.035)
kpPlotRegions(kp, data=peaks3q0, data.panel = 1, r0=0.22,r1=0.24,avoid.overlapping=FALSE,col="grey",border=NA)
kpPlotRegions(kp, data=peaks3q30, data.panel = 1, r0=0.22,r1=0.24,avoid.overlapping=FALSE,col="black",border=NA)

kp<- kpPlotBigWig(kp, data.panel = 1,  r0=0, r1=0.19, data="Coverage/N25_IgG_bwa_q0_RPM.bw", ymax=3500, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data.panel = 1,  r0=0, r1=0.19, data="Coverage/N25_IgG_bwa_q30_RPM.bw", ymax=3500, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=1000, r0=0, r1=0.19, numticks = 2,cex=0.6)
kpAddLabels(kp, labels = "Cenp-A enrichment (RPM)", srt=90, pos=3, cex=0.4, r0=0, r1=0.19, label.margin = 0.035)

dev.off()



#########  Cen Y ##########
chr=c("Y_Contig26")


pdf("Mel_Y_Contig26_otherReplicate.pdf")
pp <- getDefaultPlotParams(plot.type=1)
pp$leftmargin <- 0.15
pp$topmargin <- 10
pp$bottommargin <- 15
pp$ideogramheight <- 5
pp$data1inmargin <- 10


kp <- plotKaryotype(genome = custom.genome, cytobands = custom.cytobands, ideogram.plotter = NULL,label=NULL,cex=0.5, plot.params = pp, chromosomes=chr)

kpAddBaseNumbers(kp,tick.dist=10000,tick.len=6, cex=0.8, digits=3, units="kb" )
kpAddCytobandsAsLine(kp,color.table=vv,lwd=15)

YY=800

kp<- kpPlotBigWig(kp, data.panel = 1, r0=0.80, r1=1 , data="Coverage/N25_CID1_bwa_q0_RPM.bw", ymax=YY, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data.panel = 1, r0=0.80, r1=1, data="Coverage/N25_CID1_bwa_q30_RPM.bw", ymax=YY, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=YY, r0=0.80, r1=1 , numticks = 2,cex=0.6)
kpAddLabels(kp, labels = "Cenp-A enrichment (RPM)", r0=0.80, r1=1, srt=90, pos=3, cex=0.4, label.margin = 0.035)
kpPlotRegions(kp, data=peaks1q0, data.panel = 1,  r0=0.76,r1=0.78 ,avoid.overlapping=FALSE,col="grey",border=NA)
kpPlotRegions(kp, data=peaks1q30, data.panel = 1, r0=0.76,r1=0.78,avoid.overlapping=FALSE,col="black",border=NA)

kp<- kpPlotBigWig(kp, data.panel = 1, r0=0.53, r1=0.73, data="Coverage/N25_CID2_bwa_q0_RPM.bw", ymax=YY, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data.panel = 1, r0=0.53, r1=0.73 , data="Coverage/N25_CID2_bwa_q30_RPM.bw", ymax=YY, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=YY, r0=0.53, r1=0.73 , numticks = 2,cex=0.6)
kpAddLabels(kp, labels = "Cenp-A enrichment (RPM)", r0=0.53, r1=0.73, srt=90, pos=3, cex=0.4, label.margin = 0.035)
kpPlotRegions(kp, data=peaks2q0, data.panel = 1, r0=0.49,r1=0.51 ,avoid.overlapping=FALSE,col="grey",border=NA)
kpPlotRegions(kp, data=peaks2q30, data.panel = 1, r0=0.49,r1=0.51,avoid.overlapping=FALSE,col="black",border=NA)

kp<- kpPlotBigWig(kp, data.panel = 1, r0=0.26, r1=0.46 , data="Coverage/N25_CID3_bwa_q0_RPM.bw", ymax=YY, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data.panel = 1, r0=0.26, r1=0.46 , data="Coverage/N25_CID3_bwa_q30_RPM.bw", ymax=YY, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=YY, r0=0.26, r1=0.46 , numticks = 2,cex=0.6)
kpAddLabels(kp, labels = "Cenp-A enrichment (RPM)", r0=0.26, r1=0.46, srt=90, pos=3, cex=0.4, label.margin = 0.035)
kpPlotRegions(kp, data=peaks3q0, data.panel = 1, r0=0.22,r1=0.24,avoid.overlapping=FALSE,col="grey",border=NA)
kpPlotRegions(kp, data=peaks3q30, data.panel = 1, r0=0.22,r1=0.24,avoid.overlapping=FALSE,col="black",border=NA)

kp<- kpPlotBigWig(kp, data.panel = 1,  r0=0, r1=0.19, data="Coverage/N25_IgG_bwa_q0_RPM.bw", ymax=3500, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data.panel = 1,  r0=0, r1=0.19, data="Coverage/N25_IgG_bwa_q30_RPM.bw", ymax=3500, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=800, r0=0, r1=0.19, numticks = 2,cex=0.6)
kpAddLabels(kp, labels = "Cenp-A enrichment (RPM)", srt=90, pos=3, cex=0.4, r0=0, r1=0.19, label.margin = 0.035)

dev.off()
