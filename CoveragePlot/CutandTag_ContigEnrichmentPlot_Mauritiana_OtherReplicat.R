library(karyoploteR)
library(regioneR)
library(GenomicRanges) 
library(rtracklayer) 
library(IRanges) 
library(devtools)
library(stringr)


# Building a Genome
mygenome= read.table("Reference/dmau_scaffold2_V2.chrom.sizes",comment.char = "")
mygenome2 = data.frame(mygenome[,1],rep(1,dim(mygenome)[1]),mygenome[,2])
custom.genome <- toGRanges(mygenome2)

# Building a Cytoband
cytobands = read.csv("Reference/dmau_scaffold2_V2.Mau.042623.cytoband.mod", sep="\t", header=T)
custom.cytobands <- toGRanges(cytobands)

# Colors Cytoband
color_legend=read.csv("Reference/name-colorv3.csv", sep=";", header=F)
vv=as.character(color_legend$V2)
names(vv)=color_legend$V1

# Peaks
peaks2q0<-read.table("CUT&Tag/Macs2_Peaks/Maur_CID2_bwa_q0_macs2_peaks.narrowPeak",header=F, sep="\t", comment.char = "")
peaks2q0<- toGRanges(peaks2q0)
peaks2q30<-read.table("CUT&Tag/Macs2_Peaks/Maur_CID2_bwa_q30_macs2_peaks.narrowPeak",header=F, sep="\t", comment.char = "")
peaks2q30<- toGRanges(peaks2q30)

peaks3q0<-read.table("CUT&Tag/Macs2_Peaks/Maur_CID3_bwa_q0_macs2_peaks.narrowPeak",header=F, sep="\t", comment.char = "")
peaks3q0<- toGRanges(peaks3q0)
peaks3q30<-read.table("CUT&Tag/Macs2_Peaks/Maur_CID3_bwa_q30_macs2_peaks.narrowPeak",header=F, sep="\t", comment.char = "")
peaks3q30<- toGRanges(peaks3q30)

data.bg <- import("Dmau_k150.bedgraph", format="bedGraph")

#########  Cen 2/3 (Contig6) ##########
chr=c("Contig6")
detail.region <- toGRanges(data.frame("Contig6",250000,450000))

pdf("CUT&Tag/Figures/Mau_Contig6_otherReplicateV2.pdf")
pp <- getDefaultPlotParams(plot.type=1)
pp$leftmargin <- 0.15
pp$topmargin <- 5
pp$bottommargin <- 15
pp$ideogramheight <- 5
pp$data1inmargin <- 10

kp <- plotKaryotype(genome = custom.genome, cytobands = custom.cytobands, ideogram.plotter = NULL,label=NULL,cex=0.5, plot.params = pp, chromosomes=chr, zoom=detail.region)

kpAddBaseNumbers(kp,tick.dist=25000,tick.len=6, cex=0.7, digits=3, units="kb" )
kpAddCytobandsAsLine(kp,color.table=vv,lwd=20)

kp<- kpPlotBigWig(kp, data.panel = 1, r0=0.80, r1=1 , data="CUT&Tag/Coverage/Maur_CID2_bwa_q0_RPM.bw", ymax=200, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data.panel = 1, r0=0.80, r1=1, data="CUT&Tag/Coverage/Maur_CID2_bwa_q30_RPM.bw", ymax=200, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=200, r0=0.80, r1=1 , numticks = 2,cex=0.7)
kpAddLabels(kp, labels = "Cenp-A enrichment (RPM)", r0=0.80, r1=1, srt=90, pos=3, cex=0.5, label.margin = 0.035)
kpPlotRegions(kp, data=peaks2q0, data.panel = 1, r0=0.76,r1=0.78, avoid.overlapping=FALSE,col="grey",border=NA)
kpPlotRegions(kp, data=peaks2q30, data.panel = 1, r0=0.76,r1=0.78, avoid.overlapping=FALSE,col="black",border=NA)

kp<- kpPlotBigWig(kp, data.panel = 1,  r0=0.53, r1=0.73 , data="CUT&Tag/Coverage/Maur_CID3_bwa_q0_RPM.bw", ymax=200, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data.panel = 1,  r0=0.53, r1=0.73 , data="CUT&Tag/Coverage/Maur_CID3_bwa_q30_RPM.bw", ymax=200, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=200,  r0=0.53, r1=0.73 , numticks = 2,cex=0.7)
kpAddLabels(kp, labels = "Cenp-A enrichment (RPM)",  r0=0.53, r1=0.73, srt=90, pos=3, cex=0.5, label.margin = 0.035)
kpPlotRegions(kp, data=peaks3q0, data.panel = 1, r0=0.49,r1=0.51, avoid.overlapping=FALSE,col="grey",border=NA)
kpPlotRegions(kp, data=peaks3q30, data.panel = 1,  r0=0.49,r1=0.51, avoid.overlapping=FALSE,col="black",border=NA)

kp<- kpPlotBigWig(kp, data.panel = 1,  r0=0.26, r1=0.46, data="CUT&Tag/Coverage/Maur_IgG_bwa_q0_RPM.bw", ymax=200, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data.panel = 1,  r0=0.26, r1=0.46, data="CUT&Tag/Coverage/Maur_IgG_bwa_q30_RPM.bw", ymax=200, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=200, r0=0.26, r1=0.46, numticks = 2,cex=0.7)
kpAddLabels(kp, labels = "IgG enrichment (RPM)", srt=90, pos=3, cex=0.5,r0=0.26, r1=0.46, label.margin = 0.035)

kpArea(kp,data.panel = 1, data=data.bg, y=data.bg$score, ymin=0, ymax=1, r0=0, r1=0.19,col="black", border=NA)
kpAddLabels(kp, labels = "Mappability", r0=0, r1=0.19, srt=90, pos=3, cex=0.5, label.margin = 0.035)
kpAxis(kp, ymin=0 , ymax=1, r0=0, r1=0.19 , numticks = 2,cex=0.7)

dev.off()




#########  Cen 2/3 (Contig51) ##########
chr=c("Contig51")
detail.region <- toGRanges(data.frame("Contig51",1,200000))


pdf("CUT&Tag/Figures/Mau_Contig51_otherReplicateV2.pdf")
pp <- getDefaultPlotParams(plot.type=1)
pp$leftmargin <- 0.15
pp$topmargin <- 5
pp$bottommargin <- 15
pp$ideogramheight <- 5
pp$data1inmargin <- 10



kp <- plotKaryotype(genome = custom.genome, cytobands = custom.cytobands, ideogram.plotter = NULL,label=NULL,cex=0.5, plot.params = pp, chromosomes=chr, zoom=detail.region)

kpAddBaseNumbers(kp,tick.dist=20000,tick.len=6, cex=0.7, digits=3, units="kb" )
kpAddCytobandsAsLine(kp,color.table=vv,lwd=20)

kp<- kpPlotBigWig(kp, data.panel = 1, r0=0.80, r1=1 , data="CUT&Tag/Coverage/Maur_CID2_bwa_q0_RPM.bw", ymax=300, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data.panel = 1, r0=0.80, r1=1 , data="CUT&Tag/Coverage/Maur_CID2_bwa_q30_RPM.bw", ymax=300, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=300, r0=0.80, r1=1 , numticks = 2,cex=0.7)
kpAddLabels(kp, labels = "Cenp-A enrichment (RPM)", r0=0.80, r1=1, srt=90, pos=3, cex=0.5, label.margin = 0.035)
kpPlotRegions(kp, data=peaks2q0, data.panel = 1, r0=0.76,r1=0.78 ,avoid.overlapping=FALSE,col="grey",border=NA)
kpPlotRegions(kp, data=peaks2q30, data.panel = 1,  r0=0.76,r1=0.78,avoid.overlapping=FALSE,col="black",border=NA)

kp<- kpPlotBigWig(kp, data.panel = 1, r0=0.53, r1=0.73 , data="CUT&Tag/Coverage/Maur_CID3_bwa_q0_RPM.bw", ymax=300, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data.panel = 1, r0=0.53, r1=0.73 , data="CUT&Tag/Coverage/Maur_CID3_bwa_q30_RPM.bw", ymax=300, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=300, r0=0.53, r1=0.73 , numticks = 2,cex=0.7)
kpAddLabels(kp, labels = "Cenp-A enrichment (RPM)", r0=0.53, r1=0.73, srt=90, pos=3, cex=0.5, label.margin = 0.035)
kpPlotRegions(kp, data=peaks3q0, data.panel = 1, r0=0.49,r1=0.51,avoid.overlapping=FALSE,col="grey",border=NA)
kpPlotRegions(kp, data=peaks3q30, data.panel = 1, r0=0.49,r1=0.51,avoid.overlapping=FALSE,col="black",border=NA)

kp<- kpPlotBigWig(kp, data.panel = 1,  r0=0.26, r1=0.46, data="CUT&Tag/Coverage/Maur_IgG_bwa_q0_RPM.bw", ymax=300, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data.panel = 1,  r0=0.26, r1=0.46, data="CUT&Tag/Coverage/Maur_IgG_bwa_q30_RPM.bw", ymax=300, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=300,r0=0.26, r1=0.46, numticks = 2,cex=0.7)
kpAddLabels(kp, labels = "IgG enrichment (RPM)", srt=90, pos=3, cex=0.5,r0=0.26, r1=0.46, label.margin = 0.035)

kpArea(kp,data.panel = 1, data=data.bg, y=data.bg$score, ymin=0, ymax=1, r0=0, r1=0.19,col="black", border=NA)
kpAddLabels(kp, labels = "Mappability", r0=0, r1=0.19, srt=90, pos=3, cex=0.5, label.margin = 0.035)
kpAxis(kp, ymin=0 , ymax=1, r0=0, r1=0.19 , numticks = 2,cex=0.7)

dev.off()



#########  Cen X (Contig147) ##########
chr=c("Contig147")
detail.region <- toGRanges(data.frame("Contig147",1,150000))


pdf("CUT&Tag/Figures/Mau_Contig147_otherReplicateV2.pdf")
pp <- getDefaultPlotParams(plot.type=1)
pp$leftmargin <- 0.15
pp$topmargin <- 5
pp$bottommargin <- 15
pp$ideogramheight <- 5
pp$data1inmargin <- 10



kp <- plotKaryotype(genome = custom.genome, cytobands = custom.cytobands, ideogram.plotter = NULL,label=NULL,cex=0.5, plot.params = pp, chromosomes=chr, zoom=detail.region)

kpAddBaseNumbers(kp,tick.dist=10000,tick.len=6, cex=0.7, digits=3, units="kb" )
kpAddCytobandsAsLine(kp,color.table=vv,lwd=20)

kp<- kpPlotBigWig(kp, data.panel = 1, r0=0.80, r1=1 , data="CUT&Tag/Coverage/Maur_CID2_bwa_q0_RPM.bw", ymax=100, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data.panel = 1, r0=0.80, r1=1 , data="CUT&Tag/Coverage/Maur_CID2_bwa_q30_RPM.bw", ymax=100, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=100, r0=0.80, r1=1 , numticks = 2,cex=0.7)
kpAddLabels(kp, labels = "Cenp-A enrichment (RPM)", r0=0.80, r1=1, srt=90, pos=3, cex=0.5, label.margin = 0.035)
kpPlotRegions(kp, data=peaks2q0, data.panel = 1, r0=0.76,r1=0.78 ,avoid.overlapping=FALSE,col="grey",border=NA)
kpPlotRegions(kp, data=peaks2q30, data.panel = 1,r0=0.76,r1=0.78,avoid.overlapping=FALSE,col="black",border=NA)

kp<- kpPlotBigWig(kp, data.panel = 1, r0=0.53, r1=0.73 , data="CUT&Tag/Coverage/Maur_CID3_bwa_q0_RPM.bw", ymax=100, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data.panel = 1, r0=0.53, r1=0.73 , data="CUT&Tag/Coverage/Maur_CID3_bwa_q30_RPM.bw", ymax=100, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=100, r0=0.53, r1=0.73 , numticks = 2,cex=0.7)
kpAddLabels(kp, labels = "Cenp-A enrichment (RPM)", r0=0.53, r1=0.73, srt=90, pos=3, cex=0.5, label.margin = 0.035)
kpPlotRegions(kp, data=peaks3q0, data.panel = 1, r0=0.49,r1=0.51,avoid.overlapping=FALSE,col="grey",border=NA)
kpPlotRegions(kp, data=peaks3q30, data.panel = 1, r0=0.49,r1=0.51,avoid.overlapping=FALSE,col="black",border=NA)

kp<- kpPlotBigWig(kp, data.panel = 1,  r0=0.26, r1=0.46, data="CUT&Tag/Coverage/Maur_IgG_bwa_q0_RPM.bw", ymax=100, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data.panel = 1,  r0=0.26, r1=0.46, data="CUT&Tag/Coverage/Maur_IgG_bwa_q30_RPM.bw", ymax=100, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=100,r0=0.26, r1=0.46, numticks = 2,cex=0.7)
kpAddLabels(kp, labels = "IgG enrichment (RPM)", srt=90, pos=3, cex=0.5, r0=0.26, r1=0.46, label.margin = 0.035)

kpArea(kp,data.panel = 1, data=data.bg, y=data.bg$score, ymin=0, ymax=1, r0=0, r1=0.19,col="black", border=NA)
kpAddLabels(kp, labels = "Mappability", r0=0, r1=0.19, srt=90, pos=3, cex=0.5, label.margin = 0.035)
kpAxis(kp, ymin=0 , ymax=1, r0=0, r1=0.19 , numticks = 2,cex=0.7)

dev.off()



#########  Cen 4 (Contig134) ##########
chr=c("Contig134")

pdf("CUT&Tag/Figures/Mau_Contig134_otherReplicateV2.pdf")
pp <- getDefaultPlotParams(plot.type=1)
pp$leftmargin <- 0.15
pp$topmargin <- 5
pp$bottommargin <- 15
pp$ideogramheight <- 5
pp$data1inmargin <- 10


kp <- plotKaryotype(genome = custom.genome, cytobands = custom.cytobands, ideogram.plotter = NULL,label=NULL,cex=0.5, plot.params = pp, chromosomes=chr)

kpAddBaseNumbers(kp,tick.dist=10000,tick.len=6, cex=0.7, digits=3, units="kb" )
kpAddCytobandsAsLine(kp,color.table=vv,lwd=20)

kp<- kpPlotBigWig(kp, data.panel = 1, r0=0.80, r1=1 , data="CUT&Tag/Coverage/Maur_CID2_bwa_q0_RPM.bw", ymax=200, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data.panel = 1, r0=0.80, r1=1 , data="CUT&Tag/Coverage/Maur_CID2_bwa_q30_RPM.bw", ymax=200, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=200, r0=0.80, r1=1, numticks = 2,cex=0.7)
kpAddLabels(kp, labels = "Cenp-A enrichment (RPM)", r0=0.80, r1=1, srt=90, pos=3, cex=0.5, label.margin = 0.035)
kpPlotRegions(kp, data=peaks2q0, data.panel = 1, r0=0.76,r1=0.78 ,avoid.overlapping=FALSE,col="grey",border=NA)
kpPlotRegions(kp, data=peaks2q30, data.panel = 1,r0=0.76,r1=0.78,avoid.overlapping=FALSE,col="black",border=NA)

kp<- kpPlotBigWig(kp, data.panel = 1, r0=0.53, r1=0.73 , data="CUT&Tag/Coverage/Maur_CID3_bwa_q0_RPM.bw", ymax=200, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data.panel = 1, r0=0.53, r1=0.73 , data="CUT&Tag/Coverage/Maur_CID3_bwa_q30_RPM.bw", ymax=200, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=200, r0=0.53, r1=0.73, numticks = 2,cex=0.7)
kpAddLabels(kp, labels = "Cenp-A enrichment (RPM)", r0=0.53, r1=0.73, srt=90, pos=3, cex=0.5, label.margin = 0.035)
kpPlotRegions(kp, data=peaks3q0, data.panel = 1, r0=0.49,r1=0.51,avoid.overlapping=FALSE,col="grey",border=NA)
kpPlotRegions(kp, data=peaks3q30, data.panel = 1, r0=0.49,r1=0.51,avoid.overlapping=FALSE,col="black",border=NA)

kp<- kpPlotBigWig(kp, data.panel = 1,  r0=0.26, r1=0.46, data="CUT&Tag/Coverage/Maur_IgG_bwa_q0_RPM.bw", ymax=200, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data.panel = 1,  r0=0.26, r1=0.46, data="CUT&Tag/Coverage/Maur_IgG_bwa_q30_RPM.bw", ymax=200, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=200, r0=0.26, r1=0.46, numticks = 2,cex=0.7)
kpAddLabels(kp, labels = "IgG enrichment (RPM)", srt=90, pos=3, cex=0.5, r0=0.26, r1=0.46, label.margin = 0.035)

kpArea(kp,data.panel = 1, data=data.bg, y=data.bg$score, ymin=0, ymax=1, r0=0, r1=0.19,col="black", border=NA)
kpAddLabels(kp, labels = "Mappability", r0=0, r1=0.19, srt=90, pos=3, cex=0.5, label.margin = 0.035)
kpAxis(kp, ymin=0 , ymax=1, r0=0, r1=0.19 , numticks = 2,cex=0.7)

dev.off()


#########  Cen Y (Y_scaffold1) ##########
chr=c("Y_scaffold3")
detail.region <- toGRanges(data.frame("Y_scaffold3",2200000,2700000))

pdf("CUT&Tag/Figures/Mau_Y_scaffold3_otherReplicateV2.pdf")
pp <- getDefaultPlotParams(plot.type=1)
pp$leftmargin <- 0.15
pp$topmargin <- 5
pp$bottommargin <- 15
pp$ideogramheight <- 5
pp$data1inmargin <- 10


kp <- plotKaryotype(genome = custom.genome, cytobands = custom.cytobands, ideogram.plotter = NULL,label=NULL,cex=0.5, plot.params = pp, chromosomes=chr, zoom=detail.region)

kpAddBaseNumbers(kp,tick.dist=100000,tick.len=6, cex=0.7, digits=3, units="mb" )
kpAddCytobandsAsLine(kp,color.table=vv,lwd=20)

kp<- kpPlotBigWig(kp, data.panel = 1, r0=0.80, r1=1 , data="CUT&Tag/Coverage/Maur_CID2_bwa_q0_RPM.bw", ymax=100, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data.panel = 1, r0=0.80, r1=1 , data="CUT&Tag/Coverage/Maur_CID2_bwa_q30_RPM.bw", ymax=100, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=100, r0=0.80, r1=1 , numticks = 2,cex=0.7)
kpAddLabels(kp, labels = "Cenp-A enrichment (RPM)", r0=0.80, r1=1, srt=90, pos=3, cex=0.5, label.margin = 0.035)
kpPlotRegions(kp, data=peaks2q0, data.panel = 1, r0=0.76,r1=0.78 ,avoid.overlapping=FALSE,col="grey",border=NA)
kpPlotRegions(kp, data=peaks2q30, data.panel = 1, r0=0.76,r1=0.78,avoid.overlapping=FALSE,col="black",border=NA)

kp<- kpPlotBigWig(kp, data.panel = 1, r0=0.53, r1=0.73 , data="CUT&Tag/Coverage/Maur_CID3_bwa_q0_RPM.bw", ymax=100, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data.panel = 1, r0=0.53, r1=0.73 , data="CUT&Tag/Coverage/Maur_CID3_bwa_q30_RPM.bw", ymax=100, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=100, r0=0.53, r1=0.73 , numticks = 2,cex=0.7)
kpAddLabels(kp, labels = "Cenp-A enrichment (RPM)", r0=0.53, r1=0.73, srt=90, pos=3, cex=0.5, label.margin = 0.035)
kpPlotRegions(kp, data=peaks3q0, data.panel = 1, r0=0.49,r1=0.51, avoid.overlapping=FALSE,col="grey",border=NA)
kpPlotRegions(kp, data=peaks3q30, data.panel = 1, r0=0.49,r1=0.51, avoid.overlapping=FALSE,col="black",border=NA)

kp<- kpPlotBigWig(kp, data.panel = 1, r0=0.26, r1=0.46, data="CUT&Tag/Coverage/Maur_IgG_bwa_q0_RPM.bw", ymax=100, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data.panel = 1, r0=0.26, r1=0.46, data="CUT&Tag/Coverage/Maur_IgG_bwa_q30_RPM.bw", ymax=100, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=100,  r0=0.26, r1=0.46, numticks = 2,cex=0.7)
kpAddLabels(kp, labels = "IgG enrichment (RPM)", srt=90, pos=3, cex=0.5,  r0=0.26, r1=0.46, label.margin = 0.035)

kpArea(kp,data.panel = 1, data=data.bg, y=data.bg$score, ymin=0, ymax=1, r0=0, r1=0.19,col="black", border=NA)
kpAddLabels(kp, labels = "Mappability", r0=0, r1=0.19, srt=90, pos=3, cex=0.5, label.margin = 0.035)
kpAxis(kp, ymin=0 , ymax=1, r0=0, r1=0.19 , numticks = 2,cex=0.7)

dev.off()

