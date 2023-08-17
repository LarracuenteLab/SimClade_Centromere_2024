library(karyoploteR)
library(regioneR)
library(GenomicRanges) 
library(rtracklayer) 
library(IRanges) 
library(devtools)
library(stringr)


# Building a Genome
mygenome= read.table("Reference/dsec_scaffold2_2019.chrom.sizes",comment.char = "")
mygenome2 = data.frame(mygenome[,1],rep(1,dim(mygenome)[1]),mygenome[,2])
custom.genome <- toGRanges(mygenome2)

# Building a Cytoband
cytobands = read.csv("Reference/dsec_scaffold2_2019.Sech.042623.cytoband.mod", sep="\t", header=T)
custom.cytobands <- toGRanges(cytobands)

# Colors Cytoband
color_legend=read.csv("Reference/name-colorv3.csv", sep=";", header=F)
vv=as.character(color_legend$V2)
names(vv)=color_legend$V1

# Peaks
peaks2q0<-read.table("Macs2_Peaks/Sech_CID2_bwa_q0_macs2_peaks.narrowPeak",header=F, sep="\t", comment.char = "")
peaks2q0<- toGRanges(peaks2q0)
peaks2q30<-read.table("Macs2_Peaks/Sech_CID2_bwa_q30_macs2_peaks.narrowPeak",header=F, sep="\t", comment.char = "")
peaks2q30<- toGRanges(peaks2q30)

peaks3q0<-read.table("Macs2_Peaks/Sech_CID3_bwa_q0_macs2_peaks.narrowPeak",header=F, sep="\t", comment.char = "")
peaks3q0<- toGRanges(peaks3q0)
peaks3q30<-read.table("Macs2_Peaks/Sech_CID3_bwa_q30_macs2_peaks.narrowPeak",header=F, sep="\t", comment.char = "")
peaks3q30<- toGRanges(peaks3q30)


#########  Cen 2/3 (Contig2) ##########
chr=c("Contig46")
detail.region <- toGRanges(data.frame("Contig46",1,100000))


pdf("Sech_Contig_46_otherReplicate.pdf")
pp <- getDefaultPlotParams(plot.type=1)
pp$leftmargin <- 0.15
pp$topmargin <- 5
pp$bottommargin <- 15
pp$ideogramheight <- 5
pp$data1inmargin <- 10

kp <- plotKaryotype(genome = custom.genome, cytobands = custom.cytobands, ideogram.plotter = NULL,label=NULL,cex=0.5, plot.params = pp, chromosomes=chr, zoom=detail.region)

kpAddBaseNumbers(kp,tick.dist=10000,tick.len=6, cex=1, digits=3, units="kb" )
kpAddCytobandsAsLine(kp,color.table=vv,lwd=20)

kp<- kpPlotBigWig(kp, data.panel = 1, r0=0.64, r1=0.89 , data="Coverage/Sech_CID2_bwa_q0_RPM.bw", ymax=1700, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data.panel = 1, r0=0.64, r1=0.89 , data="Coverage/Sech_CID2_bwa_q30_RPM.bw", ymax=1700, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=1700, r0=0.64, r1=0.89 , numticks = 2,cex=1)
kpAddLabels(kp, labels = "Cenp-A enrichment (RPM)", r0=0.64, r1=0.89, srt=90, pos=3, cex=0.5, label.margin = 0.035)
kpPlotRegions(kp, data=peaks2q0, data.panel = 1, r0=0.6,r1=0.62 ,avoid.overlapping=FALSE,col="grey",border=NA)
kpPlotRegions(kp, data=peaks2q30, data.panel = 1, r0=0.6,r1=0.62,avoid.overlapping=FALSE,col="black",border=NA)

kp<- kpPlotBigWig(kp, data.panel = 1, r0=0.32, r1=0.57 , data="Coverage/Sech_CID3_bwa_q0_RPM.bw", ymax=1700, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data.panel = 1, r0=0.32, r1=0.57 , data="Coverage/Sech_CID3_bwa_q30_RPM.bw", ymax=1700, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=1700, r0=0.32, r1=0.57 , numticks = 2,cex=1)
kpAddLabels(kp, labels = "Cenp-A enrichment (RPM)", r0=0.32, r1=0.57, srt=90, pos=3, cex=0.5, label.margin = 0.035)
kpPlotRegions(kp, data=peaks3q0, data.panel = 1, r0=0.28,r1=0.30,avoid.overlapping=FALSE,col="grey",border=NA)
kpPlotRegions(kp, data=peaks3q30, data.panel = 1, r0=0.28,r1=0.30,avoid.overlapping=FALSE,col="black",border=NA)

kp<- kpPlotBigWig(kp, data.panel = 1,  r0=0, r1=0.25, data="Coverage/Sech_IgG_bwa_q0_RPM.bw", ymax=1700, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data.panel = 1,  r0=0, r1=0.25, data="Coverage/Sech_IgG_bwa_q30_RPM.bw", ymax=1700, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=1700, r0=0, r1=0.25, numticks = 2,cex=1)
kpAddLabels(kp, labels = "Cenp-A enrichment (RPM)", srt=90, pos=3, cex=0.5, r0=0, r1=0.25, label.margin = 0.035)

dev.off()



#########  Cen 2/3 (Contig188) ##########
chr=c("Contig188")
detail.region <- toGRanges(data.frame("Contig188",150000,500000))


pdf("Sech_Contig188_otherReplicate.pdf")
pp <- getDefaultPlotParams(plot.type=1)
pp$leftmargin <- 0.15
pp$topmargin <- 5
pp$bottommargin <- 15
pp$ideogramheight <- 5
pp$data1inmargin <- 10



kp <- plotKaryotype(genome = custom.genome, cytobands = custom.cytobands, ideogram.plotter = NULL,label=NULL,cex=0.5, plot.params = pp, chromosomes=chr, zoom=detail.region)

kpAddBaseNumbers(kp,tick.dist=25000,tick.len=6, cex=1, digits=3, units="kb" )
kpAddCytobandsAsLine(kp,color.table=vv,lwd=20)

kp<- kpPlotBigWig(kp, data.panel = 1, r0=0.64, r1=0.89 , data="Coverage/Sech_CID2_bwa_q0_RPM.bw", ymax=150, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data.panel = 1, r0=0.64, r1=0.89 , data="Coverage/Sech_CID2_bwa_q30_RPM.bw", ymax=150, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=150, r0=0.64, r1=0.89 , numticks = 2,cex=1)
kpAddLabels(kp, labels = "Cenp-A enrichment (RPM)", r0=0.64, r1=0.89, srt=90, pos=3, cex=0.5, label.margin = 0.035)
kpPlotRegions(kp, data=peaks2q0, data.panel = 1, r0=0.6,r1=0.62 ,avoid.overlapping=FALSE,col="grey",border=NA)
kpPlotRegions(kp, data=peaks2q30, data.panel = 1, r0=0.6,r1=0.62,avoid.overlapping=FALSE,col="black",border=NA)

kp<- kpPlotBigWig(kp, data.panel = 1, r0=0.32, r1=0.57 , data="Coverage/Sech_CID3_bwa_q0_RPM.bw", ymax=150, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data.panel = 1, r0=0.32, r1=0.57 , data="Coverage/Sech_CID3_bwa_q30_RPM.bw", ymax=150, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=150, r0=0.32, r1=0.57 , numticks = 2,cex=1)
kpAddLabels(kp, labels = "Cenp-A enrichment (RPM)", r0=0.32, r1=0.57, srt=90, pos=3, cex=0.5, label.margin = 0.035)
kpPlotRegions(kp, data=peaks3q0, data.panel = 1, r0=0.28,r1=0.30,avoid.overlapping=FALSE,col="grey",border=NA)
kpPlotRegions(kp, data=peaks3q30, data.panel = 1, r0=0.28,r1=0.30,avoid.overlapping=FALSE,col="black",border=NA)

kp<- kpPlotBigWig(kp, data.panel = 1,  r0=0, r1=0.25, data="Coverage/Sech_IgG_bwa_q0_RPM.bw", ymax=150, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data.panel = 1,  r0=0, r1=0.25, data="Coverage/Sech_IgG_bwa_q30_RPM.bw", ymax=150, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=150, r0=0, r1=0.25, numticks = 2,cex=1)
kpAddLabels(kp, labels = "Cenp-A enrichment (RPM)", srt=90, pos=3, cex=0.5, r0=0, r1=0.25, label.margin = 0.035)


dev.off()




#########  Cen X (Contig183) ##########
chr=c("Contig183")
detail.region <- toGRanges(data.frame("Contig183",550000,850000))


pdf("Sech_Contig183_otherReplicate.pdf")
pp <- getDefaultPlotParams(plot.type=1)
pp$leftmargin <- 0.15
pp$topmargin <- 5
pp$bottommargin <- 15
pp$ideogramheight <- 5
pp$data1inmargin <- 10



kp <- plotKaryotype(genome = custom.genome, cytobands = custom.cytobands, ideogram.plotter = NULL,label=NULL,cex=0.5, plot.params = pp, chromosomes=chr, zoom=detail.region)

kpAddBaseNumbers(kp,tick.dist=20000,tick.len=6, cex=1, digits=3, units="kb" )
kpAddCytobandsAsLine(kp,color.table=vv,lwd=20)

kp<- kpPlotBigWig(kp, data.panel = 1, r0=0.64, r1=0.89 , data="Coverage/Sech_CID2_bwa_q0_RPM.bw", ymax=400, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data.panel = 1, r0=0.64, r1=0.89 , data="Coverage/Sech_CID2_bwa_q30_RPM.bw", ymax=400, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=400, r0=0.64, r1=0.89 , numticks = 2,cex=1)
kpAddLabels(kp, labels = "Cenp-A enrichment (RPM)", r0=0.64, r1=0.89, srt=90, pos=3, cex=0.5, label.margin = 0.035)
kpPlotRegions(kp, data=peaks2q0, data.panel = 1, r0=0.6,r1=0.62 ,avoid.overlapping=FALSE,col="grey",border=NA)
kpPlotRegions(kp, data=peaks2q30, data.panel = 1, r0=0.6,r1=0.62,avoid.overlapping=FALSE,col="black",border=NA)

kp<- kpPlotBigWig(kp, data.panel = 1, r0=0.32, r1=0.57 , data="Coverage/Sech_CID3_bwa_q0_RPM.bw", ymax=400, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data.panel = 1, r0=0.32, r1=0.57 , data="Coverage/Sech_CID3_bwa_q30_RPM.bw", ymax=400, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=400, r0=0.32, r1=0.57 , numticks = 2,cex=1)
kpAddLabels(kp, labels = "Cenp-A enrichment (RPM)", r0=0.32, r1=0.57, srt=90, pos=3, cex=0.5, label.margin = 0.035)
kpPlotRegions(kp, data=peaks3q0, data.panel = 1, r0=0.28,r1=0.30,avoid.overlapping=FALSE,col="grey",border=NA)
kpPlotRegions(kp, data=peaks3q30, data.panel = 1, r0=0.28,r1=0.30,avoid.overlapping=FALSE,col="black",border=NA)

kp<- kpPlotBigWig(kp, data.panel = 1,  r0=0, r1=0.25, data="Coverage/Sech_IgG_bwa_q0_RPM.bw", ymax=400, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data.panel = 1,  r0=0, r1=0.25, data="Coverage/Sech_IgG_bwa_q30_RPM.bw", ymax=400, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=400, r0=0, r1=0.25, numticks = 2,cex=1)
kpAddLabels(kp, labels = "Cenp-A enrichment (RPM)", srt=90, pos=3, cex=0.5, r0=0, r1=0.25, label.margin = 0.035)

dev.off()



#########  Cen 4 (Contig40) ##########
chr=c("Contig40")

pdf("Sech_Contig40_otherReplicate.pdf")
pp <- getDefaultPlotParams(plot.type=1)
pp$leftmargin <- 0.15
pp$topmargin <- 5
pp$bottommargin <- 15
pp$ideogramheight <- 5
pp$data1inmargin <- 10


kp <- plotKaryotype(genome = custom.genome, cytobands = custom.cytobands, ideogram.plotter = NULL,label=NULL,cex=0.5, plot.params = pp, chromosomes=chr)

kpAddBaseNumbers(kp,tick.dist=10000,tick.len=6, cex=1, digits=3, units="kb" )
kpAddCytobandsAsLine(kp,color.table=vv,lwd=20)

kp<- kpPlotBigWig(kp, data.panel = 1, r0=0.64, r1=0.89 , data="Coverage/Sech_CID2_bwa_q0_RPM.bw", ymax=800, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data.panel = 1, r0=0.64, r1=0.89 , data="Coverage/Sech_CID2_bwa_q30_RPM.bw", ymax=800, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=800, r0=0.64, r1=0.89 , numticks = 2,cex=1)
kpAddLabels(kp, labels = "Cenp-A enrichment (RPM)", r0=0.64, r1=0.89, srt=90, pos=3, cex=0.5, label.margin = 0.035)
kpPlotRegions(kp, data=peaks2q0, data.panel = 1, r0=0.6,r1=0.62 ,avoid.overlapping=FALSE,col="grey",border=NA)
kpPlotRegions(kp, data=peaks2q30, data.panel = 1, r0=0.6,r1=0.62,avoid.overlapping=FALSE,col="black",border=NA)

kp<- kpPlotBigWig(kp, data.panel = 1, r0=0.32, r1=0.57 , data="Coverage/Sech_CID3_bwa_q0_RPM.bw", ymax=800, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data.panel = 1, r0=0.32, r1=0.57 , data="Coverage/Sech_CID3_bwa_q30_RPM.bw", ymax=800, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=800, r0=0.32, r1=0.57 , numticks = 2,cex=1)
kpAddLabels(kp, labels = "Cenp-A enrichment (RPM)", r0=0.32, r1=0.57, srt=90, pos=3, cex=0.5, label.margin = 0.035)
kpPlotRegions(kp, data=peaks3q0, data.panel = 1, r0=0.28,r1=0.30,avoid.overlapping=FALSE,col="grey",border=NA)
kpPlotRegions(kp, data=peaks3q30, data.panel = 1, r0=0.28,r1=0.30,avoid.overlapping=FALSE,col="black",border=NA)

kp<- kpPlotBigWig(kp, data.panel = 1,  r0=0, r1=0.25, data="Coverage/Sech_IgG_bwa_q0_RPM.bw", ymax=800, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data.panel = 1,  r0=0, r1=0.25, data="Coverage/Sech_IgG_bwa_q30_RPM.bw", ymax=800, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=800, r0=0, r1=0.25, numticks = 2,cex=1)
kpAddLabels(kp, labels = "Cenp-A enrichment (RPM)", srt=90, pos=3, cex=0.5, r0=0, r1=0.25, label.margin = 0.035)

dev.off()


#########  Cen Y (Y_scaffold1) ##########
chr=c("Y_scaffold1")
detail.region <- toGRanges(data.frame("Y_scaffold1",750000,1100000))

pdf("Sech_Y_scaffold1_otherReplicate.pdf")
pp <- getDefaultPlotParams(plot.type=1)
pp$leftmargin <- 0.15
pp$topmargin <- 5
pp$bottommargin <- 15
pp$ideogramheight <- 5
pp$data1inmargin <- 10


kp <- plotKaryotype(genome = custom.genome, cytobands = custom.cytobands, ideogram.plotter = NULL,label=NULL,cex=0.5, plot.params = pp, chromosomes=chr, zoom=detail.region)

kpAddBaseNumbers(kp,tick.dist=50000,tick.len=6, cex=1, digits=3, units="kb" )
kpAddCytobandsAsLine(kp,color.table=vv,lwd=20)

kp<- kpPlotBigWig(kp, data.panel = 1, r0=0.64, r1=0.89 , data="Coverage/Sech_CID2_bwa_q0_RPM.bw", ymax=800, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data.panel = 1, r0=0.64, r1=0.89 , data="Coverage/Sech_CID2_bwa_q30_RPM.bw", ymax=800, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=800, r0=0.64, r1=0.89 , numticks = 2,cex=1)
kpAddLabels(kp, labels = "Cenp-A enrichment (RPM)", r0=0.64, r1=0.89, srt=90, pos=3, cex=0.5, label.margin = 0.035)
kpPlotRegions(kp, data=peaks2q0, data.panel = 1, r0=0.6,r1=0.62 ,avoid.overlapping=FALSE,col="grey",border=NA)
kpPlotRegions(kp, data=peaks2q30, data.panel = 1, r0=0.6,r1=0.62,avoid.overlapping=FALSE,col="black",border=NA)

kp<- kpPlotBigWig(kp, data.panel = 1, r0=0.32, r1=0.57 , data="Coverage/Sech_CID3_bwa_q0_RPM.bw", ymax=800, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data.panel = 1, r0=0.32, r1=0.57 , data="Coverage/Sech_CID3_bwa_q30_RPM.bw", ymax=800, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=800, r0=0.32, r1=0.57 , numticks = 2,cex=1)
kpAddLabels(kp, labels = "Cenp-A enrichment (RPM)", r0=0.32, r1=0.57, srt=90, pos=3, cex=0.5, label.margin = 0.035)
kpPlotRegions(kp, data=peaks3q0, data.panel = 1, r0=0.28,r1=0.30,avoid.overlapping=FALSE,col="grey",border=NA)
kpPlotRegions(kp, data=peaks3q30, data.panel = 1, r0=0.28,r1=0.30,avoid.overlapping=FALSE,col="black",border=NA)

kp<- kpPlotBigWig(kp, data.panel = 1,  r0=0, r1=0.25, data="Coverage/Sech_IgG_bwa_q0_RPM.bw", ymax=800, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data.panel = 1,  r0=0, r1=0.25, data="Coverage/Sech_IgG_bwa_q30_RPM.bw", ymax=800, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=800, r0=0, r1=0.25, numticks = 2,cex=1)
kpAddLabels(kp, labels = "Cenp-A enrichment (RPM)", srt=90, pos=3, cex=0.5, r0=0, r1=0.25, label.margin = 0.035)

dev.off()

