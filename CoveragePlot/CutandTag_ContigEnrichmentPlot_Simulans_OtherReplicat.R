library(karyoploteR)
library(regioneR)
library(GenomicRanges) 
library(rtracklayer) 
library(IRanges) 
library(devtools)
library(stringr)


# Building a Genome
mygenome= read.table("Reference/dsim_scaffold2_2019.chrom.sizes",comment.char = "")
mygenome2 = data.frame(mygenome[,1],rep(1,dim(mygenome)[1]),mygenome[,2])
custom.genome <- toGRanges(mygenome2)

# Building a Cytoband
cytobands = read.csv("Reference/dsim_scaffold2_2019_Sim_011523.cytoband.mod", sep="\t", header=T)
custom.cytobands <- toGRanges(cytobands)

# Colors Cytoband
color_legend=read.csv("Reference/name-colorv3.csv", sep=";", header=F)
vv=as.character(color_legend$V2)
names(vv)=color_legend$V1

# Peaks
peaks2q0<-read.table("Macs2_Peaks/XD1_CID2_bwa_q0_macs2_peaks.narrowPeak",header=F, sep="\t", comment.char = "")
peaks2q0<- toGRanges(peaks2q0)
peaks2q30<-read.table("Macs2_Peaks/XD1_CID2_bwa_q30_macs2_peaks.narrowPeak",header=F, sep="\t", comment.char = "")
peaks2q30<- toGRanges(peaks2q30)

peaks3q0<-read.table("Macs2_Peaks/XD1_CID3_bwa_q0_macs2_peaks.narrowPeak",header=F, sep="\t", comment.char = "")
peaks3q0<- toGRanges(peaks3q0)
peaks3q30<-read.table("Macs2_Peaks/XD1_CID3_bwa_q30_macs2_peaks.narrowPeak",header=F, sep="\t", comment.char = "")
peaks3q30<- toGRanges(peaks3q30)


#########  Cen 2/3 (Contig2) ##########
chr=c("Contig_2")
detail.region <- toGRanges(data.frame("Contig_2",950000,1150000))


pdf("Sim_Contig_2_otherReplicate.pdf")
pp <- getDefaultPlotParams(plot.type=1)
pp$leftmargin <- 0.15
pp$topmargin <- 5
pp$bottommargin <- 15
pp$ideogramheight <- 5
pp$data1inmargin <- 10



kp <- plotKaryotype(genome = custom.genome, cytobands = custom.cytobands, ideogram.plotter = NULL,label=NULL,cex=0.5, plot.params = pp, chromosomes=chr, zoom=detail.region)

kpAddBaseNumbers(kp,tick.dist=25000,tick.len=6, cex=1, digits=3, units="mb" )
kpAddCytobandsAsLine(kp,color.table=vv,lwd=20)

kp<- kpPlotBigWig(kp, data.panel = 1, r0=0.64, r1=0.89 , data="Coverage/XD1_CID2_bwa_q0_RPM.bw", ymax=800, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data.panel = 1, r0=0.64, r1=0.89 , data="Coverage/XD1_CID2_bwa_q30_RPM.bw", ymax=800, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=800, r0=0.64, r1=0.89 , numticks = 2,cex=1)
kpAddLabels(kp, labels = "Cenp-A enrichment (RPM)", r0=0.64, r1=0.89, srt=90, pos=3, cex=0.5, label.margin = 0.035)
kpPlotRegions(kp, data=peaks2q0, data.panel = 1, r0=0.6,r1=0.62 ,avoid.overlapping=FALSE,col="grey",border=NA)
kpPlotRegions(kp, data=peaks2q30, data.panel = 1, r0=0.6,r1=0.62,avoid.overlapping=FALSE,col="black",border=NA)

kp<- kpPlotBigWig(kp, data.panel = 1, r0=0.32, r1=0.57 , data="Coverage/XD1_CID3_bwa_q0_RPM.bw", ymax=800, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data.panel = 1, r0=0.32, r1=0.57 , data="Coverage/XD1_CID3_bwa_q30_RPM.bw", ymax=800, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=800, r0=0.32, r1=0.57 , numticks = 2,cex=1)
kpAddLabels(kp, labels = "Cenp-A enrichment (RPM)", r0=0.32, r1=0.57, srt=90, pos=3, cex=0.5, label.margin = 0.035)
kpPlotRegions(kp, data=peaks3q0, data.panel = 1, r0=0.28,r1=0.30,avoid.overlapping=FALSE,col="grey",border=NA)
kpPlotRegions(kp, data=peaks3q30, data.panel = 1, r0=0.28,r1=0.30,avoid.overlapping=FALSE,col="black",border=NA)

kp<- kpPlotBigWig(kp, data.panel = 1,  r0=0, r1=0.25, data="Coverage/XD1_IgG_bwa_q0_RPM.bw", ymax=800, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data.panel = 1,  r0=0, r1=0.25, data="Coverage/XD1_IgG_bwa_q30_RPM.bw", ymax=800, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=800, r0=0, r1=0.25, numticks = 2,cex=1)
kpAddLabels(kp, labels = "Cenp-A enrichment (RPM)", srt=90, pos=3, cex=0.5, r0=0, r1=0.25, label.margin = 0.035)

dev.off()



#########  Cen 2/3 (Contig14) ##########
chr=c("Contig_14")
detail.region <- toGRanges(data.frame("Contig_14",50000,250000))


pdf("Sim_Contig_14_otherReplicate.pdf")
pp <- getDefaultPlotParams(plot.type=1)
pp$leftmargin <- 0.15
pp$topmargin <- 5
pp$bottommargin <- 15
pp$ideogramheight <- 5
pp$data1inmargin <- 10



kp <- plotKaryotype(genome = custom.genome, cytobands = custom.cytobands, ideogram.plotter = NULL,label=NULL,cex=0.5, plot.params = pp, chromosomes=chr, zoom=detail.region)

kpAddBaseNumbers(kp,tick.dist=25000,tick.len=6, cex=1, digits=3, units="mb" )
kpAddCytobandsAsLine(kp,color.table=vv,lwd=20)

kp<- kpPlotBigWig(kp, data.panel = 1, r0=0.64, r1=0.89 , data="Coverage/XD1_CID2_bwa_q0_RPM.bw", ymax=800, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data.panel = 1, r0=0.64, r1=0.89 , data="Coverage/XD1_CID2_bwa_q30_RPM.bw", ymax=800, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=800, r0=0.64, r1=0.89 , numticks = 2,cex=1)
kpAddLabels(kp, labels = "Cenp-A enrichment (RPM)", r0=0.64, r1=0.89, srt=90, pos=3, cex=0.5, label.margin = 0.035)
kpPlotRegions(kp, data=peaks2q0, data.panel = 1, r0=0.6,r1=0.62 ,avoid.overlapping=FALSE,col="grey",border=NA)
kpPlotRegions(kp, data=peaks2q30, data.panel = 1, r0=0.6,r1=0.62,avoid.overlapping=FALSE,col="black",border=NA)

kp<- kpPlotBigWig(kp, data.panel = 1, r0=0.32, r1=0.57 , data="Coverage/XD1_CID3_bwa_q0_RPM.bw", ymax=800, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data.panel = 1, r0=0.32, r1=0.57 , data="Coverage/XD1_CID3_bwa_q30_RPM.bw", ymax=800, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=800, r0=0.32, r1=0.57 , numticks = 2,cex=1)
kpAddLabels(kp, labels = "Cenp-A enrichment (RPM)", r0=0.32, r1=0.57, srt=90, pos=3, cex=0.5, label.margin = 0.035)
kpPlotRegions(kp, data=peaks3q0, data.panel = 1, r0=0.28,r1=0.30,avoid.overlapping=FALSE,col="grey",border=NA)
kpPlotRegions(kp, data=peaks3q30, data.panel = 1, r0=0.28,r1=0.30,avoid.overlapping=FALSE,col="black",border=NA)

kp<- kpPlotBigWig(kp, data.panel = 1,  r0=0, r1=0.25, data="Coverage/XD1_IgG_bwa_q0_RPM.bw", ymax=800, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data.panel = 1,  r0=0, r1=0.25, data="Coverage/XD1_IgG_bwa_q30_RPM.bw", ymax=800, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=800, r0=0, r1=0.25, numticks = 2,cex=1)
kpAddLabels(kp, labels = "Cenp-A enrichment (RPM)", srt=90, pos=3, cex=0.5, r0=0, r1=0.25, label.margin = 0.035)

dev.off()



#########  Cen X (Contig137) ##########
chr=c("Contig_137")
detail.region <- toGRanges(data.frame("Contig_137",1,150000))


pdf("Sim_Contig_137_otherReplicate.pdf")
pp <- getDefaultPlotParams(plot.type=1)
pp$leftmargin <- 0.15
pp$topmargin <- 5
pp$bottommargin <- 15
pp$ideogramheight <- 5
pp$data1inmargin <- 10



kp <- plotKaryotype(genome = custom.genome, cytobands = custom.cytobands, ideogram.plotter = NULL,label=NULL,cex=0.5, plot.params = pp, chromosomes=chr, zoom=detail.region)

kpAddBaseNumbers(kp,tick.dist=20000,tick.len=6, cex=1, digits=3, units="kb" )
kpAddCytobandsAsLine(kp,color.table=vv,lwd=20)

kp<- kpPlotBigWig(kp, data.panel = 1, r0=0.64, r1=0.89 , data="Coverage/XD1_CID2_bwa_q0_RPM.bw", ymax=350, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data.panel = 1, r0=0.64, r1=0.89 , data="Coverage/XD1_CID2_bwa_q30_RPM.bw", ymax=350, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=350, r0=0.64, r1=0.89 , numticks = 2,cex=1)
kpAddLabels(kp, labels = "Cenp-A enrichment (RPM)", r0=0.64, r1=0.89, srt=90, pos=3, cex=0.5, label.margin = 0.035)
kpPlotRegions(kp, data=peaks2q0, data.panel = 1, r0=0.6,r1=0.62 ,avoid.overlapping=FALSE,col="grey",border=NA)
kpPlotRegions(kp, data=peaks2q30, data.panel = 1, r0=0.6,r1=0.62,avoid.overlapping=FALSE,col="black",border=NA)

kp<- kpPlotBigWig(kp, data.panel = 1, r0=0.32, r1=0.57 , data="Coverage/XD1_CID3_bwa_q0_RPM.bw", ymax=350, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data.panel = 1, r0=0.32, r1=0.57 , data="Coverage/XD1_CID3_bwa_q30_RPM.bw", ymax=350, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=350, r0=0.32, r1=0.57 , numticks = 2,cex=1)
kpAddLabels(kp, labels = "Cenp-A enrichment (RPM)", r0=0.32, r1=0.57, srt=90, pos=3, cex=0.5, label.margin = 0.035)
kpPlotRegions(kp, data=peaks3q0, data.panel = 1, r0=0.28,r1=0.30,avoid.overlapping=FALSE,col="grey",border=NA)
kpPlotRegions(kp, data=peaks3q30, data.panel = 1, r0=0.28,r1=0.30,avoid.overlapping=FALSE,col="black",border=NA)

kp<- kpPlotBigWig(kp, data.panel = 1,  r0=0, r1=0.25, data="Coverage/XD1_IgG_bwa_q0_RPM.bw", ymax=350, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data.panel = 1,  r0=0, r1=0.25, data="Coverage/XD1_IgG_bwa_q30_RPM.bw", ymax=350, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=350, r0=0, r1=0.25, numticks = 2,cex=1)
kpAddLabels(kp, labels = "Cenp-A enrichment (RPM)", srt=90, pos=3, cex=0.5, r0=0, r1=0.25, label.margin = 0.035)

dev.off()



#########  Cen 4 (Contig124) ##########
chr=c("Contig_124")


pdf("Sim_Contig_124_otherReplicate.pdf")
pp <- getDefaultPlotParams(plot.type=1)
pp$leftmargin <- 0.15
pp$topmargin <- 5
pp$bottommargin <- 15
pp$ideogramheight <- 5
pp$data1inmargin <- 10



kp <- plotKaryotype(genome = custom.genome, cytobands = custom.cytobands, ideogram.plotter = NULL,label=NULL,cex=0.5, plot.params = pp, chromosomes=chr)

kpAddBaseNumbers(kp,tick.dist=20000,tick.len=6, cex=1, digits=3, units="kb" )
kpAddCytobandsAsLine(kp,color.table=vv,lwd=20)

kp<- kpPlotBigWig(kp, data.panel = 1, r0=0.64, r1=0.89 , data="Coverage/XD1_CID2_bwa_q0_RPM.bw", ymax=600, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data.panel = 1, r0=0.64, r1=0.89 , data="Coverage/XD1_CID2_bwa_q30_RPM.bw", ymax=600, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=600, r0=0.64, r1=0.89 , numticks = 2,cex=1)
kpAddLabels(kp, labels = "Cenp-A enrichment (RPM)", r0=0.64, r1=0.89, srt=90, pos=3, cex=0.5, label.margin = 0.035)
kpPlotRegions(kp, data=peaks2q0, data.panel = 1, r0=0.6,r1=0.62 ,avoid.overlapping=FALSE,col="grey",border=NA)
kpPlotRegions(kp, data=peaks2q30, data.panel = 1, r0=0.6,r1=0.62,avoid.overlapping=FALSE,col="black",border=NA)

kp<- kpPlotBigWig(kp, data.panel = 1, r0=0.32, r1=0.57 , data="Coverage/XD1_CID3_bwa_q0_RPM.bw", ymax=600, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data.panel = 1, r0=0.32, r1=0.57 , data="Coverage/XD1_CID3_bwa_q30_RPM.bw", ymax=600, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=600, r0=0.32, r1=0.57 , numticks = 2,cex=1)
kpAddLabels(kp, labels = "Cenp-A enrichment (RPM)", r0=0.32, r1=0.57, srt=90, pos=3, cex=0.5, label.margin = 0.035)
kpPlotRegions(kp, data=peaks3q0, data.panel = 1, r0=0.28,r1=0.30,avoid.overlapping=FALSE,col="grey",border=NA)
kpPlotRegions(kp, data=peaks3q30, data.panel = 1, r0=0.28,r1=0.30,avoid.overlapping=FALSE,col="black",border=NA)

kp<- kpPlotBigWig(kp, data.panel = 1,  r0=0, r1=0.25, data="Coverage/XD1_IgG_bwa_q0_RPM.bw", ymax=600, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data.panel = 1,  r0=0, r1=0.25, data="Coverage/XD1_IgG_bwa_q30_RPM.bw", ymax=600, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=600, r0=0, r1=0.25, numticks = 2,cex=1)
kpAddLabels(kp, labels = "Cenp-A enrichment (RPM)", srt=90, pos=3, cex=0.5, r0=0, r1=0.25, label.margin = 0.035)

dev.off()




#########  Cen Y (Y_Contig135) ##########
chr=c("Y_Contig_135")
detail.region <- toGRanges(data.frame("Y_Contig_135",1800000,2200000))

pdf("Sim_YContig_135_otherReplicate.pdf")
pp <- getDefaultPlotParams(plot.type=1)
pp$leftmargin <- 0.15
pp$topmargin <- 5
pp$bottommargin <- 15
pp$ideogramheight <- 5
pp$data1inmargin <- 10



kp <- plotKaryotype(genome = custom.genome, cytobands = custom.cytobands, ideogram.plotter = NULL,label=NULL,cex=0.5, plot.params = pp, chromosomes=chr, zoom=detail.region)

kpAddBaseNumbers(kp,tick.dist=50000,tick.len=6, cex=1, digits=3, units="mb" )
kpAddCytobandsAsLine(kp,color.table=vv,lwd=20)

kp<- kpPlotBigWig(kp, data.panel = 1, r0=0.64, r1=0.89 , data="Coverage/XD1_CID2_bwa_q0_RPM.bw", ymax=250, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data.panel = 1, r0=0.64, r1=0.89 , data="Coverage/XD1_CID2_bwa_q30_RPM.bw", ymax=250, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=250, r0=0.64, r1=0.89 , numticks = 2,cex=1)
kpAddLabels(kp, labels = "Cenp-A enrichment (RPM)", r0=0.64, r1=0.89, srt=90, pos=3, cex=0.5, label.margin = 0.035)
kpPlotRegions(kp, data=peaks2q0, data.panel = 1, r0=0.6,r1=0.62 ,avoid.overlapping=FALSE,col="grey",border=NA)
kpPlotRegions(kp, data=peaks2q30, data.panel = 1, r0=0.6,r1=0.62,avoid.overlapping=FALSE,col="black",border=NA)

kp<- kpPlotBigWig(kp, data.panel = 1, r0=0.32, r1=0.57 , data="Coverage/XD1_CID3_bwa_q0_RPM.bw", ymax=250, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data.panel = 1, r0=0.32, r1=0.57 , data="Coverage/XD1_CID3_bwa_q30_RPM.bw", ymax=250, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=250, r0=0.32, r1=0.57 , numticks = 2,cex=1)
kpAddLabels(kp, labels = "Cenp-A enrichment (RPM)", r0=0.32, r1=0.57, srt=90, pos=3, cex=0.5, label.margin = 0.035)
kpPlotRegions(kp, data=peaks3q0, data.panel = 1, r0=0.28,r1=0.30,avoid.overlapping=FALSE,col="grey",border=NA)
kpPlotRegions(kp, data=peaks3q30, data.panel = 1, r0=0.28,r1=0.30,avoid.overlapping=FALSE,col="black",border=NA)

kp<- kpPlotBigWig(kp, data.panel = 1,  r0=0, r1=0.25, data="Coverage/XD1_IgG_bwa_q0_RPM.bw", ymax=250, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data.panel = 1,  r0=0, r1=0.25, data="Coverage/XD1_IgG_bwa_q30_RPM.bw", ymax=250, col="black",border=NA)
kpAxis(kp, ymin=0 , ymax=250, r0=0, r1=0.25, numticks = 2,cex=1)
kpAddLabels(kp, labels = "Cenp-A enrichment (RPM)", srt=90, pos=3, cex=0.5, r0=0, r1=0.25, label.margin = 0.035)

dev.off()

