library(karyoploteR)
library(regioneR)
library(GenomicRanges) 
library(rtracklayer) 
library(IRanges) 
library(devtools)
library(stringr)

# Colors Cytoband

color_legend=read.csv("Reference/name-colorv3.csv", sep=";", header=F)
vv=as.character(color_legend$V2)
names(vv)=color_legend$V1



#######################################################
########### Simulans Cen4 ######################
############################################

# Building a Genome
mygenome= read.table("Reference/dsim_scaffold2_2019.chrom.sizes",comment.char = "")
mygenome2 = data.frame(mygenome[,1],rep(1,dim(mygenome)[1]),mygenome[,2])
custom.genome <- toGRanges(mygenome2)

# Building a Cytoband
cytobands = read.csv("Reference/dsim_scaffold2_2019_Sim_011523.cytoband.mod", sep="\t", header=T)
custom.cytobands <- toGRanges(cytobands)

#peaks
peaksq0<-read.table("Macs2_Peaks/XD1_CID1_bwa_q0_macs2_peaks.narrowPeak",header=F, sep="\t", comment.char = "")
peaksq0<- toGRanges(peaksq0)

peaksq30<-read.table("Macs2_Peaks/XD1_CID1_bwa_q30_macs2_peaks.narrowPeak",header=F, sep="\t", comment.char = "")
peaksq30<- toGRanges(peaksq30)


pdf("SimCen4.pdf",width = 15, height = 7)
pp <- getDefaultPlotParams(plot.type=1)
pp$leftmargin <- 0.15
pp$topmargin <- 20
pp$bottommargin <- 20
pp$ideogramheight <- 20
pp$data1inmargin <- 1

# Karyoplot
chr=c("Contig_124")
detail.region <- toGRanges(data.frame("Contig_124",1,120000))

kp <- plotKaryotype(genome = custom.genome, cytobands = custom.cytobands, ideogram.plotter = NULL,label=NULL,cex=0.7, plot.params = pp, chromosomes=chr, zoom=detail.region)

kpAddBaseNumbers(kp,tick.dist=10000,tick.len=6, cex=1.5, digits=3, units="kb")
kpAddCytobandsAsLine(kp,color.table=vv,lwd=25)

kp<- kpPlotBigWig(kp, data.panel = 1, r0=0.06, r1=1, data="Coverage/XD1_CID1_bwa_q0_RPM.bw", ymax=600, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data.panel = 1, r0=0.06, r1=1, data="Coverage/XD1_CID1_bwa_q30_RPM.bw", ymax=600, col="black",border=NA)

kpAxis(kp, ymin=0 , ymax=600, r0=0.06, r1=1, numticks = 2,cex=1.8)
kpAddLabels(kp, labels = "Cenp-A enrichment (RPM)", srt=90, pos=3, cex=2, label.margin = 0.035)
kpPlotRegions(kp, data=peaksq0, data.panel = 1, r0=0,r1=0.025,avoid.overlapping=FALSE,col="grey",border=NA)
kpPlotRegions(kp, data=peaksq30, data.panel = 1, r0=0,r1=0.025,avoid.overlapping=FALSE,col="black",border=NA)
dev.off()

###########



#######################################################
########### Mauritiana Cen4 ######################
############################################

# Building a Genome
mygenome= read.table("Reference/dmau_scaffold2_V2.chrom.sizes",comment.char = "")
mygenome2 = data.frame(mygenome[,1],rep(1,dim(mygenome)[1]),mygenome[,2])
custom.genome <- toGRanges(mygenome2)

# Building a Cytoband
cytobands = read.csv("Reference/dmau_scaffold2_V2.Mau.042623.cytoband.mod", sep="\t", header=T)
custom.cytobands <- toGRanges(cytobands)

#peaks
peaksq0<-read.table("Macs2_Peaks/Maur_CID1_bwa_q0_macs2_peaks.narrowPeak",header=F, sep="\t", comment.char = "")
peaksq0<- toGRanges(peaksq0)

peaksq30<-read.table("Macs2_Peaks/Maur_CID1_bwa_q30_macs2_peaks.narrowPeak",header=F, sep="\t", comment.char = "")
peaksq30<- toGRanges(peaksq30)


pdf("MauCen4.pdf",width = 15, height = 7)
pp <- getDefaultPlotParams(plot.type=1)
pp$leftmargin <- 0.15
pp$topmargin <- 20
pp$bottommargin <- 20
pp$ideogramheight <- 20
pp$data1inmargin <- 1

# Karyoplot
chr=c("Contig134")
#detail.region <- toGRanges(data.frame("Contig6",250000,450000))

kp <- plotKaryotype(genome = custom.genome, cytobands = custom.cytobands, ideogram.plotter = NULL,cex=0.7,label=NULL, plot.params = pp, chromosomes=chr)

kpAddBaseNumbers(kp,tick.dist=10000,tick.len=6, cex=1.5, digits=3, units="kb" )
kpAddCytobandsAsLine(kp,color.table=vv,lwd=25)

kp<- kpPlotBigWig(kp, data.panel = 1, r0=0.06, r1=1, data="Coverage/Maur_CID1_bwa_q0_RPM.bw", ymax=200, col="grey",border=NA)
kp<- kpPlotBigWig(kp, data.panel = 1, r0=0.06, r1=1, data="Coverage/Maur_CID1_bwa_q30_RPM.bw", ymax=200, col="black",border=NA)

kpAxis(kp, ymin=0 , ymax=200, r0=0.06, r1=1, numticks = 2,cex=1.8)
kpAddLabels(kp, labels = "Cenp-A enrichment (RPM)", srt=90, pos=3, cex=2, label.margin = 0.035)
kpPlotRegions(kp, data=peaksq0, data.panel = 1, r0=0,r1=0.025,avoid.overlapping=FALSE,col="grey",border=NA)
kpPlotRegions(kp, data=peaksq30, data.panel = 1, r0=0,r1=0.025,avoid.overlapping=FALSE,col="black",border=NA)
dev.off()

###########


