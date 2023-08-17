

########################
###### Mauritiana ######
########################
data=read.table("Maur.q0.Top-idr")
chrom=read.table("dmau_scaffold2_V2.chrom.sizes")

# Number of common peaks
data2=as.data.frame(cbind(data$V1,data$V3-data$V2))

contig=unique(data2$V1)
cc=NULL
for(i in 1:length(contig)){
 aa=sum(as.numeric(data2[which(data2$V1==contig[i]),2]))
  bb=c(contig[i], aa)
  cc=rbind(cc, bb)
}

data3=as.data.frame(cc)
data4=data3[order(as.numeric(data3$V2), decreasing=T),]

pdf("Contig_IDRpeaks_Mauritiana.pdf",width = 15, height = 7)
par(mar = c(8, 5,5, 5))
barplot(as.numeric(data4$V2), names.arg = data4$V1, main="Mauritiana q0", las=2)
dev.off()

########################
###### Simulans ######
########################
data=read.table("XD1.q0.Top-idr")
chrom=read.table("dsim_scaffold2_2019.chrom.sizes")

# Number of common peaks
data2=as.data.frame(cbind(data$V1,data$V3-data$V2))

contig=unique(data2$V1)
cc=NULL
for(i in 1:length(contig)){
  aa=sum(as.numeric(data2[which(data2$V1==contig[i]),2]))
  bb=c(contig[i], aa)
  cc=rbind(cc, bb)
}

data3=as.data.frame(cc)
data4=data3[order(as.numeric(data3$V2), decreasing=T),]

pdf("Contig_IDRpeaks_Simulans.pdf",width = 15, height = 7)
par(mar = c(8, 5,5, 5))
barplot(as.numeric(data4$V2), names.arg = data4$V1, main="Simulans q0", las=2)
dev.off()


########################
###### Sechellia ######
########################
data=read.table("Sech.q0.Top-idr")

chrom=read.table("dsec_scaffold2_2019.chrom.sizes")


# Number of common peaks
data2=as.data.frame(cbind(data$V1,data$V3-data$V2))

contig=unique(data2$V1)
cc=NULL
for(i in 1:length(contig)){
  aa=sum(as.numeric(data2[which(data2$V1==contig[i]),2]))
  bb=c(contig[i], aa)
  cc=rbind(cc, bb)
}

data3=as.data.frame(cc)
data4=data3[order(as.numeric(data3$V2), decreasing=T),]

pdf("Contig_IDRpeaks_Sechellia.pdf",width = 15, height = 7)
par(mar = c(8, 5,5, 5))
barplot(as.numeric(data4$V2), names.arg = data4$V1, main="Sechellia q30", las=2)
dev.off()


