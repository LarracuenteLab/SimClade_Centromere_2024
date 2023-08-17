Species="Maur"

pdf(paste0(Species,".pdf"), width = 10, height = 6)
par(mar=c(9,5,5,3))
### RPM 
data=read.table(paste0(Species, "_CID_RPM_Top_newname.out"))
row.names(data)=data[,1]
data=data[,-1]
data_ord=data[order(data[,1], decreasing=T),]

aa=max(data_ord)
bb=aa/15
plot(data_ord[c(1:20),1],  xaxt="n", ylab= "Normalized RPM count", xlab="",ylim=c(0,aa), cex=1.5, pch=16)
points(data_ord[c(1:20),2], pch=16, cex=1.5)
points(data_ord[c(1:20),3], pch=16, cex=1.5)
axis(1,1:20,labels=FALSE)
text(1:20,rep(-bb,20),labels=rownames(data_ord)[1:20],srt = 45,xpd=NA,adj=c(1,1), cex=1)

dev.off()



