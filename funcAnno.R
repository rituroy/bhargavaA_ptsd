####################################################################
####################################################################
## ----------------------------------------------

library(PAPi)


####################################################################
####################################################################
ann1=read.table(paste(datadir,"compound_lipid_worksheet.txt",sep=""), sep="\t", h=T, quote="", comment.char="",as.is=T,fill=T)
out=t(sapply(ann1[,1],function(x) {
    y=strsplit(x," ")[[1]]
    z=paste(y[1],y[2],sep=" ")
    if (z=="Resolving power") y=c(z,y[3:length(y)])
    y=c(y[1],paste(y[2:length(y)],collapse=" "))
    y
},USE.NAMES=F))
ann1=data.frame(id=out[,1],description=out[,2],stringsAsFactors=F)

####################################################################
####################################################################
x=unique(ann$annotation)
x=x[x!=""]
tbl=x
tbl=data.frame(id=x,stringsAsFactors=F)
tbl$id2=sapply(tbl$id,function(x) {
    y=x
    y=strsplit(y," or ")[[1]][1]
    y=strsplit(y,"; ")[[1]][1]
    y=sub(" +$","",y)
    y
},USE.NAMES=F)
#write.table(tbl$id,file=paste("compound.txt",sep=""), sep="\t", col.names=F, row.names=F, quote=F)
write.table(tbl$id2,file=paste("compound.txt",sep=""), sep="\t", col.names=F, row.names=F, quote=F)

####################################################################
####################################################################
library(MetaboAnalystR)

mSet<-InitDataObjects("conc", "pathora", FALSE)
cmpd.vec=tbl$id2
cmpd.vec=ann1$description
mSet<-Setup.MapData(mSet, cmpd.vec);
mSet<-CrossReferencing(mSet, "name");
mSet<-CreateMappingResultTable(mSet)
out=as.data.frame(mSet$dataSet$map.table,stringsAsFactors=F)
out[out$KEGG!="NA",]



####################################################################
####################################################################

datadir=""
ann2=read.table(paste(datadir,"compound_out.csv",sep=""), sep=",", h=T, quote="", comment.char="",as.is=T,fill=T)
names(ann2)=gsub(".","",sub("X","",names(ann2)),fixed=T)
for (k in 1:ncol(ann2)) {
    if (is.character(ann2[,k])) {
        ann2[,k]=gsub("\"","",ann2[,k])
        ann2[ann2[,k]=="NA",k]=NA
    }
}


x=cbind(tbl[1:nrow(ann2),],ann2,stringsAsFactors=F)
i=which(x$id!=x$Query)
#i=which(x$id!=x$Query)
x[i,][1:4,]
i1=i[1]
x$id[i1]; x$Query[i1]
