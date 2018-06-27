###########################################################
library(marray)
#source(paste(dirSrc,"functions/heatmap.5.6.R",sep=""))
#source(paste(dirSrc,"functions/heatmapAcgh.7.3.R",sep=""))
source(paste(dirSrc,"functions/heatmap.5.7.R",sep=""))
source(paste(dirSrc,"functions/heatmapAcgh.7.4.R",sep=""))

clId=rep(T,nrow(datObj$metabImp))

outFormat="pdf"
outFormat="png"

datadir=""

centrFlag="_noCentering"
centrFlag=""

subsetFlag=""

numPr=500
numPr=2000

pThres=10^-8
pThres=10^-6
pThres=0.05

metabList=unique(ann$type); metabList=metabList[which(!metabList%in%c("internal standard"))]

compList=paste("_rnd",numPr,sep="")
compList=paste("_topVar",numPr,sep="")
compList=paste("_",metabList,sep="")
compList=paste("_25percMostVarMetab_",metabList,sep="")

datFlag=""
limSl=c(330,480)

colGeneId="id"; colIdPV="pv"; colNamePV="PV"

fName1=""

tblCC=NULL
for (compFlag in compList) {
    compFName=sapply(compFlag,function(x) {y= strsplit(x," ")[[1]]; y=paste(y[1:min(c(length(y),3))],collapse="_"); gsub("-|/","_",y)}, USE.NAMES=F)
    if (length(grep("_rnd",compFlag))==1) {
        rndVec=paste("_rnd",1:4,sep="")
        #rndVec=paste("_rnd",1:20,sep="")
    } else {
        rndVec=""
    }
    
    for (rndId in rndVec) {
        limFCmmu=c(-6,6)
        if (compFlag%in%c("_topSignif")) {
            load(paste(datadir,"dge_Homo_sapiens_",compFlag,".RData",sep=""))
            switch(compFlag,
            "_topSignif"={
                stat_1=stat1_4
            }
            )
            timeThis=as.integer(sub("hrs","",strsplit(compFlag,"_")[[1]][2]))
            compName1=paste(timeThis,"hr: TGFbeta vs. untreated",sep="")
        } else {
            if (substr(compFlag,1,nchar("_topVar"))=="_topVar") {
                compName1=paste("Top ",sub("_topVar","",compFlag)," most variable probesets",sep="")
            }
            if (length(grep("percMostVarMetab",compFlag))==1) {
                x1=strsplit(compFlag,"_")[[1]]; x1=x1[2:length(x1)]; if (length(x1)==1) x1=c(x1,"")
                compName1=paste(as.integer(sub("percMostVarMetab","",x1[1])),"% most variable metabolites",sep="")
                if (x1[2]%in%unique(datObj$ann$type)) {
                    compName1=paste(capWords(x1[2]),": ",compName1,sep="")
                }
            }
            if (substr(compFlag,1,nchar("_rnd"))=="_rnd") {
                compName1=paste("Random ",sub("_rnd","",compFlag)," probesets",sep="")
            }
            if (compFlag%in%paste("_",unique(datObj$ann$type),sep="")) {
                compName1=capWords(sub("_","",compFlag))
            }
        }
        for (subsetFFlag in c("","_metabClust2")) {
            for (subsetFlag in c("")) {
                cat("\n\n================ ",compFName,", ",subsetFFlag,"\n")
                compName2=compName1
                if (subsetFFlag=="_metabClust2") {
                    subsetFName=" middle cluster"
                    datadir=paste(sub("_","",compFName),"/",sep="")
                    fName=paste(fName1,compFName,subsetFlag,centrFlag,rndId,datFlag,sep="")
                    tbl=read.table(paste(datadir,"clusterInfoFeature",fName,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
                    prId=which(datObj$ann$id%in%tbl$id[which(tbl$clustId=="cluster2")])
                    cat("No. of metabs: ",length(prId),"\n")
                } else {
                    subsetFName=""
                    prId=1:nrow(datObj$ann)
                }
                if (length(prId)<2) {cat("Need at least 2 metabolites!!!\n");  next}
                if (subsetFlag=="") {
                    subsetName=""
                    #prId=NULL
                    samId=1:nrow(datObj$phen)
                    sampleBar="cluster"
                    geneBar="clusterPr"
                    nClust=c(2,4)
                    nClust=c(2,2)
                    nClust=c(3,2)
                } else {
                    samId=NULL
                    if (F) {
                        grpUniq=sub("_","",subsetFlag)
                        subsetName=paste(", ",grpUniq,sep="")
                        if (grpUniq=="gamma") grpUniq="_"
                        fName2=paste(datadir,"clusterInfoFeature",fName1,compFlag,centrFlag,rndId,".txt",sep="")
                        fName2=paste(datadir,"clusterInfoFeature",fName1,centrFlag,rndId,".txt",sep="")
                        prId=read.table(file=fName2, header=T, sep="\t", quote="", comment.char="", as.is=T)
                        prId=prId[,"probesetid"]
                        samId=which(tolower(datObj$phen[,varList])==grpUniq)
                        sampleBar=""
                        sampleBar="cluster"
                        geneBar=""
                        nClust=c(NA,NA)
                    }
                }
                if (compFlag%in%c("_topSignif")) {
                    i1=which(stat_1[,colIdPV]<pThres)
                    if (length(i1)==0) next
                    fNameOut=paste(fName1,compFName,subsetFFlag,subsetFlag,centrFlag,"_",colNamePV,pThres,datFlag,sep="")
                    header=paste(compName2,subsetName,", ",colNamePV,"<",pThres,sep="")
                    dat0=eset$expr
                    dat0=dat0[match(stat_1[,colGeneId][i1],rownames(dat0)),]
                } else {
                    #fNameOut=paste(fName1,compFName,subsetFlag,centrFlag,rndId,datFlag,sep="")
                    fNameOut=paste(fName1,compFName,subsetFFlag,subsetFlag,centrFlag,rndId,datFlag,sep="")
                    header=paste(compName2,ifelse(subsetFName=="","",paste(", ",subsetFName,sep="")),subsetName,sep="")
                }
                #fNameOut=paste(fNameOut,varFlag,sep="")
                fNameOut=gsub(" ","_",fNameOut)
                expr=datObj$metabImp[prId,samId]
                #expr=datObj$metabRaw[,samId]
                annRow=datObj$ann[prId,]
                phen=datObj$phen[samId,]
                #phen$id2=sapply(phen$id,function(x) {strsplit(x,"_")[[1]][2]},USE.NAMES=F)
                
                if (length(grep("_primary|_steroid",compFlag))) {
                    annRow=annRow[,c("id","type")]
                }
                
                annRowAll=annRow

                i2=1:nrow(expr)
                if (length(grep("_rnd",compFlag))==1) {
                    if (compFlag==paste("_rnd",numPr,sep="")) {
                        i1=1:numPr
                    }
                    header=paste(header,": ",length(i1)," random probesets",sep="")
                    geneBar="clusterPr"
                    #set.seed(5453)
                    i2=sample(1:nrow(expr),length(i1),replace=F)
                    expr=expr[i2,]
                    annRow=annRow[i2,]
                    i=1:nrow(expr)
                } else if (length(grep("_topVar",compFlag))==1) {
                    if (is.null(prId)) {
                        varGene=apply(expr,1,var,na.rm=T)
                        i2=order(varGene,decreasing=T)[1:numPr]
                    } else {
                        compName2=paste(compName2," based on all samples",sep="")
                        i2=match(prId,annRow[,"probesetid"])
                    }
                    header=paste(compName2,subsetName,sep="")
                    expr=expr[i2,]
                    annRow=annRow[i2,]
                    i=1:nrow(expr)
                } else if (length(grep("percMostVarMetab",compFlag))==1) {
                    x1=strsplit(compFlag,"_")[[1]]; x1=x1[2:length(x1)]; if (length(x1)==1) x1=c(x1,"")
                    if (x1[2]%in%unique(annRow$type)) {
                        #i2=which(annRow$type==x1[2])
                        i2=which(annRow$type==x1[2] & apply(expr,1,function(x) mean(!is.na(x))>=.5))
                        expr=expr[i2,]
                        annRow=annRow[i2,]
                    }
                    x=apply(expr,1,sd,na.rm=T)
                    i2=order(x,decreasing=T)[1:round(quantile(1:nrow(expr),probs=as.integer(sub("percMostVarMetab","",x1[1]))/100))]
                    expr=expr[i2,]
                    annRow=annRow[i2,]
                    #header=paste(header,", n=",nrow(expr),sep="")
                    geneBar="clusterPr"
                    i=1:nrow(expr)
                } else if (length(grep("_top",compFlag))==1) {
                    header=paste(header,", n=",nrow(expr),sep="")
                    geneBar=""
                    geneBar="clusterPr"
                    #i=order(annRow$logFC)
                    i=1:nrow(expr)
                } else if (compFlag%in%paste("_",unique(annRow$type),sep="")) {
                    #i2=which(paste("_",annRow$type,sep="")==compFlag)
                    i2=which(paste("_",annRow$type,sep="")==compFlag & apply(expr,1,function(x) mean(!is.na(x))>=.5))
                    expr=expr[i2,]
                    annRow=annRow[i2,]
                    geneBar="clusterPr"
                    i=1:nrow(expr)
                } else {
                    geneBar=""
                    expr=expr[i2,]
                    annRow=cbind(annRow[i2,],logFC=stat_1$logFC[match(annRow[i2,colGeneId],stat_1[,colGeneId])])
                    i=order(annRow$logFC)
                }
                if (!is.na(nClust[1])) {
                    nClust[1]=min(c(nrow(expr)-1,nClust[1]))
                    if (nrow(expr)<5) nClust[1]=NA
                }
                
                arrayData=expr[i,j]
                annRow=annRow[i,]
                annCol=phen[j,]
                
                annColAll=datObj$phen
                annColAll$id2=annColAll$id
                
                if (centrFlag=="") {
                    centr=apply(arrayData,1,median,na.rm=T)
                    arrayData=arrayData-centr
                }
                
                if (F) {
                    varFList=names(annRow)[grep("fdrBY_",names(annRow))]
                    varFName=paste(sub("_fdrBY","",varFList)," ",sep="")
                    k=which(varFList%in%names(annRow))
                    varFList=varFList[k]
                    varFName=varFName[k]
                    varFListAll=varFList
                    varFNameAll=varFName
                }
                varFList=varFName=NULL
                
                #varList=c("ptsd","sex")
                #varName=paste(c("PTSD","Sex")," ",sep="")
                varList=c("ptsd","sex","diary_tst","psg_tst","ln_delta_nrem")
                varName=paste(c("ptsd","sex","diary_tst","psg_tst","ln_delta_nrem")," ",sep="")
                k=which(varList%in%names(annCol))
                varListAll=varList
                varNameAll=varName
                varList=varList[k]
                varName=varName[k]
                
                colList=c("skyblue","blue","yellow","purple","red")
                colList=c("brown","red","orange","yellow","green","cyan","skyblue","blue","pink","magenta","purple","darkgreen")
                colList2=c("skyblue","blue")
                colHM=c("red","blue","grey")
                
                distMethod="pearson"
                linkMethod="ward.D2"
                
                #rowName=annRow$geneSymbol
                rowName=rep("",nrow(annRow))
                if (is.null(varFList)) {
                    rowCol=NULL
                } else {
                    rowCol=matrix(nrow=length(varFList),ncol=nrow(annRow))
                    for (varId in 1:length(varFList)) {
                        x=annRowAll[,varFList[varId]]
                        x[x==""]=NA; x=as.integer(as.factor(x))
                        grpUniq=sort(unique(x))
                        x=x[match(annRow$affyId,annRowAll$affyId)]
                        if (length(grpUniq)<=length(colList2)) {
                            rowCol[varId,]=colList2[x]
                        } else if (length(grpUniq)<=length(colList)) {
                            rowCol[varId,]=colList[x]
                        } else {
                            rowCol[varId,]=rainbow(length(grpUniq))[x]
                        }
                    }
                    rownames(rowCol)=varFName
                }
                
                if (F) {
                    if (subsetFlag=="") {
                        colName=rep("",ncol(arrayData))
                    } else {
                        colName=annCol$id2
                    }
                }
                #colName=annCol$id2
                colName=rep("",nrow(annCol))
                colCol=NULL
                colCol=matrix(nrow=length(varList),ncol=nrow(annCol))
                for (varId in 1:length(varList)) {
                    if (varList[varId]%in%c("diary_tst","psg_tst","ln_delta_nrem")) {
                        j=match(annCol$id,annColAll$id)
                        if (varList[varId]%in%c("diary_tst","psg_tst")) {
                            x=round(annColAll[,varList[varId]])
                            lim=limSl
                        } else if (varList[varId]%in%c("ln_delta_nrem")) {
                            x=round(100*annColAll[,varList[varId]])
                            lim=range(x,na.rm=T)
                        } else {
                            x=round(annColAll[,varList[varId]])
                            #lim=range(x,na.rm=T)
                            lim=quantile(x,probs=c(.1,.9),na.rm=T)
                        }
                        x[x<lim[1]]=lim[1]; x[x>lim[2]]=lim[2]
                        x=x-min(x,na.rm=T)+1
                        grpUniq=lim[1]:lim[2]
                        colColUniq=maPalette(high=colList2[2],low=colList2[1],k=length(grpUniq))
                        #colColUniq=gray(0:(length(grpUniq)-1)/length(grpUniq))
                        colCol[varId,]=colColUniq[x[j]]
                    } else {
                        if (varList[varId]%in%c("time")) {
                            x=annColAll[,varList[varId]]
                        } else {
                            x=as.character(annColAll[,varList[varId]])
                        }
                        x[x==""]=NA; x=as.integer(as.factor(x))
                        grpUniq=sort(unique(x))
                        x=x[match(annCol$id,annColAll$id)]
                        if (length(grpUniq)<=length(colList2)) {
                            colCol[varId,]=colList2[x]
                        } else if (length(grpUniq)<=length(colList)) {
                            colCol[varId,]=colList[x]
                        } else {
                            colCol[varId,]=rainbow(length(grpUniq))[x]
                        }
                    }
                }
                rownames(colCol)=varName
                
                print("summary(range(c(arrayData),na.rm=T))")
                print(summary(range(c(arrayData),na.rm=T)))
                if (centrFlag=="") {
                    limit=c(-120000,120000)
                    limit=c(-10000,10000)
                    limit=c(-8,8)
                    limit=c(-1,1)
                    limit=c(-3,3)
                } else {
                    limit=c(8,13)
                }
                limit=c(-1,1)
                main=NULL
                main=header
                
                switch(distMethod,
                "pearson"={distMat=as.dist(1 - cor(arrayData,method=distMethod,use="complete.obs"))
                    if (sampleBar=="cluster") {
                        clustC=hclust(distMat, method=linkMethod)
                    } else {
                        clustC=NA
                        nClust[2]=NA
                    }
                    if (geneBar=="clusterPr") {
                        distMat=as.dist(1 - cor(t(arrayData),method=distMethod,use="complete.obs"))
                        clustR=hclust(distMat, method=linkMethod)
                    } else {
                        clustR=NA
                        nClust[1]=NA
                    }
                },
                "spearman"={distMat=as.dist(1 - cor(arrayData,method=distMethod,use="complete.obs"))
                    distMat=as.dist(1 - cor(t(arrayData),method=distMethod,use="complete.obs"))
                },
                "euclidean"={distMat=dist(t(arrayData), method=distMethod)
                    distMat=dist(arrayData, method=distMethod)
                }
                )
                #subDir=""
                #subDir=paste(sub("_","",compFName),"/",sep="")
                subDir=paste(sub("_","",compFName),subsetFFlag,"/",sep="")
                if (subDir!="" & !file.exists(subDir)) dir.create(file.path(subDir))
                if (outFormat=="png") {
                    margins=c(6,1)
                    margins=c(10,20)
                    png(paste(subDir,"heatmap",fNameOut,".png",sep=""),width=480*2,height=480*2)
                } else {
                    margins=c(12,5)
                    pdf(paste(subDir,"heatmap",fNameOut,".pdf",sep=""))
                }
                totalC=ncol(arrayData)
                #hcc=heatmap3(x=arrayData, Rowv=as.dendrogram(clustR), Colv=as.dendrogram(clustC), distfun=distMethod, hclustfun=hclust, symm=F, ColSideColors=colCol, RowSideColors=rowCol, labCol=colName, labRow=rowName, ncr=nClust[1], ncc=nClust[2], scale="none", na.rm=F, margins=margins, main=main, xlab=NULL, ylab=NULL, zlm=limit,cexCol=2, , high=colHM[1], low=colHM[2], mid=colHM[3])
                hcc=heatmap3(x=arrayData, Rowv=clustR, Colv=clustC, distfun=distMethod, hclustfun=hclust, symm=F, ColSideColors=colCol, RowSideColors=rowCol, labCol=colName, labRow=rowName, ncr=nClust[1], ncc=nClust[2], scale="none", na.rm=F, margins=margins, main=main, xlab=NULL, ylab=NULL, zlm=limit,cexCol=2, , high=colHM[1], low=colHM[2], mid=colHM[3],totalC=totalC)
                dev.off()
                
                if (is.na(nClust[1])) {
                    clustId=paste("cluster",1,sep="")
                    tbl=data.frame(variable=clustR$labels[clustR$order],clustId,order=1:nrow(arrayData),stringsAsFactors=F)
                    tbl=cbind(annRow[clustR$order,],clustId,order=1:nrow(annRow))
                    write.table(tbl, paste(subDir,"clusterInfoFeature",fNameOut,".txt",sep=""), sep="\t", col.names=T, row.names=F, quote=F)
                } else {
                    clustId=cutree(clustR,k=nClust[1])[clustR$order]
                    k1=which(!duplicated(clustId))
                    for (k in 1:length(k1)) {
                        clustId[which(clustId==clustId[k1[k]])]=paste("cluster",k,sep="")
                    }
                    #tbl=as.data.frame(as.matrix(arrayData[clustR$order,]),stringsAsFactors=F)
                    tbl=data.frame(variable=clustR$labels[clustR$order],clustId,order=1:nrow(arrayData),stringsAsFactors=F)
                    tbl=cbind(annRow[clustR$order,],clustId,order=1:nrow(annRow))
                    write.table(tbl, paste(subDir,"clusterInfoFeature",fNameOut,".txt",sep=""), sep="\t", col.names=T, row.names=F, quote=F)
                }
                
                if (is.na(nClust[2])) {
                    clustId=paste("cluster",1,sep="")
                    tbl=cbind(annCol[clustC$order,which(names(annCol)%in%c("order"))],clustId,order=1:nrow(annCol))
                    write.table(tbl, paste(subDir,"clusterInfoSample",fNameOut,".txt",sep=""), sep="\t", col.names=T, row.names=F, quote=F)
                } else {
                    if (F) {
                        #png(paste(subDir,"clusterSamples",fNameOut,".png",sep=""))
                        pdf(paste(subDir,"clusterSamples",fNameOut,".pdf",sep=""))
                        plot(clustC,main=paste("Sample clusters with ",nClust[2]," main clusters marked in red",sep=""),xlab="",sub="",ylab=NULL,axes=F, cex=.2); rect.hclust(clustC,k=nClust[2])
                        dev.off()
                    }
                    
                    clustId=cutree(clustC,k=nClust[2])[clustC$order]
                    k1=which(!duplicated(clustId))
                    for (k in 1:length(k1)) {
                        clustId[which(clustId==clustId[k1[k]])]=paste("cluster",k,sep="")
                    }
                    
                    tbl=cbind(annCol[clustC$order,which(!names(annCol)%in%c("order"))],clustId,order=1:nrow(annCol))
                    write.table(tbl, paste(subDir,"clusterInfoSample",fNameOut,".txt",sep=""), sep="\t", col.names=T, row.names=F, quote=F)
                }
            }
        }
    }
}
subDir=""
subDir="legend/"
if (subDir!="" & !file.exists(subDir)){
    dir.create(file.path(subDir))
}
if (!is.null(rowCol)) {
    for (varId in 1:length(varFListAll)) {
        if (length(grep("signif_",varFListAll[varId]))==1) {
            nm=""
            header=""
        } else {
            nm=varFListAll[varId]
            header=varFNameAll[varId]
        }
        if (outFormat=="png") {
            png(paste(subDir,"heatmapGeneColorBarLegend_",nm,".png",sep=""))
        } else {
            pdf(paste(subDir,"heatmapGeneColorBarLegend_",nm,".pdf",sep=""))
        }
        x=annRowAll[,varFListAll[varId]]
        grpUniq=table(x)
        grpUniq=names(grpUniq)
        k=1:length(grpUniq)
        if (length(grpUniq)<=length(colList2)) {
            sampleColorLegend(tls=grpUniq[k],col=colList2,legendTitle=header)
        } else if (length(grpUniq)<=length(colList)) {
            sampleColorLegend(tls=grpUniq[k],col=colList,legendTitle=header)
        } else {
            sampleColorLegend(tls=grpUniq[k],col=rainbow(length(grpUniq)),legendTitle=header)
        }
        dev.off()
    }
}
if (!is.null(colCol)) {
    for (varId in 1:length(varListAll)) {
        if (varListAll[varId]%in%c("diary_tst","psg_tst","ln_delta_nrem")) {
            if (outFormat=="png") {
                png(paste(subDir,"heatmapSampleColorBarLegend_",varListAll[varId],".png",sep=""),width=480,height=0.3*480)
            } else if (outFormat=="pdf") {
                pdf(paste(subDir,"heatmapSampleColorBarLegend_",varListAll[varId],".pdf",sep=""))
            }
            if (varList[varId]%in%c("diary_tst","psg_tst")) {
                x=round(annColAll[,varListAll[varId]])
                lim=limSl
            } else if (varList[varId]%in%c("ln_delta_nrem")) {
                x=round(100*annColAll[,varListAll[varId]])
                lim=range(x,na.rm=T)
                lim=lim/100
            } else {
                x=round(annColAll[,varListAll[varId]])
                #lim=range(x,na.rm=T)
                lim=quantile(x,probs=c(.1,.9),na.rm=T)
            }
            heatmapColorBar(limit=lim,cols=rev(colList2),main=varNameAll[varId])
            if (F) {
                grpUniq=lim[1]:lim[2]
                colColUniq=gray(0:(length(grpUniq)-1)/length(grpUniq))
                heatmapColorBar(limit=lim,cols=c(colColUniq[c(length(colColUniq),1)],median(1:length(colColUniq))),main=varNameAll[varId])
            }
        } else {
            if (outFormat=="png") {
                png(paste(subDir,"heatmapSampleColorBarLegend_",varListAll[varId],".png",sep=""))
            } else {
                pdf(paste(subDir,"heatmapSampleColorBarLegend_",varListAll[varId],".pdf",sep=""))
            }
            if (varList[varId]%in%c("time")) {
                x=annColAll[,varListAll[varId]]
            } else {
                x=as.character(annColAll[,varListAll[varId]]); x[x==""]=NA
            }
            grpUniq=table(x)
            #grpUniq=paste(names(grpUniq)," (",grpUniq,")",sep="")
            grpUniq=names(grpUniq)
            k=1:length(grpUniq)
            if (length(grpUniq)<=length(colList2)) {
                sampleColorLegend(tls=grpUniq[k],col=colList2,legendTitle=varNameAll[varId])
            } else if (length(grpUniq)<=length(colList)) {
                sampleColorLegend(tls=grpUniq[k],col=colList,legendTitle=varNameAll[varId])
            } else {
                sampleColorLegend(tls=grpUniq[k],col=rainbow(length(grpUniq)),legendTitle=varNameAll[varId])
            }
        }
        dev.off()
    }
}
if (outFormat=="png") {
    png(paste(subDir,"heatmapColorRange.png",sep=""),width=480,height=140)
} else {
    pdf(paste(subDir,"heatmapColorRange.pdf",sep=""))
}
heatmapColorBar(limit=limit,cols=colHM,main="Heatmap color range")
dev.off()

###########################################################
###########################################################


###########################################################
###########################################################
