####################################################################
####################################################################
## ----------------------------------------------

plotFlag=""
plotFlag=c("","manhattanPlot")
plotFlag="_volcanoPlot"
plotFlag="_manhattanPlot"
plotFlag="_avgLog2ExprVsWilcoxVoomPV"
plotFlag="_log2FCplot"
plotFlag="_onePlot"

parList=list(ylimM=c(0,8.5))
parList=list()

reportFlag="signif"
reportFlag=""

geneSumFlag=T
geneSumFlag=F

outlierFlag=F
outlierFlag=T

respFlag="metab"

datFlag=""
datFlag="orig"
datFlag="imputed"
datFlag="raw_allMetab"

datFlag="voom"
datFlag="raw"

adjPFlag="qv"
adjPFlag="BY"

pThres2=0.2
pThres2=99

## -------------------
metabList=unique(ann$type)
compList=c("diary_tst","psg_tst","ln_delta_nrem","diary_tst+sex","psg_tst+sex","ln_delta_nrem+sex","diary_tst*sex","psg_tst*sex","ln_delta_nrem*sex")
compList=c("ptsd+sex","ptsd*sex")
compList=c("ptsd","sex","ptsd+sex","ptsd*sex")
compList=c("diary_tst","psg_tst","ln_delta_nrem")
compList=c("ptsd","sex")

metabFlag=metabList[1]

for (testFlag in c("_t","_wilcox")) {
    if (testFlag=="_t" & !datFlag%in%c("","voom")) {cat('testFlag=="_t" & !datFlag%in%c("","voom")!!!\n'); next}
    if (testFlag=="_wilcox" & !datFlag%in%c("orig","imputed","raw","raw_allMetab","imputed_allMetabImp")) {cat('testFlag=="_wilcox" & !datFlag%in%c("orig","imputed","raw","raw_allMetab","imputed_allMetabImp")!!!\n'); next}
    if (plotFlag=="_avgLog2ExprVsWilcoxVoomPV" & testFlag=="_t") next
    if (plotFlag=="_log2FCplot" & testFlag=="_t") next
    if (testFlag=="_t") {
        #datadir=paste("results/metabO/voom/pv",ifelse(pThres2>=1,1,pThres2),"/",sep="")
        datadir=paste("results/final/voom/table/",sep="")
    } else {
        #datadir=paste("results/metabO/",datFlag,"/",sep="")
        datadir=paste("results/final/wilcoxon/table/",sep="")
        datadir=paste("raw/",sep="")
    }
    for (metabFlag in metabList) {
        if (plotFlag[1]%in%c("_qqPlot","_histogram","_volcanoPlot","_manhattanPlot","_avgLog2ExprVsWilcoxVoomPV","_log2FCplot")) {
            png(paste(sub("_","",plotFlag[1]),"_",metabFlag,"_%1d.png",sep=""),width=2*length(compList)*240, height=2*240)
            par(mfcol=c(1,length(compList)))
        }
        for (compId in compList) {
            if (testFlag=="_wilcox" & compId%in%c("ptsd+sex","ptsd*sex","diary_tst+sex","psg_tst+sex","ln_delta_nrem+sex")) next
            if (compId%in%c("ptsd+sex","ptsd*sex","diary_tst+sex","psg_tst+sex","ln_delta_nrem+sex") & pThres2<1) {
                fName3=paste(ifelse(testFlag=="_t","","_wilcox"),"_",respFlag,"Resp_",sub("*","X",sub("+","_",compId,fixed=T),fixed=T),"_pv",pThres2,"_",metabFlag,sep="")
            } else {
                fName3=paste(ifelse(testFlag=="_t","","_wilcox"),"_",respFlag,"Resp_",sub("*","X",sub("+","_",compId,fixed=T),fixed=T),"_",metabFlag,sep="")
            }
            stat2=read.table(paste(datadir,"stat",fName3,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
            #log2fc=apply(expr[,j],1,function(x,grp) {mean(x[grp==grpUniq[2]],na.rm=T)-mean(x[grp==grpUniq[1]],na.rm=T)},grp=as.integer(grp)[j])
            if (testFlag=="_wilcox") {
                fName3=paste("_",respFlag,"Resp_",sub("*","X",sub("+","_",compId,fixed=T),fixed=T),"_",metabFlag,sep="")
                stat1=read.table(paste("results/voom/table/stat",fName3,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
                #stat1=read.table(paste("results/metabO/voom/pv1/stat",fName3,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
                #stat1=read.table(paste("results/final/voom/table/stat",fName3,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
                colId=grep("log2fc_|pv_|BY_",names(stat1))
                nm=c(names(stat2),paste(names(stat1)[colId],"_tstat",sep=""))
                i=match(stat1$metabolite,stat2$metabolite); i1=which(!is.na(i)); i2=i[i1]
                for (k in colId) {
                    x=rep(NA,nrow(stat2)); x[i2]=stat1[i1,k]
                    stat2=cbind(stat2,x)
                }
                names(stat2)=nm
            }
            colListPV=names(stat2)[grep("pv_",names(stat2))]
            k=grep("_tstat",colListPV)
            if (length(k)!=0) colListPV=colListPV[-k]
            if (reportFlag=="signif") {
                x=F
                for (cId in 1:length(colListPV)) {
                    colIdPV=c(colListPV[cId],sub("pv",adjPFlag,colListPV[cId])); pThres=0.05
                    if (any(stat2[,colIdPV[2]]<pThres)) {
                        x=T
                        break
                    }
                }
                if (!x) next
            }
            for (cId in 1:length(colListPV)) {
                colIdEst=names(stat2)[grep("log2fc",names(stat2))][cId]
                colIdPV=c(colListPV[cId],sub("pv",adjPFlag,colListPV[cId])); pThres=0.05
                fName1=paste(ifelse(testFlag=="_t","","_wilcox"),"_",respFlag,"Resp_",sub("*","X",sub("+","_",compId,fixed=T),fixed=T),sub("pv","_coef",colIdPV[1]),"_",metabFlag,sep="")
                compName=paste(capWords(metabFlag),": metabolite~",compId,ifelse(length(colListPV)==1,"",paste(". Term - ",sub("log2fc_","",colIdEst),sep="")),sep="")
                #compName=paste(ifelse(testFlag=="_t","t","Wilcoxon")," test, ",compName,sep="")
                #names(stat2)=sapply(names(stat2),function(x) {strsplit(x,"_")[[1]][1]},USE.NAMES=F)
                i=match(stat2$id,ann$IlmnID)
                fName=fName1
                iA2=match(stat2$metabolite,ann$id)
                ann2=ann[iA2,]
                
                ####################################################################
                
                if (all(is.na(stat2[,colIdEst]))) {next} #{cat("No stats computed !!!\n"); next}
                
                #cat("\n\n============= ",compName,"\n",sep="")
                cat("\n",compName,"\n",sep="")
                i=1:nrow(stat2)
                stat=stat2
                if (F) {
                    if (length(colListPV)!=1) {
                        fName=paste(fName1,sub("pv_","_",colIdPV[1][grep("pv_",colIdPV[1])]),sep="")
                    }
                }
                
                if (plotFlag=="_avgLog2ExprVsWilcoxVoomPV") {
                    x1=matrix(0,nrow=2,ncol=2,dimnames=list(paste("Wilcox: ",colIdPV[2],c(">=","<"),pThres,sep=""),paste("t: ",colIdPV[2],c(">=","<"),pThres,sep="")))
                    ii=i[which(sign(stat[i,colIdEst])==sign(stat[i,paste(colIdEst,"_t",sep="")]))]
                    x2=table(stat[ii,paste(colIdPV[2],"_t",sep="")]<pThres,stat[ii,colIdPV[2]]<pThres)
                    x1[match(rownames(x2),c("FALSE","TRUE")),match(colnames(x2),c("FALSE","TRUE"))]=x2
                    cat("\nSame direction for Wilcoxon & t tests:\n")
                    print(x1)
                    x1=matrix(0,nrow=2,ncol=2,dimnames=list(paste("Wilcox: ",colIdPV[2],c(">=","<"),pThres,sep=""),paste("t: ",colIdPV[2],c(">=","<"),pThres,sep="")))
                    ii=i[which(sign(stat[i,colIdEst])!=sign(stat[i,paste(colIdEst,"_t",sep="")]))]
                    x2=table(stat[ii,paste(colIdPV[2],"_t",sep="")]<pThres,stat[ii,colIdPV[2]]<pThres)
                    x1[match(rownames(x2),c("FALSE","TRUE")),match(colnames(x2),c("FALSE","TRUE"))]=x2
                    cat("\nOpposite directions for Wilcoxon & t tests:\n")
                    print(x1)
                } else {
                    x1=matrix(0,nrow=2,ncol=2,dimnames=list(c("dn","up"),paste(colIdPV[2],c(">=","<"),pThres,sep="")))
                    ii=i[which(stat[i,colIdEst]!=0)]
                    x2=table(stat[ii,colIdEst]>0,stat[ii,colIdPV[2]]<pThres)
                    x1[match(rownames(x2),c("FALSE","TRUE")),match(colnames(x2),c("FALSE","TRUE"))]=x2
                    colnames(x1)=c("not significant","significant")
                    #print(x1)
                    if (T) {
                        cat("------------------------------\n")
                        print(x1)
                        cat("------------------------------\n")
                    }
                    if (F) {
                        cat("------------------------------\n")
                        cat("",colnames(x1),"\n")
                        cat("------------------------------\n")
                        for (k in 1:nrow(x1)) {
                            cat(rownames(x1)[k],x1[k,],"\n")
                        }
                        cat("------------------------------\n")
                    }
                    if (F) {
                        cat("------------------------------\n")
                        cat(c("",colnames(x1),"\n"),collapse="\t")
                        cat("------------------------------\n")
                        for (k in 1:nrow(x1)) {
                            cat(c(rownames(x1)[k],x1[k,],"\n"),collapse="\t")
                        }
                        cat("------------------------------\n")
                    }
                }

                if (plotFlag[1]=="_onePlot") {
                    png(paste("plots",fName,".png",sep=""),width=2*480, height=1*480)
                    par(mfcol=c(1,2))
                }
                
                if (plotFlag[1]%in%c("_qqPlot")) {
                    if (plotFlag[1]=="") {
                        png(paste("qqPlot",fName,".png",sep=""))
                        header=compName
                    } else if (plotFlag[1]=="_qqPlot") {
                        header=compName
                    } else {
                        header=""
                    }
                    pvs <- sort(na.exclude(stat[i,colIdPV[1]]))
                    qqplot(-log10(runif(length(pvs),0,1)),-log10(pvs),xlab="Expected -log10(p-values) by random",ylab="Observed -log10(p-values)",pch=19,cex.axis=1.5,cex.lab=1.5,main=header)
                    abline(0,1)
                    if (plotFlag[1]=="") {
                        dev.off()
                    }
                }
                
                if (plotFlag[1]%in%c("","_onePlot","_histogram")) {
                    if (plotFlag[1]=="") {
                        png(paste("histogram",fName,".png",sep=""))
                        header=compName
                    } else if (plotFlag[1]=="_histogram") {
                        header=compName
                    } else {
                        header=""
                    }
                    hist(stat[i,colIdPV[1]],main=header,xlab="P-value",pch=19,cex.axis=1.5,cex.lab=1.5)
                    if (plotFlag[1]=="") {
                        dev.off()
                    } else if (plotFlag[1]=="_onePlot") {
                        title(main=compName)
                    }
                }
                
                if (plotFlag[1]%in%c("","_onePlot","_volcanoPlot")) {
                    if (plotFlag[1]=="") {
                        png(paste("volcanoPlot",fName,"_",colIdPV[2],pThres,".png",sep=""))
                        header=compName
                    } else if (plotFlag[1]=="_volcanoPlot") {
                        header=compName
                    } else {
                        header=""
                        header=compName
                    }
                    iThis=i
                    if (!outlierFlag) {
                        x=quantile(abs(stat[iThis,colIdEst]),probs=.95,na.rm=T)
                        iThis=iThis[which(abs(stat[iThis,colIdEst])<=x)]
                    }
                    if ("ylimM"%in%names(parList)) yLim=parList$ylimM else yLim=NULL
                    plot(stat[iThis,colIdEst],-log10(stat[iThis,colIdPV[1]]),ylim=yLim,xlab="Estimate",ylab="-log10(p-value)",pch=19,cex.axis=1.5,cex.lab=1.5,main=header,col="grey")
                    ii=iThis[which(stat[iThis,colIdPV[2]]<pThres)]
                    points(stat[ii,colIdEst],-log10(stat[ii,colIdPV[1]]),pch=19,col="red")
                    if (plotFlag[1]=="") {
                        dev.off()
                    }
                }
                
                if (plotFlag[1]%in%c("_manhattanPlot")) {
                    if (plotFlag[1]=="") {
                        png(paste("manhattanPlot",fName,"_",colIdPV[2],pThres,".png",sep=""))
                        header=compName
                    } else if (plotFlag[1]=="_manhattanPlot") {
                        header=compName
                    } else {
                        header=""
                    }
                    iThis=i
                    if ("ylimM"%in%names(parList)) yLim=parList$ylimM else yLim=NULLL
                    plot(ann[iA2[iThis],"MAPINFO"],-log10(stat[iThis,colIdPV[1]]),xlab=paste("Chromosome ",ann[iA2[iThis][1],"CHR"],sep=""),ylim=yLim,ylab="-log10(p-value)",pch=19,cex.axis=1.5,cex.lab=1.5,main=header,col="grey")
                    ii=iThis[which(stat[iThis,colIdPV[2]]<pThres)]
                    points(ann[iA2[ii],"MAPINFO"],-log10(stat[ii,colIdPV[1]]),pch=19,col="red")
                    if (plotFlag[1]=="") {
                        dev.off()
                    }
                }
                
                if (plotFlag[1]%in%c("_pvAdjpvPlot")) {
                    if (plotFlag[1]=="") {
                        png(paste("pvAdjpvPlot",fName,"_",colIdPV[2],pThres,".png",sep=""))
                        header=compName
                    } else if (plotFlag[1]=="_pvAdjpvPlot") {
                        header=compName
                    } else {
                        header=""
                    }
                    iThis=i
                    plot(stat[iThis,"pv"],stat[iThis,adjPFlag],xlab="P-value",ylab=adjPFlag,pch=19,cex.axis=1.5,cex.lab=1.5,main=header,col="black")
                    if (plotFlag[1]=="") {
                        dev.off()
                    }
                }
                
                if (plotFlag[1]%in%c("_avgLog2ExprVsWilcoxVoomPV")) {
                    if (plotFlag[1]=="") {
                        png(paste("avgLogCPMvsVoomWilcoxPVplot",fName,"_",colIdPV[2],pThres,".png",sep=""))
                        header=compName
                    } else if (plotFlag[1]=="_avgLog2ExprVsWilcoxVoomPV") {
                        header=compName
                    } else {
                        header=""
                    }
                    iThis=i
                    lim=c(-3,3)
                    lim=NULL
                    plot(meanMetabI[iThis],log10(stat[iThis,colIdPV[1]])-log10(stat[iThis,paste(colIdPV[1],"_t",sep="")]),ylim=lim,main=paste(header,"\nN = ",sum(!is.na(stat[iThis,paste(colIdPV[1],"_t",sep="")])),sep=""),xlab="Avg log2(expression)",ylab="log10(PV-wilcox/PV-limma)",cex=.2,pch=20,cex.main=1,cex.axis=1.5,cex.lab=1.5)
                    abline(h=0,col="red")
                    if (plotFlag[1]=="") {
                        dev.off()
                    }
                }
                
                if (plotFlag[1]%in%c("_log2FCplot")) {
                    if (plotFlag[1]=="") {
                        png(paste("log2FCplot",fName,".png",sep=""))
                        header=compName
                    } else if (plotFlag[1]=="_log2FCplot") {
                        header=compName
                    } else {
                        header=""
                    }
                    iThis=i
                    lim=range(c(stat[iThis,colIdEst],stat[iThis,paste(colIdEst,"_t",sep="")]),na.rm=T)
                    plot(stat[iThis,colIdEst],stat[iThis,paste(colIdEst,"_t",sep="")],xlim=lim,ylim=lim,xlab="Observed: log2FC",ylab="Estimated: log2FC",pch=19,cex.axis=1.5,cex.lab=1.5,main=header,col="grey")
                    abline(c(0,1))
                    if (plotFlag[1]=="") {
                        dev.off()
                    }
                }
                
                if (plotFlag[1]=="_onePlot") {
                    dev.off()
                }
                
                if (F) {
                    png(paste("plots",fName,".png",sep=""),width=3*240, height=1*240)
                    par(mfcol=c(1,3))
                    plot(1:6)
                    plot(1:6)
                    title(main=sub(" ReFACTor","\nReFACTor",compName))
                    plot(1:6)
                    dev.off()
                }
                
                if (F) {
                    ii=i[order(stat[i,colIdPV[2]])]
                    ii=ii[stat[ii,colIdPV[2]]<pThres]
                    tbl=cbind(ann[iA2,][ii,],stat[ii,c(colIdEst,unique(colIdPV))])
                    write.table(tbl, file=paste("stat",fName1,"_",colIdPV[2],pThres,".txt",sep=""), append=F,col.names=T,row.names=F, sep="\t",quote=F)
                }
            }
        }
        if (plotFlag[1]%in%c("_qqPlot","_histogram","_volcanoPlot","_manhattanPlot","_avgLog2ExprVsWilcoxVoomPV","_log2FCplot")) dev.off()
    }
}

####################################################################
####################################################################
## Scatter plot

metabFlag="lipid"
modelFlag="ptsd"

pThres=0.05
colIdEst=paste("log2fc_",modelFlag,sep="")
colIdPV=paste(c("pv","BY"),"_",modelFlag,sep="")

datFlag=""


metabList=unique(ann$type)
metabList="lipid"

compList=c("diary_tst","psg_tst","ln_delta_nrem","diary_tst+sex","psg_tst+sex","ln_delta_nrem+sex","diary_tst*sex","psg_tst*sex","ln_delta_nrem*sex")
compList=c("ptsd+sex","ptsd*sex")
compList=c("ptsd","sex","ptsd+sex","ptsd*sex")
compList=c("diary_tst*sex")

#metabolite	identifier	annotation	inChIkey	species	count	esiMode	mz	rt	type	log2fc_diary_tst	pv_diary_tst	BY_diary_tst	log2fc_maleVfemale	pv_maleVfemale	BY_maleVfemale	log2fc_diary_tstGenderInteraction	pv_diary_tstGenderInteraction	BY_diary_tstGenderInteraction
#X9.88_922.01	9.88_922.01			[M+H]+	45	(+) ESI	922.0112	9.88	lipid	-0.02397346	1.07E-06	0.010205564	-22.95897124	4.89E-21	9.36E-17	0.053173639	1.24E-10	2.38E-06

colIdNum=c("diary_tst","psg_tst","ln_delta_nrem")

phen2=phen
nm=names(phen2)
for (k in colIdNum) {
    x=median(phen2[,k],na.rm=T)
    grp=rep("",nrow(phen2))
    grp[which(phen2[,k]<=x)]="low"
    grp[which(phen2[,k]>x)]="high"
    grp=paste(as.integer(phen2[,k]>x)," ",k,sep="")
    phen2=cbind(phen2,grp)
    nm=c(nm,paste(k,"_status",sep=""))
}
names(phen2)=nm

for (metabFlag in metabList) {
    for (modelFlag in compList) {
        varList=strsplit(sub("*","-",sub("+","-",modelFlag,fixed=T),fixed=T),"-")[[1]]
        if (modelFlag%in%c("ptsd","sex")) {
            #datadir="results/wilcox_metabImp/"
            datadir="results/final/wilcoxon/table/"
        } else {
            datadir="results/final/voom/table/"
        }
        if (modelFlag%in%c("ptsd+sex","ptsd*sex") & pThres2<1) {
            fName3=paste("_",respFlag,"Resp_",sub("*","X",sub("+","_",modelFlag,fixed=T),fixed=T),"_pv",pThres2,"_",metabFlag,sep="")
        } else {
            fName3=paste("_",respFlag,"Resp_",sub("*","X",sub("+","_",modelFlag,fixed=T),fixed=T),"_",metabFlag,sep="")
        }
        stat2=read.table(paste(datadir,"stat",fName3,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
        if (length(varList)==2) {
            colId=grep("log2fc_|pv_|BY_",names(stat2))
            colId=colId[grep("Interaction",names(stat2)[colId])]
            colIdEst=names(stat2)[colId][1]
            colIdPV=names(stat2)[colId][2:3]
        } else {
            colIdEst=paste("log2fc_",modelFlag,sep="")
            colIdPV=paste(c("pv","BY"),"_",modelFlag,sep="")
        }
        for (vId in 1:length(varList)) if (varList[vId]%in%colIdNum) varList[vId]=paste(varList[vId],"_status",sep="")

        ann2=ann[match(stat2$metabolite,ann$id),]
        if (modelFlag%in%c("ptsd","sex")) {
            if (F) {
                expr=metabTR[match(stat2$metabolite,ann$id),]
                expr=metabOR2[match(stat2$metabolite,ann$id),]
                expr=metabOI2[match(stat2$metabolite,ann$id),]
                expr=metabTI2[match(stat2$metabolite,ann$id),]
                expr=metabOR[match(stat2$metabolite,ann$id),]
                expr=metabTR2[match(stat2$metabolite,ann$id),]
            }
            expr=metabR[match(stat2$metabolite,ann$id),]
        } else {
            expr=metabI[match(stat2$metabolite,ann$id),]
        }

        #grp=as.integer(as.factor(phen2[,modelFlag]))
        #j=which(!is.na(phen2[,modelFlag]))
        grp=rep("",nrow(phen2))
        j=1:nrow(phen2)
        for (vId in 1:length(varList)) {
            grp=paste(grp,phen2[,varList[vId]])
            j=j[which(!is.na(phen2[j,varList[vId]]))]
        }
        grp=grp[j]
        grpUniq=unique(sort(grp))
        if (length(grpUniq)==2) {
            if (datFlag=="orig") {
                x=c(expr)
                dat=log2(expr[,j]+min(x[which(x!=0)])/10)
                log2fc=apply(dat,1,function(x,grp) {mean(x[grp==grpUniq[2]],na.rm=T)-mean(x[grp==grpUniq[1]],na.rm=T)},grp=as.integer(grp)[j])
                log2fcm=apply(dat,1,function(x,grp) {median(x[grp==grpUniq[2]],na.rm=T)-median(x[grp==grpUniq[1]],na.rm=T)},grp=as.integer(grp)[j])
            } else {
                log2fc=apply(expr[,j],1,function(x,grp) {mean(x[grp==grpUniq[2]],na.rm=T)-mean(x[grp==grpUniq[1]],na.rm=T)},grp=as.integer(grp)[j])
                log2fcm=apply(expr[,j],1,function(x,grp) {median(x[grp==grpUniq[2]],na.rm=T)-median(x[grp==grpUniq[1]],na.rm=T)},grp=as.integer(grp)[j])
            }
        } else {
            log2fc=log2fcm=stat2[,colIdEst]
        }
        i=which(stat2[,colIdPV[2]]<pThres)
        table(sign(log2fcm[i]),sign(log2fc[i]))

        subDir=""
        subDir=paste(modelFlag,"_",metabFlag,sep="")
        if (subDir!="" & !file.exists(subDir)) dir.create(file.path(subDir))
        datadir=paste(subDir,"/",sep="")
        featVec=stat2$metabolite[which(stat2[,colIdPV[2]]<pThres)]
        i=which(stat2[,colIdPV[2]]<pThres)
        #i=i[which(sign(log2fcm[i])!=sign(log2fc[i]))]
        featVec=stat2$metabolite[i]
        #featVec=featVec[1]

        for (featId in featVec) {
            i=which(ann2$id==featId)
            grp=rep("",nrow(phen2))
            j=1:nrow(phen2)
            for (vId in 1:length(varList)) {
                grp=paste(grp,phen2[,varList[vId]])
                j=j[which(!is.na(phen2[j,varList[vId]]))]
            }
            grp=grp[j]
            grpUniq=unique(sort(grp))
            #grpUniqAll=paste(rep(0:1,each=2)," ",0:1,sep="")
            grpUniqAll=grpUniq
            grpUniq=unique(sort(grp))
            ttl=grpUniq
            ttl=ttl[match(grpUniq,grpUniqAll)]
            ttl=sub("0","low",ttl)
            ttl=sub("1","high",ttl)
            x=table(grp)
            ttl=paste(ttl," (",x,")",sep="")
            header=paste(ann2$type[i],": ",ann2$id[i],ifelse(metabFlag=="lipid",ifelse(ann2$annotation[i]=="","",paste(", ",ann2$annotation[i],sep="")),""),", ",ifelse(length(grep("Interaction",colIdEst)),"Interaction - "),"log2FC ",round(log2fc[i],2),", BY-adjusted-pvalue ",signif(stat2[i,colIdPV[2]],2),sep="")
            png(paste(datadir,"boxplot_",featId,".png",sep=""),width=3*480,height=2*480)
            boxplot(expr[i,j]~grp,names=ttl,main=header,ylab="log2(metabolite expression)",las=0,cex=1.5,cex.main=2,cex.lab=1.5,cex.axis=1.5)
            for (gId in 1:length(grpUniq)) {
                jj=which(grp==grpUniq[gId])
                points(rep(gId,length(jj)),expr[i,j][jj],col="grey",cex=2)
                lines(gId+c(-.4,.4),rep(mean(expr[i,j][jj],na.rm=T),2),col="green")
            }
            #abline(h=1,lty="dotted")
            #legend(3.75,4,legend="mean",lty="solid",col="green")
            dev.off()
        }

        if (F) {
            wilcox_test(expr[i,j]~as.factor(2-as.integer(as.factor(grp))),distribution="exact")
            summary(lm(expr[i,j]~grp))
            stat2[i,]
            c(log2fc[i],log2fcm[i])
            mean(expr[i,which(grp=="PTSD")],na.rm=T)-mean(expr[i,which(grp=="Control")],na.rm=T)
        }
    }
}

####################################################################
####################################################################
## Comparison of stat tables

datadir="results/final/wilcoxon/table/"
stat1=read.table(paste(datadir,"stat_wilcox_metabResp_ptsd_lipid.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
stat2=read.table(paste(datadir,"stat_wilcox_metabResp_sex_lipid.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
pThres=0.05
table(stat1$BY_ptsd<pThres,stat2$BY_maleVfemale<pThres)

## ---------------
metabFlag="lipid"
modelFlag="ptsd"

for (metabFlag in unique(ann$type)) {
    for (modelFlag in c("ptsd","sex")) {
        cat("\n\n=============== ",metabFlag,", ",modelFlag,"\n",sep="")
        fName=paste("_wilcox_metabResp_",modelFlag,"_",metabFlag,sep="")
        datadir="results/metabO/imputed_allMetabImp/"
        stat1=read.table(paste(datadir,"stat",fName,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
        ttl="imputed_allMetabImp"
        datadir="results/metabT/raw/"
        datadir="results/metabO/raw/"
        datadir="results/metabO/imputed/"
        datadir="results/metabO/raw_allMetab/"
        stat2=read.table(paste(datadir,"stat",fName,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
        ttl=c(ttl,"raw_allMetab")

        varFlag=modelFlag
        if (modelFlag=="sex") varFlag="maleVfemale"
        colId=paste("pv_",varFlag,sep="")
        colId=paste("log2fc_",varFlag,sep="")
        colIdPV=paste(c("pv","pv"),"_",varFlag,sep=""); pThres=0.001
        colIdPV=paste(c("pv","BY"),"_",varFlag,sep=""); pThres=0.05
        
        cat("Signif in stat1 but not used in stat2: ",sum(!stat1$metabolite%in%stat2$metabolite & stat1[,colIdPV[2]]<pThres,na.rm=T),"\n",sep="")
        cat("Signif in stat2 but not used in stat1: ",sum(!stat2$metabolite%in%stat1$metabolite & stat2[,colIdPV[2]]<pThres,na.rm=T),"\n",sep="")
        
        i=match(stat1$metabolite,stat2$metabolite); i1=which(!is.na(i)); i2=i[i1]
        stat1=stat1[i1,]; stat2=stat2[i2,]

        png(paste("plot_",ttl[2],"V",ttl[1],fName,".png",sep=""))
        lim=range(c(stat1[,colId],stat2[,colId]),na.rm=T)
        header=paste(metabFlag,": ",modelFlag,"\nSignif ",ttl[1],"=",sum(stat1[,colIdPV[2]]<pThres,na.rm=T),", ",ttl[2],"=",sum(stat2[,colIdPV[2]]<pThres,na.rm=T),", both=",sum(stat1[,colIdPV[2]]<pThres & stat2[,colIdPV[2]]<pThres,na.rm=T),sep="")
        plot(stat1[,colId],stat2[,colId],xlim=lim,ylim=lim,main=header,xlab=paste(ttl[1],": log2(fold change)",sep=""),ylab=paste(ttl[2],": log2(fold change)",sep=""),col="grey")
        ii=which(stat1[,colIdPV[2]]<pThres)
        points(stat1[ii,colId],stat2[ii,colId],col="green")
        ii=which(stat2[,colIdPV[2]]<pThres)
        points(stat1[ii,colId],stat2[ii,colId],col="blue")
        ii=which(stat1[,colIdPV[2]]<pThres & stat2[,colIdPV[2]]<pThres)
        points(stat1[ii,colId],stat2[ii,colId],col="red")
        abline(c(0,1),lty="dotted")
        abline(h=0,lty="dotted")
        abline(v=0,lty="dotted")
        dev.off()
        colId=paste("BY_",varFlag,sep="")
        colId=paste("pv_",varFlag,sep="")
        table(stat1[,colId]<pThres,stat2[,colId]<pThres)
    }
}

i=which(stat1[,colId]<pThres & stat2[,colId]>=pThres)
i=which(stat1[,colId]>=pThres & stat2[,colId]<pThres)
i=which(stat1[,colId]<pThres & stat2[,colId]<pThres)
cbind(stat1[i,1:2],stat2[i,2])

####################################################################
####################################################################
## Are there loci with all effects significant?

pThres=0.05
pThres2=99

metabList=unique(ann$type)
compList=c("ptsd+sex","ptsd*sex","diary_tst+sex","psg_tst+sex","ln_delta_nrem+sex","diary_tst*sex","psg_tst*sex","ln_delta_nrem*sex")

#datadir="results/pv1/table/"
datadir="results/final/voom/table/"

for (metabFlag in metabList) {
    for (modelFlag in compList) {
        if (modelFlag%in%c("ptsd+sex","ptsd*sex") & pThres2<1) {
            fName3=paste("_",respFlag,"Resp_",sub("*","X",sub("+","_",modelFlag,fixed=T),fixed=T),"_pv",pThres2,"_",metabFlag,sep="")
        } else {
            fName3=paste("_",respFlag,"Resp_",sub("*","X",sub("+","_",modelFlag,fixed=T),fixed=T),"_",metabFlag,sep="")
        }
        stat2=read.table(paste(datadir,"stat",fName3,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
        k=grep("BY_",names(stat2))
        for (k1 in 1:(length(k)-1)) {
            for (k2 in (k1+1):length(k)) {
                if (any(stat2[,k[k1]]<pThres & stat2[,k[k2]]<pThres)) {
                    cat("\n\n============ ",metabFlag,": ",modelFlag,"\n",sep="")
                    cat(names(stat2)[k][c(k1,k2)],"\n")
                    x1=table(stat2[,k[k1]]<pThres,stat2[,k[k2]]<pThres)
                    print(x1)
                }
            }
        }
    }
}

####################################################################
####################################################################
## Heatmap
library(coin)

datadir1="results/heatmap/allMetabolites/"
datadir1="results/heatmap/25percMostVarMetabolites/"
datadir1="heatmap/allMetabolites/"
datadir1="results/final/heatmap/"
for (datadir2 in dir(datadir1)[-grep("legend",dir(datadir1))]) {
    cat("\n\n",datadir2,"\n",sep="")
    datadir=paste(datadir1,datadir2,"/",sep="")
    fileList=dir(datadir)
    fileList=fileList[grep("clusterInfoSample_",fileList)]
    tbl=read.table(paste(datadir,fileList,sep=""), sep="\t", h=T, quote="", comment.char="",as.is=T,fill=T)
    if (F) {
        grp=paste(tbl$ptsd,tbl$sex)
        pv=fisher.test(grp,tbl$clustId)$p.value
        cat("PV ",signif(pv,2))
    }
    pv=fisher.test(tbl$ptsd,tbl$clustId)$p.value
    cat("PTSD: PV ",signif(pv,2),"\n")
    pv=fisher.test(tbl$sex,tbl$clustId)$p.value
    cat("Gender: PV ",signif(pv,2),"\n")
    for (varId in c("diary_tst","psg_tst","ln_delta_nrem")) {
        pv=pvalue(wilcox_test(tbl[,varId]~as.factor(tbl$clustId),distribution="exact"))
        cat(varId,": PV ",signif(pv,2),"\n")
    }
}
"
25percMostVarMetab_lipid
PTSD: PV  0.0087
Gender: PV  1
diary_tst : PV  0.006
psg_tst : PV  0.26
ln_delta_nrem : PV  0.65


25percMostVarMetab_primary
PTSD: PV  1
Gender: PV  0.67
diary_tst : PV  0.12
psg_tst : PV  0.5
ln_delta_nrem : PV  0.73


25percMostVarMetab_steroid
PTSD: PV  0.83
Gender: PV  0.83
diary_tst : PV  0.99
psg_tst : PV  0.36
ln_delta_nrem : PV  0.48
"

####################################################################
####################################################################
"
Paper: Putative Neuroprotective and Neurotoxic Kynurenine Pathway Metabolites Are Associated with Hippocampal and Amygdalar Volumes in Subjects with Major Depressive Disorder
MSTUS normalization does not make tryptophan associated with PTSD
"

candGene=data.frame(symbol=c("IL-1RA","QA","KYN","KA","TRP","3HK","BDNF"),name=c("Interleukin-1 Receptor Antagonist","Quinolinic Acid","Kynurenine","Kynurenic Acid","Tryptophan","3-hydroxy-kynurenine","Brain-derived Neurotrophic Factor"),stringsAsFactors=F)

colId="name"
table(toupper(candGene[,colId])%in%toupper(ann$id))
table(toupper(candGene[,colId])%in%toupper(ann$annotation))
length(grep(paste(toupper(candGene[,colId]),collapse="|"),toupper(ann$id)))
length(grep(paste(toupper(candGene[,colId]),collapse="|"),toupper(ann$annotation)))

idList=c("tryptophan","creatinine")
idList=c("tryptophan")
for (id in idList) {
    i=which(ann$id==id)
    grp=paste(phen$ptsd,"/",phen$sex,sep="")
    png(paste("boxplot_",id,".png",sep=""))
    y=metabR[i,]
    y=metabI[i,]
    #y=metabR[i,]-metabR[which(ann$id=="creatinine"),]
    boxplot(y~grp,main=id,ylab="log2(Metabolite expression)")
    dev.off()
    fit=lm(y~ptsd*sex,data=phen)
    res=summary(fit)
    print(res$coef)
}


####################################################################
####################################################################
