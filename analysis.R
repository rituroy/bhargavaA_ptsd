dirSrc="/Users/royr/UCSF/"
dirSrc2=dirSrc
setwd(paste(dirSrc2,"BhargavaA",sep=""))

##############################################
source(paste(dirSrc,"functions/miscFuncs.1.3.R",sep=""))

##############################################
## Section 1

## --------------------

datadir="docs/"
phen=read.table(paste(datadir,"gender for metabolomics.txt",sep=""), sep="\t", h=T, quote="", comment.char="",as.is=T,fill=T)
tbl1=read.table(paste(datadir,"mx 315077 Aditi Bhargava_human serum_04-2017_submit TRANSPOSED.txt",sep=""), sep="\t", h=T, quote="", comment.char="",as.is=T,fill=T)
tbl2=read.table(paste(datadir,"mx 315178 Bhargava_human plasma_lipidomics_CSH-QTOF MS_05-2017_submit TRANSPOSE.txt",sep=""), sep="\t", h=T, quote="", comment.char="",as.is=T,fill=T)
tbl22=read.table(paste(datadir,"mx 315178 Bhargava_human plasma_lipidomics_CSH-QTOF MS_05-2017_submit.txt",sep=""), sep="\t", h=T, quote="", comment.char="",as.is=T,fill=T,skip=6)
tbl3=read.table(paste(datadir,"mx 315279 Bhargava bile acids_steroids_12-2017_submit TRANSPOSED.txt",sep=""), sep="\t", h=T, quote="", comment.char="",as.is=T,fill=T)

names(phen)[match(c("subjectid","gender"),names(phen))]=c("id","sex")
tbl22=tbl22[,1:8]
names(tbl22)[match(c("Identifier","Annotation","InChI.Key","Species","count","ESI.mode","m.z","RT"),names(tbl22))]=c("identifier","annotation","inChIkey","species","count","esiMode","mz","rt")
tbl22$id=tbl22$identifier
tbl22$id[which(tbl22$identifier=="5.67_821.60")[2]]="5.67_821.60.1"
tbl22$id=paste("X",tbl22$id,sep="")
colId2=which(names(tbl22)!="id")

tbl=tbl1
j=match(tbl[,1],phen$id); j1=which(!is.na(j)); j2=j[j1]; table(is.na(j))
phen=phen[j2,]
tbl=tbl[j1,]
phen$ptsd=tbl$Group
colId=13:ncol(tbl)
#ann=data.frame(id=names(tbl)[colId],type="primary",stringsAsFactors=F)
tmp=matrix("",nrow=length(colId),ncol=length(colId2)); colnames(tmp)=names(tbl22)[colId2]
ann=data.frame(id=names(tbl)[colId],tmp,type="primary",stringsAsFactors=F)
metab=t(as.matrix(tbl[,colId]))

tbl=tbl2
j=match(tbl[,1],phen$id); j1=which(!is.na(j)); j2=j[j1]; table(is.na(j))
phen=phen[j2,]
tbl=tbl[j1,]
colId=6:ncol(tbl)
#ann=rbind(ann,data.frame(id=names(tbl)[colId],type="lipid",stringsAsFactors=F))
tmp=as.matrix(tbl22[,colId2])
ann=rbind(ann,data.frame(id=names(tbl)[colId],tmp,type="lipid",stringsAsFactors=F))
metab=rbind(metab,t(as.matrix(tbl[,colId])))

tbl=tbl3
j=match(tbl[,1],phen$id); j1=which(!is.na(j)); j2=j[j1]; table(is.na(j))
phen=phen[j2,]
tbl=tbl[j1,]
colId=3:ncol(tbl)
#ann=rbind(ann,data.frame(id=names(tbl)[colId],type="steroid",stringsAsFactors=F))
tmp=matrix("",nrow=length(colId),ncol=length(colId2)); colnames(tmp)=names(tbl22)[colId2]
ann=rbind(ann,data.frame(id=names(tbl)[colId],tmp,type="steroid",stringsAsFactors=F))
metab=rbind(metab,t(as.matrix(tbl[,colId])))

ann$count=as.integer(ann$count)
ann$rt=as.numeric(ann$rt)

phen$id=paste("X",phen$id,sep="")
colnames(metab)=phen$id
metab=apply(metab,c(1,2),as.numeric)

xlim=NULL; ylim=NULL
xlim=range(metab,na.rm=T); ylim=NULL
ylim=c(0,10^-2)
grp=ann$type
grpUniq=sort(unique(grp))
for (gId in 1:length(grpUniq)) {
    i=which(grp==grpUniq[gId])
    xlim=range(metab[i,],na.rm=T)
    xlim=quantile(metab[i,],na.rm=T,probs=c(0,.975))
    png(paste("densityPlot_",grpUniq[gId],".png",sep=""))
    plot(density(metab[i[1],],na.rm=T),main=grpUniq[gId],xlim=xlim,ylim=ylim)
    for (ii in i) {
        if (sum(!is.na(metab[ii,]))<2) next
        lines(density(metab[ii,],na.rm=T))
    }
    dev.off()
}

x=apply(metab,1,function(x) mean(!is.na(x)))

metabO=metab

#library(edgeR)
#lcpmI=cpm(2^metab,log=TRUE)
#maxCntVecCpm=maxCntVec
#lcpmI=cpm(2^metab,log=TRUE)


## --------------------
"
Normalization strategies for metabonomic analysis of urine samples
Bethanne M. Warracka,b,∗, Serhiy Hnatyshyna,b, Karl-Heinz Ott a,c, Michael D. Reilya,b,
Mark Sanders a,b,2, Haiying Zhanga,b,1, Dieter M. Drexler a,d

Finally, we introduce the concept of MS
“total useful signal” (MSTUS) [26], which uses the total intensity
of components that are common to all samples, thus avoiding
xenobiotics and artifacts that would not be appropriate measures
of urine concentration. This latter approach is similar to the
common practice used in proton NMR-based metabonomics analyses
wherein each spectrum is normalized to the total integrated
proton signal, after excluding regions corresponding to xenobiotics,
internal standards and artifact-prone water and urea regions
[27,28].
"

metabT=metabO
grpUniq=unique(ann$type)
for (gId in 1:length(grpUniq)) {
    i1=which(ann$type==grpUniq[gId])
    cat(grpUniq[gId],": ",length(i1))
    i=i1[apply(metabO[i1,],1,function(x) mean(!is.na(x))==1)]
    cat(", ",length(i))
    x=apply(metabO[i,],2,sum)
    cat(", ",sum(is.na(x)),"\n")
    for (j in 1:ncol(metabT)) {
        metabT[i1,j]=metabT[i1,j]/x[j]
    }
}

## --------------------
library(impute)

normList=c("_tus","")
for (normFlag in normList) {
    if (normFlag=="_tus") {
        metab=metabT
    } else {
        metab=metabO
    }
    x=c(metab)
    metabR=log2(metab+min(x[which(x!=0)])/10)
    
    metabI=metab
    grpUniq=unique(ann$type)
    for (gId in 1:length(grpUniq)) {
        i=which(ann$type==grpUniq[gId])
        metabI[i,]=impute.knn(metabR[i,])$data
    }
    metabI=metabI
    meanMetabI=apply(metabI,1,mean,na.rm=T)
    if (normFlag=="_tus") {
        metabTR=metabR
        metabTI2=metabI
    } else {
        metabOR=metabR
        metabOI2=metabI
    }

    metabI=impute.knn(metabR)$data
    if (normFlag=="_tus") {
        metabTR=metabR
        metabTI=metabI
    } else {
        metabOR=metabR
        metabOI=metabI
    }
}

metabNoMissVec=apply(metabO,1,function(x) sum(!is.na(x)))

## --------------------
x=apply(metabO,1,function(x) mean(!is.na(x)))
i=which(x>=.5)
summary(c(metabOI2[i,])-c(metabOI[i,]))
summary(c(metabTI2[i,])-c(metabTI[i,]))
png("tmp_2.png")
par(mfrow=c(2,2))
plot(c(metabOI2[i,]),c(metabOI[i,]))
plot(c(metabTI2[i,]),c(metabTI[i,]))
dev.off()

## --------------------

datObj=list(metab=metab,metabRaw=metabR,metabImp=metabI,ann=ann,phen=phen)


##############################################
colInfo=data.frame(name1=c("(Intercept)","metabolite","ptsdPTSD","sexMale","metabolite:sexMale","ptsdPTSD:sexMale"),name2=c("intercept","metabolite","ptsd","maleVfemale","metaboliteGenderInteraction","ptsdGenderInteraction"),stringsAsFactors=F)

##############################################
## Section 2
## DE

library(limma)
library(edgeR)
library(qvalue)

adjPFlag="qv"
adjPFlag="BY"

pThres2=0.1
pThres2=0.05
pThres2=0.2
pThres2=99

minCnt=0

#colInfo=data.frame(name1=c("(Intercept)","metabolite","ptsdPTSD","sexMale","metabolite:sexMale","ptsdPTSD:sexMale"),name2=c("intercept","metabolite","ptsd","maleVfemale","metaboliteGenderInteraction","ptsdGenderInteraction"),stringsAsFactors=F)

## --------------------
for (metabFlag in unique(ann$type)) {
    for (respFlag in c("metab")) {
        for (modelFlag in c("ptsd","sex","ptsd+sex","ptsd*sex")) {
            cat("\n\n============= ",metabFlag,", ",modelFlag,"\n")
            modelThis=as.formula(paste("~",modelFlag,sep=""))
            #modelThis=as.formula(paste("~as.factor(ptsd)*as.factor(sex)",sep=""))
            i=1:10
            i=1:nrow(metab)
            i=apply(metab,1,function(x) mean(!is.na(x))>=.5)
            i=which(ann$type==metabFlag & apply(metab,1,function(x) mean(!is.na(x))>=.25))
            i=which(ann$type==metabFlag & apply(metab,1,function(x) mean(!is.na(x))>=.5))
            expr=metabI[i,]
            if (modelFlag%in%c("ptsd+sex","ptsd*sex")) {
                i=which(fitMP$p.value[,2]<pThres2 & fitMS$p.value[,2]<pThres2)
                expr=expr[i,]
            }
            design=model.matrix(modelThis,data=phen)
            colnames(design)=colInfo$name2[match(colnames(design),colInfo$name1)]
            dat <- voom(2^expr,design,save.plot=F)
            fit <- lmFit(dat,design)
            fit <- eBayes(fit)
            fit$adjP=matrix(nrow=nrow(fit$coef),ncol=ncol(fit$coef),dimnames=list(rownames(fit$coef),colnames(fit$coef)))
            for (k in 2:ncol(fit$adjP)) {
                i=which(!is.na(fit$p.value[,k]))
                if (adjPFlag=="qv") {
                    fit$adjP[i,k]=qvalue(fit$p.value[i,k])$qvalues
                } else {
                    fit$adjP[i,k]=p.adjust(fit$p.value[i,k],method=adjPFlag)
                }
            }
            switch(paste(respFlag,"~",modelFlag,sep=""),
                "metab~ptsd"={
                    fitMP=fit
                },
                "metab~sex"={
                    fitMS=fit
                },
                "metab~ptsd+sex"={
                    fitMPS=fit
                },
                "metab~ptsd*sex"={
                    fitMPxS=fit
                }
            )
            
            pThres=0.05
            for (k in 2:ncol(design)) {
                cat("\n",metabFlag,": ",colnames(design)[k],sep="")
                print(table(fit$adjP[,k]<pThres))
            }
            
            #colId=2:ncol(design)
            #top=cbind(data.frame(metabolite=rownames(fit$coef),stringsAsFactors=F),log2fc=round(fit$coef[,colId],2),pv=signif(fit$p.value[,colId],2),adjP=signif(fit$adjP[,colId],2))
            nm=c("metabolite")
            top=data.frame(metabolite=rownames(fit$coef),stringsAsFactors=F)
            if (metabFlag=="lipid") {
                colId2=which(names(ann)!="id")
                nm=c(nm,names(ann)[colId2])
                top=cbind(top,ann[match(top$metabolite,ann$id),colId2])
            }
            for (colId in 2:ncol(design)) {
                nm=c(nm,paste(c("log2fc","pv",adjPFlag),"_",colnames(design)[colId],sep=""))
                top=cbind(top,log2fc=round(fit$coef[,colId],2),pv=signif(fit$p.value[,colId],2),adjP=signif(fit$adjP[,colId],2))
            }
            names(top)=nm
            if (modelFlag%in%c("ptsd+sex","ptsd*sex") & pThres2<1) {
                fName3=paste("_",respFlag,"Resp_",sub("*","X",sub("+","_",modelFlag,fixed=T),fixed=T),"_pv",pThres2,"_",metabFlag,sep="")
            } else {
                fName3=paste("_",respFlag,"Resp_",sub("*","X",sub("+","_",modelFlag,fixed=T),fixed=T),"_",metabFlag,sep="")
            }
            write.table(top,paste("stat",fName3,".txt",sep=""), sep="\t", col.names=T, row.names=F, quote=F)
        }
    }
}

## --------------------
library(coin)

metab=metabT
metabR=metabTR
metabI=metabTI
metabI2=metabTI2

metab=metabO
metabR=metabOR
metabI=metabOI
metabI2=metabOI2

datList=c("imputed")
datList=c("raw_allMetab")
datList=c("orig","imputed","raw","raw_allMetab","imputed_allMetabImp")

for (datFlag in datList) {
    subDir=datFlag
    if (subDir!="" & !file.exists(subDir)) dir.create(file.path(subDir))
    datadir=paste(subDir,"/",sep="")
    for (metabFlag in unique(ann$type)) {
        for (respFlag in c("metab")) {
            for (modelFlag in c("ptsd","sex")) {
                cat("\n\n============= ",datFlag,", ",metabFlag,", ",modelFlag,"\n",sep="")
                modelThis=as.formula(paste("~",modelFlag,sep=""))
                #modelThis=as.formula(paste("~as.factor(ptsd)*as.factor(sex)",sep=""))
                i=1:10
                if (datFlag==c("raw_allMetab")) {
                    i=which(ann$type==metabFlag)
                } else if (datFlag=="raw") {
                    i=which(ann$type==metabFlag & apply(metab,1,function(x) mean(!is.na(x))>=.5))
                }
                if (datFlag%in%c("raw","raw_allMetab")) {
                    expr=metabR[i,]
                } else if (datFlag=="imputed") {
                    expr=metabI[i,]
                } else if (datFlag=="imputed_allMetabImp") {
                    expr=metabI2[i,]
                } else {
                    expr=metab[i,]
                }
                if (modelFlag%in%c("ptsd+sex","ptsd*sex")) {
                    i=which(fitMP$p.value[,2]<pThres2 & fitMS$p.value[,2]<pThres2)
                    expr=expr[i,]
                }
                design=model.matrix(modelThis,data=phen)
                colnames(design)=colInfo$name2[match(colnames(design),colInfo$name1)]
                grp=as.factor(phen[,modelFlag])
                j=which(!is.na(grp))
                tmp=rep(NA,nrow(expr))
                if (datFlag=="orig") {
                    x=c(expr[,j])
                    dat=log2(expr[,j]+min(x[which(x!=0)])/10)
                    log2fc=apply(dat,1,function(x,grp) {mean(x[grp==2],na.rm=T)-mean(x[grp==1],na.rm=T)},grp=as.integer(grp)[j])
                } else {
                    log2fc=apply(expr[,j],1,function(x,grp) {mean(x[grp==2],na.rm=T)-mean(x[grp==1],na.rm=T)},grp=as.integer(grp)[j])
                }
                nm=c("metabolite","numOfSamples")
                top=data.frame(metabolite=rownames(expr),numOfSamples=apply(expr,1,function(x) sum(!is.na(x))),stringsAsFactors=F)
                if (metabFlag=="lipid") {
                    colId2=which(names(ann)!="id")
                    nm=c(nm,names(ann)[colId2])
                    top=cbind(top,ann[match(top$metabolite,ann$id),colId2])
                }
                colId=2
                nm=c(nm,paste(c("log2fc","pv",adjPFlag),"_",colnames(design)[colId],sep=""))
                top=cbind(top,log2fc=log2fc,pv=tmp,adjP=tmp,stringsAsFactors=F)
                for (i in 1:nrow(expr)) {
                    fit=try(wilcox_test(expr[i,]~grp,distribution="exact"))
                    if (class(fit)=="ScalarIndependenceTest") {
                        top$pv[i]=pvalue(fit)
                    }
                }
                top$adjP=p.adjust(top$pv,method=adjPFlag)
                
                pThres=0.05
                k=2
                cat("\n",metabFlag,": ",colnames(design)[k],sep="")
                print(table(top$adjP<pThres))
                names(top)=nm
                fName3=paste("_",respFlag,"Resp_",sub("*","X",sub("+","_",modelFlag,fixed=T),fixed=T),"_",metabFlag,sep="")
                write.table(top,paste(datadir,"stat_wilcox",fName3,".txt",sep=""), sep="\t", col.names=T, row.names=F, quote=F)
            }
        }
    }
}

##############################################
## PCA

#source(paste(dirSrc,"functions/heatmapAcgh.7.3.R",sep=""))
source(paste(dirSrc,"functions/heatmapAcgh.7.4.R",sep=""))

metab=metabO
metabR=metabOR
metabI=metabOI
metabI=metabOI2

testFlag="linear"
testFlag="logistic"

datFlag="raw"
datFlag="imputed"

colorInfo=data.frame(grp=c("Control/Female","PTSD/Female","Control/Male","PTSD/Male"),col=c("orange","red","skyblue","blue"),stringsAsFactors=F)

colVec=rep("black",nrow(phen))
colVec[which(phen$ptsd=="Control" & phen$sex=="Female")]="orange"
colVec[which(phen$ptsd=="PTSD" & phen$sex=="Female")]="red"
colVec[which(phen$ptsd=="Control" & phen$sex=="Male")]="skyblue"
colVec[which(phen$ptsd=="PTSD" & phen$sex=="Male")]="blue"

if (testFlag=="linear") {
    modelList=paste("metabolite~",c("ptsd","sex","ptsd+sex","ptsd*sex"),sep="")
} else {
    modelList=c("ptsd~metabolite","sex~metabolite","ptsd~metabolite+sex","ptsd~metabolite*sex")
}

cutoff=3

statTbl=NULL
for (filtFlag in c("","_25percMostVar")) {
    for (metabFlag in unique(ann$type)) {
        cat("\n\n===================================================\n",sep="")
        cat("===================================================\n\n",sep="")
        fNameOut=paste(filtFlag,"_",metabFlag,sep="")
        header=paste(capWords(metabFlag),": ",ifelse(filtFlag=="","All",paste(as.integer(gsub("_|percMostVar","",filtFlag)),"% most variable",sep=""))," metabolites",sep="")
        if (length(grep("percMostVar",filtFlag))==1) {
            i=which(ann$type==metabFlag)
            i=which(ann$type==metabFlag & apply(metab,1,function(x) mean(!is.na(x))>=.5))
            if (datFlag=="raw") {
                arrayData=metab[i,]
            } else {
                arrayData=metabI[i,]
            }
            x=apply(arrayData,1,sd,na.rm=T)
            i=order(x,decreasing=T)[1:round(quantile(1:nrow(arrayData),probs=as.integer(gsub("_|percMostVar","",filtFlag))/100))]
            arrayData=arrayData[i,]
        } else {
            i=which(ann$type==metabFlag & apply(metab,1,function(x) mean(!is.na(x))>=.5))
            if (datFlag=="raw") {
                arrayData=metab[i,]
            } else {
                arrayData=metabI[i,]
            }
        }
        arrayData=t(arrayData)
        fit=prcomp(arrayData, center=T, scale=T)
        png(paste("screePlot_pca",fNameOut,".png",sep=""))
        screeplot(fit,npcs=length(fit$sdev),main=paste("PCA screeplot\n",header,sep=""))
        #if (plotCutoffFlag) abline(v=cutoff+1,col="red",lty="dotted")
        dev.off()
        if (F) {
            png("tmp.png")
            barplot(fit$sdev,main=paste("PCA screeplot\n",header,sep=""),axes=T,names.arg=paste("PCA ",1:length(fit$sdev)),las=3)
            dev.off()
        }
        
        if (F) {
        png(paste("biPlot_pca",fNameOut,"_%1d.png",sep=""))
        biplot(fit,choices=c(1,2),main=paste("PCA biplot: ",header,sep=""))
        biplot(fit,choices=c(1,3),main=paste("PCA biplot: ",header,sep=""))
        biplot(fit,choices=c(2,3),main=paste("PCA biplot: ",header,sep=""))
        dev.off()
        
        png(paste("rotationPlot_pca",fNameOut,".png",sep=""))
        par(mfcol=c(2,2))
        plot(fit$rotation[,"PC1"],fit$rotation[,"PC2"],main=paste("PCA: ",header,sep=""),xlab="Rotation: PC1",ylab="Rotation: PC2")
        plot(fit$rotation[,"PC1"],fit$rotation[,"PC3"],main=paste("PCA: ",header,sep=""),xlab="Rotation: PC1",ylab="Rotation: PC3")
        plot(fit$rotation[,"PC2"],fit$rotation[,"PC3"],main=paste("PCA: ",header,sep=""),xlab="Rotation: PC2",ylab="Rotation: PC3")
        dev.off()
        }
        
        png(paste("scorePlot_pca",fNameOut,".png",sep=""),width=1.5*480,height=1.5*480)
        par(mfcol=c(2,2))
        plot(fit$x[,"PC1"],fit$x[,"PC2"],main=paste("PCA\n",header,sep=""),xlab="Score: PC1",ylab="Score: PC2",col=colVec)
        plot(fit$x[,"PC1"],fit$x[,"PC3"],main=paste("PCA\n",header,sep=""),xlab="Score: PC1",ylab="Score: PC3",col=colVec)
        plot(fit$x[,"PC2"],fit$x[,"PC3"],main=paste("PCA\n",header,sep=""),xlab="Score: PC2",ylab="Score: PC3",col=colVec)
        sampleColorLegend(tls=colorInfo$grp,col=colorInfo$col,legendTitle=NULL)
        dev.off()

        tbl=cbind(feature=rownames(fit$x),fit$x[,1:cutoff])
        write.table(tbl, paste("prinComp",fNameOut,".txt",sep=""), sep="\t", col.names=T, row.names=F, quote=F)

        for (modelFlag in modelList) {
            cat("\n\n============= ",header,"\n",sep="")
            cat("============= Model: ",modelFlag,"\n\n",sep="")
            dat=cbind(metabolite=fit$x[,1],phen)
            for (k in c("ptsd","sex")) {
                dat[,k]=as.factor(dat[,k])
            }
            if (testFlag=="linear") {
                modelThis=as.formula(modelFlag)
                fit2=lm(modelThis,data=dat)
                res=summary(fit2)$coef
                rownames(res)=colInfo$name2[match(rownames(res),colInfo$name1)]
                colnames(res)=c("coef","stdErr","t","pv")
                tbl=data.frame(basedOn=rep(header,nrow(res)),model=rep(paste("",modelFlag,sep=""),nrow(res)),term=rownames(res),stringsAsFactors=F)
            } else {
                modelThis=as.formula(modelFlag)
                fit2=tryCatch(glm(modelThis, family="binomial",data=dat),error = function(e) e)
                res=summary(fit2)$coef
                rownames(res)=colInfo$name2[match(rownames(res),colInfo$name1)]
                colnames(res)=c("coef","stdErr","z","pv")
                tbl=data.frame(basedOn=rep(header,nrow(res)),model=rep(paste("",modelFlag,sep=""),nrow(res)),term=rownames(res),stringsAsFactors=F)
            }
            tbl=cbind(tbl,res)
            statTbl=rbind(statTbl,tbl[2:nrow(tbl),])
            for (k in which(colnames(res)%in%c("coef","stdErr","t","z"))) {
                res[,k]=round(res[,k],2)
            }
            for (k in c("pv")) {
                res[,k]=signif(res[,k],2)
            }
            print(res)
        }
    }
}
rownames(statTbl)=NULL

pThres=.05
tbl=statTbl
for (k in which(colnames(tbl)%in%c("coef","stdErr","t","z"))) {
    tbl[,k]=round(tbl[,k],3)
}
for (k in c("pv")) {
    tbl[,k]=signif(tbl[,k],3)
}
tbl[which(statTbl$pv<pThres),c("basedOn","model","term","coef","pv")]
"
imputed
basedOn               model       term   coef       pv
8                  Lipid: All metabolites     ptsd~metabolite metabolite  0.044 4.80e-04
10                 Lipid: All metabolites ptsd~metabolite+sex metabolite  0.044 4.79e-04
12                 Lipid: All metabolites ptsd~metabolite*sex metabolite  0.035 4.70e-02
16               Steroid: All metabolites      sex~metabolite metabolite -0.381 8.79e-05
29   Lipid: 25% most variable metabolites     ptsd~metabolite metabolite  0.087 5.31e-03
31   Lipid: 25% most variable metabolites ptsd~metabolite+sex metabolite  0.088 5.09e-03
37 Steroid: 25% most variable metabolites      sex~metabolite metabolite  1.495 8.23e-07

imputed_allMetabImp
basedOn               model       term   coef       pv
8                    Lipid: All metabolites     ptsd~metabolite metabolite  0.044 4.80e-04
10                   Lipid: All metabolites ptsd~metabolite+sex metabolite  0.044 4.79e-04
12                   Lipid: All metabolites ptsd~metabolite*sex metabolite  0.035 4.68e-02
16                 Steroid: All metabolites      sex~metabolite metabolite -0.381 8.76e-05
29     Lipid: 25% most variable metabolites     ptsd~metabolite metabolite  0.088 5.20e-03
31     Lipid: 25% most variable metabolites ptsd~metabolite+sex metabolite  0.088 5.00e-03
37   Steroid: 25% most variable metabolites      sex~metabolite metabolite  1.495 8.23e-07
"

fName=paste("stat_prinCompMetabolite.txt",sep="")
tbl2=tbl[grep("25% most variable metabolites",tbl$basedOn),]
grp=paste(tbl2$basedOn,tbl2$model)
grpUniq=unique(grp)
tbl1=NULL
write.table("Metabolite meta variable from PCA\n",fName, sep="\t", col.names=F, row.names=F, quote=F)
for (gId in 1:length(grpUniq)) {
    k=which(grp==grpUniq[gId])
    if (length(k)>1) {
        tbl2[k[2:length(k)],"basedOn"]=""
        tbl2[k[2:length(k)],"model"]=""
    }
    tbl1=tbl2[k,]
    write.table(tbl1,fName, sep="\t", col.names=T, row.names=F, quote=F,append=T)
    write.table("",fName, sep="\t", col.names=F, row.names=F, quote=F,append=T)
}

##############################################
