#  total sleep time (tst) measured by self report diary and by polysomnography, and delta sleep from polysomnography.

dirSrc="/Users/royr/UCSF/"
dirSrc2=dirSrc
setwd(paste(dirSrc2,"BhargavaA/ptsd",sep=""))

##############################################
source(paste(dirSrc,"functions/miscFuncs.1.3.R",sep=""))

##############################################
## Section 1


## --------------------
datadir="docs/"

tbl3=read.table(paste(datadir,"mx 315279 Bhargava_bile acids steroids_human serum_03-2018_submit (1) - Data Submit.txt",sep=""), sep="\t", h=T, quote="", comment.char="",as.is=T,fill=T,skip=7)
tbl32=read.table(paste(datadir,"mx 315279 Bhargava_bile acids steroids_human serum_03-2018_submit (1) - Data Submit.txt",sep=""), sep="\t", h=T, quote="", comment.char="",as.is=T,fill=T,nrow=10)
x=unlist(tbl32[3,],use.names=F)
k=grep("label",x)+1
names(tbl3)[k:ncol(tbl3)]=x[k:ncol(tbl3)]
tbl31=tbl3

tbl3=read.table(paste(datadir,"May4-ABedits- mx 315279 Bhargava_bile acids steroids_human serum_03-2018_submit (1) - Data Submit.txt",sep=""), sep="\t", h=T, quote="", comment.char="",as.is=T,fill=T,skip=7)
tbl32=read.table(paste(datadir,"May4-ABedits- mx 315279 Bhargava_bile acids steroids_human serum_03-2018_submit (1) - Data Submit.txt",sep=""), sep="\t", h=T, quote="", comment.char="",as.is=T,fill=T,nrow=10)
x=unlist(tbl32[3,],use.names=F)
k=grep("label",x)+1
names(tbl3)[k:ncol(tbl3)]=x[k:ncol(tbl3)]
tbl32=tbl3

k=match(names(tbl31),names(tbl32)); k1=which(!is.na(k)); k2=k[k1]
for (k in 1:length(k1)) {
    if (any(tbl31[,k1[k]]!=tbl32[,k2[k]],na.rm=T) | any(is.na(tbl31[,k1[k]])!=is.na(tbl32[,k2[k]]))) print(k)
}

## --------------------

datadir="docs/"
phen=read.table(paste(datadir,"gender for metabolomics.txt",sep=""), sep="\t", h=T, quote="", comment.char="",as.is=T,fill=T)
ann0=read.table(paste(datadir,"NIH West Coast Metabolomics Center_UC Davis_ 1083 compounds platforms list 11-2012.txt",sep=""), sep="\t", h=T, quote="", comment.char="",as.is=T,fill=T,skip=1)
tbl1=read.table(paste(datadir,"mx 315077 Aditi Bhargava_human serum_04-2017_submit TRANSPOSED.txt",sep=""), sep="\t", h=T, quote="", comment.char="",as.is=T,fill=T)
tbl2=read.table(paste(datadir,"mx 315178 Bhargava_human plasma_lipidomics_CSH-QTOF MS_05-2017_submit TRANSPOSE.txt",sep=""), sep="\t", h=T, quote="", comment.char="",as.is=T,fill=T)
tbl22=read.table(paste(datadir,"mx 315178 Bhargava_human plasma_lipidomics_CSH-QTOF MS_05-2017_submit.txt",sep=""), sep="\t", h=T, quote="", comment.char="",as.is=T,fill=T,skip=6)
##tbl3=read.table(paste(datadir,"mx 315279 Bhargava bile acids_steroids_12-2017_submit TRANSPOSED.txt",sep=""), sep="\t", h=T, quote="", comment.char="",as.is=T,fill=T)
#tbl3=read.table(paste(datadir,"mx 315279 Bhargava_bile acids steroids_human serum_03-2018_submit (1) - Data Submit.txt",sep=""), sep="\t", h=T, quote="", comment.char="",as.is=T,fill=T,skip=7)
#tbl32=read.table(paste(datadir,"mx 315279 Bhargava_bile acids steroids_human serum_03-2018_submit (1) - Data Submit.txt",sep=""), sep="\t", h=T, quote="", comment.char="",as.is=T,fill=T,nrow=10)
tbl3=read.table(paste(datadir,"May4-ABedits- mx 315279 Bhargava_bile acids steroids_human serum_03-2018_submit (1) - Data Submit.txt",sep=""), sep="\t", h=T, quote="", comment.char="",as.is=T,fill=T,skip=7)
tbl32=read.table(paste(datadir,"May4-ABedits- mx 315279 Bhargava_bile acids steroids_human serum_03-2018_submit (1) - Data Submit.txt",sep=""), sep="\t", h=T, quote="", comment.char="",as.is=T,fill=T,nrow=10)
tbl3$type[which(tbl3$type=="internal statndard")]="internal standard"
x=unlist(tbl32[3,],use.names=F)
k=grep("label",x)+1
names(tbl3)[k:ncol(tbl3)]=x[k:ncol(tbl3)]

if (F) {
    phen2=read.table(paste(datadir,"Sleep R01 Metabolism Paper Database.txt",sep=""), sep="\t", h=T, quote="", comment.char="",as.is=T,fill=T)
    phen3=read.table(paste(datadir,"GTT.txt",sep=""), sep="\t", h=T, quote="", comment.char="",as.is=T,fill=T)
    phen4=read.table(paste(datadir,"delta_acth_means_wide.txt",sep=""), sep="\t", h=T, quote="", comment.char="",as.is=T,fill=T)
}
phen2=read.table(paste(datadir,"sleep_metab.txt",sep=""), sep="\t", h=T, quote="", comment.char="",as.is=T,fill=T)
names(phen2)[match(c("subjectid","groupnum"),names(phen2))]=c("id","ptsd")

names(phen)[match(c("subjectid","gender"),names(phen))]=c("id","sex")
tbl22=tbl22[,1:8]
#ann2=tbl22
names(tbl22)[match(c("Identifier","Annotation","InChI.Key","Species","count","ESI.mode","m.z","RT"),names(tbl22))]=c("identifier","annotation","inChIkey","species","count","esiMode","mz","rt")
tbl22$id=tbl22$identifier
tbl22$id[which(tbl22$identifier=="5.67_821.60")[2]]="5.67_821.60.1"
tbl22$id=paste("X",tbl22$id,sep="")
colId2=which(names(tbl22)!="id")

tbl=tbl1
j=match(tbl[,1],phen$id); j1=which(!is.na(j)); j2=j[j1]; table(is.na(j))
phen=phen[j2,]
tbl=tbl[j1,]
#phen$ptsd=tbl$Group
colId=13:ncol(tbl)
#ann=data.frame(id=names(tbl)[colId],type="primary",stringsAsFactors=F)
tmp=matrix("",nrow=length(colId),ncol=length(colId2)); colnames(tmp)=names(tbl22)[colId2]
ann=data.frame(id=names(tbl)[colId],tmp,type="primary",stringsAsFactors=F)
metab=t(as.matrix(tbl[,colId]))

phen=cbind(phen,phen2[match(phen$id,phen2$id),!names(phen2)%in%names(phen)])
phen$ptsd[which(phen$ptsd=="PTSD+")]="PTSD"

tbl=tbl2
j=match(tbl[,1],phen$id); j1=which(!is.na(j)); j2=j[j1]; table(is.na(j))
phen=phen[j2,]
tbl=tbl[j1,]
colId=6:ncol(tbl)
#ann=rbind(ann,data.frame(id=names(tbl)[colId],type="lipid",stringsAsFactors=F))
tmp=as.matrix(tbl22[,colId2])
ann=rbind(ann,data.frame(id=names(tbl)[colId],tmp,type="lipid",stringsAsFactors=F))
metab=rbind(metab,t(as.matrix(tbl[,colId])))

if (F) {
    tbl=tbl3
    j=match(tbl[,1],phen$id); j1=which(!is.na(j)); j2=j[j1]; table(is.na(j))
    phen=phen[j2,]
    tbl=tbl[j1,]
    colId=3:ncol(tbl)
    #ann=rbind(ann,data.frame(id=names(tbl)[colId],type="steroid",stringsAsFactors=F))
    tmp=matrix("",nrow=length(colId),ncol=length(colId2)); colnames(tmp)=names(tbl22)[colId2]
    ann=rbind(ann,data.frame(id=names(tbl)[colId],tmp,type="steroid",stringsAsFactors=F))
    metab=rbind(metab,t(as.matrix(tbl[,colId])))
}
tbl=tbl3
j=match(names(tbl),phen$id); j1=which(!is.na(j)); j2=j[j1]; table(is.na(j))
ann3=tbl[,is.na(j)]
names(ann3)[match(c("name","InChIKey","Human.Metabolome.DB","KEGG","LipidMAPS","PubChem.CID","short.key","acc.mass..neutral.","formula","acc.mass..M.H...resp...M.H..","ESI.mode","MRM","ret.time","DP","CE","EP.Entrance.Potential","CXP.Collision.Cell.Exil.Potential","LOD..nM.","LOQ..nM."),names(ann3))]=
c("identifier","inChIkey","Human.Metabolome.DB","KEGG","LipidMAPS","PubChem.CID","shortKey","acc.mass..neutral.","formula","acc.mass..M.H...resp...M.H..","ESI.mode","MRM","ret.time","DP","CE","EP.Entrance.Potential","CXP.Collision.Cell.Exil.Potential","LOD..nM.","LOQ..nM.")
ann3$id=sapply(ann3$identifier,function(x) {
    y=tolower(gsub("_+",".",gsub("-| ","_",sub(" +^","",sub(" +$","",x)))))
    if (length(grep(pattern="[[:digit:]]", x=substr(y,1,1)))==1) y=paste("X",y,sep="")
    y
},USE.NAMES=F)
grpUniq=unique(ann3$id[duplicated(ann3$id)])
for (gId in 1:length(grpUniq)) {
    i=which(ann3$id==grpUniq[gId])
    ann3$id[i]=paste(ann3$id[i],"_",1:length(i),sep="")
}
#ann3$id
rownames(tbl)=ann3$id
phen=phen[j2,]
tbl=tbl[,j1]
tmp=as.data.frame(matrix("",nrow=nrow(ann3),ncol=ncol(ann)),stringsAsFactors=F); colnames(tmp)=names(ann)
k=match(names(ann3),names(tmp)); k1=which(!is.na(k)); k2=k[k1]
tmp[,k2]=ann3[,k1]
#tmp$type="steroid"
ann=rbind(ann,tmp)
rm(ann3,tmp)
metab=rbind(metab,as.matrix(tbl))

ann$count=as.integer(ann$count)
ann$rt=as.numeric(ann$rt)

phen$id=paste("X",phen$id,sep="")
rownames(phen)=phen$id
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
    if (any(is.na(xlim))) next
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

if (F) {
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

    for (gId in 1:length(grpUniq)) {
        i1=which(ann$type==grpUniq[gId])
        cat(grpUniq[gId],": ",length(i1))
        x=c(metab[i1,])
        cat(", ",mean(x==0,na.rm=T),", ",min(x[which(x!=0)]),"\n")
        print(summary(x))

    }
}

## --------------------
library(impute)

#normList=c("_tus","")
normList=c("")
for (normFlag in normList) {
    if (normFlag=="_tus") {
        metab=metabT
    } else {
        metab=metabO
    }
    x=c(metab)
    metabR=log2(metab+min(x[which(x!=0)])/10)
    metabR1=metabR
    metabR=metab
    metabR[metab==0]=min(x[which(x!=0)])/10
    metabR=log2(metabR)
    metabR2=metabR
    
    if (F) {
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

    x=c(metab)
    metabI=metab
    metabI[is.na(metab) | metab==0]=min(x[which(x!=0)])/10
    metabI=log2(metabI)
}

metabNoMissVec=apply(metabO,1,function(x) sum(!is.na(x)))

## --------------------
if (F) {
    x=apply(metabO,1,function(x) mean(!is.na(x)))
    i=which(x>=.5)
    summary(c(metabOI2[i,])-c(metabOI[i,]))
    summary(c(metabTI2[i,])-c(metabTI[i,]))
    png("tmp_2.png")
    par(mfrow=c(2,2))
    plot(c(metabOI2[i,]),c(metabOI[i,]))
    plot(c(metabTI2[i,]),c(metabTI[i,]))
    dev.off()
}

## --------------------

datObj=list(metab=metab,metabRaw=metabR,metabImp=metabI,ann=ann,phen=phen)


##############################################
colInfo=data.frame(name1=c("(Intercept)","metabolite","ptsdPTSD","sexMale","metabolite:sexMale","ptsdPTSD:sexMale","diary_tst","psg_tst","ln_delta_nrem","ptsdPTSD:diary_tst","ptsdPTSD:psg_tst","ptsdPTSD:ln_delta_nrem","diary_tst:sexMale","psg_tst:sexMale","ln_delta_nrem:sexMale","metabolite:diary_tst","metabolite:psg_tst","metabolite:ln_delta_nrem"),name2=c("intercept","metabolite","ptsd","maleVfemale","metaboliteSexInteraction","ptsdSexInteraction","diary_tst","psg_tst","ln_delta_nrem","ptsdDiary_tstInteraction","ptsdPsg_tstInteraction","ptsdLn_delta_nremInteraction","diary_tstSexInteraction","psg_tstSexInteraction","ln_delta_nremSexInteraction","metaboliteDiary_tstInteraction","metabolitePsg_tstInteraction","metaboliteLn_delta_nremInteraction"),stringsAsFactors=F)

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

metabI=datObj$metabImp

modelList=c("ptsd","sex","ptsd+sex","ptsd*sex","diary_tst","psg_tst","ln_delta_nrem")
modelList=c("diary_tst+sex","psg_tst+sex","ln_delta_nrem+sex")
modelList=c("ptsd","sex","ptsd+sex","ptsd*sex","diary_tst","psg_tst","ln_delta_nrem","diary_tst+sex","psg_tst+sex","ln_delta_nrem+sex","diary_tst*sex","psg_tst*sex","ln_delta_nrem*sex")
modelList=c("ptsd+sex","ptsd*sex","diary_tst","psg_tst","ln_delta_nrem","ptsd+diary_tst","ptsd+psg_tst","ptsd+ln_delta_nrem","ptsd*diary_tst","ptsd*psg_tst","ptsd*ln_delta_nrem","diary_tst+sex","psg_tst+sex","ln_delta_nrem+sex","diary_tst*sex","psg_tst*sex","ln_delta_nrem*sex")
modelList=c("ptsd","sex","ptsd+sex","ptsd*sex","diary_tst","psg_tst","ln_delta_nrem","ptsd+diary_tst","ptsd+psg_tst","ptsd+ln_delta_nrem","ptsd*diary_tst","ptsd*psg_tst","ptsd*ln_delta_nrem","diary_tst+sex","psg_tst+sex","ln_delta_nrem+sex","diary_tst*sex","psg_tst*sex","ln_delta_nrem*sex")

metabList=unique(ann$type); metabList=metabList[which(!metabList%in%c("internal standard"))]

modelList=c("sex","diary_tst","diary_tst+sex","diary_tst*sex")
metabList="lipid"

## --------------------
for (metabFlag in metabList) {
    for (respFlag in c("metab")) {
        for (modelFlag in modelList) {
            cat("\n\n============= ",metabFlag,", ",modelFlag,"\n")
            modelThis=as.formula(paste("~",modelFlag,sep=""))
            i=1:10
            i=1:nrow(metab)
            i=apply(metab,1,function(x) mean(!is.na(x))>=.5)
            i=which(ann$type==metabFlag & apply(metab,1,function(x) mean(!is.na(x))>=.25))
            i=which(ann$type==metabFlag & apply(metab,1,function(x) mean(!is.na(x))>=.5))
            expr=metabI[i,]
            if (modelFlag%in%c("ptsd+sex","ptsd*sex")) {
                i=which(fitMP$p.value[,2]<pThres2 & fitMS$p.value[,2]<pThres2)
                expr=expr[i,]
            } else if (modelFlag%in%c("ptsd+diary_tst","ptsd*diary_tst")) {
                i=which(fitMDT$p.value[,2]<pThres2 & fitMP$p.value[,2]<pThres2)
                expr=expr[i,]
            } else if (modelFlag%in%c("ptsd+psg_tst","ptsd*psg_tst")) {
                i=which(fitMPT$p.value[,2]<pThres2 & fitMP$p.value[,2]<pThres2)
                expr=expr[i,]
            } else if (modelFlag%in%c("ptsd+ln_delta_nrem","ptsd*ln_delta_nrem")) {
                i=which(fitMDR$p.value[,2]<pThres2 & fitMP$p.value[,2]<pThres2)
                expr=expr[i,]
            } else if (modelFlag%in%c("diary_tst+sex","diary_tst*sex")) {
                i=which(fitMDT$p.value[,2]<pThres2 & fitMS$p.value[,2]<pThres2)
                expr=expr[i,]
            } else if (modelFlag%in%c("psg_tst+sex","psg_tst*sex")) {
                i=which(fitMPT$p.value[,2]<pThres2 & fitMS$p.value[,2]<pThres2)
                expr=expr[i,]
            } else if (modelFlag%in%c("ln_delta_nrem+sex","ln_delta_nrem*sex")) {
                i=which(fitMDR$p.value[,2]<pThres2 & fitMS$p.value[,2]<pThres2)
                expr=expr[i,]
            }
            design=model.matrix(modelThis,data=phen)
            colnames(design)=colInfo$name2[match(colnames(design),colInfo$name1)]
            expr=expr[,match(rownames(design),colnames(expr))]
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
                },
                "metab~diary_tst"={
                    fitMDT=fit
                },
                "metab~psg_tst"={
                    fitMPT=fit
                },
                "metab~ln_delta_nrem"={
                    fitMDR=fit
                }
            )
            
            pThres=0.05
            for (k in 2:ncol(design)) {
                cat("\n",metabFlag,": ",colnames(design)[k],sep="")
                print(table(fit$adjP[,k]<pThres))
            }
            
            #colId=2:ncol(design)
            #top=cbind(data.frame(metabolite=rownames(fit$coef),stringsAsFactors=F),log2fc=round(fit$coef[,colId],5),pv=signif(fit$p.value[,colId],2),adjP=signif(fit$adjP[,colId],2))
            nm=c("metabolite")
            top=data.frame(metabolite=rownames(fit$coef),stringsAsFactors=F)
            if (metabFlag=="lipid") {
                colId2=which(names(ann)!="id")
                nm=c(nm,names(ann)[colId2])
                top=cbind(top,ann[match(top$metabolite,ann$id),colId2])
            }
            for (colId in 2:ncol(design)) {
                nm=c(nm,paste(c("log2fc","pv",adjPFlag),"_",colnames(design)[colId],sep=""))
                #top=cbind(top,log2fc=round(fit$coef[,colId],5),pv=signif(fit$p.value[,colId],2),adjP=signif(fit$adjP[,colId],2))
                top=cbind(top,log2fc=fit$coef[,colId],pv=fit$p.value[,colId],adjP=fit$adjP[,colId])
            }
            names(top)=nm
            if (modelFlag%in%c("ptsd+sex","ptsd*sex","diary_tst+sex","psg_tst+sex","ln_delta_nrem+sex") & pThres2<1) {
                fName3=paste("_",respFlag,"Resp_",sub("*","X",sub("+","_",modelFlag,fixed=T),fixed=T),"_pv",pThres2,"_",metabFlag,sep="")
            } else {
                fName3=paste("_",respFlag,"Resp_",sub("*","X",sub("+","_",modelFlag,fixed=T),fixed=T),"_",metabFlag,sep="")
            }
            write.table(top,paste("stat",fName3,".txt",sep=""), sep="\t", col.names=T, row.names=F, quote=F)
        }
    }
}

## --------------------
## Logistic regression

if (F) {
    metabI=datObj$metabImp

    metabList=unique(ann$type); metabList=metabList[which(!metabList%in%c("internal standard"))]

    for (datFlag in datList) {
        subDir=datFlag
        if (subDir!="" & !file.exists(subDir)) dir.create(file.path(subDir))
        datadir=paste(subDir,"/",sep="")
        for (metabFlag in metabList) {
            for (respFlag in c("metab")) {
                for (modelFlag in c("ptsd","sex")) {
                    cat("\n\n============= ",datFlag,", ",metabFlag,", ",modelFlag,"\n",sep="")
                    modelThis=as.formula(paste("~",modelFlag,sep=""))
                    i=1:10
                    if (datFlag==c("raw_allMetab")) {
                        i=which(ann$type==metabFlag)
                    } else if (datFlag=="raw") {
                        #i=which(ann$type==metabFlag & apply(metab,1,function(x) mean(!is.na(x))>=.5))
                        i=which(ann$type==metabFlag)
                    } else if (datFlag=="imputed") {
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
                    } else if (modelFlag%in%c("diary_tst+sex")) {
                        i=which(fitMDT$p.value[,2]<pThres2 & fitMS$p.value[,2]<pThres2)
                        expr=expr[i,]
                    } else if (modelFlag%in%c("psg_tst+sex")) {
                        i=which(fitMPT$p.value[,2]<pThres2 & fitMS$p.value[,2]<pThres2)
                        expr=expr[i,]
                    } else if (modelFlag%in%c("ln_delta_nrem+sex")) {
                        i=which(fitMDR$p.value[,2]<pThres2 & fitMS$p.value[,2]<pThres2)
                        expr=expr[i,]
                    }
                    design=model.matrix(modelThis,data=phen)
                    colnames(design)=colInfo$name2[match(colnames(design),colInfo$name1)]
                    expr=expr[,match(rownames(design),colnames(expr))]
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
                        modelThis=as.formula(modelFlag)
                        fit2=tryCatch(glm(modelThis, family="binomial",data=dat),error = function(e) e)
                        res=summary(fit2)$coef
                        rownames(res)=colInfo$name2[match(rownames(res),colInfo$name1)]
                        colnames(res)=c("coef","stdErr","z","pv")
                        tbl=data.frame(basedOn=rep(header,nrow(res)),model=rep(paste("",modelFlag,sep=""),nrow(res)),term=rownames(res),stringsAsFactors=F)
                        
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
}

## --------------------
library(coin)

if (F) {
    metab=metabT
    metabR=metabTR
    metabI=metabTI
    metabI2=metabTI2

    metab=metabO
    metabR=metabOR
    metabI=metabOI
    metabI2=metabOI2
}

metabI=datObj$metabImp

datList=c("raw_allMetab")
datList=c("orig","imputed","raw","raw_allMetab","imputed_allMetabImp")
datList=c("imputed")
datList=c("raw")

metabList=unique(ann$type); metabList=metabList[which(!metabList%in%c("internal standard"))]

for (datFlag in datList) {
    subDir=datFlag
    if (subDir!="" & !file.exists(subDir)) dir.create(file.path(subDir))
    datadir=paste(subDir,"/",sep="")
    for (metabFlag in metabList) {
        for (respFlag in c("metab")) {
            for (modelFlag in c("ptsd","sex")) {
                cat("\n\n============= ",datFlag,", ",metabFlag,", ",modelFlag,"\n",sep="")
                modelThis=as.formula(paste("~",modelFlag,sep=""))
                i=1:10
                if (datFlag==c("raw_allMetab")) {
                    i=which(ann$type==metabFlag)
                } else if (datFlag=="raw") {
                    #i=which(ann$type==metabFlag & apply(metab,1,function(x) mean(!is.na(x))>=.5))
                    i=which(ann$type==metabFlag)
                } else if (datFlag=="imputed") {
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
                } else if (modelFlag%in%c("diary_tst+sex")) {
                    i=which(fitMDT$p.value[,2]<pThres2 & fitMS$p.value[,2]<pThres2)
                    expr=expr[i,]
                } else if (modelFlag%in%c("psg_tst+sex")) {
                    i=which(fitMPT$p.value[,2]<pThres2 & fitMS$p.value[,2]<pThres2)
                    expr=expr[i,]
                } else if (modelFlag%in%c("ln_delta_nrem+sex")) {
                    i=which(fitMDR$p.value[,2]<pThres2 & fitMS$p.value[,2]<pThres2)
                    expr=expr[i,]
                }
                design=model.matrix(modelThis,data=phen)
                colnames(design)=colInfo$name2[match(colnames(design),colInfo$name1)]
                expr=expr[,match(rownames(design),colnames(expr))]
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

if (F) {
    metab=metabO
    metabR=metabOR
    metabI=metabOI
    metabI=metabOI2
}

metabI=datObj$metabImp

verbose=F

filtList=c("","_25percMostVar")
filtList=c("_25percMostVar")

testFlag="linear"
testFlag="logistic"

datFlag="raw"
datFlag="imputed"

metabList=unique(ann$type); metabList=metabList[which(!metabList%in%c("internal standard"))]

colorInfo=data.frame(grp=c("Control/Female","PTSD/Female","Control/Male","PTSD/Male"),col=c("orange","red","skyblue","blue"),stringsAsFactors=F)

colVec=rep("black",nrow(phen))
colVec[which(phen$ptsd=="Control" & phen$sex=="Female")]="orange"
colVec[which(phen$ptsd=="PTSD" & phen$sex=="Female")]="red"
colVec[which(phen$ptsd=="Control" & phen$sex=="Male")]="skyblue"
colVec[which(phen$ptsd=="PTSD" & phen$sex=="Male")]="blue"

cutoff=3
for (testFlag in c("linear","logistic")) {
    #fName2=paste("stat_prinCompMetabolite.txt",sep="")
    fName2=paste("stat_prinCompMetabolite_",testFlag,".txt",sep="")
    write.table("Metabolite meta variable from PCA\n",fName2, sep="\t", col.names=F, row.names=F, quote=F)
    if (testFlag=="linear") {
        #modelList=paste("metabolite~",c("ptsd","sex","ptsd+sex","ptsd*sex","diary_tst","psg_tst","ln_delta_nrem"),sep="")
        modelList=paste("metabolite~",c("diary_tst","psg_tst","ln_delta_nrem"),sep="")
        modelList=paste("metabolite~",c("diary_tst","psg_tst","ln_delta_nrem","diary_tst+sex","psg_tst+sex","ln_delta_nrem+sex"),sep="")
        modelList=paste("metabolite~",c("diary_tst","psg_tst","ln_delta_nrem","diary_tst+sex","psg_tst+sex","ln_delta_nrem+sex","diary_tst*sex","psg_tst*sex","ln_delta_nrem*sex"),sep="")
    } else {
        modelList=c("ptsd~metabolite","sex~metabolite","ptsd~metabolite+sex","ptsd~metabolite*sex")
        modelList=c("ptsd~metabolite","sex~metabolite","ptsd~metabolite+sex","ptsd~metabolite*sex","ptsd~metabolite*diary_tst","ptsd~metabolite*psg_tst","ptsd~metabolite*ln_delta_nrem")
    }

    statTbl=NULL
    for (filtFlag in filtList) {
        for (metabFlag in metabList) {
            if (verbose) {
                cat("\n\n===================================================\n",sep="")
                cat("===================================================\n\n",sep="")
            }
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
                if (verbose) {
                    cat("\n\n============= ",header,"\n",sep="")
                    cat("============= Model: ",modelFlag,"\n\n",sep="")
                }
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
                if (T) {
                    for (k in which(colnames(res)%in%c("coef","stdErr","t","z"))) {
                        res[,k]=round(res[,k],5)
                    }
                    for (k in c("pv")) {
                        res[,k]=signif(res[,k],2)
                    }
                }
                if (verbose) {
                    print(res)
                }
            }
        }
    }
    rownames(statTbl)=NULL
    
    pThres=.05
    tbl=statTbl
    if (T) {
        for (k in which(colnames(tbl)%in%c("coef","stdErr","t","z"))) {
            tbl[,k]=round(tbl[,k],3)
        }
        for (k in c("pv")) {
            tbl[,k]=signif(tbl[,k],3)
        }
    }
    cat('tbl[which(statTbl$pv<pThres),c("basedOn","model","term","coef","pv")]\n')
    print(tbl[which(statTbl$pv<pThres),c("basedOn","model","term","coef","pv")])
    
    tbl2=tbl[grep("25% most variable metabolites",tbl$basedOn),]
    #tbl2=statTbl[grep("25% most variable metabolites",tbl$basedOn),]
    grp=paste(tbl2$basedOn,tbl2$model)
    grpUniq=unique(grp)
    tbl1=NULL
    for (gId in 1:length(grpUniq)) {
        k=which(grp==grpUniq[gId])
        if (length(k)>1) {
            tbl2[k[2:length(k)],"basedOn"]=""
            tbl2[k[2:length(k)],"model"]=""
        }
        tbl1=tbl2[k,]
        write.table(tbl1,fName2, sep="\t", col.names=T, row.names=F, quote=F,append=T)
        write.table("",fName2, sep="\t", col.names=F, row.names=F, quote=F,append=T)
    }
}

##############################################
## Association among clinical variables

library(coin)
phen=datObj$phen

varList=c("ptsd","sex")
varList=c("diary_tst","psg_tst","ln_delta_nrem")
varList=c("ptsd","sex","diary_tst","psg_tst","ln_delta_nrem")
n=length(varList)*(length(varList)-1)/2
tmp=rep(NA,n); tmpC=rep("",n)
out=data.frame(var1=tmpC,var2=tmpC,testType=tmpC,pv=tmp,corrPears=tmp,corrSpear=tmp,stringsAsFactors=F)
k=1
for (vId1 in 1:(length(varList)-1)) {
    x1=sum(!duplicated(phen[!is.na(phen[,varList[vId1]]),varList[vId1]]))
    for (vId2 in (vId1+1):length(varList)) {
        x2=sum(!duplicated(phen[!is.na(phen[,varList[vId2]]),varList[vId2]]))
        out$var1[k]=varList[vId1]
        out$var2[k]=varList[vId2]
        if (x1<5 & x2<5) {
            out$testType[k]="Fisher's exact"
            out$pv[k]=fisher.test(phen[,varList[vId1]],phen[,varList[vId2]])$p.value
        } else if (x1==2) {
            out$testType[k]="Wilcoxon's rank sum"
            out$pv[k]=pvalue(wilcox_test(phen[,varList[vId2]]~as.factor(phen[,varList[vId1]]),distribution="exact"))
        } else if (x2==2) {
            out$testType[k]="Wilcoxon's rank sum"
            out$pv[k]=pvalue(wilcox_test(phen[,varList[vId1]]~as.factor(phen[,varList[vId2]]),distribution="exact"))
        } else if (x1<5) {
            out$testType[k]="Kruskal-Wallis"
            out$pv[k]=pvalue(kruskal_test(phen[,varList[vId2]]~as.factor(phen[,varList[vId1]]),distribution="exact"))
        } else if (x2<5) {
            out$testType[k]="Kruskal-Wallis"
            out$pv[k]=pvalue(kruskal_test(phen[,varList[vId1]]~as.factor(phen[,varList[vId2]]),distribution="exact"))
        } else {
            out$corrPears[k]=cor(phen[,varList[vId1]],phen[,varList[vId2]],use="complete.obs",method="pearson")
            out$corrSpear[k]=cor(phen[,varList[vId1]],phen[,varList[vId2]],use="complete.obs",method="spearman")
        }
        k=k+1
    }
}
tbl=out
for (k in c("pv")) {
    tbl[,k]=signif(tbl[,k],2)
}
for (k in c("corrPears","corrSpear")) {
    tbl[,k]=round(tbl[,k],2)
}
names(tbl)=c("variable1","variable2","testType","pValue","corrPearson","corrSpearman")
fName="stat_variables.txt"
write.table(tbl,fName, sep="\t", col.names=T, row.names=F, quote=F)

##############################################
