#R script implmenting GSEA (Gene Set Enrichment Analysis)
#Author: Minghui Wang
#GSEAfunc 
#Input
#(1) x, a data.frame with at least two columns:
#	strength: strength of correlation (eg, log[fold change] of differential expression) between gene and trait
#	tag: a vector of TRUE/FALSE indicating whether or not the genes are present in a gene set of interest
#(2) alpha, power to scale the weights: 0 (unweighted = Kolmogorov-Smirnov), 1 (weighted), and 2 or larger (over-weighted)
#(3) isOrdered, whether the input has been ordered (only set its value when you known what you are doing)
#
#Output
#a list with three elements: ES (enrichment score), indicator (a vector of integers indicating the location of the genes of interest), and RES (a vector of running enrichment scores)
GSEAfunc=function(x,alpha=1,isOrdered=FALSE){
	if(!isOrdered) x=x[order(- x$strength),]
	Nh=sum(x$tag)
	n=nrow(x)
	xid=which(x$tag)
	if(length(xid)<2) return(list(ES=0,indicator=xid,RES=rep(0,n)))
	f1=iy=rep(0,n)
	iy[xid]=1
	f1[xid]=abs(x$strength[xid])^alpha
	f1=cumsum(f1)
	f1=f1/f1[n]
	f2=cumsum(1-iy)/(n-Nh)
	RES=f1-f2
	a=max(RES)
	b=min(RES)
	ES=ifelse(a > -b, a, b)
	return(list(ES=ES,indicator=xid,RES=RES))
}
GSEAfunc2=function(x,alpha=1,isOrdered=FALSE){
#compute ES as GSEAfunc but much faster (7x)
	if(!isOrdered) x=x[order(- x$strength),]
	Nh=sum(x$tag)
	n=nrow(x)
	Nm=n-Nh
	xid=which(x$tag)
	if(length(xid)<2) return(0)
	xNum=xid-c(0,xid[-Nh])-1
	if(alpha == 0){
		ys=rep(1,Nh)
	}else if (alpha == 1){
		ys=abs(x$strength[xid])
	}else if (alpha == 2){
		ys=abs(x$strength[xid])
		ys=ys*ys
	}else{
		ys=abs(x$strength[xid])^alpha
	}
	ys=ys/sum(ys)
	f1=cumsum(ys)
	f2=cumsum(xNum)/Nm
	f3=f1-f2
	a=max(f3)
	b=min(f3-ys)
	ifelse(a>-b,a,b)
}
GSEA=function(x, go, alpha=1, permutations=1000,ncores=1,iseed=12345){
#x, a data.frame with first column gene ids and second column phenotype association measures
#go, a vector of gene ids
	x1=data.frame(strength=x[,2],tag=x[,1] %in% go,stringsAsFactors=FALSE)
	ntag=sum(x1$tag)
	if(ntag<2){
		Res=list(ES=0,indicator=which(x1$tag),RES=rep(0,nrow(x1)),NES=0,Pvalue=1)
		class(Res)='gsea'
		return(Res)
	}
	x1=x1[order(- x1$strength),] #we will keep the order from now on
	isOrdered=TRUE
	x1$strength=abs(x1$strength) ^ alpha
	Res=GSEAfunc(x1, alpha=alpha, isOrdered=isOrdered)
	x1$tag=FALSE #reset for repeated sampling in permutations
	func1=function(i,x1,ntag,alpha,isOrdered=FALSE){
		x1[sample.int(nrow(x1),ntag,replace=FALSE),2]=TRUE #we will keep the order of column "strength" such that there is no need to sort it again
		GSEAfunc2(x1,alpha=alpha,isOrdered=isOrdered)
	}
	if(ncores>1){
		cl=makeCluster(ncores)
		cat('Recruited',ncores,'cores\n')
		clusterSetRNGStream(cl, iseed)
		clusterExport(cl,varlist=c('GSEAfunc2'),envir=environment())
		ESnull=parSapply(cl,1:permutations,func1,x1=x1,ntag=ntag,alpha=alpha,isOrdered=isOrdered)
		stopCluster(cl)
	}else{
		set.seed(iseed)
		ESnull=sapply(1:permutations,func1,x1=x1,ntag=ntag,alpha=alpha,isOrdered=isOrdered)
	}
	if(Res$ES>=0){
		ESnull=ESnull[ESnull>=0]
		Res$NES=Res$ES/mean(ESnull)
		Res$Pvalue=sum(ESnull >= Res$ES)/length(ESnull)
	}else{
		ESnull=ESnull[ESnull < 0]
		Res$NES=Res$ES/abs(mean(ESnull))
		Res$Pvalue=sum(ESnull <= Res$ES)/length(ESnull)
	}
	class(Res)='gsea'
	Res
}
#perform GSEA analysis for a set of GO categoris in phenotype association signatures
callGSEAfunc=function(Set1, Set2, alpha=1, permutations=1000,ncores=1,iseed = 12345){
#Set1, usually Gene Ontology
#Set2, data.frame of phenotype association signatures
	stopifnot(is.data.frame(Set2))
	colnames(Set2)=c('tag','phenotype','strength')
	Res=NULL
	if(ncores>1){
		cl=makeCluster(ncores)
		cat('Recruit',ncores,'cores\n')
		clusterSetRNGStream(cl, iseed)
	}else{
		set.seed(iseed)
	}
	fun1=function(k,dat1,gos,alpha,isOrdered=FALSE){
		dat1$tag=FALSE #method 2
		sapply(gos,function(iGO,dat1,alpha,isOrdered=FALSE){
			dat1$tag[sample(nrow(dat1),length(iGO))]=TRUE #method 2
			GSEAfunc2(dat1,alpha=alpha,isOrdered=isOrdered)
		},dat1=dat1,alpha=alpha,isOrdered=isOrdered)
	}
	for(dat in split(Set2,Set2$phenotype)){ #for each phenotype
		Set11=Set1[Set1[,1] %in% dat$tag,]
		if(nrow(Set11)==0) stop('No common elements between sets\n')
		Category.Size=table(Set1[,2])
		Category.Size=setNames(as.vector(Category.Size),names(Category.Size))
		sNoOverlap=unique(Set1[,2])
		sNoOverlap=sNoOverlap[!sNoOverlap %in% unique(Set11[,2])]
		ResNoOverlap=NULL
		phenotype=dat$phenotype[1]
		if(length(sNoOverlap)>0){
			ResNoOverlap=data.frame(Term=sNoOverlap,Input=phenotype,Overlap.Size=0,Input.Size=nrow(dat),Category.Size=as.vector(Category.Size[sNoOverlap]),ES=0,NES=0,Pvalue=1,P.adj=1,Elements='',stringsAsFactors=FALSE)
		}
		gos=split(Set11[,1],Set11[,2])
		dat=dat[order(- dat$strength),colnames(dat)!='phenotype'] #we will sort the data by the strength and keep the order from now on
		isOrdered=TRUE
		tab2=lapply(gos,function(iGO,dat,alpha,isOrdered=FALSE){
			Overlap=dat$tag
			dat$tag = dat$tag %in% iGO
			ESs=0
			if(any(dat$tag)) ESs=GSEAfunc2(dat,alpha=alpha,isOrdered=isOrdered)
			Overlap=sort(Overlap[dat$tag])
			data.frame(Overlap.Size=length(Overlap),Input.Size=nrow(dat),Category.Size=length(iGO),ES=ESs,NES=NA,Pvalue=NA,P.adj=NA,Elements=paste0(Overlap,collapse=';'),stringsAsFactors=FALSE)
		},dat=dat,alpha=alpha,isOrdered=isOrdered)
		tab2=do.call(rbind,tab2)
		#Compute significance and control for multiple tests by permutations
		if(ncores>1){
			ESnull=parLapply(cl,1:permutations,fun1,dat1=dat,gos=gos,alpha=alpha,isOrdered=isOrdered)
		}else{
			ESnull=lapply(1:permutations,fun1,dat1=dat,gos=gos,alpha=alpha,isOrdered=isOrdered)
		}
		ESnull=do.call(cbind,ESnull) #number of GO terms x number of permutations
		NESnull=vector('list',nrow(tab2))
		for(i in 1:nrow(tab2)){
			if(tab2$ES[i]==0){
				NESnull[[i]]=0
				tab2$NES[i]=0
				#tab2$Pvalue[i]=1
				next
			}else if(tab2$ES[i]<0){
				s1=ESnull[i,ESnull[i,]<0]
				#tab2$Pvalue[i]=(sum(s1 <= tab2$ES[i])+1)/ (length(s1)+1)
			}else{
				s1=ESnull[i,ESnull[i,]>0]
				#tab2$Pvalue[i]=(sum(s1 >= tab2$ES[i])+1) / (length(s1)+1)
			}
			if(length(s1)==0) next #this happens when there is insufficient number of permutations
			m1=abs(mean(s1))
			NESnull[[i]]=s1/m1
			tab2$NES[i]=tab2$ES[i]/m1
		}
		NESnull=unlist(NESnull)
		#Method 1: separate positive and negative scores (no longer used because it leads to different significance P values for the same absolute NES with opposite signs)
		#NESnull=sort(NESnull)
		#wpid = tab2$NES >= 0
		#if(nrow(tab2)>1){
			#if(any(wpid)) tab2$P.adj[wpid]=getFDR.0(tab2$NES[wpid],NESnull[NESnull >= 0],truncated.size=0)
			#if(any(! wpid)) tab2$P.adj[ ! wpid ]=getFDR.0(abs(tab2$NES[ ! wpid]),abs(NESnull[NESnull < 0]),truncated.size=0)
		#}
		#Method 2: combine positive and negative
		NESnull=abs(NESnull)
		pq=getFDR.0(abs(no.na(tab2$NES)),NESnull,truncated.size=0,return.p.value=TRUE)
		tab2$Pvalue[!is.na(tab2$NES)]=pq$P.value
		tab2$P.adj[!is.na(tab2$NES)]=pq$P.adj
		Res=rbind(Res,ResNoOverlap,data.frame(Term=names(gos),Input=phenotype,tab2,stringsAsFactors=FALSE))
	}
	if(ncores>1) stopCluster(cl)
	colnames(Res)=c('Category','Input','Overlap.Size','Input.Size','Category.Size','ES','NES','Pvalue','P.adj','Elements')
	rownames(Res)=NULL
	Res
}
