#' \code{tree.gen} simulate trees with representative sampling
#' 
#' @param pars a vector consists of speciation rate of state 0, speciation rate of state 1, extinction rate of state 0, extinction rate of state 1, transition rate of state 0 to state 1, transition rate of state 1 to state 0
#' @param max.taxa the number of total extant taxa of the tree,
#' @param x0 the root state of the tree,
#' @param sampling.c percentage of groups being sampled.
#' @param sampling.f the maximum number of representatives sampled for each group.
#' @param define.clade how a group is identified in the tree
#'        "state": a clade is a group if all its extant members are in the same state
#'        "time": a clade is a group if its most recent common ancestor existed after half of the evolutionary time of the tree
#' @return rawtree a fully sampled tree
#' @return bonetree the tree with one representative per group
#' @return ftree the tree that is collpased to apply unresolved clade correction
#' @return xtree the tree with only the sampled tips
#' @return unresolved an input matrix to apply unresolved clade correction
#' @return sfraction an input matrix to apply our correction
#' @return ssize an input list for the tlistcal function
#' 
tree.gen <- function(pars,max.taxa,x0,sampling.c,sampling.f,define.clade=c("state","time")) {
#generating tree
rawtree <- tree.bisse(pars=pars,max.taxa=max.taxa,x0=x0)
while(sum(rawtree$tip.state)>length(rawtree$tip.state) || sum(rawtree$tip.state)==0) {
	rawtree <- tree.bisse(pars=pars,max.taxa=max.taxa,x0=x0)
}
while(is.null(rawtree)) {rawtree <- tree.bisse(pars=pars,max.taxa=max.taxa,x0=x0)}
subtree <- subtrees(rawtree,wait=F)
if (define.clade=="state") {
a <- sapply(subtree,function (i) ifelse(sum(rawtree$tip.state[i$tip.label])==0,1,0))
subtree0 <- subtree[a==1]
a <- sapply(subtree,function (i) ifelse(sum(rawtree$tip.state[i$tip.label])==length(i$tip.label),1,0))
subtree1 <- subtree[a==1]
b <- sapply(1:length(subtree0),function (i) sapply(c(1:length(subtree0))[-i],function (j) ifelse(length(intersect(subtree0[[i]]$tip.label,subtree0[[j]]$tip.label))==length(subtree0[[i]]$tip.label),j,0)))
subtree0 <- subtree0[which(colSums(b)==0)]
b <- sapply(1:length(subtree1),function (i) sapply(c(1:length(subtree1))[-i],function (j) ifelse(length(intersect(subtree1[[i]]$tip.label,subtree1[[j]]$tip.label))==length(subtree1[[i]]$tip.label),j,0)))
subtree1 <- subtree1[which(colSums(b)==0)]
}
if (define.clade=="time") {
  time.limit <- max(branching.times(rawtree))/2
  a <- sapply(subtree,function (i) ifelse(max(branching.times(i))<=time.limit,1,0))
  subtree0 <- subtree[a==1]
  b <- sapply(1:length(subtree0),function (i) sapply(c(1:length(subtree0))[-i],function (j) ifelse(length(intersect(subtree0[[i]]$tip.label,subtree0[[j]]$tip.label))==length(subtree0[[i]]$tip.label),j,0)))
  subtree0 <- subtree0[which(colSums(b)==0)]
  subtree1 <- NULL
}
singletip <- rawtree$tip.label
for (i in c(subtree0,subtree1)) {
  singletip <- setdiff(singletip,i$tip.label)
}
#generating bonetree
bonetree <- rawtree
tip.label <- numeric()
Nc <- numeric()
n0 <- numeric()
n1 <- numeric()
for (i in subtree0) {
    bonetree <- drop.tip.fixed(bonetree,i$tip.label[-1])
    tip.label <- c(tip.label,i$tip.label[1])
    Nc <- c(Nc,length(i$tip.label))
    n0 <- c(n0,sum(rawtree$tip.state[i$tip.label]==0))
    n1 <- c(n1,sum(rawtree$tip.state[i$tip.label]==1))
}
if (!is.null(subtree1)) {
for (i in subtree1) {
  bonetree <- drop.tip.fixed(bonetree,i$tip.label[-1])
  tip.label <- c(tip.label,i$tip.label[1])
  Nc <- c(Nc,length(i$tip.label))
  n0 <- c(n0,sum(rawtree$tip.state[i$tip.label]==0))
  n1 <- c(n1,sum(rawtree$tip.state[i$tip.label]==1))
}
}
bonetree$tip.state <- rawtree$tip.state[bonetree$tip.label]
tip.label <- c(tip.label,singletip)
Nc <- c(Nc,rep(1,length(singletip)))
n0 <- c(n0,as.numeric(rawtree$tip.state[singletip]==0))
n1 <- c(n1,as.numeric(rawtree$tip.state[singletip]==1))
names(Nc) <- tip.label
unresolved <- data.frame(Nc,n0,n1)
#generating xtree to run with our correction; ftree to run with the unresolved clade correction
bNtips <- length(bonetree$tip.label)
bname <- sapply(bonetree$edge[,2],function (j) ifelse(j>bNtips,bonetree$node.label[j-bNtips],bonetree$tip.label[j]))
rownames(bonetree$edge) <- bname
xtree <- rawtree
ftree <- bonetree
unsampled.tip <- numeric()
included.tip0 <- vector("list",length(subtree0))
included.tip1 <- vector("list",length(subtree1))
for (i in 1:length(subtree0)) {
    n <- min(length(subtree0[[i]]$tip.label),sample.int(sampling.f,1))
    if (n==1) {
      tmp <- 1
    } else {
      tmp <- c(1,ifelse(length(subtree0[[i]]$tip.label)>2,sample(2:length(subtree0[[i]]$tip.label),size=n-1),2))
    }
    if (runif(1)<=sampling.c) {
      if (length(tmp)<length(subtree0[[i]]$tip.label)) {xtree <- drop.tip.fixed(xtree,subtree0[[i]]$tip.label[-tmp])}
      included.tip0[[i]] <- subtree0[[i]]$tip.label[tmp]
    } else {
      xtree <- drop.tip.fixed(xtree,subtree0[[i]]$tip.label)
      ftree <- drop.tip.fixed(ftree,subtree0[[i]]$tip.label[1])
      unsampled.tip <- c(unsampled.tip,subtree0[[i]]$tip.label[1])
    }
}
if (!is.null(subtree1)) {
for (i in 1:length(subtree1)) {
  n <- min(length(subtree1[[i]]$tip.label),sample.int(sampling.f,1))
  if (n==1) {
    tmp <- 1
  } else {
    tmp <- c(1,ifelse(length(subtree1[[i]]$tip.label)>2,sample(2:length(subtree1[[i]]$tip.label),size=n-1),2))
  }
  if (runif(1)<=sampling.c) {
  	 if (length(tmp)<length(subtree1[[i]]$tip.label)) {xtree <- drop.tip.fixed(xtree,subtree1[[i]]$tip.label[-tmp])}
 	 included.tip1[[i]] <- subtree1[[i]]$tip.label[tmp]
  } else {
      xtree <- drop.tip.fixed(xtree,subtree1[[i]]$tip.label)
      ftree <- drop.tip.fixed(ftree,subtree1[[i]]$tip.label[1])
      unsampled.tip <- c(unsampled.tip,subtree1[[i]]$tip.label[1])
    }
}
}
if (length(singletip)>0) {
for (i in singletip) {
  if (bonetree$tip.state[i]==0) {
  	if (runif(1)>sampling.c) {
    		xtree <- drop.tip.fixed(xtree,i)
   		  ftree <- drop.tip.fixed(ftree,i)
    		unsampled.tip <- c(unsampled.tip,i)
    	} else {
    		included.tip0 <- c(included.tip0,i)
    	}
  } else {
  	if (runif(1)>sampling.c) {
    	xtree <- drop.tip.fixed(xtree,i)
   		ftree <- drop.tip.fixed(ftree,i)
    	unsampled.tip <- c(unsampled.tip,i)
    } else {
		included.tip1 <- c(included.tip1,i)
	}
  }
}
}
xtree$tip.state <- rawtree$tip.state[xtree$tip.label]
#calculate unresolved table for ftree
fNtips <- length(ftree$tip.label)
anclist <- numeric(length(unsampled.tip))
names(anclist) <- unsampled.tip
unresolved2 <- matrix(NA,fNtips,3)
colnames(unresolved2) <- c("Nc","n0","n1")
rownames(unresolved2) <- ftree$tip.label
for (i in unsampled.tip) {
	anc <- bonetree$node.label[bonetree$edge[i,1]-bNtips]
	if (as.numeric(substring(anc,3))>=as.numeric(substring(ftree$node.label[1],3))) {
		while(!is.element(anc,ftree$node.label)) {
			anc <- bonetree$node.label[bonetree$edge[anc,1]-bNtips]
			if (as.numeric(substring(anc,3))<as.numeric(substring(ftree$node.label[1],3))) {
				anc <- 0
				break
			}
		}
	} else {
		anc <- 0
	}
	anclist[i] <- anc
}
anclist <- anclist[anclist!=0]
anclist <- anclist[order(as.numeric(sapply(anclist,substring,3)))]
if (length(anclist)>0) {
for (j in 1:length(anclist)) {
	i <- anclist[j]
	if (is.element(i,ftree$node.label)) {
	tmp <- which(bonetree$node.label==i)+bNtips
	tmp <- bonetree$edge[which(bonetree$edge[,1]==tmp),2]
	des <- bonetree$node.label[tmp[tmp>bNtips]-bNtips]
	tips1 <- extract.clade(bonetree,des[1])$tip.label
	if (is.element(names(i),tips1)) {
		tips <- intersect(tips1,ftree$tip.label)
		unresolved2[tips[1],] <- as.numeric(colSums(unresolved[tips1,]))
		if (length(tips)>1) {
			tmp <- sapply(tips[-1],function (z) which(ftree$tip.label==z))
			unresolved2 <- unresolved2[-tmp,]
			ftree <- drop.tip.fixed(ftree,tips[-1])
		}
	}
	if (length(des)>1) {
	tips2 <- extract.clade(bonetree,des[2])$tip.label
	if (is.element(names(i),tips2)) {
		tips <- intersect(tips2,ftree$tip.label)
		unresolved2[tips[1],] <- as.numeric(colSums(unresolved[tips2,]))
		if (length(tips)>1) {
			tmp <- sapply(tips[-1],function (z) which(ftree$tip.label==z))
			unresolved2 <- unresolved2[-tmp,]
			ftree <- drop.tip.fixed(ftree,tips[-1])
		}
	}
	}
	}
}
}
tmp <- names(which(is.na(unresolved2[,1])))
unresolved2[tmp,] <- as.matrix(unresolved[tmp,],length(tmp),3)
ftree$tip.state <- rawtree$tip.state[ftree$tip.label]
unresolved2 <- data.frame(tip.label=ftree$tip.label,unresolved2)
#calculate sfraction table and ssize list for xtree
xNtips <- length(xtree$tip.label)
sfraction <- matrix(0,dim(xtree$edge)[1],5)
xname <- sapply(xtree$edge[,2],function (j) ifelse(j>xNtips,xtree$node.label[j-xNtips],xtree$tip.label[j]))
rownames(sfraction) <- xname
ssize.new <- vector("list",1)
if (length(included.tip0)>0) {
for (i in included.tip0) {
  if (!is.null(i)) {
  if (define.clade=="state") {
    sfraction[i,1:3] <- c(rep(0,length(i)),rep(length(i),length(i)),rep(unresolved[i[[1]],"Nc"],length(i)))
  }
  if (define.clade=="time") {
    tmp <- rawtree$tip.state[i]
    sfraction[i,] <- c(tmp,rep(sum(tmp==0),length(i)),rep(unresolved[i[[1]],"n0"],length(i)),rep(sum(tmp==1),length(i)),rep(unresolved[i[[1]],"n1"],length(i)))
    if (sfraction[i[1],2]==0 && sfraction[i[1],3]>0) {
    	ssize.tmp <- vector("list",1)
    	if (length(i)>1) {
    		idx <- sapply(i,function (x) which(xname==x))
    		ssize.tmp[[1]]$edge <- xname[(min(idx)-1):max(idx)]
    	} else {
    		ssize.tmp[[1]]$edge <- i
    	}
    	ssize.tmp[[1]]$state <- c(n01=1/as.numeric(sfraction[i[1],3]),n11=0,e01=1-1/as.numeric(sfraction[i[1],3]),e11=1-as.numeric(sfraction[i[1],4]/sfraction[i[1],5]))
    	ssize.tmp[[1]]$ordered <- FALSE
    	sfraction[i,2] <- rep(1,length(i))
    	ssize.new <- c(ssize.new,ssize.tmp)
    }
    if (sfraction[i[1],4]==0 && sfraction[i[1],5]>0) {
    	ssize.tmp <- vector("list",1)
    	if (length(i)>1) {
    		idx <- sapply(i,function (x) which(xname==x))
    		ssize.tmp[[1]]$edge <- xname[(min(idx)-1):max(idx)]
    	} else {
    		ssize.tmp[[1]]$edge <- i
    	}
    	ssize.tmp[[1]]$state <- c(n01=0,n11=1/as.numeric(sfraction[i[1],5]),e01=1-as.numeric(sfraction[i[1],2])/as.numeric(sfraction[i[1],3]),e11=1-1/as.numeric(sfraction[i[1],5]))
    	ssize.tmp[[1]]$ordered <- FALSE
    	sfraction[i,4] <- rep(1,length(i))
    	ssize.new <- c(ssize.new,ssize.tmp)
    }
  }
  }
}
}
if (length(included.tip1)>0) {
for (i in included.tip1) {
  if (!is.null(i)) {
    sfraction[i,c(1,4,5)] <- c(rep(1,length(i)),rep(length(i),length(i)),rep(unresolved[i[[1]],"Nc"],length(i)))
  }
}
}
anclist <- numeric(length(unsampled.tip))
for (i in unsampled.tip) {
	anc <- bonetree$node.label[bonetree$edge[i,1]-bNtips]
	if (as.numeric(substring(anc,3))>=as.numeric(substring(xtree$node.label[1],3))) {
		while(!is.element(anc,xtree$node.label)) {
			anc <- bonetree$node.label[bonetree$edge[anc,1]-bNtips]
			if (as.numeric(substring(anc,3))<as.numeric(substring(xtree$node.label[1],3))) {
				anc <- 0
				break
			}
		}
	} else {
		anc <- 0
	}
	anclist[i] <- anc
}
ssize <- numeric()
anclist <- anclist[anclist!=0]
anclist <- anclist[order(as.numeric(sapply(anclist,substring,3)))]
for (i in unique(anclist)) {
	tmp <- names(anclist)[which(anclist==i)]
	des <- bname[bonetree$edge[,1]==(which(bonetree$node.label==i)+bNtips)]
	xdes <- xname[xtree$edge[,1]==(which(xtree$node.label==i)+xNtips)]
	if (substring(des[1],1,2)=="nd") {
		clade1 <- extract.clade(bonetree,des[1])$tip.label
	} else {
		clade1 <- des[1]
	}
	if (substring(des[2],1,2)=="nd") {
		clade2 <- extract.clade(bonetree,des[2])$tip.label
	} else {
		clade2 <- des[2]
	}
	if (substring(xdes[1],1,2)=="nd") {
		tmp3 <- extract.clade(xtree,xdes[1])
		ll <- tmp3$tip.label
	} else {
		tmp3 <- xdes[1]
		ll <- tmp3
	}
	if (length(intersect(clade1,ll))>0) {
		xclade1 <- tmp3
		if (substring(xdes[2],1,2)=="nd") {
			xclade2 <- extract.clade(xtree,xdes[2])
		} else {
			xclade2 <- xdes[2]
		}
	} else {
		if (substring(xdes[2],1,2)=="nd") {
			xclade1 <- extract.clade(xtree,xdes[2])
		} else {
			xclade1 <- xdes[2]
		}
		xclade2 <- tmp3
	}
	n <- length(tmp)	
	res <- numeric()
		j=n
		while(j>1) {
			com <- combn(tmp,j)
			tmp2 <- com[,sapply(1:dim(com)[2], function (z) is.monophyletic(bonetree,com[,z]))]
			if (length(tmp2)>0) {
				if (!is.matrix(tmp2)) {
					j <- n-j+1
					tmp <- tmp[!tmp%in%tmp2]
					n <- length(tmp)
					if (length(intersect(clade1,tmp2))>0) {
						res <- rbind(res,c(tmp2[1],ifelse(substring(xdes[1],1,2)=="nd",xclade1$node.label[1],xclade1),colSums(unresolved[tmp2,])))
					} else {
						res <- rbind(res,c(tmp2[1],ifelse(substring(xdes[2],1,2)=="nd",xclade2$node.label[1],xclade2),colSums(unresolved[tmp2,])))
					}
				} else {
					m <- dim(tmp2)[2]
					j <- n-j*m+1
					tmp <- tmp[!tmp%in%tmp2]
					n <- length(tmp)
					for (z in 1:m) {
						if (length(intersect(clade1,tmp2[,z]))>0) {
							res <- rbind(res,c(tmp2[1,z],ifelse(substring(xdes[1],1,2)=="nd",xclade1$node.label[1],xclade1),colSums(unresolved[tmp2[,z],])))
						} else {
							res <- rbind(res,c(tmp2[1,z],ifelse(substring(xdes[2],1,2)=="nd",xclade2$node.label[1],xclade2),colSums(unresolved[tmp2[,z],])))
						}
					}
				}
			}
			j <- j-1
		}
		if (length(tmp)>0) {
			for (j in 1:length(tmp)) {
				if (length(intersect(clade1,tmp[j]))>0) {
					res <- rbind(res,c(tmp[j],ifelse(substring(xdes[1],1,2)=="nd",xclade1$node.label[1],xclade1),unresolved[tmp[j],]))
				} else {
					res <- rbind(res,c(tmp[j],ifelse(substring(xdes[2],1,2)=="nd",xclade2$node.label[1],xclade2),unresolved[tmp[j],]))
				}
			}
		}
		if (length(res)>5) {
		m <- dim(res)[1]
		pos <- numeric()
		for (z in 1:m) {
			if (substring(res[z,2],1,2)=="nd") {
				tips <- extract.clade(rawtree,res[z,2][[1]])$tip.label
			} else {
				tips <- res[z,2][[1]]
			}
			pos <- c(pos,getMRCA(phy=rawtree,tip=c(res[z,1][[1]],tips)))
		}
		pos <- order(pos,decreasing=T)
		res <- res[pos,]
		ssize <- rbind(ssize,res[,-1])
		} else {
		ssize <- rbind(ssize,res[-1])
		}
}
if (length(ssize)>0) {
tmp <- unlist(unique(ssize[,1]))
for (i in 1:length(tmp)) {
	idx <- which(ssize[,1]==tmp[i])
	ssize.tmp <- vector("list",1)
	if (length(idx)==1) {
		ssize.tmp[[1]] <- list(edge=tmp[i],state=c(n01=as.numeric(ssize[idx,3])/as.numeric(ssize[idx,2])^2,n11=as.numeric(ssize[idx,4])/as.numeric(ssize[idx,2])^2,e01=ifelse(as.numeric(ssize[idx,3])>0,1-as.numeric(ssize[idx,3])/as.numeric(ssize[idx,2])^2,0),e11=ifelse(as.numeric(ssize[idx,4])>0,1-as.numeric(ssize[idx,4])/as.numeric(ssize[idx,2])^2,0)),ordered=TRUE)
	} else {
		ssize.tmp[[1]] <- list(edge=tmp[i],state=cbind(n01=as.numeric(ssize[idx,3])/as.numeric(ssize[idx,2])^2,n11=as.numeric(ssize[idx,4])/as.numeric(ssize[idx,2])^2,e01=ifelse(as.numeric(ssize[idx,3])>0,as.numeric(ssize[idx,3])/as.numeric(ssize[idx,2])^2,0),e11=ifelse(as.numeric(ssize[idx,4])>0,as.numeric(ssize[idx,4])/as.numeric(ssize[idx,2])^2,0)),ordered=TRUE)
	}
	ssize.new <- c(ssize.new,ssize.tmp)
}
}
ssize.new <- ssize.new[-1]

list(rawtree=rawtree,bonetree=bonetree,ftree=ftree,xtree=xtree,unresolved=unresolved2,sfraction=sfraction,ssize=ssize.new)
}