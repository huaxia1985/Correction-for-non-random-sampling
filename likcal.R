#' \code{likcal} calculate the likelihood under the BiSSE framework when sampling is not random
#' 
#' @param p the parameter values of the BiSSE model, a vector consists of diversification rate of state 0, diversification rate of state 1, extinction rate of state 0, extinction rate of state 1, transition rate from state 0 to state 1, transition rate from state 1 to state 0
#' @param tree a phylogeny (of class phylo)
#' @param sfraction a matrix consists of the sampling fraction for each tip of the tree. Rownames are the tipnames.
#'        The first column is the state of the tip, 0 for state 0 and 1 for state 1.
#'        The second column is the number of sampled tips with state 0 in the group that the tip is sampled as a representative.
#'        The third column is the total number of extant members with state 0 in the group.
#'        The fourth column is the number of sampled tips with state 1 in the group that the top is sampled as a representative.
#'        The fifth column is the total number of extant members with state 1 in the group.
#' @param root.state the root.state of the root. If root.state is not given, than the solution of in the appendix of FitzJohn et al. 2009 is used to calculate the overall likelihood at the root
#' @param tlist the output from the tlistcal function, which list all possible locations where an un-sampled group or state can be attached to the tree
#' @param constrained whether the BiSSE model is constrained to have equal diversification and extinction rates for different states.
#'        The default is FALSE.
#'        If constrained = TRUE, then input p has four values diversification rate, extinction rate, transition rate from state 0 to state 1, transition rate from state 1 to state 0
#' @return the overall likelihood function of a BiSSE model for the tree
#' @examples \dontrun{
#' library(deSolve)
#' p <- c(0.1,0.1,0.03,0.03,0.01,0.001)
#' px <- c(p[1]-p[3],p[2]-p[4],p[3],p[4],p[5],p[6])
#' treelist <- tree.gen(pars=p,max.taxa=500,x0=0,sampling.c=1,sampling.f=5,define.clade="state") #simulate a tree with representative sampling
#' tlist <- tlistcal(tree=treelist$xtree,ssize=treelist$ssize,dt=1) #calculate tlist
#' xtrue <- likcal(p=px,tree=treelist$xtree,sfraction=treelist$sfraction,ssize=treelist$ssize,root.state=NULL,tlist=tlist) #calculate the likelihood of the model that simulates the tree

likcal <- function (p,tree,sfraction,root.state,tlist,constrained=FALSE) {
  if (constrained==TRUE) {
    pars <- as.numeric(c(p[1]+p[2],p[1]+p[2],p[2],p[2],p[3],p[4]))
  } else {
    pars <- as.numeric(c(p[1]+p[3],p[2]+p[4],p[3],p[4],p[5],p[6]))
  }
  calD <- function (time,state,par) {
    with(as.list(c(state,par)), {
      s0 <- par[1]
      s1 <- par[2]
      x0 <- par[3]
      x1 <- par[4]
      d01 <- par[5]
      d10 <- par[6]
      de0 <- -(s0+d01+x0)*e0+x0+d01*e1+s0*e0*e0
      de1 <- -(s1+d10+x1)*e1+x1+d10*e0+s1*e1*e1
      dn0 <- -(s0+d01+x0)*n0+d01*n1+2*s0*n0*e0
      dn1 <- -(s1+d10+x1)*n1+d10*n0+2*s1*n1*e1
      return(list(c(dn0,dn1,de0,de1)))
    })
  } 
  calE <- function (time,state,par) {
    with(as.list(c(state,par)), {
      s0 <- par[1]
      s1 <- par[2]
      x0 <- par[3]
      x1 <- par[4]
      d01 <- par[5]
      d10 <- par[6]
      de0 <- -(s0+d01+x0)*e0+x0+d01*e1+s0*e0*e0
      de1 <- -(s1+d10+x1)*e1+x1+d10*e0+s1*e1*e1
      return(list(c(de0,de1)))
    })
  } 
  intD <- function (pars,tlist,summ,edge,initial=NULL) {
  	integrand <- function (t,initial,state,Tt2,Tt1,pars) {
  		t <- cumsum(t)+Tt1
  		t <- as.numeric(c(t,Tt2))
  		D <- as.numeric(ode(c(n0=as.numeric(initial[1]),n1=as.numeric(initial[2]),e0=as.numeric(initial[3]),e1=as.numeric(initial[4])),c(Tt1,t[1]),calD,pars)[2,2:3])
  		D[D<0] <- 10^-6
  		if (!is.matrix(state)) {
  			D1 <- as.numeric(ode(c(n0=as.numeric(state[1]),n1=as.numeric(state[2]),e0=as.numeric(state[3]),e1=as.numeric(state[4])),c(0,t[1]),calD,pars)[2,2:3])
  			D1[D1<0] <- 10^-6
  			E <- as.numeric(ode(c(e0=0,e1=0),c(0,t[1]),calE,pars)[2,2:3])
  			E[E<0] <- 0
  			D <- as.numeric(ode(c(n0=D1[1]*D[1]*pars[1],n1=D1[2]*D[2]*pars[2],e0=E[1],e1=E[2]),c(t[1],t[1+1]),calD,pars)[2,2:3])
  			D[D<0] <- 10^-6
  		} else {
		  for (i in 1:dim(state)[1]) {
			  D1 <- as.numeric(ode(c(n0=as.numeric(state[i,1]),n1=as.numeric(state[i,2]),e0=as.numeric(state[i,3]),e1=as.numeric(state[i,4])),c(0,t[i]),calD,pars)[2,2:3])
  			D1[D1<0] <- 10^-6
  			E <- as.numeric(ode(c(e0=0,e1=0),c(0,t[i]),calE,pars)[2,2:3])
  			E[E<0] <- 0
  			D <- as.numeric(ode(c(n0=D1[1]*D[1]*pars[1],n1=D1[2]*D[2]*pars[2],e0=E[1],e1=E[2]),c(t[i],t[i+1]),calD,pars)[2,2:3])
  			D[D<0] <- 10^-6
		  }
		}
		D
	  }
  	integrand2 <- function (t,initial,state,Tt2,Tt1,pars) {
  		t <- cumsum(t)+Tt1
  		t <- as.numeric(c(t,Tt2))
  		tmp <- ode(c(n0=as.numeric(initial[1]),n1=as.numeric(initial[2]),e0=as.numeric(initial[3]),e1=as.numeric(initial[4])),c(Tt1,t[1]),calD,pars)
  		D <- tmp[2,2:3]
  		D[D<0] <- 10^-6
  		E <- tmp[2,4:5]
  		E[E<0] <- 0
  		if (!is.matrix(state)) {
  			D1 <- as.numeric(ode(c(n0=as.numeric(state[1]),n1=as.numeric(state[2]),e0=as.numeric(state[3]),e1=as.numeric(state[4])),c(0,t[1]),calD,pars)[2,2:3])
  			D1[D1<0] <- 10^-6
  			D <- as.numeric(ode(c(n0=as.numeric(D1[1]*D[1]*pars[1]),n1=as.numeric(D1[2]*D[2]*pars[2]),e0=as.numeric(E[1]),e1=as.numeric(E[2])),c(t[1],t[1+1]),calD,pars)[2,2:3])
  			D[D<0] <- 10^-6
  		} else {
  			D1 <- as.numeric(ode(c(n0=as.numeric(state[1,1]),n1=as.numeric(state[1,2]),e0=as.numeric(state[1,3]),e1=as.numeric(state[1,4])),c(0,t[1]),calD,pars)[2,2:3])
  			D1[D1<0] <- 10^-6
  			D <- as.numeric(ode(c(n0=as.numeric(D1[1]*D[1]*pars[1]),n1=as.numeric(D1[2]*D[2]*pars[2]),e0=as.numeric(E[1]),e1=as.numeric(E[2])),c(t[1],t[1+1]),calD,pars)[2,2:3])
  			D[D<0] <- 10^-6
			  for (i in 2:dim(state)[1]) {
				  D1 <- as.numeric(ode(c(n0=as.numeric(state[i,1]),n1=as.numeric(state[i,2]),e0=as.numeric(state[i,3]),e1=as.numeric(state[i,4])),c(0,t[i]),calD,pars)[2,2:3])
  				D1[D1<0] <- 10^-6
  				E <- as.numeric(ode(c(e0=0,e1=0),c(0,t[i]),calE,pars)[2,2:3])
  				E[E<0] <- 0
  				D <- as.numeric(ode(c(n0=D1[1]*D[1]*pars[1],n1=D1[2]*D[2]*pars[2],e0=E[1],e1=E[2]),c(t[i],t[i+1]),calD,pars)[2,2:3])
  				D[D<0] <- 10^-6
			}
		}
		D
	  }
	if (tlist$type==1) {
		if (is.null(initial)) {
			initial <- c(n0=as.numeric(summ[tlist$edge,"d0"]),n1=as.numeric(summ[tlist$edge,"d1"]),e0=as.numeric(summ[tlist$edge,"e0"]),e1=as.numeric(summ[tlist$edge,"e1"]))
		}
		if (length(tlist$time)==1 || !is.matrix(tlist$time)) {
			res <- integrand(t=tlist$time,initial,tlist$state,Tt2=summ[tlist$edge,1],Tt1=summ[tlist$edge,2],pars)
		} else {
  	  res <- sapply(1:dim(tlist$time)[1], function (z) integrand(t=tlist$time[z,],initial,tlist$state,Tt2=summ[tlist$edge,1],Tt1=summ[tlist$edge,2],pars))
  		res1 <- res[1,]^2/colSums(res)+res[2,]^2/colSums(res)
  		res <- c(sum(res[1,]*res1),sum(res[2,]*res1))/sum(res1)
  	}
  }
  if (tlist$type==2) {
  	tips <- tlist$edge[which(substr(tlist$edge,1,2)=="sp")]
  	tmpD <- matrix(NA,length(tlist$edge),2)
  	tmpE <- matrix(NA,length(tlist$edge),2)
  	D <- matrix(NA,length(tlist$edge),2)
  	rownames(tmpD) <- tlist$edge
  	rownames(tmpE) <- tlist$edge
  	rownames(D) <- tlist$edge
  	anc <- NULL
  	for (i in tips) {
  		initial <- c(n0=as.numeric(summ[i,"d0"]),n1=as.numeric(summ[i,"d1"]),e0=as.numeric(summ[i,"e0"]),e1=as.numeric(summ[i,"e1"]))
  		tmp <- ode(initial,as.numeric(c(summ[i,"t1"],summ[i,"t2"])),calD,pars)
    	tmpD[i,] <- tmp[2,2:3]
  		tmpE[i,] <- tmp[2,4:5]
			tmpD[i,tmpD[i,]<0] <- 10^-6
			tmpE[i,tmpE[i,]<0] <- 0
			tmpt <- tlist$time[[which(tlist$edge==i)]]
			if (length(tmpt)==1) {
				D[i,] <- integrand2(t=tmpt,initial,tlist$state,Tt2=summ[i,1],Tt1=summ[i,2],pars)
			} else {
				res <- sapply(1:dim(tmpt)[1], function (z) integrand2(t=tmpt[z,],initial,tlist$state,Tt2=summ[i,1],Tt1=summ[i,2],pars))
				res1 <- res[1,]^2/colSums(res)+res[2,]^2/colSums(res)
  			D[i,] <- c(sum(res[1,]*res1),sum(res[2,]*res1))/sum(res1)
			}
			anc <- c(anc,rownames(edge)[which(edge[,2]==edge[i,1])])
  	}
  	anc1 <- anc[duplicated(anc)]
  	while (length(anc1)>0) {
  		for (i in anc1) {
  			des <- rownames(edge)[which(edge[,1]==edge[i,2])]
  			initial <- c(n0=as.numeric(prod(tmpD[des,1])*pars[1]),n1=as.numeric(prod(tmpD[des,2])*pars[2]),e0=as.numeric(tmpE[des[1],1]),e1=as.numeric(tmpE[des[1],2]))
  			if (i!=tlist$edge[1]) {
      		tmp <- ode(initial,as.numeric(c(summ[i,"t1"],summ[i,"t2"])),calD,pars)
      		tmpD[i,] <- tmp[2,2:3]
  	  		tmpE[i,] <- tmp[2,4:5]
      		tmpD[i,tmpD[i,]<0] <- 10^-6
	  			tmpE[i,tmpE[i,]<0] <- 0
	  		}
	  		tmpt <- tlist$time[[which(tlist$edge==i)]]
			  if (length(tmpt)==1) {
				  D[i,] <- integrand2(t=tmpt,initial,tlist$state,Tt2=summ[i,1],Tt1=summ[i,2],pars)
			  } else {
				  res <- sapply(1:dim(tmpt)[1], function (z) integrand2(t=tmpt[z,],initial,tlist$state,Tt2=summ[i,1],Tt1=summ[i,2],pars))
				  res1 <- res[1,]^2/colSums(res)+res[2,]^2/colSums(res)
  				D[i,] <- c(sum(res[1,]*res1),sum(res[2,]*res1))/sum(res1)
  			}
		  	anc <- c(anc[-which(anc==i)],rownames(edge)[which(edge[,2]==edge[i,1])])
  	  }
  		anc1 <- anc[duplicated(anc)]
  	}
  	for (i in tlist$edge[-1]) {
  		anc <- rownames(edge)[which(edge[,2]==edge[i,1])]
  		sis <- rownames(edge)[which(edge[,1]==edge[i,1])]
  		sis <- sis[sis!=i]
  		initial <- c(n0=as.numeric(tmpD[sis,1]*D[i,1]*pars[1]),n1=as.numeric(tmpD[sis,2]*D[i,2]*pars[2]),e0=as.numeric(tmpE[sis,1]),e1=as.numeric(tmpE[sis,2]))
  		tmp <- ode(initial,as.numeric(c(summ[anc,"t1"],summ[anc,"t2"])),calD,pars)
  		D[i,] <- tmp[2,2:3]
    	D[i,D[i,]<0] <- 10^-6
    	while(anc!=tlist$edge[1]) {
    		anc1 <- anc
    		anc <- rownames(edge)[which(edge[,2]==edge[anc1,1])]
  			sis <- rownames(edge)[which(edge[,1]==edge[anc,1])]
  			sis <- sis[sis!=anc1]
  			initial <- c(n0=as.numeric(tmpD[sis,1]*D[i,1]*pars[1]),n1=as.numeric(tmpD[sis,2]*D[i,2]*pars[2]),e0=as.numeric(tmpE[sis,1]),e1=as.numeric(tmpE[sis,2]))
  			tmp <- ode(initial,as.numeric(c(summ[anc,"t1"],summ[anc,"t2"])),calD,pars)
    		D[i,] <- tmp[2,2:3]
    		D[i,D[i,]<0] <- 10^-6
      }
  	}
  	res1 <- D[,1]^2/rowSums(D)+D[,2]^2/rowSums(D)
  	res <- c(sum(D[,1]*res1),sum(D[,2]*res1))/sum(res1)
  }
  if (tlist$type==3) {
  	tips <- tlist$edge[which(substr(tlist$edge,1,2)=="sp")]
  	tmpD <- matrix(NA,length(tlist$edge),2)
  	tmpE <- matrix(NA,length(tlist$edge),2)
  	tmpD2 <- matrix(NA,length(tlist$edge),2)
  	tmpE2 <- matrix(NA,length(tlist$edge),2)
  	D <- matrix(NA,length(tlist$edge),2)
  	rownames(tmpD) <- tlist$edge
  	rownames(tmpE) <- tlist$edge
  	rownames(tmpD2) <- tlist$edge
  	rownames(tmpE2) <- tlist$edge
  	rownames(D) <- tlist$edge
  	anc <- NULL
  	for (i in tips) {
  		initial <- c(n0=as.numeric(summ[i,"d0"]),n1=as.numeric(summ[i,"d1"]),e0=as.numeric(summ[i,"e0"]),e1=as.numeric(summ[i,"e1"]))
  		tmp <- ode(initial,as.numeric(c(summ[i,"t1"],summ[i,"t2"])),calD,pars)
    	tmpD[i,] <- tmp[2,2:3]
  		tmpE[i,] <- tmp[2,4:5]
			tmpD[i,tmpD[i,]<0] <- 10^-6
			tmpE[i,tmpE[i,]<0] <- 0
			if (length(tlist$edge)==1) {
				tmpt <- tlist$time1
				if (length(tmpt)==1 || !is.matrix(tmpt)) {
					D[i,] <- integrand2(t=tmpt,initial,tlist$state,Tt2=summ[i,1],Tt1=summ[i,2],pars)
				} else {
					res <- sapply(1:dim(tmpt)[1], function (z) integrand2(t=tmpt[z,],initial,tlist$state,Tt2=summ[i,1],Tt1=summ[i,2],pars))
					res1 <- res[1,]^2/colSums(res)+res[2,]^2/colSums(res)
  				D[i,] <- c(sum(res[1,]*res1),sum(res[2,]*res1))/sum(res1)
				}
			} else {
				tmpt <- tlist$time[[which(tlist$edge==i)-1]]
				if (length(tmpt)==1) {
					D[i,] <- integrand2(t=tmpt,initial,tlist$state[1,],Tt2=summ[i,1],Tt1=summ[i,2],pars)
				} else {
					res <- sapply(1:dim(tmpt)[1], function (z) integrand2(t=tmpt[z,],initial,tlist$state[1,],Tt2=summ[i,1],Tt1=summ[i,2],pars))
					res1 <- res[1,]^2/colSums(res)+res[2,]^2/colSums(res)
  				D[i,] <- c(sum(res[1,]*res1),sum(res[2,]*res1))/sum(res1)
				}
			}
			anc <- c(anc,rownames(edge)[which(edge[,2]==edge[i,1])])
  	}
  	anc1 <- anc[duplicated(anc)]
  	while (length(anc1)>0) {
  		for (i in anc1) {
  			des <- rownames(edge)[which(edge[,1]==edge[i,2])]
  			initial <- c(n0=as.numeric(prod(tmpD[des,1])*pars[1]),n1=as.numeric(prod(tmpD[des,2])*pars[2]),e0=as.numeric(tmpE[des[1],1]),e1=as.numeric(tmpE[des[1],2]))
  			if (i!=tlist$edge[1]) {
      		tmp <- ode(initial,as.numeric(c(summ[i,"t1"],summ[i,"t2"])),calD,pars)
      		tmpD[i,] <- tmp[2,2:3]
  	  		tmpE[i,] <- tmp[2,4:5]
      		tmpD[i,tmpD[i,]<0] <- 10^-6
	  			tmpE[i,tmpE[i,]<0] <- 0
	  			tmpt <- tlist$time[[which(tlist$edge==i)-1]]
				  if (length(tmpt)==1) {
					  D[i,] <- integrand2(t=tmpt,initial,tlist$state[1,],Tt2=summ[i,1],Tt1=summ[i,2],pars)
				  } else {
					  res <- sapply(1:dim(tmpt)[1], function (z) integrand2(t=tmpt[z,],initial,tlist$state[1,],Tt2=summ[i,1],Tt1=summ[i,2],pars))
					  res1 <- res[1,]^2/colSums(res)+res[2,]^2/colSums(res)
  					D[i,] <- c(sum(res[1,]*res1),sum(res[2,]*res1))/sum(res1)
  				}
	  		} else {
	  			tmpt <- tlist$time1
	  			if (length(tmpt)==1 || !is.matrix(tmpt)) {
					  D[i,] <- integrand2(t=tmpt,initial,tlist$state,Tt2=summ[i,1],Tt1=summ[i,2],pars)
				  } else {
					  res <- sapply(1:dim(tmpt)[1], function (z) integrand2(t=tmpt[z,],initial,tlist$state,Tt2=summ[i,1],Tt1=summ[i,2],pars))
					  res1 <- res[1,]^2/colSums(res)+res[2,]^2/colSums(res)
  					D[i,] <- c(sum(res[1,]*res1),sum(res[2,]*res1))/sum(res1)
  				}
			  }
			  anc <- c(anc[-which(anc==i)],rownames(edge)[which(edge[,2]==edge[i,1])])
  		}
  		anc1 <- anc[duplicated(anc)]
  	}
  	for (i in tlist$edge[-1]) {
  		anc <- rownames(edge)[which(edge[,2]==edge[i,1])]
  		sis <- rownames(edge)[which(edge[,1]==edge[i,1])]
  		sis <- sis[sis!=i]
  		initial <- c(n0=as.numeric(tmpD[sis,1]*D[i,1]*pars[1]),n1=as.numeric(tmpD[sis,2]*D[i,2]*pars[2]),e0=as.numeric(tmpE[sis,1]),e1=as.numeric(tmpE[sis,2]))
  		tmp <- ode(initial,as.numeric(c(summ[anc,"t1"],summ[anc,"t2"])),calD,pars)
      D[i,] <- tmp[2,2:3]
      D[i,D[i,]<0] <- 10^-6
      while(anc!=tlist$edge[1]) {
      	anc1 <- anc
      	anc <- rownames(edge)[which(edge[,2]==edge[anc1,1])]
  			sis <- rownames(edge)[which(edge[,1]==edge[anc,1])]
  			sis <- sis[sis!=anc1]
  			initial <- c(n0=as.numeric(tmpD[sis,1]*D[i,1]*pars[1]),n1=as.numeric(tmpD[sis,2]*D[i,2]*pars[2]),e0=as.numeric(tmpE[sis,1]),e1=as.numeric(tmpE[sis,2]))
  			if (anc==tlist$edge[1]) {
  				tmpt <- tlist$time2
  				if (length(tmpt)==1) {
						D[i,] <- integrand2(t=tmpt,initial,tlist$state[-1,],Tt2=summ[i,1],Tt1=summ[i,2],pars)
					} else {
  					res <- sapply(1:dim(tmpt)[1], function (z) integrand2(t=tmpt[z,],initial,tlist$state[-1,],Tt2=summ[i,1],Tt1=summ[i,2],pars))
  					res1 <- res[1,]^2/colSums(res)+res[2,]^2/colSums(res)
  					D[i,] <- c(sum(res[1,]*res1),sum(res[2,]*res1))/sum(res1)
  				}
  			} else {
  				tmp <- ode(initial,as.numeric(c(summ[anc,"t1"],summ[anc,"t2"])),calD,pars)
      		D[i,] <- tmp[2,2:3]
      		D[i,D[i,]<0] <- 10^-6
      	}
      }
  	}
  	res1 <- D[,1]^2/rowSums(D)+D[,2]^2/rowSums(D)
  	res <- c(sum(D[,1]*res1),sum(D[,2]*res1))/sum(res1)
  }
  if (is.na(res[1])||is.infinite(res[1])) {res[1] <- 10^-6}
  if (is.na(res[2])||is.infinite(res[2])) {res[2] <- 10^-6}
  res
  }
  res <- 0
  ntips <- length(tree$tip.label)
  nedge <- length(tree$edge[,1])
  rname <- sapply(tree$edge[,2],function (j) ifelse(j>ntips,tree$node.label[j-ntips],tree$tip.label[j]))
  rownames(tree$edge) <- rname
  node.height <- tree$edge+NA
  anc <- which(tree$edge[,2]<=ntips)
  node.height[anc,2] <- 0
  node.height[anc,1] <- tree$edge.length[anc]
  anc2 <- as.numeric(sapply(anc,function (i) which(tree$edge[,2]==tree$edge[i,1])))
  while(is.na(sum(node.height[anc2[!is.na(anc2)],2]))) {
    node.height[anc2[!is.na(anc2)],2] <- node.height[anc[!is.na(anc2)],1]
    node.height[anc2[!is.na(anc2)],1] <- node.height[anc2[!is.na(anc2)],2]+tree$edge.length[anc2[!is.na(anc2)]]
    anc <- anc2[!is.na(anc2)]
    anc2 <- as.numeric(sapply(anc[!is.na(anc)],function (i) which(tree$edge[,2]==tree$edge[i,1])))
  }
  if (!is.null(tlist)) {
  	nam <- unique(unlist(sapply(tlist,function (x) x$edge)))
  } else {
  	nam <- NULL
  }
  if (dim(sfraction)[1] < nedge) {
    tmp <- matrix(0,nedge,5)
    rownames(tmp) <- rname
    tmp[rownames(sfraction),] <- sfraction
    sfraction <- tmp
  }
  summ <- cbind(node.height,sapply(1:nedge, function (i) as.numeric(ifelse(sfraction[i,1]==0,sfraction[i,2]/sfraction[i,3],0))),sapply(1:nedge, function (i) as.numeric(ifelse(sfraction[i,1]==1,sfraction[i,4]/sfraction[i,5],0))),sapply(1:nedge, function (i) as.numeric(ifelse(sfraction[i,3]>0,1-sfraction[i,2]/sfraction[i,3],0))),sapply(1:nedge, function (i) as.numeric(ifelse(sfraction[i,5]>0,1-sfraction[i,4]/sfraction[i,5],0))),sapply(1:nedge, function (i) is.element(rname[i],nam)))
  colnames(summ) <- c("t2","t1","d0","d1","e0","e1","m")
  D <- tree$edge+NA
  E <- tree$edge+NA
  tips <- intersect(rname[summ[,"m"]==0],tree$tip.label)
  for (i in tips) {
    tmp <- ode(c(n0=as.numeric(summ[i,"d0"]),n1=as.numeric(summ[i,"d1"]),e0=as.numeric(summ[i,"e0"]),e1=as.numeric(summ[i,"e1"])),as.numeric(c(summ[i,"t1"],summ[i,"t2"])),calD,pars)
    D[i,] <- tmp[2,2:3]
  	E[i,] <- tmp[2,4:5]
	  D[i,D[i,]<0] <- 10^-6
	  E[i,E[i,]<0] <- 0
	  res <- res+log(sum(D[i,]))
	  D[i,] <- D[i,]/sum(D[i,])
  }
  anc <- tree$node.label[tree$edge[tips,1]-ntips]
  idx <- anc[duplicated(anc)]
  while(length(idx)>0&&idx!=tree$node.label[1]) {
    for (i in idx) {
    	if (!is.element(i,nam)) {
      	des <- rname[which(tree$edge[,1]==tree$edge[i,2])]
      	if (as.numeric(E[des[1],1])==as.numeric(E[des[2],1]) && as.numeric(E[des[1],2])==as.numeric(E[des[2],2])) {
      		tmp <- ode(c(n0=as.numeric(prod(D[des,1])*pars[1]),n1=as.numeric(prod(D[des,2])*pars[2]),e0=as.numeric(E[des[1],1]),e1=as.numeric(E[des[1],2])),as.numeric(c(node.height[i,2],node.height[i,1])),calD,pars)
      	} else {
      		Etmp <- as.numeric(ode(c(e0=0,e1=0),c(0,summ[i,"t1"]),calE,pars)[2,2:3])
    		  Etmp[Etmp<0] <- 0
			    tmp <- ode(c(n0=as.numeric(prod(D[des,1])*pars[1]),n1=as.numeric(prod(D[des,2])*pars[2]),e0=as.numeric(Etmp[1]),e1=as.numeric(Etmp[2])),as.numeric(c(node.height[i,2],node.height[i,1])),calD,pars)
      	}
      	D[i,] <- tmp[2,2:3]
  	  	E[i,] <- tmp[2,4:5]
      	D[i,D[i,]<0] <- 10^-6
	  	  E[i,E[i,]<0] <- 0
	  	  res <- res+log(sum(D[i,]))
	  	  D[i,] <- D[i,]/sum(D[i,])
	    } else {
	  	  des <- rname[which(tree$edge[,1]==tree$edge[i,2])]
    	  Etmp <- as.numeric(ode(c(e0=0,e1=0),c(0,summ[i,"t1"]),calE,pars)[2,2:3])
    	  Etmp[Etmp<0] <- 0
    	  initial <- c(n0=as.numeric(prod(D[des,1])*pars[1]),n1=as.numeric(prod(D[des,2])*pars[2]),e0=as.numeric(Etmp[1]),e1=as.numeric(Etmp[2]))
    	  idx1 <- which(unlist(sapply(tlist, function (x) is.element(i,x$edge))))
    	  D[i,] <- intD(pars,tlist[[idx1]],summ,tree$edge,initial)
     	  D[i,D[i,]<0] <- 10^-6
    	  res <- res+log(sum(D[i,]))
    	  D[i,] <- D[i,]/sum(D[i,])
    	  E[i,] <- as.numeric(ode(c(e0=0,e1=0),c(0,summ[i,"t2"]),calE,pars)[2,2:3])
    	  E[i,E[i,]<0] <- 0  	
	    }
      anc <- anc[-which(anc==i)]
      anc <- c(anc,tree$node.label[tree$edge[i,1]-ntips])
    }
    idx <-  anc[duplicated(anc)]
  }
  if (!is.null(tlist)) {
  tmpD <- matrix(NA,length(tlist),2)
  tmpE <- matrix(NA,length(tlist),2)
  nodes <- NULL
  for (i in 1:length(tlist)) {
  	is.tips <- which(substr(tlist[[i]]$edge,1,2)=="sp")
  	if (length(is.tips)>0) {
  	  tmp <- intD(pars,tlist[[i]],summ,tree$edge)
  	  tmpD[i,] <- tmp
  	  res <- res+log(sum(tmpD[i,]))
  	  tmpD[i,] <- tmpD[i,]/sum(tmpD[i,])
  	  nodes <- c(nodes,tlist[[i]]$edge[1])
  	  tmp <- as.numeric(ode(c(e0=0,e1=0),c(0,summ[tlist[[i]]$edge[1],"t2"]),calE,pars)[2,2:3])
  	  tmp[tmp<0] <- 0
  	  tmpE[i,] <- tmp
  	}
  }
  idx <- which(!is.na(tmpD[,1]))
  D[nodes,] <- tmpD[idx,]
  E[nodes,] <- tmpE[idx,]
  anc <- c(anc,tree$node.label[tree$edge[nodes,1]-ntips])
  idx <- anc[duplicated(anc)]
  if (is.element(tree$node.label[1],idx)) {idx <- idx[-which(idx==tree$node.label[1])]}
  while(length(idx)>0) {
    for (i in idx) {
    	des <- rname[which(tree$edge[,1]==tree$edge[i,2])]
    	E <- as.numeric(ode(c(e0=0,e1=0),c(0,summ[i,"t1"]),calE,pars)[2,2:3])
    	E[E<0] <- 0
    	initial <- c(n0=as.numeric(prod(D[des,1])*pars[1]),n1=as.numeric(prod(D[des,2])*pars[2]),e0=as.numeric(E[1]),e1=as.numeric(E[2]))
    	if (is.element(i,nam)) {
    		idx1 <- which(unlist(sapply(tlist, function (x) is.element(i,x$edge))))
    		D[i,] <- intD(pars,tlist[[idx1]],summ,tree$edge,initial)
    	} else {
    		D[i,] <- as.numeric(ode(initial,as.numeric(c(node.height[i,2],node.height[i,1])),calD,pars)[2,2:3])
    	}
    	D[i,D[i,]<0] <- 10^-6
    	res <- res+log(sum(D[i,]))
    	D[i,] <- D[i,]/sum(D[i,])   	
	  	anc <- anc[-which(anc==i)]
  		anc <- c(anc,tree$node.label[tree$edge[i,1]-ntips])
	  }
	  idx <-  anc[duplicated(anc)]
	  if (is.element(tree$node.label[1],idx)) {idx <- idx[-which(idx==tree$node.label[1])]}
  }
  }	
  res1 <- D[which(tree$edge[,1]==ntips+1),]
  if (!is.null(root.state)) {
    res1 <- sum(c(res1[1,1]*res1[2,1],res1[1,2]*res1[2,2])*root.state)
  } else {
    res1 <- c(res1[1,1]*res1[2,1],res1[1,2]*res1[2,2])
    res1 <- res1[1]^2/sum(res1)+res1[2]^2/sum(res1)
  }
  res <- -log(res1)-res
  if (is.na(res)||is.infinite(res)||res<0) {res <- 10000} 
  res
}
