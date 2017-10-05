#' \code{tlistcal} lists all the possible locations where an un-sampled group or state can be attached to a tree
#' 
#' @param tree a phylogeny (of class phylo)
#' @param ssize a list consists of the sampling pattern in the tree. List unsampled state first, then unsampled groups.
#'        $edge gives the node or tip name where the unsampled group or state can be attached to the tree. Edges that are closer to the root are list first.
#'        $state gives the initial D and E values for the unsampled group or state, with names "n01" for the initial value of D0, "n11" for the initial value of D1, "e01" for the initial value of E0, "e11" for the initial value of E1
#'        If there are more than one unsampled groups, $state is a matrix, where each row corresponds to each unsampled group.
#'        $ordered gives whether the unsampled groups are in ordered location on the $edge, that is, if one unsampled group is clustered with other tips before another unsampled group
#'        If $ordered = TRUE, then first row in $state is the first unsampled group to be clustered with other tips in the tree
#' @param dt the time interval between two possible location where an unsampled group or state can be attached to the tree. The defaul is 1 unit of branch length
#' @return tlist input to likcal
#' @note the current version only accept one edge on which unsampled groups can be attached, that is the full topology of unsampled groups are known. Future version will allow incomplete topology
#' 
tlistcal <- function (tree,ssize,dt=1) {
	sub <- function(n, r, v, N) {
            if (r == 1) 
                matrix(v, n, 1)
            else if (n == 1) 
                matrix(v, 1, r)
            else {
                inner <- Recall(n, r - 1, v, N)
                out <- cbind(rep(v, rep(nrow(inner), n)), matrix(t(inner), 
                  ncol = ncol(inner), nrow = nrow(inner) * n, 
                  byrow = TRUE))
                out[rowSums(out)<=N,]
            }
        }
  if (length(ssize)>0) {
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
    nam <- unlist(sapply(ssize,function (x) x$edge))
    dupedge <- nam[which(duplicated(nam))]
	  tlist <- vector("list",length(ssize))
	  for (i in length(ssize):1) {
		if (length(ssize[[i]])>1) {
		  tmp <- is.element(ssize[[i]]$edge,dupedge)
		  if (sum(tmp)==0 && ssize[[i]]$ordered) {
			  tlist[[i]]$type <- 1 #1=groups are state, 2=one state is not sampled in a group, 3=groups are state and one state is not sampled
			  tlist[[i]]$edge <- ssize[[i]]$edge
  	    tlist[[i]]$state <- ssize[[i]]$state
			  Tt2 <- as.numeric(node.height[tlist[[i]]$edge,1])
			  Tt1 <- as.numeric(node.height[tlist[[i]]$edge,2])
			  z <- length(ssize[[i]]$state)/4
  		  N <- floor ((Tt2-Tt1)/dt)
  			if (N<=(z+1)) {
  				tlist[[i]]$time <- rep((Tt2-Tt1)/(z+1),z)
  			} else {
  				if (z==1) {
  					tlist[[i]]$time <- as.matrix(c(1:(N-1))*dt,1,(N-1))
				} else {
  					tmp <- sub(N-z,z,c(0:(N-z-1)),N-z-1)
  					tmp <- tmp + 1
  					tlist[[i]]$time <- tmp*dt
  					if (is.matrix(tlist[[i]]$time)) {
  						while (dim(tlist[[i]]$time)[1]>1000) {
  							N <- N-1
  							dt <- (Tt2-Tt1)/N
  							tmp <- sub(N-z,z,c(0:(N-z-1)),N-z-1)
  							tmp <- tmp + 1
  							tlist[[i]]$time <- tmp*dt
  						}
  					}
  				}
  			}
  		} else if (sum(tmp)==0 && !ssize[[i]]$ordered) {
  			tlist[[i]]$type <- 2
  			tlist[[i]]$edge <- ssize[[i]]$edge
 			  tlist[[i]]$state <- ssize[[i]]$state
			  tlist[[i]]$time <- vector("list",length(tlist[[i]]$edge))
			  for (j in 1:length(tlist[[i]]$edge)) {
				  Tt2 <- as.numeric(node.height[tlist[[i]]$edge[j],1])
				  Tt1 <- as.numeric(node.height[tlist[[i]]$edge[j],2])
  				N <- floor ((Tt2-Tt1)/dt)
  				if (N<=2) {
  					tlist[[i]]$time[[j]] <- (Tt2-Tt1)/2
  				} else {
  					tlist[[i]]$time[[j]] <- as.matrix(c(1:(N-1))*dt,1,(N-1))
  				}
  			}
  		} else {
  			tlist[[i]]$type <- 3
  			idx <- sapply(ssize, function (x) is.element(ssize[[i]]$edge[1],x[[1]])*1)
			  idx <- which(idx==1)[1]
  			tlist[[i]]$edge <- ssize[[idx]]$edge # the first edge is always the root of a group
  			tlist[[i]]$state <- rbind(ssize[[idx]]$state,ssize[[i]]$state) #state state first, then state group
			  #when all the state lineages are in the root of the group
			  Tt2 <- as.numeric(node.height[tlist[[i]]$edge[1],1])
			  Tt1 <- as.numeric(node.height[tlist[[i]]$edge[1],2])
			  z <- length(ssize[[i]]$state)/4+1
  			N <- floor ((Tt2-Tt1)/dt)
  			if (N<=(z+1)) {
  				tlist[[i]]$time1 <- rep((Tt2-Tt1)/(z+1),z)
  			} else {
  				if (z==1) {
  					tlist[[i]]$time1 <- as.matrix(c(1:(N-1))*dt,1,(N-1))
				} else {
  					tmp <- sub(N-z,z,c(0:(N-z-1)),N-z-1)
  					tmp <- tmp + 1
  					tlist[[i]]$time1 <- tmp*dt
  					if (is.matrix(tlist[[i]]$time1)) {
  						while (dim(tlist[[i]]$time1)[1]>1000) {
  							N <- N-1
  							dt <- (Tt2-Tt1)/N
  							tmp <- sub(N-z,z,c(0:(N-z-1)),N-z-1)
  							tmp <- tmp + 1
  							tlist[[i]]$time1 <- tmp*dt
  						}
  					}
  				}
  			}
  			if (length(tlist[[i]]$edge)>1) {
  			#when state state is not in the root of the group
			  z <- length(ssize[[i]]$state)/4
  			N <- floor ((Tt2-Tt1)/dt)
  			if (N<=(z+1)) {
  				tlist[[i]]$time2 <- rep((Tt2-Tt1)/(z+1),z)
  			} else {
  				if (z==1) {
  					tlist[[i]]$time2 <- as.matrix(c(1:(N-1))*dt,1,(N-1))
				} else {
  					tmp <- sub(N-z,z,c(0:(N-z-1)),N-z-1)
  					tmp <- tmp + 1
  					tlist[[i]]$time2 <- tmp*dt
  					if (is.matrix(tlist[[i]]$time2)) {
  						while (dim(tlist[[i]]$time2)[1]>1000) {
  							N <- N-1
  							dt <- (Tt2-Tt1)/N
  							tmp <- sub(N-z,z,c(0:(N-z-1)),N-z-1)
  							tmp <- tmp + 1
  							tlist[[i]]$time2 <- tmp*dt
  						}
  					}
  				}
  			}
  			#now calcualte the possible location for the state state when not in the root
  			tlist[[i]]$time <- vector("list",length(tlist[[i]]$edge)-1)
			  for (j in 2:length(tlist[[i]]$edge)) {
				  Tt2 <- as.numeric(node.height[tlist[[i]]$edge[j],1])
				  Tt1 <- as.numeric(node.height[tlist[[i]]$edge[j],2])
  				N <- floor ((Tt2-Tt1)/dt)
  				if (N<=2) {
  					tlist[[i]]$time[[j-1]] <- (Tt2-Tt1)/2
  				} else {
  					tlist[[i]]$time[[j-1]] <- as.matrix(c(1:(N-1))*dt,1,(N-1))
  				}
  			}
  			}
  			ssize[[idx]] <- 1
  		}
  	}
  	}
  	idx <- which(unlist(sapply(tlist,function (x) is.null(x))))
  	if (length(idx)>0) {tlist <- tlist[-idx]}
	tlist
  }
}
