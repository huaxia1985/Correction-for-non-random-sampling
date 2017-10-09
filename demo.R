#This is the demo on how to use our correction for the impact of non-random sampling on the performance of BiSSE framework
#The first part is used to generate the results in our paper "using simulated tree "The influence of non-random species-sampling on macroevolutionary and macroecological inference from phylogenies"
#The second part gives instructions on how to use our codes to account for non-random sampling using users' own tree and sampling pattern.

######## PART I ##########
library(diversitree)
library(deSolve)
library(nloptr)

#Here lists all the parameter sets used in our study (Table 1).
p <- c(0.1,0.1,0.03,0.03,0.01,0.001)
#p <- c(0.3,0.1,0.03,0.03,0.01,0.001)
#p <- c(0.1,0.1,0.01,0.03,0.01,0.001)
#p <- c(0.03,0.1,0.03,0.03,0.01,0.001)
#p <- c(0.1,0.1,0.09,0.03,0.01,0.001)
#p <- c(0.1,0.1,0.03,0.03,0.01,0.01)
#p <- c(0.3,0.1,0.03,0.03,0.01,0.01)
#p <- c(0.1,0.1,0.01,0.03,0.01,0.01)
#p <- c(0.03,0.1,0.03,0.03,0.01,0.01)
#p <- c(0.1,0.1,0.09,0.03,0.01,0.01)
px <- c(p[1]-p[3],p[2]-p[4],p[3],p[4],p[5],p[6])

#simulate a tree, here the tree has root state 0 and 500 extant tips. A group is defined by the "Time" scheme. All the groups are sampled in the tree. A random integer between 1 and 5 representatives are sampled for each group.
treelist <- try(tree.gen(pars=p,max.taxa=500,x0=0,sampling.c=1,sampling.f=5,define.clade="time"),silent=T)
#the while loop is to make sure that the tree does not lead to extinction
while(inherits(treelist,"try-error")) {
  treelist <- try(tree.gen(pars=p,max.taxa=500,x0=0,sampling.c=1,sampling.f=5,define.clade="time"),silent=T)
}

#fit BiSSE model to the tree using our correction
tlist <- tlistcal(tree=treelist$xtree,ssize=treelist$ssize,dt=1)
xfit <- sbplx(px,likcal,tree=treelist$xtree,sfraction=treelist$sfraction,root.state=NULL,tlist=tlist,lower=rep(0,length(p)),upper=rep(Inf,length(p)))

#fit BiSSE model to the tree using the unresolved clade correction
lik <- make.bisse(treelist$ftree,treelist$ftree$tip.state,unresolved=as.data.frame(treelist$unresolved))
ffit <- find.mle(lik,p)

pct <- proc.time()
#calculate the likelihood of the true model using our correction
xtrue <- likcal(p=px,tree=treelist$xtree,sfraction=treelist$sfraction,root.state=NULL,tlist=tlist)
#calculate the computing time of our correction
xpct <- pct-proc.time()

pct <- proc.time()
#calculate the likelihood of the true model using the unresolved clade correction
ftrue <- lik(p)
#calculate the computing time of the unresolved clade correction
fpct <- pct-proc.time()

#fit constrained BiSSE model to the tree using our correction
xfix <- sbplx(c(0.07,0.03,p[5],p[6]),likcal,tree=treelist$xtree,sfraction=treelist$sfraction,root.state=NULL,tlist=tlist,constrained=TRUE,lower=rep(0,4),upper=rep(Inf,4))

#fit constrainted BiSSE model to the tree using the unresolved clade correction
lik2 <- constrain(lik,lambda0~lambda1,mu0~mu1)
ffix <- find.mle(lik2,c(0.1,0.1,0.03,0.03,p[5],p[6]))

######## PART II ##########
library(deSolve)
library(nloptr)

#This is the demo on how to formulate user's own sampling pattern to apply our correction for non-random sampling
#To do this, users need to prepare:
#1) tree: a tree of class phylo that consists of the sampled tips in the tree
#2) sfraction: a matrix that consists of the sampling fraction for each tip of the tree. Make sure that the rownames of the matrix are the tipnames in the tree
#   The first column is the state of the tip, 0 for state 0 and 1 for state 1.
#   The second column is the number of sampled tips with state 0 in the group that the tip is sampled as a representative.
#   The third column is the total number of extant members with state 0 in the group.
#   The fourth column is the number of sampled tips with state 1 in the group that the top is sampled as a representative.
#   The fifth column is the total number of extant members with state 1 in the group.
#   For example:
#   If your tree consists of 10 tips, than the matrix should contain 10 rows, each row has the name corresponding to each tip.
#   If "sp 1" and "sp 2" are representatives for the same group and if sp 1 is in state 0 and sp 2 in state 1, so the number of sampled tips in that group for state 0 is 1 and the number of sampled tips for state 1 is 1.
#   If the total number of members in that group include 5 members with state 0 and 6 members with state 1, then the row for sp 1 and sp 2 looks like:
#   "sp 1" 0 1 5 1 6
#   "sp 2" 1 1 5 1 6
#3) ssize: a list that consists of all the information for unsampled state or group
#   Please list unsampled state first, then unsampled group
#   For example:
#   If your tree has a group with 2 sampled tips "sp 1" and "sp 2", both in state 1 and joins at "nd 1". The total number of members in that group include 5 members with state 0 and 6 members with state 1
#   So state 0 is missing in that group, then create a cell in the list that includes the following elments:
#   $edge: a vector of all the edges that the unsampled state can be attached. In this example, the vector is c("nd 1", "sp 1", "sp 2").
#   For internal branch, edge name is the node.label of its descendent node in the tree.
#   $state: a vector of the initial D and E values for the unsampled state. In this example, the vector is c(n01=1/5,n11=0,e01=1-1/5,e11=1-2/6). Here n01 is the initial D value for state 0, n11 is the initial D value for state 1, e01 is the initial E value for state 0, e11 is the initial E value for state 1
#   $ordered: set to FALSE for unsampled state
#   Now if you have listed all the unsampled state, let's start listing the unsampled groups
#   For the current version, we assume full topology is known for the relationship among groups. Later version will generalize the code to account for incomplete topological information
#   For example:
#   If your tree has two unsampled groups a and b and group a is clustered with some other sampled group first, which is then clustered with group b, and if both groups are attached to edge "nd 1"
#   $edge = "nd 1"
#   If group a has 5 extant members with state 0 and 6 members with state 1, and group b has 0 extant members with state 0 and 3 members with state 1
#   $state is a matrix where the first row is group a, the group that clustered with other sampled group first; the second row is group b.
#   $state looks like:
#   n01        n11        e01          e11
#   5/(5+6)^2  6/(5+6)^2  1-5/(5+6)^2  1-6/(5+6)^2
#   0          1/3        0            1-1/3
#   $ordered = TRUE, this is to state that the location of group a and group b should be ordered on edge "nd 1" because group a is clustered with other sampled groups first

# Now use tlistcal to list all possible locations that unsampled state or group can be attached to the tree
# tlist <- tlistcal(tree=your tree,ssize=your ssize list,dt=1)
# Here, dt is the time interval between two locations and default is 1 unit of branch length. The smaller the inveral is, the more accurate the likelihod value, but takes longer to compute

# Now fit BiSSE model to your tree using our correction for your sampling pattern
# xfit <- sbplx(p,likcal,tree=your tree,sfraction=your sfraction matrix,ssize=your ssize list,root.state=NULL,tlist=your tlistcal output,constrained=FALSE,lower=rep(0,length(p)),upper=rep(Inf,length(p)))
# Here, p is the starting parameter values to search for the best-fit model
#       root.state is the state of the root. If root.state=NULL, the solution in appendix 1 in FitzJohn et al. 2009 is used to calcualte the overall likelihood of the tree
#       constrained is whether you want to contrain speciation rate and extinction rate to be equal between different states. The default is FALSE.
#       sbplx is the function that applies subplex method to search for the maximum likelihood, with lower and higher defines the boundary of the parameter space to search
#       users can use other method by calling likcal function if the subplex method does not converge well.

