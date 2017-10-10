# Correction-for-non-random-sampling
Correction for non-random sampling under the BiSSE framework

This contains the codes used in paper:
"The influence of non-random species-sampling on macroevolutionary and macroecological inference from phylogenies"
by Xia Hua and Robert Lanfear

Please cite the paper if you use the codes.

The codes implement our proposed correction for non-random sampling under the BiSSE framework in R.
Each function is documented in the beginning of the file, explaining each input and output of the function.

The following instructions are on:
1) use the code to generate the results in the paper
2) use the code to account for non-random sampling in user's own tree and sampling pattern.

######## Use the code to generate the results in the paper ##########

Step 1: copy all the functions in the folder to R and load the related library
library(diversitree)
library(deSolve)
library(nloptr)

Step 2: specify a set of parameter values that include six parameters in order of speciation rate of state 0, speciation rate of state 1, extinction rate of state 0, extinction rate of state 1, transition rate from state 0 to state 1, and transition rate from state 1 to state 0. Here lists all sets of parameter values used in the paper. If you want to generate results under the first set of parameter values, then remove the hash sign in the begining of the first line.
#p <- c(0.1,0.1,0.03,0.03,0.01,0.001)
#p <- c(0.3,0.1,0.03,0.03,0.01,0.001)
#p <- c(0.1,0.1,0.01,0.03,0.01,0.001)
#p <- c(0.03,0.1,0.03,0.03,0.01,0.001)
#p <- c(0.1,0.1,0.09,0.03,0.01,0.001)
#p <- c(0.1,0.1,0.03,0.03,0.01,0.01)
#p <- c(0.3,0.1,0.03,0.03,0.01,0.01)
#p <- c(0.1,0.1,0.01,0.03,0.01,0.01)
#p <- c(0.03,0.1,0.03,0.03,0.01,0.01)
#p <- c(0.1,0.1,0.09,0.03,0.01,0.01)

Step 3: simulate a tree, here the tree has root state 0 and 500 extant tips. A group is defined by the "Time" scheme. All the groups are sampled in the tree. A random integer between 1 and 5 representatives are sampled for each group. Please read the definition of each input and output of the tree.gen function in tree.gen.R.
treelist <- try(tree.gen(pars=p,max.taxa=500,x0=0,sampling.c=1,sampling.f=5,define.clade="time"),silent=T)
#the while loop is to make sure that the tree does not lead to extinction
while(inherits(treelist,"try-error")) {
  treelist <- try(tree.gen(pars=p,max.taxa=500,x0=0,sampling.c=1,sampling.f=5,define.clade="time"),silent=T)
}

Step 4: fit BiSSE model to the tree using our correction. The "tlistcal" function lists all the possible locations that unsampled states and/or groups can be attached to the simulated tree. The "likcal" function calculates the overall likelihood of the tree using our correction. xfit returns the best-fit BiSEE model in xfit$par and the maximum likelihood xfit$value.
tlist <- tlistcal(tree=treelist$xtree,ssize=treelist$ssize,dt=1)
px <- c(p[1]-p[3],p[2]-p[4],p[3],p[4],p[5],p[6])
xfit <- sbplx(px,likcal,tree=treelist$xtree,sfraction=treelist$sfraction,root.state=NULL,tlist=tlist,lower=rep(0,length(p)),upper=rep(Inf,length(p)))

Step 5: fit BiSSE model to the tree using the unresolved clade correction
lik <- make.bisse(treelist$ftree,treelist$ftree$tip.state,unresolved=as.data.frame(treelist$unresolved))
ffit <- find.mle(lik,p)

Step 6: calculate the likelihood of the true model using our correction stored in variable xtrue and record the computing time in variable xpct.
pct <- proc.time()
xtrue <- likcal(p=px,tree=treelist$xtree,sfraction=treelist$sfraction,root.state=NULL,tlist=tlist)
xpct <- pct-proc.time()

Step 7: calculate the likelihood of the true model using the unresolved clade correction stored in variable ftrue and record the computing time in variable fpct.
pct <- proc.time()
ftrue <- lik(p)
fpct <- pct-proc.time()

Step 8: fit constrained BiSSE model to the tree using our correction. The BiSSE model is constrained to have equal speciation and extinction rates for both states. xfix returns the best-fit constrained BiSEE model in xfix$par and the maximum likelihood xfix$value.
xfix <- sbplx(c(0.07,0.03,p[5],p[6]),likcal,tree=treelist$xtree,sfraction=treelist$sfraction,root.state=NULL,tlist=tlist,constrained=TRUE,lower=rep(0,4),upper=rep(Inf,4))

Step 9: fit constrainted BiSSE model to the tree using the unresolved clade correction
lik2 <- constrain(lik,lambda0~lambda1,mu0~mu1)
ffix <- find.mle(lik2,c(0.1,0.1,0.03,0.03,p[5],p[6]))

######## use the code to account for non-random sampling in user's own tree and sampling pattern ##########

Step 1. copy all the functions in the folder to R and load the related library
library(deSolve)
library(nloptr)

Step 2. prepare a tree of class phylo that consists of all the sampled tips

Step 3. prepare a matrix that consists of the sampling fraction of each tip of the tree. Make sure that the rownames of the matrix are the tipnames in the tree.
The matrix should includes:
1) The first column is the state of the tip, 0 for state 0 and 1 for state 1.
2) The second column is the number of sampled tips with state 0 in the group that the tip is sampled as a representative.
3) The third column is the total number of extant members with state 0 in the group.
4) The fourth column is the number of sampled tips with state 1 in the group that the top is sampled as a representative.
5) The fifth column is the total number of extant members with state 1 in the group.
For example:
If your tree consists of 10 tips, than the matrix should contain 10 rows, each row has the name corresponding to each tip,
if "sp 1" and "sp 2" are representatives for the same group and if sp 1 is in state 0 and sp 2 in state 1, so the number of sampled tips in that group for state 0 is 1 and the number of sampled tips for state 1 is 1,
if the total number of members in that group include 5 members with state 0 and 6 members with state 1,
then the row for sp 1 and sp 2 looks like:
"sp 1" 0 1 5 1 6
"sp 2" 1 1 5 1 6

Step 4: prepare a list that consists of all the information for unsampled state or group, with unsampled state listed first, then unsampled group.
For example:
If your tree has a group with 2 sampled tips "sp 1" and "sp 2", both in state 1 and joins at "nd 1". The total number of members in that group include 5 members with state 0 and 6 members with state 1,
then state 0 is missing in that group, so create a cell in the list that includes the following elments:
1) $edge: a vector of all the edges that the unsampled state can be attached. In this example, the vector is c("nd 1", "sp 1", "sp 2"). For internal branch, edge name is the node.label of its descendent node in the tree.
2) $state: a vector of the initial D and E values for the unsampled state. In this example, the vector is c(n01=1/5,n11=0,e01=1-1/5,e11=1-2/6). Here n01 is the initial D value for state 0, n11 is the initial D value for state 1, e01 is the initial E value for state 0, e11 is the initial E value for state 1.
3) $ordered: set to FALSE for unsampled state

Now if you have listed all the unsampled state, let's start listing the unsampled groups. For the current version, we assume full topology is known for the relationship among groups. Later version will generalize the code to account for incomplete topological information.
For example:
If your tree has two unsampled groups a and b and group a is clustered with some other sampled group first, which is then clustered with group b, 
if both groups are attached to edge "nd 1",
if group a has 5 extant members with state 0 and 6 members with state 1,
if group b has 0 extant members with state 0 and 3 members with state 1,
then create a cell in the list that includes the following elments:
1) $edge = "nd 1"
2) $state: a matrix where the first row are the initial D and E values for group a, and the second row are the initial D and E values for group b. $state looks like:
n01        n11        e01          e11
5/(5+6)^2  6/(5+6)^2  1-5/(5+6)^2  1-6/(5+6)^2
0          1/3        0            1-1/3
3) $ordered = TRUE, this is to state that the location of group a and group b should be ordered on edge "nd 1" because group a is clustered with other sampled groups first

Step 5: use tlistcal to list all possible locations that unsampled states and/or groups can be attached to the tree
tlist <- tlistcal(tree=your tree in step 2,ssize=your list in step 4,dt=1)
Here, dt is the time interval between two locations and default is 1 unit of branch length. The smaller the inveral is, the more accurate the likelihod value, but takes longer to compute

Step 6: fit BiSSE model to your tree using our correction for your sampling pattern
xfit <- sbplx(p,likcal,tree=your tree in step 2,sfraction=your matrix in step 3,ssize=your list in step 4,root.state=NULL,tlist=your tlistcal output,constrained=FALSE,lower=rep(0,length(p)),upper=rep(Inf,length(p)))
Here, p is the starting parameter values to search for the best-fit model, in order of diversification rate in state 0, diversification rate in state 1, extinction rate in state 0, extinction rate in state 1, transition rate from state 0 to state 1, transition rate from state 1 to state 0. Please find the definition of each input and output of function likcal in likcal.R.

sbplx is the function that applies subplex method to search for the maximum likelihood, with lower and higher defines the boundary of the parameter space to search. Users can use other method by calling likcal function if the subplex method does not converge well.
