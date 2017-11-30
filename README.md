# Correction-for-non-random-sampling
Correction for non-random sampling under the BiSSE framework

This contains the codes used in paper:
"The influence of non-random species-sampling on macroevolutionary and macroecological inference from phylogenies"
by Xia Hua and Robert Lanfear

Please cite the paper if you use the code.

The code implements our proposed correction for non-random sampling under the BiSSE framework in R.
Each function is documented at the beginning of the file.

This repo gives information on:

1. replicating the results in the paper
2. using the method no your own data

## Replicating the results in the paper

#### Step 1

Download this github repository to your computer. From here on in, we'll assume that this is now located at `~/Desktop/Correction-for-non-random-sampling/`

#### Step 2

Source all the functions to R and load the related libraries

```r
  library(diversitree)
  library(deSolve)
  library(nloptr)

  file.sources = list.files(path = "~/Desktop/Correction-for-non-random-sampling/", pattern="*.R", full.names = T)
  sapply(file.sources,source,.GlobalEnv)
```

#### Step 3

specify a set of parameter values that include six parameters in the following order: 

* speciation rate of state 0
* speciation rate of state 1
* extinction rate of state 0
* extinction rate of state 1
* transition rate from state 0 to state 1
* transition rate from state 1 to state 0

Here are lists of all sets of parameter values used in the paper. To use a parameter set, just choose the relevant line.

```r
  p <- c(0.1,0.1,0.03,0.03,0.01,0.001)
  p <- c(0.3,0.1,0.03,0.03,0.01,0.001)
  p <- c(0.1,0.1,0.01,0.03,0.01,0.001)
  p <- c(0.03,0.1,0.03,0.03,0.01,0.001)
  p <- c(0.1,0.1,0.09,0.03,0.01,0.001)
  p <- c(0.1,0.1,0.03,0.03,0.01,0.01)
  p <- c(0.3,0.1,0.03,0.03,0.01,0.01)
  p <- c(0.1,0.1,0.01,0.03,0.01,0.01)
  p <- c(0.03,0.1,0.03,0.03,0.01,0.01)
  p <- c(0.1,0.1,0.09,0.03,0.01,0.01)
```

#### Step 4

Simulate a tree. Here the tree has root state 0 and 500 extant tips, and higher taxa are defined by the "Time" scheme described in the paper. All higher groups in the tree are sampled with between 1 and 5 representatives sampled for each group. Please read the definition of each input and output of the tree.gen function in tree.gen.R.

```r

  treelist <- try(tree.gen(pars=p,max.taxa=500,x0=0,sampling.c=1,sampling.f=5,define.clade="time"),silent=T)

  # the while loop is to make sure that the tree does not lead to extinction
  while(inherits(treelist,"try-error")) {
    treelist <- try(tree.gen(pars=p,max.taxa=500,x0=0,sampling.c=1,sampling.f=5,define.clade="time"),silent=T)
  }  
```

#### Step 5

Fit the BiSSE model to the tree using our correction. The `tlistcal()` function lists all the possible locations that unsampled states and/or groups can be attached to the simulated tree. The `likcal()` function calculates the overall likelihood of the tree using our correction. `xfit` contains the best-fit BiSEE model in `xfit$par` and the maximum likelihood value in `xfit$value`.

```r
  tlist <- tlistcal(tree=treelist$xtree,ssize=treelist$ssize,dt=1)

  px <- c(p[1]-p[3],p[2]-p[4],p[3],p[4],p[5],p[6])

  xfit <- sbplx(px,likcal,tree=treelist$xtree,sfraction=treelist$sfraction,root.state=NULL,tlist=tlist,lower=rep(0,length(p)),upper=rep(Inf,length(p)))
```

#### Step 6
Fit BiSSE model to the tree using the unresolved clade correction

```r
  lik <- make.bisse(treelist$ftree,treelist$ftree$tip.state,unresolved=as.data.frame(treelist$unresolved))

  ffit <- find.mle(lik,p)
```

#### Step 7
Calculate the likelihood of the true model using our correction. The likelihood is stored in `xtrue`, and the computing time in `xpct`.

```r
  pct <- proc.time()

  xtrue <- likcal(p=px,tree=treelist$xtree,sfraction=treelist$sfraction,root.state=NULL,tlist=tlist)

  xpct <- pct-proc.time()
```

#### Step 8
Calculate the likelihood of the true model using the unresolved clade correction. The likelihood is stored in `ftrue` and  the computing time in `fpct`.

```r
  pct <- proc.time()

  ftrue <- lik(p)

  fpct <- pct-proc.time()
```

#### Step 9

Fit the constrained BiSSE model to the tree using our correction. The BiSSE model is constrained to have equal speciation and extinction rates for both states. `xfix` contains the best-fit constrained BiSEE model in `xfix$par` and the maximum likelihood in `xfix$value`.

```r
  xfix <- sbplx(c(0.07,0.03,p[5],p[6]),likcal,tree=treelist$xtree,sfraction=treelist$sfraction,root.state=NULL,tlist=tlist,constrained=TRUE,lower=rep(0,4),upper=rep(Inf,4))
```

#### Step 10

Fit the constrainted BiSSE model to the tree using the unresolved clade correction

```r
  lik2 <- constrain(lik,lambda0 ~ lambda1,mu0 ~ mu1)
  ffix <- find.mle(lik2,c(0.1,0.1,0.03,0.03,p[5],p[6]))
```

Now we have all of the information required to compare these models for one iteration of one set of parameter values. 


### Use the method on your own data

#### Step 1

Download this github repository to your computer. From here on in, we'll assume that this is now located at `~/Desktop/Correction-for-non-random-sampling/`

#### Step 2

Source all the functions to R and load the related libraries

```r
  library(deSolve)
  library(nloptr)
  library(ape)

  file.sources = list.files(path = "~/Desktop/Correction-for-non-random-sampling/", pattern="*.R", full.names = T)
  sapply(file.sources,source,.GlobalEnv)
```

#### Step 3

Load your tree (with all sampled tips) into R:

```r
  t = read.tree(path_to_your_tree)
```

#### Step 4

Build a matrix that consists of the sampling fraction for each tip of the tree.

The rownames of the matrix should be the tips of the tree, and the columns shoul be

1. The first column is the state of the tip, 0 for state 0 and 1 for state 1.
2. The second column is the number of sampled tips with state 0 in the group that the tip is sampled as a representative.
3. The third column is the total number of extant members with state 0 in the group.
4. The fourth column is the number of sampled tips with state 1 in the group that the top is sampled as a representative.
5. The fifth column is the total number of extant members with state 1 in the group.

For example:
If your tree consists of 10 tips, than the matrix should contain 10 rows, where each row is named to match a tip. If "sp 1" and "sp 2" are representatives for the same group and if sp 1 is in state 0 and sp 2 in state 1, then the number of sampled tips in that group for state 0 is 1 and the number of sampled tips for state 1 is 1. Let's also assume that the total number of members in that group include 5 members with state 0 and 6 members with state 1. In this case the first two rows of the matrix will look like this:

```r
my.matrix = matrix(c("sp 1", 0, 1, 5, 1, 6,
                     "sp 2", 1, 1, 5, 1, 6), nrow = 2, byrow = T)
```

so entering `my.matrix` gives

```r
     [,1]   [,2] [,3] [,4] [,5] [,6]
[1,] "sp 1" "0"  "1"  "5"  "1"  "6" 
[2,] "sp 2" "1"  "1"  "5"  "1"  "6" 
```

#### Step 4

Build a list that consists of all the information for unsampled states or groups, with unsampled states listed first.

For example:
If your tree has a group with 2 sampled tips "sp 1" and "sp 2", both in state 1, and that join at "nd 1" in the tree. In this example (as above) the total number of members in that group include 5 members with state 0 and 6 members with state 1, then state 0 is missing in that group, so create an entry in the list that includes the following elments:

1. `$edge`: a vector of all the edges that the unsampled state can be attached. In this example, the vector is `c("nd 1", "sp 1", "sp 2")`. For internal branches, edge name is the node.label of its descendent node in the tree.
2. `$state`: a vector of the initial D and E values for the unsampled state. In this example, the vector is `c(n01=1/5,n11=0,e01=1-1/5,e11=1-2/6)`. Here n01 is the initial D value for state 0, n11 is the initial D value for state 1, e01 is the initial E value for state 0, e11 is the initial E value for state 1.
3. `$ordered`: set to FALSE for unsampled state.

Now if you have listed all the unsampled states, let's start listing the unsampled groups. For the current version, we assume full topology is known for the relationship among groups. Later versions will generalize the code to account for incomplete topological information.

For example:
Let's have a tree with two unsampled groups a and b, where and group a is clustered with some other sampled group first, which is then clustered with group b. Assume that 
* both groups are attached to edge "nd 1",
* group a has 5 extant members with state 0 and 6 members with state 1,
* group b has 0 extant members with state 0 and 3 members with state 1,
then create a cell in the list that includes the following elments:
1. `$edge = "nd 1"`
2. `$state`: a matrix where the first row is the initial D and E values for group a, and the second row is the initial D and E values for group b. `$state` looks like:
```
n01        n11        e01          e11
5/(5+6)^2  6/(5+6)^2  1-5/(5+6)^2  1-6/(5+6)^2
0          1/3        0            1-1/3
```
3. `$ordered = TRUE`: this is to state that the location of group a and group b should be ordered on edge "nd 1" because group a is clustered with other sampled groups first

#### Step 5 
Use tlistcal to list all possible locations that unsampled states and/or groups can be attached to the tree
```r
  tlist <- tlistcal(tree=your tree in step 2,ssize=your list in step 4,dt=1)
```
Here, dt is the time interval between two locations and default is 1 unit of branch length. The smaller the inveral is, the more accurate the likelihod value, but takes longer to compute

#### Step 6
Fit the BiSSE model to your tree using our correction for your sampling pattern
```r
  xfit <- sbplx(p,likcal,tree=your tree in step 2,sfraction=your matrix in step 3,ssize=your list in step 4,root.state=NULL,tlist=your tlistcal output,constrained=FALSE,lower=rep(0,length(p)),upper=rep(Inf,length(p)))
```

Here, p is the starting parameter values to search for the best-fit model, in order of diversification rate in state 0, diversification rate in state 1, extinction rate in state 0, extinction rate in state 1, transition rate from state 0 to state 1, transition rate from state 1 to state 0. Please find the definition of each input and output of function likcal in likcal.R.

sbplx is the function that applies subplex method to search for the maximum likelihood, with lower and higher defines the boundary of the parameter space to search. Users can use other method by calling likcal function if the subplex method does not converge well.
