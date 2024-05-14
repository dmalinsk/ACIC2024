## install.packages("BiocManager")
## BiocManager::install(c('graph','RBGL','Rgraphviz'))
## install.packages("pcalg")

library(pcalg)

## true DAG
## X1 -> X2 <- L -> X3 <- X4

set.seed(123)
n <- 10000
X1 <- rnorm(n, 0 , .5)
X4 <- rnorm(n, 0, 1.4)
L <- rnorm(n, 0, 0.7)
X2 <- rnorm(n, -1 + 0.8*X1 + 0.8*L, 1)
X3 <- rnorm(n, -1 + 2.2*X4 + 1.8*L, 1)

data <- cbind(X1,X2,X3,X4) ## L is not observed in data

indepTest <- gaussCItest ## specify the independence test
suffStat <- list(C = cor(data), n = n) ## using the correlation matrix

## RUNNING FCI
fci.est <- fci(suffStat, indepTest, alpha = 0.05, p = 4, verbose=TRUE) ## estimate a PAG
plot(fci.est)

## RUNNING PC
pc.est <- pc(suffStat, indepTest, alpha = 0.05, p = 4, verbose=TRUE) ## estimate at CPDAG
plot(pc.est) ## something seems to be wrong/incompatible with latest update of 'pcalg' and Rstudio plotting
as(pc.est@graph,"matrix") ## look at the result as a matrix directly

## RUNNING GES
score <- new("GaussL0penObsScore", data) ## define a score function, this is the Gaussian BIC score

ges.est <- ges(score, verbose=TRUE) ## estimate at CPDAG with GES
plot(ges.est$essgraph)

## WARNING: the default plotting method for CPDAGs in package 'pcalg' will display undirected (--) edges
## as if they are bidirected (<->). This is just a somewhat annoying feature of the built-in plot function, they
## are not really bidirected edges, they are undirected. CPDAGs may *only* contain directed and undirected edges!

## IMPORTANT:
## pcalg uses two different internal matrix representations for graphs. 
## i) a matrix where each element is 0 or 1. 
#### Here if element (ij)==1 and (ji)==0, then i-->j. 
#### If (ij)==(ji)==1 then i--j which is (annoyingly!) displayed as i<->j
## ii) a matrix where each element is 0, 1, 2, or 3.
#### Here if element (ij)==1/2/3 there is a circle/arrowhead/tail at j from i, respectively. 0 means no edge.
#### This allows one to encode PAGs with all types of edges: o-->, o--o, <-->, -->, o--, --- etc.

## compare:
as(pc.est@graph,"matrix") ## matrix of 0s and 1s
## with
fci.est@amat ## matrix of 0/1/2/3

## It is aways good to look at the matrix representations of the graphs!

###############################################################################
## MORE EXAMPLES
###############################################################################

## True DAG: Z1 -> X <-- Z2 ; X --> Y <-- L ; possibly L -?-> X

set.seed(123)
n <- 10000
Z1 <- rnorm(n, 0, .5)
Z2 <- rnorm(n, 0 , 1.5)
L <- rnorm(n, 0, 1)
X <- rnorm(n, 0.8*Z1 + 0.8*Z2 + 1.5*L, 1) ## try to replace coefficient for L with a zero, see result change
Y <- rnorm(n, 1.6*X + 0.6*L, 1)

data <- cbind(Z1,Z2,X,Y) ## L is not observed in data

indepTest <- gaussCItest ## specify the independence test
suffStat <- list(C = cor(data), n = n) ## using the correlation matrix

fci.est <- fci(suffStat, indepTest, alpha = 0.05, labels = c("Z1","Z2","X","Y"), verbose=TRUE) ## estimate a PAG
plot(fci.est)

## now try to replace the coefficient on L with a zero (so it is no longer a common cause of X and Y) and run again!

## also try a random DAG

p <- 10
n <- 10000
myDAG <- pcalg::randomDAG(n = p, prob = 0.4)
plot(myDAG)
data <- rmvDAG(n, myDAG, errDist = "normal")
data <- data[,-c(4,5)] ## "hide" two of the variables

indepTest <- gaussCItest ## specify the independence test
suffStat <- list(C = cor(data), n = n) ## using the correlation matrix

fci.est <- fci(suffStat, indepTest, alpha = 0.01, labels = c("1","2","3","6","7","8","9","10"), verbose=TRUE) ## estimate a PAG
plot(fci.est)

pc.est <- pc(suffStat, indepTest, alpha = 0.01, labels = c("1","2","3","6","7","8","9","10"), verbose=TRUE) ## estimate at CPDAG
plot(pc.est)

score <- new("GaussL0penObsScore", data) ## define a score function

ges.est <- ges(score, verbose=TRUE) ## estimate a CPDAG with GES
plot(ges.est$essgraph)


## IMPORTANT: every algorithm has various settings/parameters you can choose. Look at the documentation!
## For example, with PC and FCI, it may be prudent to restrict the max size of conditioning set: set m.max. 
## Higher-order independence tests may be unreliable, and this can also save time in computationally-intensive problems.
## Must also choose appropriate conditional independence tests! (There are many out there.)
## May also consider "conservative" or "majority rule" orientations of colliders. (See documentation.)
## May also consider different ways to handle conflicts among tests, e.g., "solve.confl=TRUE". (See documentation.)