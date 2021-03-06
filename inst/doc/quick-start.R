## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----eval=FALSE---------------------------------------------------------------
#  install.packages("Rdimtools")

## ----eval=FALSE---------------------------------------------------------------
#  ## install.packages("devtools")
#  devtools::install_github("kisungyou/Rdimtools")

## ----setup, message=FALSE, warning=FALSE--------------------------------------
library(Rdimtools)

## ----echo=FALSE, include=FALSE------------------------------------------------
vernow = utils::packageVersion("Rdimtools")
ndo    = (sum(unlist(lapply(ls("package:Rdimtools"), startsWith, "do."))))
nest   = (sum(unlist(lapply(ls("package:Rdimtools"), startsWith, "est."))))

## ----message=FALSE, warning=FALSE, fig.align='center', fig.width=7, fig.height=3----
# load the iris data
X   = as.matrix(iris[,1:4])
lab = as.factor(iris[,5])

# we will compare 5 methods (out of 17 methods from version 1.0.0)
vecd = rep(0,5)
vecd[1] = est.Ustat(X)$estdim       # convergence rate of U-statistic on manifold
vecd[2] = est.correlation(X)$estdim # correlation dimension
vecd[3] = est.made(X)$estdim        # manifold-adaptive dimension estimation
vecd[4] = est.mle1(X)$estdim        # MLE with Poisson process
vecd[5] = est.twonn(X)$estdim       # minimal neighborhood information

# let's visualize
plot(1:5, vecd, type="b", ylim=c(1.5,3.5), 
     main="estimating dimension of iris data",
     xaxt="n",xlab="",ylab="estimated dimension")
xtick = seq(1,5,by=1)
axis(side=1, at=xtick, labels = FALSE)
text(x=xtick,  par("usr")[3], 
     labels = c("Ustat","correlation","made","mle1","twonn"), pos=1, xpd = TRUE)

## ----message=FALSE, warning=FALSE, fig.align='center', fig.width=7------------
# run 3 algorithms mentioned above
mypca = do.pca(X, ndim=2)
mylap = do.lscore(X, ndim=2)
mydfm = do.dm(X, ndim=2, bandwidth=10)

# extract embeddings from each method
Y1 = mypca$Y
Y2 = mylap$Y
Y3 = mydfm$Y
  
# visualize
par(mfrow=c(1,3))
plot(Y1, pch=19, col=lab, xlab="axis 1", ylab="axis 2", main="PCA")
plot(Y2, pch=19, col=lab, xlab="axis 1", ylab="axis 2", main="Laplacian Score")
plot(Y3, pch=19, col=lab, xlab="axis 1", ylab="axis 2", main="Diffusion Maps")

