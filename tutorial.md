# Tutorial
As described in [`README.md`](https://github.com/siriolegramanti/CUSP/blob/master/README.md), this tutorial contains guidelines and code to perform posterior inference for Gaussian factor models endowed with the CUSP prior, leveraging the source code [`cusp.R`](https://github.com/siriolegramanti/CUSP/blob/master/cusp.R). 
In particular, it provides a step-by-step guide to reproduce the analysis of the `bfi` dataset from the `R` package `psych` performed in Section 5 of the article "Bayesian cumulative shrinkage for infinite factorizations" ([Biometrika, 107(3), 745-752](https://academic.oup.com/biomet/advance-article-abstract/doi/10.1093/biomet/asaa008/5847840); also on [arXiv](http://arxiv.org/abs/1902.04349)).
We recommend to execute the code below in the order in which it is presented.

## Setup and data uploading
Before starting the analysis, set the working directory to the location where you saved the file [`cusp.R`](https://github.com/siriolegramanti/CUSP/blob/master/cusp.R).

Then clean the workspace and load the source code [`cusp.R`](https://github.com/siriolegramanti/CUSP/blob/master/cusp.R), which contains the function that will be used for posterior inference.
``` r
rm(list=ls())
source("cusp.R")
``` 

We proceed by loading the `bfi` dataset from the `R` package `psych` and removing the rows with missing values.
Then, as described in the paper, we focus on the n = 126 individuals
older than 50 years and on the first p = 25 variables, consisting in six-point-scale answers about five personality traits.
As suggested in the [`R` documentation of the `bfi` dataset](https://www.rdocumentation.org/packages/psych/versions/1.9.12.31/topics/bfi), we change the sign of variables 1, 9, 10, 11, 12, 22 and 25 to have coherent answers within each personality trait.
Then, recalling common implementations of factor models, we center the 25 items and compute the data correlation matrix S.
``` r
library(psych)
data(bfi)
bfi<-na.omit(bfi)
bfi<-bfi[which(bfi$age>50),1:25]
bfi[,1]<- -bfi[,1]
bfi[,9]<- -bfi[,9]
bfi[,10]<- -bfi[,10]
bfi[,11]<- -bfi[,11]
bfi[,12]<- -bfi[,12]
bfi[,22]<- -bfi[,22]
bfi[,25]<- -bfi[,25]
for (j in 1:25){
  bfi[,j]<-(bfi[,j]-mean(bfi[,j]))
}
y<-as.matrix(bfi)
S<-cor(y)
``` 

## Posterior computation and post-processing
Posterior computation is carried out through the function `cusp_factor_adapt()` from [`cusp.R`](https://github.com/siriolegramanti/CUSP/blob/master/cusp.R), with parameters set as reported in the article.
``` r
gibbs<-cusp_factor_adapt(y,my_seed=6784,N_sampl=15000,
                  alpha=5,a_sig=1,b_sig=0.3,a_theta=2,b_theta=2,theta_inf=0.05,
                  start_adapt=500,Hmax=dim(y)[2]+1,alpha0=-1,alpha1=-5*10^(-4))
```

Post-processing starts by removing the first 5000 burn-in samples and  thinning every five. Then the posterior for the covariance matrix <img src="https://render.githubusercontent.com/render/math?math=\Omega"> is derived from those of the loadings matrix <img src="https://render.githubusercontent.com/render/math?math=\Lambda"> and of the noise covariance matrix <img src="https://render.githubusercontent.com/render/math?math=\Sigma"> through the well-known decomposition <img src="https://render.githubusercontent.com/render/math?math=\Omega = \Lambda \Lambda^T %2B \Sigma"> and, by averaging over posterior samples, we compute the Mean Square Error matrix between the data correlation matrix S and the posterior over S derived from that on <img src="https://render.githubusercontent.com/render/math?math=\Omega">.
``` r
burn_in=5000
thin=5
thinning <- seq(burn_in+1,gibbs$N_sampl,by=thin)
p <- dim(gibbs$y)[2]
Omega_post <- array(NA,c(p,p,length(thinning)))
SE <- array(dim=dim(Omega_post))
for (i in 1:length(thinning)){
  Omega_post[,,i]<-gibbs$Lambda_post[[thinning[i]]]%*%t(gibbs$Lambda_post[[thinning[i]]])+
    diag(1/gibbs$inv_sigma_sq_post[,thinning[i]])
  S_post<-diag(1/sqrt(diag(Omega_post[,,i])))%*%Omega_post[,,i]%*%diag(1/sqrt(diag(Omega_post[,,i])))
  SE[,,i] <- (S_post-S)^2
}
MSE <- apply(SE,c(1,2),mean)
``` 

Finally we compute:

* a scalar MSE, obtained by averaging over the upper-diagonal elements (diagonal included) of the MSE matrix:
``` r
mean(MSE[upper.tri(MSE,diag=TRUE)])
``` 
* the expected number of active factors:
``` r
mean(gibbs$Hstar[thinning])
``` 

* the runtime:
``` r
gibbs$runtime[1]
``` 

Please refer to Section 5 of the article "Bayesian cumulative shrinkage for infinite factorizations" ([Biometrika, 107(3), 745-752](https://academic.oup.com/biomet/advance-article-abstract/doi/10.1093/biomet/asaa008/5847840); also on [arXiv](http://arxiv.org/abs/1902.04349)) for detailed comments on these results in comparison to the performance of the *multiplicative gamma process* proposed in 
> Bhattacharya, A. and Dunson, D. B. (2011). Sparse Bayesian infinite factor models. Biometrika 98, 291–306.
