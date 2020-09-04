##########################################################################
# Adaptive Gibbs sampler for Gaussian factor models under the CUSP prior
##########################################################################

library(MASS) # for mvrnorm (multivariate Normal)
library(LaplacesDemon) # for dmvt (multivariate Student's t)

cusp_factor_adapt <- function(y,my_seed,N_sampl,alpha,a_sig,b_sig,a_theta,b_theta,
                              theta_inf,start_adapt,Hmax,alpha0,alpha1){
  set.seed(my_seed)
  # statistical units
  n<-dim(y)[1]
  # number of variables
  p<-dim(y)[2]
  # uniform rvs for adaptation
  u<-runif(N_sampl)
  ########################################################################
  # VARIABLES TO UPDATE
  ########################################################################
  # number of factors
  H<-rep(NA,N_sampl)
  # number of active factors
  Hstar<-rep(NA,N_sampl)
  # factor loadings matrix
  Lambda_post<-vector("list",N_sampl)
  # factor scores matrix
  eta_post<-vector("list",N_sampl)
  # column-specific loadings precisions
  theta_inv_post<-vector("list",N_sampl)
  # stick-breaking weights
  w_post<-vector("list",N_sampl)
  # augmented data
  z_post<-vector("list",N_sampl)
  # error variances
  inv_sigma_sq_post<-matrix(NA,p,N_sampl)
  ########################################################################
  #INITIALIZATION
  ########################################################################
  H[1]<-p+1
  Hstar[1]<-p
  Lambda_post[[1]]<-matrix(rnorm(p*H[1]),p,H[1])
  eta_post[[1]]<-matrix(rnorm(n*H[1]),n,H[1])
  theta_inv_post[[1]]<-rep(1,H[1])
  w_post[[1]]<-rep(1/H[1],H[1])
  z_post[[1]]<-rep(H[1],H[1])
  inv_sigma_sq_post[,1]<-1
  ########################################################################
  # GIBBS SAMPLING #######################################################
  ########################################################################
  t0 <- proc.time()
  for (t in 2:N_sampl){
    Lambda_post[[t]]<-matrix(NA,p,H[t-1])
    eta_post[[t]]<-matrix(NA,n,H[t-1])
    theta_inv_post[[t]]<-rep(NA,H[t-1])
    w_post[[t]]<-rep(NA,H[t-1])
    z_post[[t]]<-rep(NA,H[t-1])
    # 1) sample the j-th row of Lambda
    for (j in 1:p){
      V_j<-chol2inv(chol(diag(theta_inv_post[[t-1]],nrow=H[t-1])+inv_sigma_sq_post[j,t-1]*(t(eta_post[[t-1]])%*%eta_post[[t-1]])))
      mu_j<-V_j%*%t(eta_post[[t-1]])%*%(y[,j])*inv_sigma_sq_post[j,t-1]
      Lambda_post[[t]][j,]<-mvrnorm(1,mu_j,V_j)
    }
    # 2) sample sigma^{-2}
    diff_y<-apply((y-(eta_post[[t-1]]%*%t(Lambda_post[[t]])))^2,2,sum)
    for (j in 1:p){
      inv_sigma_sq_post[j,t]<-rgamma(n=1,shape=a_sig+0.5*n,rate=b_sig+0.5*diff_y[j])
    }		
    # 3) sample eta
    V_eta<-chol2inv(chol(diag(H[t-1])+t(Lambda_post[[t]])%*%diag(inv_sigma_sq_post[,t])%*%Lambda_post[[t]]))
    for (i in 1:n){
      mu_eta_i<-V_eta%*%t(Lambda_post[[t]])%*%diag(inv_sigma_sq_post[,t])%*%y[i,]
      eta_post[[t]][i,]<-mvrnorm(1,mu_eta_i,V_eta)
    }	
    # 4) sample z
    lhd_spike<-rep(0,H[t-1])
    lhd_slab<-rep(0,H[t-1])
    for (h in 1:H[t-1]){
      lhd_spike[h]<-exp(sum(log(dnorm(Lambda_post[[t]][,h], mean = 0, sd = theta_inf^(1/2), log = FALSE))))
      lhd_slab[h]<-dmvt(x=Lambda_post[[t]][,h], mu=rep(0,p), S=(b_theta/a_theta)*diag(p), df=2*a_theta)
      prob_h<-w_post[[t-1]]*c(rep(lhd_spike[h],h),rep(lhd_slab[h],H[t-1]-h))
      if (sum(prob_h)==0){
        prob_h<-c(rep(0,H[t-1]-1),1)
      }
      else{
        prob_h<-prob_h/sum(prob_h)
      }
      z_post[[t]][h]<-c(1:H[t-1])%*%rmultinom(n=1, size=1, prob=prob_h)
    }
    # 5) sample v and update w
    v<-rep(NA,H[t-1])
    for (h in 1:(H[t-1]-1)){
      v[h]<-rbeta(1, shape1 = 1+sum(z_post[[t]]==h), shape2 = alpha+sum(z_post[[t]]>h))
    }
    v[H[t-1]]<-1
    w_post[[t]][1]<-v[1]
    for (h in 2:H[t-1]){
      w_post[[t]][h]<-v[h]*prod(1-v[1:(h-1)])  
    }
    # 6) sample theta^{-1}
    for (h in 1:H[t-1]){
      if (z_post[[t]][h]<=h){
        theta_inv_post[[t]][h]<-theta_inf^(-1)
      }
      else{
        theta_inv_post[[t]][h]<-rgamma(n=1,shape=a_theta+0.5*p,rate=b_theta+0.5*t(Lambda_post[[t]][,h])%*%Lambda_post[[t]][,h])
      }
    }
    # 7) update H[t]
    active <- which(z_post[[t]]>c(1:H[t-1]))
    Hstar[t] <- length(active)
    H[t]<-H[t-1]
    # by default, keep the same truncation, unless...
    if (t>=start_adapt & u[t]<=exp(alpha0+alpha1*t)){
      if (Hstar[t]<H[t-1]-1){
        # set truncation to Hstar[t] and subset all variables, keeping only active columns
        H[t]<-Hstar[t]+1
        eta_post[[t]]<-cbind(eta_post[[t]][,active],rnorm(n))
        theta_inv_post[[t]]<-c(theta_inv_post[[t]][active],theta_inf^(-1))
        w_post[[t]]<-c(w_post[[t]][active],1-sum(w_post[[t]][active]))
        Lambda_post[[t]]<-cbind(Lambda_post[[t]][,active],rnorm(p,mean=0,sd=sqrt(theta_inf)))
      } else if (H[t-1]<Hmax) {
        # increase truncation by 1 and extend all variables, sampling from the prior/model
        H[t]<-H[t-1]+1
        eta_post[[t]]<-cbind(eta_post[[t]],rnorm(n))
        v[H[t-1]]<-rbeta(1,shape1=1,shape2=alpha)
        v<-c(v,1)
        w_post[[t]]<-rep(NA,H[t])
        w_post[[t]][1]<-v[1]
        for (h in 2:H[t]){
          w_post[[t]][h]<-v[h]*prod(1-v[1:(h-1)])  
        }
        theta_inv_post[[t]]<-c(theta_inv_post[[t]],theta_inf^(-1))
        Lambda_post[[t]]<-cbind(Lambda_post[[t]],rnorm(p,mean=0,sd=sqrt(theta_inf)))
      }
    }
    print(t)
  }
  runtime <- proc.time()-t0
  output<-list("y"=y,"my_seed"=my_seed,"N_sampl"=N_sampl,
               "alpha"=alpha,"a_sig"=a_sig,"b_sig"=b_sig,"a_theta"=a_theta,"b_theta"=b_theta,
               "theta_inf"=theta_inf,"start_adapt"=start_adapt,"Hmax"=Hmax,
               "alpha0"=alpha0,"alpha1"=alpha1,
               "H"=H,"Hstar"=Hstar,"Lambda_post"=Lambda_post,"eta_post"=eta_post,
               "theta_inv_post"=theta_inv_post,"w_post"=w_post,"z_post"=z_post,
               "inv_sigma_sq_post"=inv_sigma_sq_post,"runtime"=runtime)
  return(output)
}
