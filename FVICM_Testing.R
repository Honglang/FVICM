library(nlme)
library(RLRsim)
library(MASS)
library(mvtnorm)
library(Matrix)

#########################################################################
####################### Power Calculation ###############################
#########################################################################

all_tau = seq(0, 0.2,0.025)
power=rep(0,length(all_tau))

for (kkk in 1:length(all_tau)){
  
  tau = all_tau[kkk]
  
  
  ########## p-spline functions ################################
  
  spline_fix = function(xx, d){
    
    S = matrix(0, nrow=length(xx), ncol=d+1)	 
    
    for(i in 1:(d+1)){
      S[,i]=xx^(i-1)
    }
    
    S
  }
  
  
  spline_ran=function(xx,d,K){
    knots=quantile(xx, probs = (1:K)/(K+1))
    Z = matrix(0, nrow=length(xx), ncol=K)
    for(j in 1:K){
      
      Z[,j]=((xx-knots[j])^d)*(xx>=knots[j]) 
    } 
    Z
    
  }
  
  #####################################################################
  
  #set.seed(56789)
  
  powerL=0
  
  step1=0
  while(step1<20){
    
    #genarate data
    set.seed(100*(step1+1))    
    N = 500
    m = 10
    mvec = rep(m,N)
    Nobs = sum(mvec)
    mvec0 = cumsum(c(1, mvec[-N])) #start location for each id
    mvec1 = cumsum(mvec) #end location for each id
    id = rep(c(1:N),mvec)
    
    x1 = runif(Nobs)
    x2 = runif(Nobs)
    x3 = runif(Nobs)
    x = cbind(x1,x2,x3)
    p = ncol(x)
    
    pA = 0.5
    G0 = sample(c(0,1,2), size = N, replace = T, prob=c((1-pA)^2, 2*pA*(1-pA), pA^2)) 
    G = rep(G0,mvec) 
    
    
    ### exchangeable
    rho = 0.8
    R = matrix(0,m,m)
    for (k in 1:m){
      for (l in 1:m){
        if(k==l){
          R[k,l]=1} else{
            R[k,l]=rho        }         }    ## exchangeable
    }
    sigma2 = 0.1*R
    eps = c(t(rmvnorm(N, mean=rep(0,m), sigma=sigma2)))
    
    
    ####true values 
    A1 = sqrt(3)/2-1.645/sqrt(12)
    A2 = sqrt(3)/2+1.645/sqrt(12)
    beta0_true = c(sqrt(5/13),sqrt(4/13),sqrt(4/13))  
    beta1_true = c(sqrt(1/3),sqrt(1/3),sqrt(1/3))
    beta_true=c(beta0_true, beta1_true)
    
    m1_H0 = 0.3 + x%*%beta1_true
    m1_full = sin(pi*(x%*%beta1_true-A1)/(A2-A1))
    m1_H1 = m1_H0 + tau*(m1_full-m1_H0)
    
    Y = cos(pi*x%*%beta0_true) + m1_H1*G + eps
    
    beta0_est=c(0.6201687, 0.5544204, 0.5548358)  
    beta1_est=c(0.5765780, 0.5773831, 0.5779934) # estimators from the estimation step
    
    u0=x%*%beta0_est
    u1=x%*%beta1_est
    d=4
    K0=K1=1
    X1=spline_fix(u0,d)
    X20=spline_fix(u1,d)
    G_ext1 = matrix(rep(G,d+1), nrow=Nobs, ncol=d+1)
    X2=X20*G_ext1
    X=cbind(X1,X2)
    Z1=spline_ran(u0,d,K0)
    Z20=spline_ran(u1,d,K1)
    G_ext2 = matrix(rep(G,K1), nrow=Nobs, ncol=K1)
    Z2=Z20*G_ext2
    
    X_0 = X[,1:(d+1+2)] #design matrix under H0
    q=ncol(X_0)
    p=ncol(X)
    group1=group=rep(1,Nobs)
    fitt=try(lme(Y~-1+X, random=list(group=pdIdent(~-1+Z1), group=pdIdent(~-1+Z2), id=~1),method='ML'))
    if(class(fitt)!="try-error"){
      
      var1 = as.numeric(VarCorr(fitt)[2,1])
      var2 = as.numeric(VarCorr(fitt)[3+3,1])
      var_a = as.numeric(VarCorr(fitt)[7,1])
      var_e = fitt$sigma^2
      Z3 = kronecker(diag(N),rep(1,m))
      Z3_2 = kronecker(diag(N),matrix(1,m,m))
      Z1_2 = Z1%*%t(Z1)
      Var = var1*Z1_2+var2*Z2%*%t(Z2)+var_a*Z3_2+var_e*diag(Nobs)
      coff = coef(fitt,level=1)[1:10]
      
      temp = chol2inv(chol(Var))
      #temp = solve(Var)
      pred = (var1*Z1_2+var_a*Z3_2)%*%(temp%*%(Y-X%*%t(coff)))
      
      Y_s = Y-pred   ##pseudo-response
      
    }
    # critical value for LRT 
    sim.LRT = LRTSim(X=cbind(X1,X2), Z=Z2, q=d+1-2, sqrt.Sigma=diag(ncol(Z2)), nsim=5000)
    thresL=quantile(sim.LRT, prob=0.95)
    
    step2=0
    while(step2<50){
      
      ### exchangeable
      rho = 0.8
      R = matrix(0,m,m)
      for (k in 1:m){
        for (l in 1:m){
          if(k==l){
            R[k,l]=1} else{
              R[k,l]=rho        }         }    ## exchangeable
      }
      sigma2 = 0.1*R
      eps = c(t(rmvnorm(N, mean=rep(0,m), sigma=sigma2)))   
      
      Y = cos(pi*x%*%beta0_true) + m1_H1*G + eps  
      Y_s = Y-pred
      group1=group=rep(1,Nobs)
      X_0 = X[,1:(d+1+2)] #design matrix under H0
      q=ncol(X_0)
      p=ncol(X)
      fit=try(lme(Y_s~-1+X, random=list(group=pdIdent(~-1+Z2)),method='ML'))
      fit0=try(lm(Y_s~-1+X_0))
      if(all(c(class(fit), class(fit0))!="try-error")){
        
        # compute the LRT test statistics
        LRTobs=anova(fit,fit0)$L.Ratio[2]
        # power of the LRT
        powerL=powerL+(LRTobs>thresL)*1
      }
      #}
      step2=step2+1
      
    }
    step1=step1+1
  }
  
  power[kkk]=powerL/2000
}
power/2


