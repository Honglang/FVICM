library("MASS")
library("mvtnorm")
library("Matrix")

#########################################################################
####################### Generate Demo Data  #############################
#########################################################################
####true values
beta0_true = c(sqrt(5/13),sqrt(4/13),sqrt(4/13))  
beta1_true = c(sqrt(1/3),sqrt(1/3),sqrt(1/3))
beta_true=c(beta0_true, beta1_true)


#genarate data
N = 200
m = 10 
mvec = rep(0,N)
indx = vector("list", N)
for (i in 1:N){
  visit = c(1,rbinom((m-1),1,1)) # 30% missing except 1st time
  mvec[i] = sum(visit)  
  indx[[i]] = which(visit==1)
}
Nobs = sum(mvec)
mvec0 = cumsum(c(1, mvec[-N])) #start loc for each id
mvec1 = cumsum(mvec) #end loc for each id


x1 = runif(Nobs)
x2 = runif(Nobs)
x3 = runif(Nobs)
x = cbind(x1,x2,x3)
p = ncol(x)

pr = 0.5
G0 = sample(c(0,1,2), size = N, replace = T, prob=c(pr^2, 2*pr*(1-pr), (1-pr)^2))
G = rep(G0,mvec) 

rho = 0.5
R = matrix(0,m,m)
for (k in 1:m){
  for (l in 1:m){
    if(k==l){
      R[k,l]=1} else{
        R[k,l]=rho        }         }    ## exchangeable
}
sigma2 = 0.1*R

eps = vector("list", N)
for(i in 1:N){
  eps.old = c(rmvnorm(1, mean=rep(0,m), sigma=sigma2))
  eps[[i]] = eps.old[indx[[i]]]
}

eps = unlist(eps)


A1 = sqrt(3)/2-1.645/sqrt(12)
A2 = sqrt(3)/2+1.645/sqrt(12)

y = cos(pi*x%*%beta0_true) + sin(pi*(x%*%beta1_true-A1)/(A2-A1))*G + eps 

#####################some matrix needed later################

### corr matrix: exchangable #####
m0 = list()
m1 = list()
for (i in 1:N){
  ni = mvec[i]
  ## corr matrix: exchangable  
  m0[[i]] = diag(ni)           
  m1[[i]] = matrix(rep(1,ni*ni),ni) - m0[[i]]
}
##############################################################

#########################################################################
####################### Functions  ######################################
#########################################################################


########## p-spline functions ################################

spline = function(xx, K, d){
  Jn = d+K+1
  
  S = matrix(0, nrow=length(xx), ncol=Jn)	 
  if (K == 0){ 
    for(i in 1:(d+1)){
      S[,i]=xx^(i-1)
    }
  }  else {
    for(i in 1:(d+1)){
      S[,i]=xx^(i-1)
    }
    
    knots=quantile(xx, probs = (1:K)/(K+1))
    
    for(j in 1:K){
      
      S[,d+1+j]=((xx-knots[j])^d)*(xx>=knots[j])
    } 
  }
  S
}

### first derivative of splines

spline_dev = function(xx, K, d){
  
  Jn = d+K+1
  
  S_dev = matrix(0, nrow=length(xx), ncol=Jn)
  if (K==0){
    for(i in 2:(d+1)){
      S_dev[,i]=(i-1)*(xx^(i-2))
    }
  }  else{ for(i in 2:(d+1)){
    S_dev[,i]=(i-1)*(xx^(i-2))
  }
    
    knots=quantile(xx, probs = (1:K)/(K+1))
    
    for(j in 1:K){
      
      S_dev[,d+1+j]=d*((xx-knots[j])^(d-1))*(xx>=knots[j]) 
    }
  } 
  S_dev
}




#################### estimation function ###############

#### estimation function ##########
estimation = function(x, y, G, K, d, lambda){
  
  D = diag(rep(c(rep(0,d+1),rep(1,K)),2)) # K: num of knots, d: order
  Jn = d+K+1
  p = ncol(x)
  
  ##step 0: initial values
  beta.old = c(sqrt(1/3),sqrt(1/3),sqrt(1/3), sqrt(5/13),sqrt(4/13),sqrt(4/13))
  beta.old0 = beta.old[1:p]
  beta.old1 = beta.old[(p+1):(2*p)]
  u0 = x%*%beta.old0
  u1 = x%*%beta.old1 
  u = cbind(u0, u1)  
  
  B0 = spline(u0, K, d)
  B1 = spline(u1, K, d)
  
  G_ext = matrix(rep(G,K+d+1), nrow=Nobs, ncol=K+d+1)
  B = cbind(B0, B1*G_ext)
  
  gamma.old = ginv(t(B)%*%B)%*%t(B)%*%y  ##simple LS, initial value 
  
  maxstep = 0   ##loop for gamma and beta
  while(maxstep <= 50){
    maxstep = maxstep+1
    
    ##################
    #step 1: b-spline approximation
    beta.old0 = beta.old[1:p]
    beta.old1 = beta.old[(p+1):(2*p)]
    u0 = x%*%beta.old0
    u1 = x%*%beta.old1 
    u = cbind(u0, u1)  
    
    B0 = spline(u0, K, d)
    B1 = spline(u1, K, d)
    
    B = cbind(B0, B1*G_ext)
    
    ###########################
    ## step 2: gamma est
    ### Newton-Raphson for estimating gamma
    
    run = 0  
    while(run <= 50){
      
      run = run+1
      
      mu = B %*% gamma.old    
      mudot_gamma = B
      ngamma = 2*(K+d+1)
      arsumg = matrix(rep(0,2*ngamma),nrow=2*ngamma)
      arsumc = matrix(rep(0,2*ngamma*2*ngamma),nrow=2*ngamma)
      #gi = matrix(rep(0,2*ngamma),nrow=2*ngamma)
      arsumgfirstdev1 = matrix(rep(0,2*ngamma*ngamma),nrow=2*ngamma)
      #firstdev = matrix(rep(0,2*ngamma*ngamma),nrow=2*ngamma)
      
      for(i in 1:N){
        seq = c(mvec0[i]:mvec1[i])
        ni = mvec[i]
        yi = y[seq]
        xi = x[seq,]
        mui = mu[seq]
        mudoti = mudot_gamma[seq,]
        fmui = mui  #link function
        fmui_dev = diag(ni)  # first dev of link function
        vmui = diag(ni) #variance of mui
        
        wi = t(mudoti) %*% fmui_dev %*% vmui %*% m0[[i]] %*% vmui
        zi = t(mudoti) %*% fmui_dev %*% vmui %*% m1[[i]] %*% vmui
        
        gi0 = (1/N)* wi %*% (yi-mui)
        gi1 = (1/N)* zi %*% (yi-mui)
        
        gi = c(gi0, gi1)
        
        arsumc = arsumc + gi %*% t(gi) 
        arsumg = arsumg + gi 
        
        di0 = - (1/N)* wi %*% fmui_dev %*% mudoti #first dev of gi0
        di1 = - (1/N)* zi %*% fmui_dev %*% mudoti #first dev of gi1
        
        firstdev = rbind(di0,di1)
        arsumgfirstdev1 = arsumgfirstdev1 + firstdev
        
      }
      
      arcinv1 = ginv(arsumc)
      
      Q1 = t(arsumg) %*% arcinv1 %*% arsumg 
      
      arqif1dev1 = (2*t(arsumgfirstdev1) %*% arcinv1 %*% arsumg)/N + 2*lambda*D%*%gamma.old
      arqif2dev1 = (2*t(arsumgfirstdev1) %*% arcinv1 %*% arsumgfirstdev1)/N + 2*lambda*D 	
      
      invarqif2dev1 = ginv(arqif2dev1)
      
      gamma.new = gamma.old - invarqif2dev1 %*% arqif1dev1
      gammadiff = max(abs(gamma.new - gamma.old))
      if(gammadiff<1e-6){break}
      
      gamma.old = gamma.new
      #print(gamma.old)
    } #loop for gamma
    
    
    #########################               
    # step 3: est for beta
    ### Newton-Raphson for estimating beta
    
    NRbeta.old=beta.old
    run2 = 0  
    while(run2 <= 50){
      
      run2 = run2+1
      
      NRbeta.old0 = NRbeta.old[1:p]
      NRbeta.old1 = NRbeta.old[(p+1):(2*p)]
      u0 = x%*%NRbeta.old0
      u1 = x%*%NRbeta.old1 
      u = cbind(u0, u1)  
      
      B0 = spline(u0, K, d)
      B1 = spline(u1, K, d)
      B = cbind(B0, B1*G_ext)
      mu = B %*% gamma.old
      
      B0d = spline_dev(u0, K, d)
      B1d = spline_dev(u1, K, d)
      
      gamma0.old = gamma.old[1:(K+d+1)]
      gamma1.old = gamma.old[(K+d+1+1):(2*(K+d+1))]
      
      # first derivative of mu 
      mudot_beta = cbind(matrix(rep(B0d%*%gamma0.old,p),nrow=Nobs, ncol=p)*x, matrix(rep(B1d%*%gamma1.old*G,p),nrow=Nobs, ncol=p)*x)        
      
      nbeta = 2*p
      arsumg = matrix(rep(0,2*nbeta),nrow=2*nbeta)
      arsumc = matrix(rep(0,2*nbeta*2*nbeta),nrow=2*nbeta)
      #gi = matrix(rep(0,2*nbeta),nrow=2*nbeta)
      arsumgfirstdev = matrix(rep(0,2*nbeta*nbeta),nrow=2*nbeta)
      #firstdev = matrix(rep(0,2*nbeta*nbeta),nrow=2*nbeta)
      
      for(i in 1:N){
        seq = c(mvec0[i]:mvec1[i])
        ni = mvec[i]
        yi = y[seq]
        xi = x[seq,]
        mui = mu[seq]
        mudoti = mudot_beta[seq,]
        fmui = mui  #link function
        fmui_dev = diag(ni)  # first dev of link function
        vmui = diag(ni) #variance of mui
        
        wi = t(mudoti) %*% fmui_dev %*% vmui %*% m0[[i]] %*% vmui
        zi = t(mudoti) %*% fmui_dev %*% vmui %*% m1[[i]] %*% vmui
        
        gi0 = (1/N)* wi %*% (yi-mui)
        gi1 = (1/N)* zi %*% (yi-mui)
        
        gi = c(gi0, gi1)
        
        arsumc = arsumc + gi %*% t(gi) 
        arsumg = arsumg + gi 
        
        di0 = - (1/N)* wi %*% fmui_dev %*% mudoti #first dev of gi0
        di1 = - (1/N)* zi %*% fmui_dev %*% mudoti #first dev of gi1
        
        firstdev = rbind(di0,di1)
        arsumgfirstdev = arsumgfirstdev + firstdev
        
      }
      
      arcinv=ginv(arsumc)
      
      Q2 = t(arsumg) %*% arcinv %*% arsumg
      
      arqif1dev2 = 2*t(arsumgfirstdev) %*% arcinv %*% arsumg
      arqif2dev2 = 2*t(arsumgfirstdev) %*% arcinv %*% arsumgfirstdev   	
      
      invarqif2dev2 = ginv(arqif2dev2)
      
      NRbeta.new = NRbeta.old - invarqif2dev2 %*% arqif1dev2
      
      NRbeta.new1 = sign(NRbeta.new[1])*NRbeta.new[1:p]/sqrt(sum(NRbeta.new[1:p]^2))
      NRbeta.new2 = sign(NRbeta.new[1+p])*NRbeta.new[(p+1):(2*p)]/sqrt(sum(NRbeta.new[(p+1):(2*p)]^2))
      NRbeta.new_norm = c(NRbeta.new1, NRbeta.new2)
      
      betadiff = max(abs(NRbeta.new_norm - NRbeta.old))
      if(betadiff<1e-6){break}
      
      NRbeta.old = NRbeta.new_norm
      #print(NRbeta.old)
    } #loop for beta
    
    beta.new = NRbeta.old
    
    dif = max(abs(beta.new-beta.old))
    #dif = sqrt(sum((beta.new-beta.old)^2))
    if(dif<1e-4){break}
    #print(dif)
    #print(gamma.new)
    #print(beta.new)
    beta.old = beta.new
    
  }#outer iteration for beta and gamma
  
  
  # based on the estimation, calculate Q, covariance
  QNdotdot = 2*t(arsumgfirstdev1) %*% arcinv1 %*% arsumgfirstdev1
  QN = Q1
  df = sum(diag(ginv(QNdotdot + 2*N*lambda*D)%*%QNdotdot))
  GCV = (QN/N)/(1-df/N)^2
  
  lam_gam = lambda*D%*%gamma.old
  arsums = matrix(rep(0,ngamma*ngamma),nrow=ngamma)
  
  for (i in 1:N){
    seq = c(mvec0[i]:mvec1[i])
    ni = mvec[i]
    yi = y[seq]
    xi = x[seq,]
    mui = mu[seq]
    mudoti = mudot_gamma[seq,]
    fmui = mui  #link function
    fmui_dev = diag(ni)  # first dev of link function
    vmui = diag(ni) #variance of mui
    
    wi = t(mudoti) %*% fmui_dev %*% vmui %*% m0[[i]] %*% vmui
    zi = t(mudoti) %*% fmui_dev %*% vmui %*% m1[[i]] %*% vmui
    
    gi0 = wi %*% (yi-mui)
    gi1 = zi %*% (yi-mui)
    
    gi = c(gi0, gi1)
    si = t(arsumgfirstdev1) %*% (arcinv1/N) %*% gi + lam_gam
    
    arsums = arsums + si%*%t(si)
    
  }
  
  cov_part = ginv(arqif2dev1*N/2)
  cov_gamma = cov_part%*%arsums%*%cov_part
  cov_gamma1 = cov_gamma[(Jn+1):(2*Jn),(Jn+1):(2*Jn)]
  cov_gamma0 = cov_gamma[1:Jn,1:Jn]
  
  
  w = seq(0.25,1.25,0.01)
  sp_matrix = spline(w,K,d)
  m1SE = sqrt(diag(sp_matrix%*%cov_gamma1%*%t(sp_matrix)))
  m0SE = sqrt(diag(sp_matrix%*%cov_gamma0%*%t(sp_matrix)))
  
  asy.cov = ginv(arqif2dev2/2)
  SE_beta = sqrt(diag(asy.cov))
  
  CIlower = beta.old - 1.96*SE_beta   #confidence interval
  CIupper = beta.old + 1.96*SE_beta
  
  cp = (beta_true<=CIupper) & (beta_true>=CIlower)
  
  
  list(GCV, beta.old, gamma.old, SE_beta, m0SE, m1SE, cp, QN)
  
  
} #estimation function



###################################################################
### golden search method

#f0 = function(x){estimation_0(X, Y_0, K, d, x)[[1]]}
f1 = function(z){estimation(x, y, G, K, d, z)[[1]]}

golden = function(f, lower, upper, K, d) { 
  
  golden.ratio = 2/(sqrt(5) + 1)
  x1 = upper - golden.ratio*(upper - lower)
  x2 = lower + golden.ratio*(upper - lower)
  
  f1 = f(x1)
  f2 = f(x2)
  
  iteration = 0
  
  while (abs(upper - lower) > 0.1)
  {
    iteration = iteration + 1
    
    if (f2 > f1)
    {
      ### Set the new upper bound
      upper = x2
      ### Set the new upper test point
      ### Use the special result of the golden ratio
      x2 = x1
      f2 = f1
      
      ### Set the new lower test point
      x1 = upper - golden.ratio*(upper - lower)
      f1 = f(x1)
    } 
    else 
    {
      # the minimum is to the right of x1
      # let x1 be the new lower bound
      # let x2 be the new lower test point
      
      ### Set the new lower bound
      lower = x1
      
      ### Set the new lower test point
      x1 = x2
      
      f1 = f2
      
      ### Set the new upper test point
      x2 = lower + golden.ratio*(upper - lower)
      f2 = f(x2)
    }
  }
  
  ### Use the mid-point of the final interval as the estimate of the optimzer
  estimated.minimizer = (lower + upper)/2
  estimated.minimizer
}

#########################################################################################

#########################################################################
####################### Estimation ######################################
#########################################################################
########## after determination of knots and order
K=1
d=4

lambda_est = golden(f1, 0.01, 2, K, d)

result2 = estimation(x, y, G, K, d, lambda_est)
beta.old = result2[[2]]
gamma.old = result2[[3]]
SE_beta = result2[[4]]
m0SE = result2[[5]]
m1SE = result2[[6]]
cp = result2[[7]]














