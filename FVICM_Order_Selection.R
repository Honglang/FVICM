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
mvec = rep(m,N)
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
eps = c(t(rmvnorm(N, mean=rep(0,m), sigma=sigma2)))


A1 = sqrt(3)/2-1.645/sqrt(12)
A2 = sqrt(3)/2+1.645/sqrt(12)

y = cos(pi*x%*%beta0_true) + sin(pi*(x%*%beta1_true-A1)/(A2-A1))*G + eps 
y_0 = cos(pi*x%*%beta0_true) + eps

####################Some matrix needed later#################
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
estimation0 = function(x, y, K, d, lambda){
  
  D = diag(c(rep(0,d+1),rep(1,K))) # K: num of knots, d: order
  Jn = d+K+1
  p = ncol(x)
  
  ##step 0: initial values
  beta.old = c(sqrt(1/3),sqrt(1/3),sqrt(1/3))
  u = x%*%beta.old
  
  B = spline(u, K, d)
  
  gamma.old = ginv(t(B)%*%B)%*%t(B)%*%y  ##simple LS, initial value 
  
  maxstep = 0   ##loop for gamma and beta
  while(maxstep <= 50){
    maxstep = maxstep+1
    
    ##################
    #step 1: b-spline approximation
    u = x%*%beta.old        
    B = spline(u, K, d)
    
    ###########################
    ## step 2: gamma est
    ### Newton-Raphson for estimating gamma
    
    run = 0  
    while(run <= 50){
      
      run = run+1
      
      mu = B %*% gamma.old    
      mudot_gamma = B
      ngamma = K+d+1
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
      
      u = x%*%NRbeta.old
      
      B = spline(u, K, d)
      mu = B %*% gamma.old
      
      Bd = spline_dev(u, K, d)
      
      # first derivative of mu 
      mudot_beta = matrix(rep(Bd%*%gamma.old,p),nrow=Nobs, ncol=p)*x    
      nbeta = p
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
      
      NRbeta.new_norm = sign(NRbeta.new[1])*NRbeta.new/sqrt(sum(NRbeta.new^2))
      
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
  
  
  
  list(GCV, beta.old, gamma.old, QN)
  
  
} #estimation0 function




###################################################################
### golden search method

f0 = function(z){estimation0(x, y_0, K, d, z)[[1]]}

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



########################################################################

Qstat0 = function(x, y, K, d, beta, gamma){      
  # based on the estimation, calculate Q
  u = x%*%beta
  
  B = spline(u, K, d)
  mu = B %*% gamma
  
  Bd = spline_dev(u, K, d)
  
  #J1 = -beta[-1]/beta[1]
  #J = rbind(J1,diag(p-1))
  
  
  # first derivative of mu wrt beta   
  #mudot_beta = (matrix(rep(Bd%*%gamma,p),nrow=Nobs, ncol=p)*x)%*%J
  
  # first derivative of mu 
  mudot_gamma = B
  mudot = mudot_gamma     
  #mudot = cbind(mudot_beta, mudot_gamma)
  
  npar = length(gamma)
  #npar = length(beta)+length(gamma)-1  #number of parameters
  arsumg = matrix(rep(0,2*npar),nrow=2*npar)
  arsumc = matrix(rep(0,2*npar*2*npar),nrow=2*npar)
  #gi = matrix(rep(0,2*npar),nrow=2*npar)
  arsumgfirstdev = matrix(rep(0,2*npar*npar),nrow=2*npar)
  #firstdev = matrix(rep(0,2*npar*npar),nrow=2*npar)
  
  for(i in 1:N){
    seq = c(mvec0[i]:mvec1[i])
    ni = mvec[i]
    yi = y[seq]
    xi = x[seq,]
    mui = mu[seq]
    mudoti = mudot[seq,]
    fmui = mui  #link function
    fmui_dev = diag(ni)  # first dev of link function
    vmui = diag(ni) #variance of mui
    
    wi = t(mudoti) %*% fmui_dev %*% vmui %*% m0[[i]] %*% vmui
    zi = t(mudoti) %*% fmui_dev %*% vmui %*% m1[[i]] %*% vmui
    
    gi0 = (1/N)* wi %*% (yi-mui)
    gi1 = (1/N)* zi %*% (yi-mui)
    
    gi = c(gi0, gi1)
    
    arsumc = arsumc + gi %*% t(gi) #CN_bar =  1/N^2*sum(gi*t(gi))
    arsumg = arsumg + gi  #gN_bar
    
    di0 = - (1/N)* wi %*% fmui_dev %*% mudoti #first dev of gi0
    di1 = - (1/N)* zi %*% fmui_dev %*% mudoti #first dev of gi1
    
    firstdev = rbind(di0,di1)
    arsumgfirstdev = arsumgfirstdev + firstdev
    
  }
  
  arcinv=ginv(arsumc)
  
  Q = t(arsumg) %*% arcinv %*% arsumg
  Q
  
} #Qstat0 function



##############################################################
####### choose order and number of knots #####################
##############################################################

knot = c(1,2,3,4,5)

#comb = as.matrix(expand.grid(knot,order))

p.value1 = rep(0,length(knot))
BIC1 = rep(0,length(knot))
BIC1_1 = rep(0,length(knot))
for(k in 1:length(knot)){
  K = knot[k]
  d = 3
  
  lambda_est = golden(f0, 0.01, 2, K, d)
  
  fit = estimation0(x,y_0,K,d,lambda_est)
  
  beta = fit[[2]]
  gamma = fit[[3]]
  gamma_all = c(gamma[1:(d+1)],0,gamma[(d+2):(K+d+1)],rep(0,5-K))
  Q = Qstat0(x,y_0,5,4,beta,gamma_all)
  
  lamb = lambda_est
  D = diag(c(rep(0,d+1),rep(1,K)))
  
  npar = length(gamma) 
  p.value1[k] = pchisq(Q,npar+N*lamb*t(gamma)%*%D%*%gamma,lower.tail=F)
  
  BIC1[k] = Q+N*lamb*t(gamma)%*%D%*%gamma + npar*log(N) 
  BIC1_1[k] = Q+npar*log(N) 
}


p.value2 = rep(0,length(knot))
BIC2 = rep(0,length(knot))
BIC2_1 = rep(0,length(knot))
for(k in 1:length(knot)){
  K = knot[k]
  d = 4
  
  lambda_est = golden(f0, 0.01, 2, K, d)
  
  fit = estimation0(x,y_0,K,d,lambda_est)
  
  beta = fit[[2]]
  gamma = fit[[3]]
  gamma_all = c(gamma[1:(d+1)],gamma[(d+2):(K+d+1)],rep(0,5-K))
  Q = Qstat0(x,y_0,5,4,beta,gamma_all)
  
  lamb = lambda_est
  D = diag(c(rep(0,d+1),rep(1,K)))
  
  npar = length(gamma)+length(beta)-1 
  p.value2[k] = pchisq(Q,npar,lower.tail=F)
  
  BIC2[k] = Q+N*lamb*t(gamma)%*%D%*%gamma + npar*log(N) 
  BIC2_1[k] = Q+npar*log(N) 
}

p.value = c(p.value1,p.value2)
BIC = c(BIC1,BIC2)






