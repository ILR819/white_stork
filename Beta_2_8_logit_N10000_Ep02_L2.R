#!/bin/sh test

###############################################################
# Script for simulations:
# - Generate capture-recapture data using model A=Beta(2,8)
# - Estimate population size N using model B=logit-normal
# - True population size N=10000
# - Sampling occasions T=16
# - Number of cells C=4/5/6/7
# - Number of iterations s=100
###############################################################

#install.packages("pso")
#install.packages("logitnorm")
#library("pso")
library("logitnorm")

N=10000
T=16
N_hat_cl=rep(0,100)
N_hat_fl=rep(0,100)

#A function to compute ln N!
Fa=function(N){
  sum=0
  for (i in 1:N) {
    sum=sum+log(i)
  }
  return(sum)
}

#continuous function beta
#beta distribution(continuous) theta=c(a,b)
g=function(p,theta){
  return(dbeta(p,theta[1],theta[2]))
}
#for continuous function
F=function(x,theta){
  f=function(p){
    (p^x)*((1-p)^(T-x))*g(p,theta)
  }
  return(f)
}
pieg=function(x,theta) { 
  choose(T, x)*integrate(F(x,theta),0,1)$value 
}
piecg=function(x,theta) pieg(x,theta)/(1-pieg(0,theta))

E=function(x,theta) theta[3]*pieg(x,theta)



for (s in 1:100) {
  x=rep(0,N)
  f=rep(0,T)
  p=rbeta(N,2,8)
  CH=matrix(rep(0,N*T),nrow=T,ncol=N)
  
  #MAIN to get frequencies
  for (i in 1:T) {
    for (j in 1:N) {
      u=runif(1)
      if (p[j]>u) {
        x[j]=x[j]+1 
        CH[i,j]=1
      }
    }
  }
  
  for (i in 1:T) {f[i]=length(x[x==i])}
  
  n=sum(f)
  m1=sum(CH[1,])
  m2=0
  u2=0
  for (i in 1:N) {
    if (CH[2,i]==1) {
      if (CH[1,i]==1) m2=m2+1
      else u2=u2+1
    }
  }
  LP=m1*u2/m2+m1
  
  if (LP==Inf | LP>=n*1.5 | n>=LP) N0=n else N0=LP
  
  NK=rep(0,50)
  for (i in 1:50){
    NK[i]=N0+50*i
  }
  

  #maximum likelihood 
  flg=function(theta) { 
    sum=Fa(theta[3])-Fa(theta[3]-n)+(theta[3]-n)*log(pieg(0,theta))
    for (i in 1:T) {
      sum=sum+f[i]*log(pieg(i,theta))
    }
    return (-sum)
  }   #log likelihood
  #Sensitivity for initial values
  fixed_n_CS=rep(0,50)
  N_hat_temp=rep(0,50)
  for (i in 1:50) {
    Mfg=optim(c(2,8,NK[i]),flg,method="L-BFGS-B",lower=c(0.2,0.2,N0+50),upper=c(50,50,20000))
    fixed_n_CS[i]=Mfg$value
    N_hat_temp[i]=Mfg$par[3]
  }
  N_hat_fl[s]=N_hat_temp[which(fixed_n_CS==min(fixed_n_CS))]
  
  #conditional likelihood*-1
  clg=function(theta) { 
    sum=0
    for (i in 1:T) {
      sum=sum+f[i]*log(piecg(i,theta))
    }
    return (-sum)
  }   #log likelihood
  
  Mcg=optim(c(2,8),clg,method="L-BFGS-B",lower=c(0.2,0.2),upper=c(50,50))
  if (Mcg$convergence==0) {
    N_hat_cl[s]=n/(1-pieg(0,Mcg$par))
  } else {
    N_hat_cl[s]=-1
  }
  ####################################################################################
}

print("Conditional Likelihood")
N_hat_cl
print("Full Likelihood")
N_hat_fl