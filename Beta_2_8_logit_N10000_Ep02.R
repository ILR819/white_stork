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

# 100 iterations
alpha_hat4=rep(0,100)
beta_hat4=rep(0,100)
N_hat4=rep(0,100)

alpha_hat5=rep(0,100)
beta_hat5=rep(0,100)
N_hat5=rep(0,100)

alpha_hat6=rep(0,100)
beta_hat6=rep(0,100)
N_hat6=rep(0,100)

alpha_hat7=rep(0,100)
beta_hat7=rep(0,100)
N_hat7=rep(0,100)

for (s in 1:100) {
  N=10000
  T=16
  x=rep(0,N)
  f=rep(0,T)
  p=rbeta(N,2,8)
  H=matrix(rep(0,N*T),nrow=T,ncol=N)
  for(i in 1:T){
    for(j in 1:N){
      q=runif(1)
      if(p[j]>=q){
        H[i,j]=1
        x[j]=x[j]+1
      }else{
        H[i,j]=0
      }
    }
  }
  
  for(i in 1:T){
    f[i]=length(x[x==i])
  }
  n=sum(f)
  n
  #compute ln(N!)
  ln_factorialN=function(N){
    sum=0
    for(i in 1:N){
      sum=sum+log(i)
    }
    return(sum)
  }
  
  g=function(p,theta){
    return(dlogitnorm(p,mu=theta[1],sigma=theta[2]))
  }
  
  F=function(x,theta){
    f=function(p){
      (p^x)*((1-p)^(T-x))*g(p,theta)
    }
    return(f)
  }
  
  pi_g=function(x,theta){
    return(choose(T,x)*integrate(F(x,theta),0,1)$value)
  }
  
  pic_g=function(x,theta){
    return(pi_g(x,theta)/(1-pi_g(0,theta)))
  }
  
  ##############################################################################
  k=3
  O=rep(0,k)
  sum=0
  for (i in 1:(k-1)) {
    O[i]=f[i]
    sum=sum+O[i]
  } 
  O[k]=n-sum
  
  fixed_n=rep(0,180)
  fixed_n_CS=rep(0,180)
  param1=rep(0,180)
  param2=rep(0,180)
  for (i in 1:180){
    fixed_n[i]=5000+50*i
    #expected values
    EXP=function(x,theta){
      return(fixed_n[i]*pi_g(x,theta))
    } 
    E_combined=function(theta){
      sum=0
      for(j in k:T){
        sum=sum+EXP(j,theta)
      }
      return(sum)
    } 
    chi_square=function(theta){
      O_o=fixed_n[i]-sum(O)
      sum=(O_o-EXP(0,theta))^2/EXP(0,theta)
      for(i in 1:(k-1)){
        sum=sum+(O[i]-EXP(i,theta))^2/EXP(i,theta)
      }
      sum=sum+(O[k]-E_combined(theta))^2/E_combined(theta)
      return(sum)
    }
    CS=optim(c(-2,1),chi_square,method="L-BFGS-B",lower=c(-5,0.5),upper=c(5,2))
    fixed_n_CS[i]=chi_square(CS$par)
    param1[i]=CS$par[1]
    param2[i]=CS$par[2]
  }
  startN=5000+50*which(fixed_n_CS==min(fixed_n_CS))
  
  fixed_n=rep(0,100)
  fixed_n_CS=rep(0,100)
  param1=rep(0,100)
  param2=rep(0,100)
  for (i in 1:100){
    fixed_n[i]=startN-50+i
    #expected values
    EXP=function(x,theta){
      return(fixed_n[i]*pi_g(x,theta))
    } 
    E_combined=function(theta){
      sum=0
      for(j in k:T){
        sum=sum+EXP(j,theta)
      }
      return(sum)
    } 
    chi_square=function(theta){
      O_o=fixed_n[i]-sum(O)
      sum=(O_o-EXP(0,theta))^2/EXP(0,theta)
      for(i in 1:(k-1)){
        sum=sum+(O[i]-EXP(i,theta))^2/EXP(i,theta)
      }
      sum=sum+(O[k]-E_combined(theta))^2/E_combined(theta)
      return(sum)
    }
    CS=optim(c(-2,1),chi_square,method="L-BFGS-B",lower=c(-5,0.5),upper=c(5,2))
    fixed_n_CS[i]=chi_square(CS$par)
    param1[i]=CS$par[1]
    param2[i]=CS$par[2]
  }
  print("C=4")
  print(c(param1[which(fixed_n_CS==min(fixed_n_CS))], param2[which(fixed_n_CS==min(fixed_n_CS))], min(fixed_n_CS), startN-50+which(fixed_n_CS==min(fixed_n_CS))))
  alpha_hat4[s]=param1[which(fixed_n_CS==min(fixed_n_CS))]
  beta_hat4[s]=param2[which(fixed_n_CS==min(fixed_n_CS))]
  N_hat4[s]=startN-50+which(fixed_n_CS==min(fixed_n_CS))
  
  ##############################################################################
  k=4
  O=rep(0,k)
  sum=0
  for (i in 1:(k-1)) {
    O[i]=f[i]
    sum=sum+O[i]
  }
  O[k]=n-sum
  
  fixed_n=rep(0,180)
  fixed_n_CS=rep(0,180)
  param1=rep(0,180)
  param2=rep(0,180)
  for (i in 1:180){
    fixed_n[i]=5000+50*i
    #expected values
    EXP=function(x,theta){
      return(fixed_n[i]*pi_g(x,theta))
    } 
    E_combined=function(theta){
      sum=0
      for(j in k:T){
        sum=sum+EXP(j,theta)
      }
      return(sum)
    } 
    chi_square=function(theta){
      O_o=fixed_n[i]-sum(O)
      sum=(O_o-EXP(0,theta))^2/EXP(0,theta)
      for(i in 1:(k-1)){
        sum=sum+(O[i]-EXP(i,theta))^2/EXP(i,theta)
      }
      sum=sum+(O[k]-E_combined(theta))^2/E_combined(theta)
      return(sum)
    }
    CS=optim(c(-2,1),chi_square,method="L-BFGS-B",lower=c(-5,0.5),upper=c(5,2))
    fixed_n_CS[i]=chi_square(CS$par)
    param1[i]=CS$par[1]
    param2[i]=CS$par[2]
  }
  startN=5000+50*which(fixed_n_CS==min(fixed_n_CS))
  
  fixed_n=rep(0,100)
  fixed_n_CS=rep(0,100)
  param1=rep(0,100)
  param2=rep(0,100)
  for (i in 1:100){
    fixed_n[i]=startN-50+i
    #expected values
    EXP=function(x,theta){
      return(fixed_n[i]*pi_g(x,theta))
    } 
    E_combined=function(theta){
      sum=0
      for(j in k:T){
        sum=sum+EXP(j,theta)
      }
      return(sum)
    } 
    chi_square=function(theta){
      O_o=fixed_n[i]-sum(O)
      sum=(O_o-EXP(0,theta))^2/EXP(0,theta)
      for(i in 1:(k-1)){
        sum=sum+(O[i]-EXP(i,theta))^2/EXP(i,theta)
      }
      sum=sum+(O[k]-E_combined(theta))^2/E_combined(theta)
      return(sum)
    }
    CS=optim(c(-2,1),chi_square,method="L-BFGS-B",lower=c(-5,0.5),upper=c(5,2))
    fixed_n_CS[i]=chi_square(CS$par)
    param1[i]=CS$par[1]
    param2[i]=CS$par[2]
  }
  print("C=5")
  print(c(param1[which(fixed_n_CS==min(fixed_n_CS))],param2[which(fixed_n_CS==min(fixed_n_CS))],min(fixed_n_CS),startN-50+which(fixed_n_CS==min(fixed_n_CS))))
  
  alpha_hat5[s]=param1[which(fixed_n_CS==min(fixed_n_CS))]
  beta_hat5[s]=param2[which(fixed_n_CS==min(fixed_n_CS))]
  N_hat5[s]=startN-50+which(fixed_n_CS==min(fixed_n_CS))
  
  ##############################################################################
  k=5
  O=rep(0,k)
  sum=0
  for (i in 1:(k-1)) {
    O[i]=f[i]
    sum=sum+O[i]
  }
  O[k]=n-sum
  
  fixed_n=rep(0,180)
  fixed_n_CS=rep(0,180)
  param1=rep(0,180)
  param2=rep(0,180)
  for (i in 1:180){
    fixed_n[i]=5000+50*i
    #expected values
    EXP=function(x,theta){
      return(fixed_n[i]*pi_g(x,theta))
    } 
    E_combined=function(theta){
      sum=0
      for(j in k:T){
        sum=sum+EXP(j,theta)
      }
      return(sum)
    } 
    chi_square=function(theta){
      O_o=fixed_n[i]-sum(O)
      sum=(O_o-EXP(0,theta))^2/EXP(0,theta)
      for(i in 1:(k-1)){
        sum=sum+(O[i]-EXP(i,theta))^2/EXP(i,theta)
      }
      sum=sum+(O[k]-E_combined(theta))^2/E_combined(theta)
      return(sum)
    }
    CS=optim(c(-2,1),chi_square,method="L-BFGS-B",lower=c(-5,0.5),upper=c(5,2))
    fixed_n_CS[i]=chi_square(CS$par)
    param1[i]=CS$par[1]
    param2[i]=CS$par[2]
  }
  startN=5000+50*which(fixed_n_CS==min(fixed_n_CS))
  
  fixed_n=rep(0,100)
  fixed_n_CS=rep(0,100)
  param1=rep(0,100)
  param2=rep(0,100)
  for (i in 1:100){
    fixed_n[i]=startN-50+i
    #expected values
    EXP=function(x,theta){
      return(fixed_n[i]*pi_g(x,theta))
    } 
    E_combined=function(theta){
      sum=0
      for(j in k:T){
        sum=sum+EXP(j,theta)
      }
      return(sum)
    } 
    chi_square=function(theta){
      O_o=fixed_n[i]-sum(O)
      sum=(O_o-EXP(0,theta))^2/EXP(0,theta)
      for(i in 1:(k-1)){
        sum=sum+(O[i]-EXP(i,theta))^2/EXP(i,theta)
      }
      sum=sum+(O[k]-E_combined(theta))^2/E_combined(theta)
      return(sum)
    }
    CS=optim(c(-2,1),chi_square,method="L-BFGS-B",lower=c(-5,0.5),upper=c(5,2))
    fixed_n_CS[i]=chi_square(CS$par)
    param1[i]=CS$par[1]
    param2[i]=CS$par[2]
  }
  print("C=6")
  print(c(param1[which(fixed_n_CS==min(fixed_n_CS))],param2[which(fixed_n_CS==min(fixed_n_CS))],min(fixed_n_CS),startN-50+which(fixed_n_CS==min(fixed_n_CS))))
  
  alpha_hat6[s]=param1[which(fixed_n_CS==min(fixed_n_CS))]
  beta_hat6[s]=param2[which(fixed_n_CS==min(fixed_n_CS))]
  N_hat6[s]=startN-50+which(fixed_n_CS==min(fixed_n_CS))
  
  ##############################################################################
  k=6
  O=rep(0,k)
  sum=0
  for (i in 1:(k-1)) {
    O[i]=f[i]
    sum=sum+O[i]
  }
  O[k]=n-sum
  
  fixed_n=rep(0,180)
  fixed_n_CS=rep(0,180)
  param1=rep(0,180)
  param2=rep(0,180)
  for (i in 1:180){
    fixed_n[i]=5000+50*i
    #expected values
    EXP=function(x,theta){
      return(fixed_n[i]*pi_g(x,theta))
    } 
    E_combined=function(theta){
      sum=0
      for(j in k:T){
        sum=sum+EXP(j,theta)
      }
      return(sum)
    } 
    chi_square=function(theta){
      O_o=fixed_n[i]-sum(O)
      sum=(O_o-EXP(0,theta))^2/EXP(0,theta)
      for(i in 1:(k-1)){
        sum=sum+(O[i]-EXP(i,theta))^2/EXP(i,theta)
      }
      sum=sum+(O[k]-E_combined(theta))^2/E_combined(theta)
      return(sum)
    }
    CS=optim(c(-2,1),chi_square,method="L-BFGS-B",lower=c(-5,0.5),upper=c(5,2))
    fixed_n_CS[i]=chi_square(CS$par)
    param1[i]=CS$par[1]
    param2[i]=CS$par[2]
  }
  startN=5000+50*which(fixed_n_CS==min(fixed_n_CS))
  
  fixed_n=rep(0,100)
  fixed_n_CS=rep(0,100)
  param1=rep(0,100)
  param2=rep(0,100)
  for (i in 1:100){
    fixed_n[i]=startN-50+i
    #expected values
    EXP=function(x,theta){
      return(fixed_n[i]*pi_g(x,theta))
    } 
    E_combined=function(theta){
      sum=0
      for(j in k:T){
        sum=sum+EXP(j,theta)
      }
      return(sum)
    } 
    chi_square=function(theta){
      O_o=fixed_n[i]-sum(O)
      sum=(O_o-EXP(0,theta))^2/EXP(0,theta)
      for(i in 1:(k-1)){
        sum=sum+(O[i]-EXP(i,theta))^2/EXP(i,theta)
      }
      sum=sum+(O[k]-E_combined(theta))^2/E_combined(theta)
      return(sum)
    }
    CS=optim(c(-2,1),chi_square,method="L-BFGS-B",lower=c(-5,0.5),upper=c(5,2))
    fixed_n_CS[i]=chi_square(CS$par)
    param1[i]=CS$par[1]
    param2[i]=CS$par[2]
  }
  print("C=7")
  print(c(param1[which(fixed_n_CS==min(fixed_n_CS))],param2[which(fixed_n_CS==min(fixed_n_CS))],min(fixed_n_CS),startN-50+which(fixed_n_CS==min(fixed_n_CS))))
  
  alpha_hat7[s]=param1[which(fixed_n_CS==min(fixed_n_CS))]
  beta_hat7[s]=param2[which(fixed_n_CS==min(fixed_n_CS))]
  N_hat7[s]=startN-50+which(fixed_n_CS==min(fixed_n_CS))
  ##############################################################################
}

N_hat4
N_hat5
N_hat6
N_hat7
