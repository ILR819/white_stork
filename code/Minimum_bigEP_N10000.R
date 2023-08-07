
N=10000
T=16

N_hat4=rep(0,100)
N_hat5=rep(0,100)
N_hat6=rep(0,100)
N_hat7=rep(0,100)
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
g=function(p,theta) dbeta(p,theta[1],theta[2])
F=function(x,theta){
  f=function(p){
    p^x*(1-p)^(T-x)*g(p,theta)
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
  a=rnorm(N,-1.5,1)
  p=exp(a)/(1+exp(a))#probability distribution/probabilites are fixed
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
  m=length(NK)
  

  
  #maximum likelihood 
  flg=function(theta) { 
    sum=Fa(theta[3])-Fa(theta[3]-n)+(theta[3]-n)*log(pieg(0,theta))
    for (i in 1:T) {
      sum=sum+f[i]*log(pieg(i,theta))
    }
    return (-sum)
  }   #log likelihood
  #Sensitivity for initial values
  fixed_n_CS=rep(0,m)
  N_hat_temp=rep(0,m)
  for (i in 1:m) {
    Mfg=optim(c(2,8,NK[i]),flg,method="L-BFGS-B",lower=c(0.2,0.2,N0+1),upper=c(50,50,20000))
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
  
  #chi_square 
  chi_square_optim <- function(cell) {
    k=cell
    O=rep(0,k)
    sum=0
    for (i in 1:(k-1)) {
      O[i]=f[i]
      sum=sum+O[i]
    }
    O[k]=n-sum
    
    EL=function(theta) {
      sum=0
      for (j in k:T) {
        sum=sum+E(j,theta)
      }
      return (sum) 
    }
    
    chi_square=function(theta) { 
      OO=theta[3]-sum(O)
      sum=(OO-E(0,theta))^2/E(0,theta)
      for (i in 1:(k-1)) {
        sum=sum+(O[i]-E(i,theta))^2/E(i,theta)
      }
      sum=sum+(O[k]-EL(theta))^2/EL(theta)
      return (sum)
    } 
    
    fixed_n_CS=rep(0,m)
    N_hat_temp=rep(0,m) 
    for (i in 1:m){
      CS=optim(c(2,8,NK[i]),chi_square,method="L-BFGS-B",lower=c(0.2,0.2,N0+1),upper=c(50,50,20000))
      fixed_n_CS[i]=CS$value
      N_hat_temp[i]=CS$par[3]
    }
    N_hat=N_hat_temp[which(fixed_n_CS==min(fixed_n_CS))]
    return(N_hat)
  }
  N_hat4[s]=chi_square_optim(3)
  N_hat5[s]=chi_square_optim(4)
  N_hat6[s]=chi_square_optim(5)
  N_hat7[s]=chi_square_optim(6)
  
}

results=data.frame(N_hat_fl,N_hat_cl,N_hat4, N_hat5,N_hat6,N_hat7) 
write.csv(results,"results10000_bigEP.csv")


MFL=mean(N_hat_fl)
MCL=mean(N_hat_cl)
M4=mean(N_hat4)
M5=mean(N_hat5)
M6=mean(N_hat6)
M7=mean(N_hat7)

BFL=abs(MFL-N)
BCL=abs(MCL-N)
B4=abs(M4-N)
B5=abs(M5-N)
B6=abs(M6-N)
B7=abs(M7-N)


mseFL=sum((N_hat_fl-N)^2)/100
mseCL=sum((N_hat_cl-N)^2)/100
mse4=sum((N_hat4-N)^2)/100
mse5=sum((N_hat5-N)^2)/100
mse6=sum((N_hat6-N)^2)/100
mse7=sum((N_hat7-N)^2)/100

smseFL=sqrt(mseFL)
smseCL=sqrt(mseCL)
smse4=sqrt(mse4)
smse5=sqrt(mse5)
smse6=sqrt(mse6)
smse7=sqrt(mse7)



mean_N <- c(MFL,MCL,M4,M5,M6,M7)
bias_N <- c(BFL,BCL,B4,B5,B6,B7)
smse_N <- c(smseFL,smseCL,smse4,smse5,smse6,smse7)

results1=data.frame(mean_N,bias_N,smse_N) 
write.csv(results1,"results10000_bigEP_final.csv")