install.packages("MCMCpack")
require(MCMCpack)
install.packages("LaplacesDemon")
require(LaplacesDemon)
rcat(n=10, p=c(0.1,0.3,0.6))

plot(density(rbeta(1000, 1, 1)))

N=10000 # Número de individuos
t.tot=1000 # Número de steps
## Modelo Simulaciones
{
  ## ----------------------- Prior distributions -----------------------
  ## state process var-covar matrix, Sigma
  
  Sigma=riwish(3, matrix(c(1,.3,.3,1),2,2))
  #rWishart(1, 2, )
  #S=matrix(c(1,0,0,1), ncol=2, byrow=T)
  #iSigma=inverse(rWishart(1, 2, S))
  #iSigma[1:2,1:2] ~ dwish(Omega[,], 2)
  #Sigma[1:2,1:2] <- inverse(iSigma[,])
  
  ## mean turn angle: theta[1] (transient), theta[2] (ARS)
  tmp=vector(length = 2)
  tmp[1] = rbeta(1, 1, 1)
  tmp[2] = rbeta(1, 1, 1)
  theta=vector(length = 2)
  theta[1] <- (2 * tmp[1] - 1) * pi
  theta[2] <- (tmp[2] * pi * 2)
  
  ## move persistence: gamma[1] (transient), gamma[2] (ARS)
  gamma=vector(length = 2)
  gamma[1] = rbeta(1,5, 2)
  gamma[2] = rbeta(1,2, 5)
  
  ## behavioural state switching probabilities: alpha ``matrix''
  alpha=vector(length = 2)
  alpha[1] = rbeta(1, 1, 1)
  alpha[2] = rbeta(1, 1, 1)
  
  ## probabilities of initial behavioural state (transient or ARS)
  lambda=vector(length = 2)
  lambda[1] = rbeta(1, 1, 1)
  lambda[2] <- 1 - lambda[1]
  
  ## Loop over the N animals
  for(k in 1:N){
    ## ---------------------- State Process Model ----------------------
    ## randomly specify initial behavioural state, b[1] for animal k
    b[t.tot*(k-1)+1] = rcat(2,lambda)
    
    ## randomly specify location of first state, x[1,1:2] for animal k
    x[t.tot*(k-1)+1,1] ~ dt(first.loc[k,1], tau[t.tot*(k-1)+1,1]^-2, nu[t.tot*(k-1)+1,1])
    x[t.tot*(k-1)+1,2] ~ dt(first.loc[k,2], tau[t.tot*(k-1)+1,2]^-2, nu[t.tot*(k-1)+1,2])
    
    ## randomly specify location of second state, x[2,1:2] for animal k
    x[t.tot*(k-1)+2,1:2] ~ dmnorm(x[t.tot*(k-1)+1,], iSigma[,])
    
    ## Loop over the 2nd to t.tot-1 time steps for animal k
    for(t in (t.tot*(k-1)+2):(t.tot*k-1)){
      ## randomly specify the time t behavioural state, b[t] for animal k
      phi[t,1] <- alpha[b[t-1]]
      phi[t,2] <- 1 - alpha[b[t-1]]
      b[t] ~ dcat(phi[t,])
      
      ## randomly specify the time t+1 location state, x[t+1,1:2] for animal k
      x.mn[t,1] <- x[t,1] + (cos(theta[b[t]]) * (x[t,1] - x[t-1,1]) -
                               sin(theta[b[t]]) * (x[t,2] - x[t-1,2])) * gamma[b[t]]
      x.mn[t,2] <- x[t,2] + (sin(theta[b[t]]) * (x[t,1] - x[t-1,1]) +
                               cos(theta[b[t]]) * (x[t,2] - x[t-1,2])) * gamma[b[t]]
      x[t+1,1:2] ~ dmnorm(x.mn[t,], iSigma[,])
    }
    
    ## randomly specify the last behavioural state, b[t.tot] for animal k
    zeta[k,1] <- alpha[b[t.tot*k-1]]
    zeta[k,2] <- 1 - zeta[k,1]
    b[t.tot*k] ~ dcat(zeta[k,])
    
    ## -----------------------------------------------------------------
    ## ----------------------- Observation Model -----------------------
    for(t in (t.tot*(k-1)+2):(t.tot*k)){
      y[t,1] ~ dt(x[t,1], tau[t,1]^-2, nu[t,1])
      y[t,2] ~ dt(x[t,2], tau[t,2]^-2, nu[t,2])
    }
  } 
}
