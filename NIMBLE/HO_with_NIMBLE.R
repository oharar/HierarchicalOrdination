library(nimble)
library(mvabund)
data("antTraits")
Y<-antTraits$abund
X<-scale(antTraits$env)
TR<-scale(model.matrix(~0+.,antTraits$traits))
nLVs <- 2
n <- nrow(Y)
m <- ncol(Y)
t <- ncol(TR)
p = ncol(X)

HO <- nimbleCode({
  for (i in 1:n) {
    for (j in 1:m) {
      eta[i,j] <- beta0[j] + sum(gammas[j,1:d]*LVsd[1:d]*zs[i,1:d])
      log(lambda[i, j]) <- eta[i, j]
      Y[i, j] ~ dpois(lambda[i, j])
    }      
    for (q in 1:d) {
      XB[i, q] <- sum(X[i, 1:p]*B[1:p, q])
      epsilon[i,q] ~ dnorm(0,Sitesd[q]^2)#Residual
      z[i,q] <- XB[i,q] + epsilon[i,q]
    }
  }
  
  for(j in 1:m) {
    for (q in 1:d) {
      omegaTR[j, q] <- sum(TR[j, 1:t]*O[1:t, q])
      varepsilon[j,q] ~ dnorm(0,Speciesd[q]^2) # Residual
      gamma[j,q] <- omegaTR[j,q] + varepsilon[j,q]
    }
    
    beta0[j] ~ dnorm(0, sd=1)
  }
  # Constraints to 0 on upper diagonal
  # stole some code from Boral for this - thanks Francis
  for(i in 1:(d-1)) { 
    for(j in (i+1):(d)) {
      B[i,j] <- 0 
      O[i,j] <- 0
    } 
  }
  
  for(i in 1:d) { 
    ## Sign constraints on diagonal elements
    B[i,i] ~ T(dnorm(0,1),0,Inf)
    O[i,i] ~ dnorm(0,1)#T(dnorm(0,1),0,Inf)
    ## standardizing z and gamma
    zmu[i] <- mean(z[1:n,i])
    zs[1:n,i] <- (z[1:n,i]-zmu[i])/mean((z[1:n,i] - zmu[i])^2) #scale z to unit sd and center
    gammamu[i] <- mean(gamma[1:m,i])
    gammas[1:m,i] <- (gamma[1:m,i]-gammamu[i])/mean((gamma[1:m,i] - gammamu[i])^2) #scale gamma to unit sd and center
    # priors for scales
    Sitesd[i] ~ dexp(1)
    Speciesd[i] ~ dexp(1)
    LVsd[i] ~ dexp(1)
    } 
  
  ## Free lower diagonals
  for(i in 2:d) { 
    for(j in 1:(i-1)) { 
      B[i,j] ~ dnorm(0,1)
      O[i,j] ~ dnorm(0,1)
      } 
    }
  for(i in (d+1):p) { for(j in 1:(d)) { B[i,j] ~dnorm(0,1) } } ## All other elements
  for(i in (d+1):t) { for(j in 1:(d)) { O[i,j] ~dnorm(0,1) } } ## All other elements

})

consts <- list(d = nLVs,
                n = n,
               p = p,
               t = t,
               m = m)

dat <- list(Y = Y, X = X, TR = TR)

inits<-function(consts){
  
  B = matrix(rnorm(consts$d*consts$p),ncol=consts$d)
  B[upper.tri(B)] = 0
  O = matrix(rnorm(consts$d*consts$t),nrow=consts$t)
  O[upper.tri(O)] = 0
  varepsilon = mvtnorm::rmvnorm(consts$m,rep(0,consts$d),diag(consts$d))
  epsilon = mvtnorm::rmvnorm(consts$n,rep(0,consts$d),diag(consts$d))
  list(
    B = B,
    O = O,
    epsilon = epsilon,
    varepsilon = varepsilon,
    Sitesd = rexp(consts$d),
    Speciesd = rexp(consts$d),
    beta0 = rnorm(consts$m),
    LVsd = rexp(consts$d)
  )
}
mod <- nimbleModel(code = HO, name = "HO", constants = consts, inits = inits(consts),
                    data = dat)
model <- compileNimble(mod)
conf <- configureMCMC(model, monitors = c("beta0","Speciesd","Sitesd","LVsd","B","O","z","gamma","epsilon","XB","omegaTR","varepsilon","beta0"), print = TRUE)
mcmc <- buildMCMC(conf)
cmcmc <- compileNimble(mcmc, project = model)


samples <- runMCMC(cmcmc,  niter=55000, nburnin = 5000, thin=10,nchains = 4, samplesAsCodaMCMC = T)
pdf("Posts.pdf")
plot(samples,trace=T,auto.layout = F,density=F)
dev.off()


