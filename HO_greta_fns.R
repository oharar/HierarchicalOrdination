create_betas <- function(n_row,n_lat, sd, name="lat", Abs =FALSE){
  res_raw <- normal(0, 1, dim = c(n_row, n_lat))
  res <- res_raw * sd 
  if(Abs) res[1] <- abs(res[1])
  return(res)
}


#function that creates slopes for traits
create_omega <- function(n_traits,n_latent, sd){
  omega_raw <- normal(0, 1, dim = c(n_latent, n_traits))
  omega <- sweep(omega_raw, 1, sd, FUN = "*") 
  
  return(omega)
}

#function that creates traits for predictors
create_B <- function(n_covs,n_latent, sd){
  B_raw <- normal(0, 1, dim = c(n_covs, n_latent))
  B <- sweep(B_raw, 2, sd, FUN = "*") 
  
  return(B)
}
#function that creates residual for LVs
create_epsilon <- function(n_sites, n_latent, sd){
  if(length(sd)!=n_latent){
    stop("Wrong length for sd.")
  }
  
  epsilon_raw <- normal(0, 1, dim = c(n_sites, n_latent))
  epsilon <- sweep(epsilon_raw, 2, sd, FUN = "*") 
  
  return(epsilon)
}
#function that creates residual for species
create_varepsilon <- function(n_species, n_latent, sd){
  if(length(sd)!=n_latent){
    stop("Wrong length for sd.")
  }
  varepsilon_raw <- normal(0, 1, dim = c(n_latent, n_species))
  varepsilon <- sweep(varepsilon_raw, 1, sd, FUN = "*") 
  return(varepsilon)
}
#function that creates intercepts
create_int <- function(n_species){
  int <- normal(0,1,dim=c(n_species))
  return(int)
}
#function that creates SDs for residuals
create_sds <- function(n_latent){
  sds <- cauchy(0, 3, truncation = c(0, Inf),dim=n_latent)
  return(sds)
}

applyVar <- function(X, MARGIN, sd=FALSE) {
  x.c <- sweep(X, MARGIN, apply(X, MARGIN, "mean"), "-")
  x.var <- apply(x.c^2, MARGIN, "mean")
  if(sd) {
    res <- sqrt(x.var)
  } else {
    res <- x.var
  }
  return(res)
  
}

ScaleVar <- function(x, MARGIN=2) { # column as default
  x.c <- sweep(x, MARGIN, apply(x, MARGIN, "mean"), "-")
  x.sd <- applyVar(x, MARGIN, sd = TRUE)
  x.sc <- sweep(x.c, MARGIN, x.sd, "/")
  return(x.sc)
}

CalcR2 <- function(tot, eps, MARGIN=2) { # column as default
  tot.var <- applyVar(tot, MARGIN)
  eps.var <- applyVar(eps, MARGIN)
  R2 <- 1 - eps.var/tot.var
  return(R2)
}



MakeInits <- function(y, n_latent=2) {
  mod <- gllvm(y, num.lv=n_latent, family="poisson")
  modB <- lm(mod$lvs~-1+X1)
  modO <- lm(mod$params$theta~-1+TR1)
  
  #  init <- initials()
  
  #make initial values for greta
  SimVals <- function(x) {
    x.vec <- c(x)
    e <- rnorm(length(x.vec), 0, sd(x.vec)/sqrt(length(x.vec)))
    x + e
  }
  as.l.df <- function(x) as.list(as.data.frame(x))
  #need to do it a bit special, since we can only have vector-based initial values
  init_int <- list(int=SimVals(mod$params$beta0))
  init_epsilon <- as.l.df(SimVals(residuals(modB)))
  names(init_epsilon) <- paste("epsilon",1:n_latent,sep="")
  init_varepsilon <- as.l.df(SimVals(residuals(modO)))
  names(init_varepsilon) <- paste("varepsilon",1:n_latent,sep="")
  init_B <- as.l.df(SimVals(coef(modB)))
  names(init_B) <- paste("B",1:n_latent,sep="")
  init_omega <- as.l.df(SimVals(coef(modO)))
  names(init_omega) <- paste("omega", 1:n_latent, sep="")
  
  #then we use some magic to make greta accept the initial values
  init <- c(init_int, init_epsilon,init_varepsilon,init_B,init_omega)
  values <- lapply(init, greta:::as_2d_array)
  class(values) <- "initials"
  values
}
