#' Fits a GLLVM to data
#'
#' @param abund Matrix of site by species abundance data
#' @param traits Matrix of species by trait data of species traits.
#' @param env Matrix of site by environmental variable data.
#' @param n_latent The number of latent variables required. Defaults to 1
#' @param FixOmega Either NULL (the default, or the index of the value of Omega to be forced to be positive)
#' @param FixB Either NULL (the default, or the index of the value of B to be forced to be positive)
#' @param ... More arguments (none used at the moment)
#' @return A greta model object
#' @export
#' @import greta

CreateModel <- function(abund, traits, env, n_latent=1, 
                        FixOmega=NULL, FixB=NULL, ...) {
  if(nrow(abund)!=nrow(env)) stop("number of rows of abund and env should be equal")
  if(ncol(abund)!=nrow(traits)) stop("number of columns of abund should equal number of rows of env")
  # Format data
  n_sites <- nrow(abund)
  n_species <- ncol(abund)
  
  # Trait and environmental data
  Env <- scale(env)
  
  Traits <- model.matrix(as.formula(paste("~",paste(colnames(traits),collapse="+"))),traits)[,-1]
  Traits <- t(scale(Traits))
  
  n_traits <- nrow(Traits)
  n_covs <- ncol(Env)
  
  # Make Greta objects of data
  EnvCov <- as_data(Env)
  TraitCov <- as_data(Traits) # traits
  
  # The Model
  # This is defined in HO_model.R
  #  Y ~ Po(lambda)
  # log(lambda) = intB + z.scale%*%LVscales%*%gamma.scale
  # z <- X %*% B + epsilon
  # gamma <- omega %*% TR + varepsilon
  
  int <- normal(0, 1, dim = c(n_species))
  
  # Site Residuals
  epsilon_sd <- exponential(0.99, dim = n_latent)
  epsilon <- create_epsilon(n_sites, n_latent, epsilon_sd)
  epsilon <- create_epsilon(n_sites, n_latent, rep(1, n_latent))
  
  # Species residuals
  varepsilon_sd <- exponential(1, dim = n_latent)
  # varepsilon <- create_varepsilon(n_species, n_latent, varepsilon_sd)
  varepsilon <- t(create_epsilon(n_species, n_latent, varepsilon_sd))
  
  # Species and site predictor coefficients
  omega <- t(create_betas(n_traits, n_latent, sd = 1, 
                          name = "omega", Abs = FixOmega))
  B <- create_betas(n_covs, n_latent, sd = 1, name = "B", 
                    Abs = FixB)
  
  # Latent variable scales
  LVpars <- exponential(1, dim = n_latent)
  LVscales <- rev(cumsum(LVpars))
  
  # row and column scores
  z <- X %*% B + epsilon # site
  gamma <- omega %*% TR + varepsilon # species
  R2.z <- CalcR2(tot = z, eps = epsilon, MARGIN = 2)
  R2.gamma <- CalcR2(tot = gamma, eps = varepsilon, MARGIN = 1)
  
  # scale row & column scores to be mean 0, sd 1
  z.sc <- ScaleVar(z, MARGIN = 2)
  gamma.sc <- ScaleVar(gamma, MARGIN = 1)
  
  thingo <-  sweep(z.sc, 2, LVscales, "*")
  eta <- sweep(thingo %*% gamma.sc, 2, int, FUN = "+")
  
  lambda <- exp(eta)
  distribution(y) <- poisson(lambda)
  
  m <- model(int, B, omega, epsilon_sd,varepsilon_sd, epsilon, varepsilon,
             LVscales, R2.gamma, R2.z)
  return(m)
}
  


# This is for an infinite factor model. That probably needs a different function
# a1 <- gamma(2, 1, dim = 1)
# a2 <- gamma(2, 1, dim = 1)
# delta1 <- gamma(a1, 1, dim = 1)
# delta2 <- gamma(a2, 1, dim = n_latent - 1)
# delta <- c(delta1, delta2)
# phi <- gamma(1.5, 1.5, dim = n_latent)
# tau <- phi * cumprod(delta)
# LVscales <- 1 / tau
