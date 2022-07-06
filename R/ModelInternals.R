
# Create parameters for covariate effects on row/column LVs
# arguments 
#   n_row: number of rows
#   n_lat: number of LVs 
#   sd: standard deviation of betas
#   name: Name. Defaults to "lat", currently not used
#   Abs: Whether a beta should be forced to be positive. 
#      Defaults to NULL, i.e. nothing is contrained. 
#      Alternative is an integer to say which beta should be constrained
  
create_betas <- function(n_row,n_lat, sd, name="lat", Abs = NULL){
  if(!is.null(Abs)) {
    if(!is.integer(Abs)) stop("Abs should be NULL or an integer")
    if(Abs>n_row) stop("Abs cannot be larger than n_row")
  }
  res_raw <- normal(0, 1, dim = c(n_row, n_lat))
  res <- res_raw * sd 
  if(!is.null(Abs)) res[Abs,] <- abs(res[Abs,])
  return(res)
}


#function that creates slopes for traits
#  Deprecated, because create_betas() does the same job
# create_omega <- function(n_traits,n_latent, sd){
#   omega_raw <- normal(0, 1, dim = c(n_latent, n_traits))
#   omega <- sweep(omega_raw, 1, sd, FUN = "*") 
#   
#   return(omega)
# }


#function that creates traits for predictors
# Deprecated
# create_B <- function(n_covs,n_latent, sd){
#   B_raw <- normal(0, 1, dim = c(n_covs, n_latent))
#   B <- sweep(B_raw, 2, sd, FUN = "*") 
#   
#   return(B)
# }

# function that creates residual for LVs
# Arguments
# n_sites: number of sites/species, 
# n_latent: number of LVs
# sd: standard deviation
create_epsilon <- function(n_sites, n_latent, sd){
  if(length(sd)!=n_latent){
    stop("Wrong length for sd.")
  }
  
  epsilon_raw <- normal(0, 1, dim = c(n_sites, n_latent))
  epsilon <- sweep(epsilon_raw, 2, sd, FUN = "*") 
  
  return(epsilon)
}

#function that creates residual for species
# Deprecated as does same job as create_epsilon()
# create_varepsilon <- function(n_species, n_latent, sd){
#   if(length(sd)!=n_latent){
#     stop("Wrong length for sd.")
#   }
#   varepsilon_raw <- normal(0, 1, dim = c(n_latent, n_species))
#   varepsilon <- sweep(varepsilon_raw, 1, sd, FUN = "*") 
#   return(varepsilon)
# }

#function that creates intercepts
# Argument
#  n_species: number of species
create_int <- function(n_species){
  int <- normal(0,1,dim=c(n_species))
  return(int)
}
#function that creates SDs for residuals
# Argument
#  n_latent: number of LVs
create_sds <- function(n_latent){
  sds <- cauchy(0, 3, truncation = c(0, Inf),dim=n_latent)
  return(sds)
}

# Calculate a variance of rows/columns 
#  Arguments
#  X: matrix, 
#  MARGIN: Margin of X that variance will be calculated over
#  sd: Should we return the standard deviation, 
# even though the function is called applyVar(). Defaults to FALSE
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

# Function to centre and scale a matrix, by MARGIN
# Arguments
#   x: matrix 
#   MARGIN: which margin to sweep over (defaults to 2, i.e. scale each column)
ScaleVar <- function(x, MARGIN=2) { # column as default
  x.c <- sweep(x, MARGIN, apply(x, MARGIN, "mean"), "-")
  x.sd <- applyVar(x, MARGIN, sd = TRUE)
  x.sc <- sweep(x.c, MARGIN, x.sd, "/")
  return(x.sc)
}

# Function to calculates something like R^2
#  Arguments:
#    tot: matrix of total variance
#    eps: matrix of residual variance
#    MARGIN: which margin to calculare R^2 for (defaults to column)
CalcR2 <- function(tot, eps, MARGIN=2) { # column as default
  tot.var <- applyVar(tot, MARGIN)
  eps.var <- applyVar(eps, MARGIN)
  R2 <- 1 - eps.var/tot.var
  return(R2)
}

