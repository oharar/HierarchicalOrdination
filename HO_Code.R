source("HO_greta_fns.R")
library(greta)

# Create data
# Use ant data from mvabund
library(mvabund)
data(antTraits)
y <- antTraits$abund
n_sites <- nrow(antTraits$abund)
n_species <- ncol(antTraits$abund)

# Trait and environmental data
X1 <- scale(antTraits$env)

TR1 <- antTraits$traits
TR1 <- model.matrix(as.formula(paste("~",paste(colnames(antTraits$traits),collapse="+"))),TR1)[,-1]
TR <- t(scale(TR1))

n_latent <- 2# Set 5 LVs
n_traits <- nrow(TR)
n_covs <- ncol(X1)

# Make Greta objects of data
X <- as_data(X1)
TR <- as_data(TR) # traits

# The Model
# This is defined in HO_model.R
#  Y ~ Po(lambda)
# log(lambda) = intB + z.scale%*%LVscales%*%gamma.scale
# z <- X %*% B + epsilon
# gamma <- omega %*% TR + varepsilon

#The current model forces B[1]>0
source("HO_model.R")

m <- model(int, B, omega, epsilon_sd,varepsilon_sd, epsilon, varepsilon,
           LVscales, R2.gamma, R2.z)

# n_samples <- 500; n_warmup <- 500; n_thin <- 1; n_chains <- 10
n_samples <- 1e4;  n_warmup <- 1e3; n_thin <- 1e1; n_chains <- 10

# set initial values only for the delta part of the regularising prior on latent variables
init <- initials(
  delta1 = 10,
  delta2 = rep(10, n_latent - 1)
)

init <- initials(
  LVpars = rep(1,n_latent)
)

#and sample away
draws <- mcmc(m,
              chains = n_chains,
              n_samples = n_samples,
              warmup = n_warmup,
              initial_values = init,
              thin = n_thin,
              sampler = hmc(Lmin = 20, Lmax = 25)
              )
pdf("Posts.pdf")
plot(draws)
dev.off()

# save(draws,file="draws.RData")
draws2 <- extra_samples(draws, n_samples = 4.5e3, thin=10)
# save(draws2,file="draws2.RData")


neff <- coda::effectiveSize(draws)
sort(neff)
