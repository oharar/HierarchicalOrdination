int <- normal(0, 1, dim = c(n_species))

# Site Residuals
# epsilon_sd <- cauchy(0, 3, truncation = c(0, Inf),dim=n_latent)
epsilon_sd <- exponential(1, dim = n_latent)
epsilon <- create_epsilon(n_sites, n_latent, epsilon_sd)

# Species residuals
# varepsilon_sd <- cauchy(0, 3, truncation = c(0, Inf),dim=n_latent)
varepsilon_sd <- exponential(1, dim = n_latent)
varepsilon <- create_varepsilon(n_species, n_latent, varepsilon_sd)

# Species and site predictor coefficients
# omega <- create_omega(n_traits,n_latent,10)
# B <- create_B(n_covs, n_latent, 10)

omega <- t(create_betas(n_traits, n_latent, sd = 1, name = "omega"))
B <- create_betas(n_covs, n_latent, sd = 1, name = "B", Abs = TRUE)

# Latent variable scales
 LVpars <- exponential(1, dim = n_latent)
 LVscales <- rev(cumsum(LVpars))

# a1 <- gamma(2, 1, dim = 1)
# a2 <- gamma(2, 1, dim = 1)
# delta1 <- gamma(a1, 1, dim = 1)
# delta2 <- gamma(a2, 1, dim = n_latent - 1)
# delta <- c(delta1, delta2)
# phi <- gamma(1.5, 1.5, dim = n_latent)
# tau <- phi * cumprod(delta)
# LVscales <- 1 / tau

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

# thingo <-  sweep(z, 2, LVscales, "*")
# eta <- sweep(thingo %*% gamma, 2, int, FUN = "+")

lambda <- exp(eta)
distribution(y) <- poisson(lambda)
