##--- Run this script in order to check the convergence of the mcmc inference chains ---##

library(coda)
library(ggmcmc)

# I - Load the data

mcmc_chains <- readRDS("./inst/extdata/mcmc_chains.rds")

# II - Run the different diagnosis tools

mcmcobj = mcmc.list(mcmc(mcmc_chains[[1]]), mcmc(mcmc_chains[[2]]), mcmc(mcmc_chains[[3]]))
S <- ggs(mcmcobj)

ggs_density(S)

ggs_traceplot(S)

ggs_running(S)

ggs_compare_partial(S)

ggs_autocorrelation(S)

ggs_crosscorrelation(S)

ggs_Rhat(S)
