##--- Run this script in order test the pde model ---##

library(spreadR)
library(ggplot2)

# I - Load the data

grid2D <- readRDS("inst/extdata/grid2D.rds")
grid2D_release_sites <- readRDS("inst/extdata/grid2D_release_sites.rds")
mortalityParDSM <- readRDS("inst/extdata/mortalityParDSM.rds")
mortalityParWT <- readRDS("inst/extdata/mortalityParWT.rds")
diffusionParDSM <- readRDS("inst/extdata/diffusionParDSM.rds")
diffusionParWT <- readRDS("inst/extdata/diffusionParWT.rds")
data_mrr_s_path <- system.file("extdata", "mrr_swarms_data.csv", package="rRiskDSMspread")
MRR_swarm <- readr::read_csv(data_mrr_s_path)

# II - Generate random variables from prior distributions
set.seed(1)

nSim = 3

deltaL <- (par2sample(par = diffusionParWT, type = "diffusion",  n = nSim, weights = rep(1, nrow(diffusionParWT)))*1000)^2 / (3.14 * 1)
muL <- par2sample(par = mortalityParWT, type = "mortality",  n = nSim, weights = rep(1, nrow(mortalityParWT)))
alphaL = rlnorm(nSim, -2, 1)
sigmaL = rlnorm(nSim, 3.5, 4)

paramL = list()
for(j in 1:nSim){
    tmp = list(delta=deltaL[j],alpha=alphaL[j],nu=0,mu=muL[j],beta=1)
    paramL[[length(paramL)+1]] = tmp
    
}

# III - Run a simulation from a single release site

prep.covars <- prepare.covariates( grid=grid2D,
                                   covariates=list(as.data.frame(MRR_swarm[,c("X_m","Y_m")])),
                                   sigma= mean(sigmaL))


lastDay = 3
times = seq(0,lastDay,1)

set.seed(1)
simM <- simMozzies(params = paramL[[1]],releaseSite = grid2D_release_sites[1,],
                   grid=grid2D, covariates=prep.covars, times=times, 
                   releaseSize = 5000, boundaryCond = "Neumann")

# IV - Plot the results

mozImage(simM, grid2D, mfrow = c(2,2))

# V - Run a simulation from multiple release sites

paramL = releaseSiteL = releaseSizeL = list()
for(j in 1:nSim){
    tmp = list(delta=deltaL[j],alpha=alphaL[j],nu=0,mu=muL[j], beta=1)
    paramL[[length(paramL)+1]] = tmp
    releaseSiteL[[j]] = grid2D_release_sites[j,]
    releaseSizeL[j] = 5000/nSim
}
lastDay = 3
times = seq(0,lastDay,1)

multiParamSim <- function(i){
    pL = paramL[[i]]
    rSt = releaseSiteL[[i]]
    rSz = releaseSizeL[[i]]
    simM <- simMozzies(params = pL,releaseSite = rSt,
                       grid=grid2D, covariates=prep.covars, times=times, 
                       releaseSize = rSz, boundaryCond = "Neumann")
    simM
}
testParallel = parallel::mclapply(1:3, multiParamSim)

aggMozzies = abind::abind(testParallel, along=3)

# VI - Plot the results

par(mfrow = c(2,2))
for(i in 1:length(times)){
    soln <- matrix(apply(aggMozzies, MARGIN = 1:2, FUN = sum)[i,-1], nrow = grid2D$x.N, ncol = grid2D$y.N)
    image(x=grid2D$x.mid, y=grid2D$y.mid, soln)
    usr <- par( "usr")
    nAtTime <- rowSums( soln[,-1])
    text( usr[2], usr[3], paste0("n=",round( sum(soln), 3)), adj=c(1.25,-0.5))
}
