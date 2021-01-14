##--- Run this script in order to the the mcmc inference ---##

library(rRiskDSMspread)
library(spreadR)
library(ggplot2)
library(abind)

# I - Load the data

grid2D <- readRDS("inst/extdata/grid2D.rds")
grid2D_release_sites <- readRDS("inst/extdata/grid2D_release_sites.rds")
mortalityParDSM <- readRDS("inst/extdata/mortalityParDSM.rds")
mortalityParWT <- readRDS("inst/extdata/mortalityParWT.rds")
diffusionParDSM <- readRDS("inst/extdata/diffusionParDSM.rds")
diffusionParWT <- readRDS("inst/extdata/diffusionParWT.rds")
mrr_data <- readRDS("inst/extdata/mrr_data.rds")
data_mrr_s_path <- system.file("extdata", "mrr_swarms_data.csv", package="rRiskDSMspread")
MRR_swarm <- readr::read_csv(data_mrr_s_path)
data_mrr_c_path <- system.file("extdata", "mrr_compounds_data.csv", package="rRiskDSMspread")
MRR_compounds <- readr::read_csv(data_mrr_c_path)

# II - Define prior and likelihood functions

mPrior = mortalityParWT
mPrior[3,1:2] = mPrior[3, 2:1]
dPrior = diffusionParWT

prior_fun = function(p, trans_p = FALSE){
    
    if(trans_p){
        trans_fun = function(x) exp(x)
        trans_fun_prob = function(x) exp(x) / (1 + exp(x))
    } else{
        trans_fun = function(x) x
        trans_fun_prob = function(x) x
    }
    m = trans_fun(p[1])
    d = trans_fun(p[2])
    a = trans_fun(p[3])
    s = trans_fun(p[4]) 
    pp = trans_fun_prob(p[5])
    pc = trans_fun_prob(p[6])
    ps = trans_fun_prob(p[7])
    priorMu = log(sum(dbeta(exp(-m), mPrior$alpha, mPrior$beta, log = FALSE)*1/nrow(mPrior)))
    priorDelta = log(sum(dlnorm(sqrt((d^2)*(4*7))/1000, dPrior$mu, dPrior$sigma, log=FALSE)*1/nrow(dPrior)))
    priorAlpha = dlnorm(a, -2,1, log = T)
    priorSigma = dlnorm(s, 3.5,.4, log = T)
    dpPSC = dbeta(pp,1.4,1, log = T)
    dpCFR = dbeta(pc,1,1, log = T)
    dpSwarm = dbeta(ps,10,20, log = T)
    
    dprior = priorMu + priorDelta + priorAlpha + priorSigma + dpPSC + dpCFR + dpSwarm
    
    return(-2*(dprior))
}

LL_fun = function(par,df, dist = "", prior = FALSE, many = FALSE, trans_p = FALSE, print_par = FALSE, linearize = FALSE){
    
    idxCFR = which(df$Col_Meth == "Pot")
    idxPSC = which(df$Col_Meth %in% "PSC")
    idxSwarm = which(df$Col_Meth %in% "Swarm")
    
    
    if(trans_p){
        trans_fun = function(x) exp(x)
        trans_fun_prob = function(x) exp(x) / (1 + exp(x))
    } else{
        trans_fun = function(x) x
        trans_fun_prob = function(x) x
    }
    m = trans_fun(par[1])
    d = trans_fun(par[2])
    a = trans_fun(par[3])
    s = trans_fun(par[4]) 
    pp = trans_fun_prob(par[5])
    pc = trans_fun_prob(par[6])
    ps = trans_fun_prob(par[7])
    
    prep.covars <- spreadR::prepare.covariates(grid=grid2D, 
                                               covariates=list(as.data.frame(MRR_swarm[,c("X_m","Y_m")])), sigma=s)
    
    
    
    param=list(delta = d, alpha = a, beta=1, nu=0, mu = m)
    
    if(print_par){
        print(c(m, d, a, s, pp, pc, ps))
    }
    
    # Process model
    
    # In order to really micmic the MRR, we need to simulate 12 (3 x 4) PDE with same parameters except for: 
    #      -  # of released mosquitoes
    #      - Location of the release.
    
    times = seq(0,7.5,by = 0.25)
    releasedMoz = c(1146, 1103, 1158, 1878, 1734, 1655, 1665, 1684, 1673, 1807, 1653, 1813)
    
    
    
    
    if(!linearize){
        timeID = match(df$DayDiff, times)
        obsID = cbind(df$Loc_IDs, timeID)
        
        paramL = releaseSiteL = releaseSizeL = list()
        for(j in 1:12){
            paramL[[j]] = param
            releaseSiteL[[j]] = grid2D_release_sites[j,]
            releaseSizeL[j] =  releasedMoz[j]
        }
        # MRR 1-12
        if(!many){
            SimMoz=NULL
            for(rel in 1:12){
                sim <- simMozzies(params = paramL[[rel]], #param, 
                                  grid=grid2D, 
                                  releaseSite = releaseSiteL[[rel]], #grid2D_release_sites[rel,], 
                                  releaseSize = releaseSizeL[[rel]], #releasedMoz[rel],
                                  covariates=prep.covars, times=times)
                if(dim(sim)[1] < length(times)){
                    sim = matrix(-Inf, nrow = length(times), ncol = grid2D$x.N*grid2D$y.N)
                }
                SimMoz = c(SimMoz, apply(obsID, 1, function(x) sim[x[2], x[1]]))
            }
            SimMoz = matrix(SimMoz, nrow = nrow(df), byrow=F)
        }
        if(many){
            multiParamSim <- function(i){
                pL = paramL[[i]]
                rSt = releaseSiteL[[i]]
                rSz = releaseSizeL[[i]]
                simM <- simMozzies(params = pL,releaseSite = rSt,
                                   grid=grid2D, covariates=prep.covars, times=times, 
                                   releaseSize = rSz)
                simM
            }
            simMany = parallel::mclapply(1:12, multiParamSim)
            arrayMozzies = abind(simMany, along=3)
            SimMoz = sapply(1:12, function(y) apply(obsID, 1, function(x) arrayMozzies[x[2], x[1],y]))
        }
        relE = 0
        
        
        FinalSim=data.frame(site = df$relSite, MRR = df$MRR, Nmoz = NA)
        for(mrr_rel in 1:4){
            for(site in c("A", "B", "C")){
                relE = relE + 1
                idxDF = which(df$MRR == mrr_rel & df$relSite == site)
                idxFS = which(FinalSim$MRR == mrr_rel & FinalSim$site == site)
                FinalSim$Nmoz[idxFS] = SimMoz[idxDF,relE]
            }
        }
        
        # Final Obs are the expected number of moz at the different observation time and location. The actual mozzie population displays additional variation not considered in the pde model, which will be modelled here .
        
        # if n ~ P(lambda), and y | n ~ B(n, p) then y ~ P(p*lambda), properties of Poisson and Binomial dist. We use it here:
        
        LL_CFR= dpois(df$N_mos[idxCFR], pc * sapply(FinalSim$Nmoz[idxCFR], function(x) max(x,1e-10)), log=T)
        LL_PSC= dpois(df$N_mos[idxPSC], pp * sapply(FinalSim$Nmoz[idxPSC], function(x) max(x,1e-10)), log=T)
        LL_Swarm = dpois(df$N_mos[idxSwarm], ps * sapply(FinalSim$Nmoz[idxSwarm], function(x) max(x,1e-10)), log=T)
        
    }
    # We are assuming that we can approximate mutliple releases at the exact same location by one release, where the total number of mozz is the sum of the three releases.
    if(linearize){
        
        df_agg = df %>% group_by(relSite, Col_Meth, DayDiff, Loc_IDs) %>% summarise(N_mos = sum(N_mos)) %>% ungroup()
        
        timeID = match(df_agg$DayDiff, times)
        obsID = cbind(df_agg$Loc_IDs, timeID)
        
        aggregated_releases = dplyr::bind_cols(as.data.frame(grid2D_release_sites), as.data.frame(releasedMoz)) %>%
            group_by(row, col) %>% summarise(releasedMoz = sum(releasedMoz)) %>% ungroup()
        
        SimMoz=NULL
        for(relLoc in 1:3){
            sim <- simMozzies(params = param, grid=grid2D, 
                              releaseSite = as.matrix(aggregated_releases[relLoc,1:2]), 
                              releaseSize = as.numeric(aggregated_releases[relLoc,3]),
                              covariates=prep.covars, times=times)
            if(dim(sim)[1] < length(times)){
                sim = matrix(-Inf, nrow = length(times), ncol = grid2D$x.N*grid2D$y.N)
            }
            SimMoz = c(SimMoz, apply(obsID, 1, function(x) sim[x[2], x[1]]))
        }
        SimMoz = matrix(SimMoz, nrow = nrow(df_agg), byrow=F)
        
        
        relE = 0
        FinalSim=data.frame(site = df_agg$relSite, Nmoz = NA)
        for(site in c("A", "B", "C")){
            relE = relE + 1
            idxDF = which(df_agg$relSite == site)
            idxFS = which( FinalSim$site == site)
            FinalSim$Nmoz[idxFS] = SimMoz[idxDF,relE]
        }
        
        LL_CFR= dpois(df_agg$N_mos[idxCFR], pc * sapply(FinalSim$Nmoz[idxCFR], function(x) max(x,1e-10)), log=T)
        LL_PSC= dpois(df_agg$N_mos[idxPSC], pp * sapply(FinalSim$Nmoz[idxPSC], function(x) max(x,1e-10)), log=T)
        LL_Swarm = dpois(df_agg$N_mos[idxSwarm], ps * sapply(FinalSim$Nmoz[idxSwarm], function(x) max(x,1e-10)), log=T)
    }
    
    LL = sum(LL_CFR, na.rm=TRUE) + sum(LL_PSC, na.rm=TRUE) + sum(LL_Swarm, na.rm=TRUE)
    
    if(sum(FinalSim$Nmoz<0)>0){
        LL = LL - sum(FinalSim$Nmoz<0)*10
    }
    
    if(is.infinite(LL)){print(par)}
    if(is.nan(LL)){print(par)}
    ret = -2*LL
    
    if(prior){
        ret = ret + prior_fun(par)
    }
    
    return(ret)
}

# Testing some set of values:

prior_fun(p = c(-3.63, 5.34, -1.54, 2.81, -1, -1, -1), trans_p = TRUE)

LL_fun(par = c(-3.63, 5.34, -1.54, 2.81, -1, -1, -1), df = mrr_data, trans_p = TRUE, many = TRUE)

# III - Try an optimization routine

resOptM = optim(par = c(-2.63, 10, -3.98, 5.21, -1.8, -3, -2),
                fn = LL_fun, df = mrr_data, prior= FALSE, many = TRUE, trans_p = TRUE, control = list(maxit = 500), hessian = FALSE)


# IV - Run the mcmc inference 

# This takes a long time. Please consider loading the chains from the data repository instead.
mcmc_chains = list()
system.time(
    for(nc in 1:3){
        mcmc_chains[nc] <- modMCMC(f = LL_fun,
                        prior = prior_fun,
                        p = c(0.12, 130, 0.10, 16, 0.20, 0.04, 0.30),
                        df = mrr_data,
                        many = TRUE,
                        print_par = FALSE,
                        niter = 7000,#10000
                        burninlength = 1400,#2000
                        updatecov = 1400,#200
                        upper = c(1, 250, 0.5, 80, 0.5,0.5,0.5),
                        lower = c(0.01, 50, 0.01, 5, 0.01, 0.01, 0.01),
                        outputlength =3500)#5000
    }
)

# V - Plot the chain results

mcmcobj = mcmc.list(mcmc(mcmc_chains[[1]]), mcmc(mcmc_chains[[2]]), mcmc(mcmc_chains[[3]]))
mcmc_mat = as.matrix(mcmcobj)

corMargPlot(mcmc_mat[!duplicated(mcmc_mat),1:4])
corMargPlot(mcmc_mat[!duplicated(mcmc_mat),5:7])

# VI - Save the data to re-use in later scripts

saveRDS(mcmc_chains, file = "/inst/extdata/mcmc_chains.rds")