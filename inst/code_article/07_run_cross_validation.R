##--- Run this script in order to run the cross validation on MRR 5 ---##

library(rRiskDSMspread)
library(spreadR)
library(ggplot2)
library(scales)

# I - Load the data

grid2D <- readRDS("inst/extdata/grid2D.rds")
grid2D_release_sites <- readRDS("inst/extdata/grid2D_release_sites.rds")
mortalityParDSM <- readRDS("inst/extdata/mortalityParDSM.rds")
mortalityParWT <- readRDS("inst/extdata/mortalityParWT.rds")
diffusionParDSM <- readRDS("inst/extdata/diffusionParDSM.rds")
diffusionParWT <- readRDS("inst/extdata/diffusionParWT.rds")
data_mrr_s_path <- system.file("extdata", "mrr_swarms_data.csv", package="rRiskDSMspread")
MRR_swarm <- readr::read_csv(data_mrr_s_path)
data_mrr_c_path <- system.file("extdata", "mrr_compounds_data.csv", package="rRiskDSMspread")
MRR_compounds <- readr::read_csv(data_mrr_c_path)
data_mrr_l_path <- system.file("extdata", "mrr_location_data.csv", package="rRiskDSMspread")
LocMRR <- readr::read_csv(data_mrr_l_path)
mrr_data <- readRDS("inst/extdata/mrr_data.rds")
mcmc_chains <- readRDS("./inst/extdata/mcmc_chains.rds")

# II - Setting up MRR 5 data 

MRR = MRR[which(MRR$MRR == 5),]
MRR$Date = as.Date(as.character(MRR$Date), format="%Y-%m-%d")
MRR$Date = as.character(MRR$Date)
names(MRR)<- c("MRR", "Cap_Date","Long_CP", "Lat_CP", "Col_Meth","cVar", 
               "N_mos", "relSite")

DateMRR = unique(MRR[,c(1,2)])
DateMRR = DateMRR[ DateMRR$Cap_Date == ave(DateMRR$Cap_Date, DateMRR$MRR, FUN=min), ]
names(DateMRR)[2] <- c("Release date")
names(LocMRR) = c("MRR","Release site", "Longitude", "Latitude", "Number released")

MRR$DayDiff =  as.Date(MRR$Cap_Date) - as.Date(DateMRR$`Release date`)
MRR$Col_Meth = as.factor(MRR$Col_Meth)


LocMRR = LongLatToUTM(LocMRR$Longitude[13], LocMRR$Latitude[13], 30)
LocMRR$ID = c("C")
names(LocMRR)[1] <- c("relSite")

releaseSites = LocMRR

IDs_ReleaseSites <- as.numeric( class::knn1( expand.grid( grid2D$x.mid, grid2D$y.mid), releaseSites[,2:3], 1:(grid2D$x.N*grid2D$y.N)))

releasedMoz = c(5992)

releaseSitesCoded = NULL
for(rel in 1:1){
    mat = matrix(0, nrow = grid2D$x.N, ncol = grid2D$y.N)
    mat[IDs_ReleaseSites[rel]]<- releasedMoz[rel]
    releaseSitesCoded = rbind(releaseSitesCoded ,which(mat>0, arr.ind=T))
}

# III - Draw samples from mcmc samples

mcmcobj = mcmc.list(mcmc(mcmc_chains[[1]]), mcmc(mcmc_chains[[2]]), mcmc(mcmc_chains[[3]]))
mcmc_mat = as.matrix(mcmcobj)



nSim = 1000

paramL = releaseSiteL = releaseSizeL = list()
for(j in 1:nSim){
    idxJ = sample(nrow(mcmc_mat),1)
    tmp = list(delta=exp(mcmc_mat[idxJ,2]),alpha=exp(mcmc_mat[idxJ,3]),beta=1,nu=0,mu=exp(mcmc_mat[idxJ,1]))
    paramL[[length(paramL)+1]] = tmp
    releaseSiteL[[length(paramL)]] = releaseSitesCoded
    releaseSizeL[[length(paramL)]] = 5992
    
}
prep.covars <- spreadR::prepare.covariates(grid=grid2D, 
                                           covariates=list(as.data.frame(MRR_swarm[,c("X_m","Y_m")])), sigma=mean(exp(mcmc_mat[,4])))

RelS = 1
lastDaySim = as.numeric(max(MRR$DayDiff[MRR$Col_Meth == "Swarm"]))+1
timesSim = seq(0,lastDaySim,0.125)


# IV - Run the simulation 

simMRR5 <- manyReleaseMozzies(manyParams=paramL, mc.cores = 6,
                              manyReleaseSites = releaseSiteL,
                              manyReleaseSizes = releaseSizeL,
                              grid=grid2D, covariates=prep.covars, times=timesSim)

# V - Match simulation outputs and observations from MRR 5

CP_Coord_UTM = LongLatToUTM(MRR$Long_CP, MRR$Lat_CP, 30)

colnames(CP_Coord_UTM) <- c("ID", "Cap_UTM_x", "Cap_UTM_y")

MRR = cbind(MRR, CP_Coord_UTM[,2:3])

locations5 = MRR[,c("Cap_UTM_x", "Cap_UTM_y")]
IDs5 <- as.numeric( class::knn1( expand.grid( grid2D$x.mid, grid2D$y.mid), locations5, 1:(grid2D$x.N*grid2D$y.N)))
MRR$Loc_IDs = IDs5

dfObs = ddply(MRR, .(relSite,Col_Meth, DayDiff, Loc_IDs),  summarize, N_mos = sum(N_mos))
idxSwarm = which(MRR$Col_Meth %in% "Swarm")
IdsSwarm5 = unique(dfObs$Loc_IDs[intersect(idxSwarm, which(dfObs$Col_Meth %in% "Swarm" & dfObs$N_mos > 0))])


simMRR5_CV = simMRR5
simMRR5_CV$mods = simMRR5$mods[,IdsSwarm5,]
attr(simMRR5_CV$mods,"times") <- timesSim

# VI - Plot comparison

par(mfrow = c(5,3), mar = c(2,2,1,1))
prob_catch_swarm = mean(exp(mcmc_mat[,7]) / (1 + exp(mcmc_mat[,7]) ))
for(i in 1:14){
    plot(timesSim,simMRR5_CV$mods[,i,1]*prob_catch_swarm, type = "l", ylim = range(simMRR5_CV$mods[,i,]*prob_catch_swarm), 
         col = alpha(rgb(0,0,0), 0.01))
    for(j in 2:dim(simMRR5_CV$mods)[3]){
        points(timesSim, simMRR5_CV$mods[,i,j]*prob_catch_swarm, type=  "l", col = alpha(rgb(0,0,0), 0.01))
    }
    points(timesSim, rowMeans(simMRR5_CV$mods[,i,]*prob_catch_swarm), type = "l", lwd = 2)
    CV_OBS = dfObs[dfObs$Loc_IDs == IdsSwarm5[i],]
    points(CV_OBS$DayDiff+0.125, CV_OBS$N_mos, col = "red")
}
plot(MRR_compounds$X_m[MRR_compounds$Locality == "Village"], MRR_compounds$Y_m[MRR_compounds$Locality == "Village"])
points(LocMRR$X, LocMRR$Y, col = "blue")
points(expand.grid( grid2D$x.mid, grid2D$y.mid)[IdsSwarm5,], col = "red", cex=0.1)
text(expand.grid( grid2D$x.mid, grid2D$y.mid)[IdsSwarm5,], labels = 1:14, col = "red")

