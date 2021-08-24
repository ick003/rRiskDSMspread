##--- Run this script in order to run the cross validation on MRR 5 ---##

library(rRiskDSMspread)
library(spreadR)
library(ggplot2)
library(scales)
library(dplyr)
library(coda)

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

data_mrr_path <- system.file("extdata", "mrr_data.csv", package="rRiskDSMspread")

# II - Setting up MRR 5 data 

MRR <- readr::read_csv(data_mrr_path)
MRR = MRR %>% filter(.data$MRR == 5)
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



nSim = 100

paramL = releaseSiteL = releaseSizeL = list()
for(j in 1:nSim){
    idxJ = sample(nrow(mcmc_mat),1)
    tmp = list(delta=mcmc_mat[idxJ,2],alpha=mcmc_mat[idxJ,3],beta=1,nu=0,mu=mcmc_mat[idxJ,1])
    paramL[[length(paramL)+1]] = tmp
    releaseSiteL[[length(paramL)]] = releaseSitesCoded
    releaseSizeL[[length(paramL)]] = 5992
    
}
prep.covars <- spreadR::prepare.covariates(grid=grid2D, 
                                           covariates=list(as.data.frame(MRR_swarm[,c("X_m","Y_m")])), sigma=mean((mcmc_mat[,4])))

RelS = 1
lastDaySim = as.numeric(max(MRR$DayDiff[MRR$Col_Meth == "Swarm"]))+1
timesSim = seq(0,lastDaySim,0.125)


# IV - Run the simulation 

simMRR5 <- manyReleaseMozzies(manyParams=paramL, mc.cores = 6,
                              manyReleaseSites = releaseSiteL,
                              manyReleaseSizes = releaseSizeL,
                              grid=grid2D, covariates=prep.covars, times=timesSim)



# V - Match simulation outputs and observations from MRR 5

library(plyr)

CP_Coord_UTM = LongLatToUTM(MRR$Long_CP, MRR$Lat_CP, 30)

colnames(CP_Coord_UTM) <- c("ID", "Cap_UTM_x", "Cap_UTM_y")

MRR = cbind(MRR, CP_Coord_UTM[,2:3])

locations5 = MRR[,c("Cap_UTM_x", "Cap_UTM_y")]
IDs5 <- as.numeric( class::knn1( expand.grid( grid2D$x.mid, grid2D$y.mid), locations5, 1:(grid2D$x.N*grid2D$y.N)))
MRR$Loc_IDs = IDs5

dfObs = ddply(MRR, .(relSite,Col_Meth, DayDiff, Loc_IDs),  summarize, N_mos = sum(N_mos))
#dfObs <- MRR %>% group_by(relSite, Col_Meth, DayDiff, Loc_IDs) %>% summarise(N_mos = sum(N_mos)) %>% ungroup()
catchType = c("PSC","Swarm")
idxCatch = which(MRR$Col_Meth %in% catchType)
IdsCatch5 = unique(dfObs$Loc_IDs[intersect(idxCatch, which(dfObs$Col_Meth %in% catchType & dfObs$N_mos > 0))])


simMRR5_CV = simMRR5
simMRR5_CV$mods = simMRR5$mods[,IdsCatch5,]
attr(simMRR5_CV$mods,"times") <- timesSim

# VI - Plot comparison

png("./inst/manuscript/Draft/Figure3.png", units="in", width=10, height=10, res=600)
par(mfrow = c(5,3), mar = c(2,4,1,1))
prob_catch_swarm = mean((mcmc_mat[,7]) / (1 + (mcmc_mat[,7]) ))
SSE = nse = nci = 0
for(i in 1:length(IdsCatch5)){
    CV_OBS = dfObs[dfObs$Loc_IDs == IdsCatch5[i],]
    SSE = SSE + sum((rowMeans(simMRR5_CV$mods[,i,]*prob_catch_swarm)[timesSim %in% CV_OBS$DayDiff] - CV_OBS$N_mos)^2)
    nse = nse + nrow(CV_OBS)
    cil = qpois(0.05, rowMeans(simMRR5_CV$mods[,i,]*prob_catch_swarm))[timesSim %in% CV_OBS$DayDiff]
    ciu = qpois(0.95, rowMeans(simMRR5_CV$mods[,i,]*prob_catch_swarm))[timesSim %in% CV_OBS$DayDiff]
    nci = nci + sum(sapply(1:nrow(CV_OBS), function(x) between(CV_OBS$N_mos[x], cil[x], ciu[x])), na.rm=TRUE)
        
    if(max(CV_OBS$N_mos)==0) next
    plot(timesSim,simMRR5_CV$mods[,i,1]*prob_catch_swarm, type = "l", 
         #ylim = range(simMRR5_CV$mods[,i,]*prob_catch_swarm) + c(0,1), 
         ylim = c(0,6),
         col = alpha(rgb(1,0.6,0), 0.1), bty = "n", 
         xlab = "", ylab = expression(lambda))
    polygon(x = c(timesSim, rev(timesSim)), 
            y = c(qpois(0.05, rowMeans(simMRR5_CV$mods[,i,]*prob_catch_swarm)), 
                  rev(qpois(0.95, rowMeans(simMRR5_CV$mods[,i,]*prob_catch_swarm)))), 
            border = NA, col = alpha(rgb(1,0.6,0),0.15))
    #for(j in 2:dim(simMRR5_CV$mods)[3]){
    #    points(timesSim, simMRR5_CV$mods[,i,j]*prob_catch_swarm, type=  "l", col = alpha(rgb(0,0,0), 0.1))
    #}
    points(timesSim, rowMeans(simMRR5_CV$mods[,i,]*prob_catch_swarm), type = "l", lwd = 2, 
           col = alpha(rgb(1,0.6,0), 0.75))
    
    
    points(CV_OBS$DayDiff+0.125, CV_OBS$N_mos, col = "red", pch = 3, cex = 0.5)
    text(x = 4, y = 6,label = paste("Loc:",IdsCatch5[i]))
    
    
}

plot.new()
legend("center", c("Expected # mosquitoes", 
                   "90% credible interval", 
                   "Observation"), 
       lty = c(1, 1, NA), cex = 1.25, col = c(alpha(rgb(1,0.6,0), 0.75), alpha(rgb(1,0.6,0), 0.15), "red"), 
       lwd = c(2, 8, NA), pch = c(NA, NA, 3), box.lwd = 0.1)
dev.off()

par(mfrow=c(1,1))
plot(MRR_compounds$X_m[MRR_compounds$Locality == "Village"], MRR_compounds$Y_m[MRR_compounds$Locality == "Village"])
points(LocMRR$X, LocMRR$Y, col = "blue")
points(expand.grid( grid2D$x.mid, grid2D$y.mid)[IdsSwarm5,], col = "red", cex=0.5)
text(expand.grid( grid2D$x.mid, grid2D$y.mid)[IdsSwarm5,], labels = 1:14, col = "red")

