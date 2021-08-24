##--- Run this script in order test the pde model ---##

library(rRiskDSMspread)
library(spreadR)
library(ggplot2)
library(stringr)

# I - Load the data

data_mrr_path <- system.file("extdata", "mrr_data.csv", package="rRiskDSMspread")
grid2D <- readRDS("inst/extdata/grid2D.rds")
grid2D_release_sites <- readRDS("inst/extdata/grid2D_release_sites.rds")
mortalityParDSM <- readRDS("inst/extdata/mortalityParDSM.rds")
mortalityParWT <- readRDS("inst/extdata/mortalityParWT.rds")
diffusionParDSM <- readRDS("inst/extdata/diffusionParDSM.rds")
diffusionParWT <- readRDS("inst/extdata/diffusionParWT.rds")
data_mrr_s_path <- system.file("extdata", "mrr_swarms_data.csv", package="rRiskDSMspread")
MRR_swarm <- readr::read_csv(data_mrr_s_path)
mcmc_chains <- readRDS("inst/extdata/mcmc_chains.rds")
parDSM <- readRDS("inst/extdata/parDSM.rds")
data_mrr_l_path <- system.file("extdata", "mrr_location_data.csv", package="rRiskDSMspread")
# II - Setup simulation environmental conditions

MRR <- readr::read_csv(data_mrr_path)
MRR = MRR[which(MRR$MRR == 1),]
MRR$Date = as.Date(as.character(MRR$Date), format="%Y-%m-%d")
MRR$Date = as.character(MRR$Date)
names(MRR)<- c("MRR", "Cap_Date","Long_CP", "Lat_CP", "Col_Meth","cVar", 
               "N_mos", "relSite")

DateMRR = unique(MRR[,c(1,2)])
DateMRR = DateMRR[ DateMRR$Cap_Date == ave(DateMRR$Cap_Date, DateMRR$MRR, FUN=min), ]
names(DateMRR)[2] <- c("Release date")

LocMRR <- readr::read_csv(data_mrr_l_path)
names(LocMRR) = c("MRR","Release site", "Longitude", "Latitude", "Number released")

MRR$DayDiff =  as.Date(MRR$Cap_Date) - as.Date(DateMRR$`Release date`)
MRR$Col_Meth = as.factor(MRR$Col_Meth)


LocMRR = LongLatToUTM(LocMRR$Longitude[1], LocMRR$Latitude[1], 30)
LocMRR$ID = c("A")
names(LocMRR)[1] <- c("relSite")

releaseSites = LocMRR

IDs_ReleaseSites <- as.numeric( class::knn1( expand.grid( grid2D$x.mid, grid2D$y.mid), releaseSites[,2:3], 1:(grid2D$x.N*grid2D$y.N)))

releasedMoz = c(5000)

releaseSitesCoded = NULL
for(rel in 1:1){
    mat = matrix(0, nrow = grid2D$x.N, ncol = grid2D$y.N)
    mat[IDs_ReleaseSites[rel]]<- releasedMoz[rel]
    releaseSitesCoded = rbind(releaseSitesCoded ,which(mat>0, arr.ind=T))
}



# III - Generate random variables from prior distributions

set.seed(1)

nSim = nrow(parDSM)/10

paramL = releaseSiteL = releaseSizeL = list()
for(j in 1:nSim){
    idxJ = sample(nrow(parDSM),1)
    tmp = list(delta=parDSM$delta[j],alpha=parDSM$alpha[j],beta=1,nu=0,mu=parDSM$mu[j])
    paramL[[length(paramL)+1]] = tmp
    releaseSiteL[[length(paramL)]] = releaseSitesCoded
    releaseSizeL[[length(paramL)]] = 5000
    
}
prep.covars <- spreadR::prepare.covariates(grid=grid2D, 
                                           covariates=list(as.data.frame(MRR_swarm[,c("X_m","Y_m")])), sigma=mean(parDSM$sigma))

RelS = 1
lastDaySim = 12
timesSim = seq(0,lastDaySim,1)

# IV - Run Ag(DSM) pde simulation

# It will take some time, especially for a large number of simulations. You can load the saved data instead.
sim_mrr_full <- manyReleaseMozzies(manyParams=paramL, mc.cores = 7,
                              manyReleaseSites = releaseSiteL,
                              manyReleaseSizes = releaseSizeL,
                              grid=grid2D, covariates=prep.covars, times=timesSim)
# To save space we only keep the slices of time we want to plot
sim_mrr <- sim_mrr_full
timesPlot = c(2,5,9,12)
sim_mrr$mods <- sim_mrr$mods[timesPlot,,]

# V - Plots

# Plot survival
png(filename = "inst/manuscript/Draft/figure/Figure4.png", width = 20, height = 15, res = 600, units = "cm")
timeSurv = seq(0,30,1)
par(mfrow=c(1,1), mar = c(4,4,1,1))
par(mar = c(4.5,4.5,1,1))
layout(matrix(c(1,2), 2, 1, byrow = TRUE),
       widths=c(1), heights=c(3,1))
TrFun = function(x){log(x)}
plot(timeSurv,TrFun(5000*exp(-mean(parDSM$mu[1:nSim])*timeSurv)), type="l", xlab = "Days since release", ylab=expression(paste("Expected # sterile males (",lambda, ")")), 
     bty='n', axes=F, col = rgb(1,0,0,0.5), pch=3, cex=0.5, ylim = TrFun(c(1e-1, 6000)))

axis(1, at = seq(0,30,length.out = 11));
axis(2, at = TrFun(c(1 %o% 10^(-1:3))), labels = c(1 %o% 10^(-1:3)), las=1, cex.axis=0.75)
for(i in 1:nSim){
    points(timeSurv,TrFun(5000*exp(-parDSM$mu[i]*timeSurv)), type="l", col=rgb(0,0,0,0.05), pch=3, cex=0.2)
    #  points(attr(simM$mods,"times"),rowSums(simM$mods[, -1,i]), col = rgb(0,0,0,0.1), pch=3, cex=0.5)
}
points(timeSurv,TrFun(5000*exp(-mean(parDSM$mu)*timeSurv)), type="l",col = rgb(1,0,0,1), pch=3, cex=0.5)
points(timeSurv,TrFun(5000*exp(-quantile(parDSM$mu,0.5)*timeSurv)), type="l",col = rgb(0,0,1,0.6), pch=3, cex=0.5)
points(timeSurv,TrFun(5000*exp(-quantile(parDSM$mu,0.05)*timeSurv)), type="l",col = rgb(0,0,1,0.6), lty=2, pch=3, cex=0.5)
points(timeSurv,TrFun(5000*exp(-quantile(parDSM$mu,0.95)*timeSurv)), type="l",col = rgb(0,0,1,0.6), lty=2, pch=3, cex=0.5)

arrows(x0=15, y0=TrFun(500), x1=8.6, y1=TrFun(1.5), col='black', length=0.1, lwd=1)
text(x = 15.5, y = TrFun(500), stringr::str_wrap("Expecting less than 1 sterile male by Day 8", width = 25), adj=c(0,0))
points(x = 8.3, y = TrFun(1), col = "red", cex=0.5)

arrows(x0=19, y0=TrFun(50), x1=10.9, y1=TrFun(0.15), col='black', length=0.1, lwd=1)
text(x = 19.5, y = TrFun(50), stringr::str_wrap("Expected 90% probability of no sterile male by Day 10", width = 27), adj=c(0,0))
points(x = 10.6, y = TrFun(0.1), col = "red", cex=0.5)

par(mar = c(1.5,0,0,0))
plot(c(0,0), xlim = c(1,10), ylim=c(8,10),type="n", bty='n', axes=F, xlab="", ylab="")
legend("center", lty=c(1,1,2,1),
       col=c(rgb(1,0,0,0.6), rgb(0,0,1,0.6), rgb(0,0,1,0.6), rgb(0,0,0,0.1)), 
       legend = c("Mean", "Median", "5th and 95th percentiles", "Simulation"), ncol=2)
dev.off()
# Plot dispersal

png(filename = "inst/manuscript/Draft/figure/Figure5bis.png", res = 600, units = "cm", width = 20, height = 20)
set.seed(3)
par(mfrow=c(2,2), mar = c(0.1,0.1,0.1,0.1))
rangeX = range(MRR_swarm$X_m)+50*c(-1,1)
rangeY = range(MRR_swarm$Y_m)+50*c(-1,1)
for(j in  1:length(timesPlot)){
    if(nSim > 1){zM = matrix(apply(sim_mrr$mods[j,,],1,mean), nrow = grid2D$x.N, ncol = grid2D$y.N)} else{
        zM = matrix(sim_mrr$mods[j,], nrow = grid2D$x.N, ncol = grid2D$y.N)}
    ticks<-10^seq(-2,3,1)
    plot(MRR_swarm$X_m, MRR_swarm$Y_m, cex=0.2, pch = 8, yaxt="n", xaxt="n", xlim = rangeX, ylim = rangeY)
    # contour(x = grid2D$x.mid, y = grid2D$y.mid, z = zM,add = TRUE, levels = ticks[3], col = "red", labels = "1",
    #         xlim = rangeX, ylim = rangeY)
    contour(x = grid2D$x.mid, y = grid2D$y.mid, z = zM,add = TRUE, levels = ticks[3], col = "red", labels =  ">=1",
            xlim = rangeX, ylim = rangeY, lwd = 2)
    contour(x = grid2D$x.mid, y = grid2D$y.mid, z = zM,add = TRUE, levels = ticks[1], col = "orange", labels = "P99",
            xlim = rangeX, ylim = rangeY, lwd = 2)
    lines(c(max(rangeX) - 50, max(rangeX)), rep(min(rangeY)+50,2))
    text(max(rangeX) - 25, min(rangeY)+50, "50m", adj=c(0.5, 1.5), cex = 0.75)
    text(mean(rangeX), max(rangeY), paste("Day", timesPlot[j]),adj=c(0.5, 1.5), cex = 0.75 )
}
dev.off()

