##--- Run this script in order to set the MRR data ---##

library(rRiskDSMspread)
library(spreadR)
library(ggplot2)

# I - Load the data

data_mrr_path <- system.file("extdata", "mrr_data.csv", package="rRiskDSMspread")
data_mrr_c_path <- system.file("extdata", "mrr_compounds_data.csv", package="rRiskDSMspread")
data_mrr_s_path <- system.file("extdata", "mrr_swarms_data.csv", package="rRiskDSMspread")
data_mrr_l_path <- system.file("extdata", "mrr_location_data.csv", package="rRiskDSMspread")

MRR <- readr::read_csv(data_mrr_path)
MRR_compounds <- readr::read_csv(data_mrr_c_path)
MRR_swarm <- readr::read_csv(data_mrr_s_path)
LocMRR <- readr::read_csv(data_mrr_l_path)

# II - Define the pde grid using spatio temporal informaton from MRRs

MRR = dplyr::filter(MRR, .data$MRR != 5)
MRR$Date = as.Date(as.character(MRR$Date), format="%Y-%m-%d")
MRR$Date = as.character(MRR$Date)
MRR = MRR[-which(MRR$Date == "2013-09-18"),]
names(MRR)<- c("MRR", "Cap_Date","Long_CP", "Lat_CP", "Col_Meth","cVar", 
               "N_mos", "relSite")

DateMRR = unique(MRR[,c(1,2)])
DateMRR = DateMRR[ DateMRR$Cap_Date == ave(DateMRR$Cap_Date, DateMRR$MRR, FUN=min), ]
names(DateMRR)[2] <- c("Release date")
names(LocMRR) = c("MRR","Release site", "Longitude", "Latitude", "Number released")

dist_per_release <- full_join(LocMRR, MRR, by  = c("MRR", "Release site" = "relSite")) %>% 
    rename("relSite" = "Release site") %>%
    filter(NumMoz > 0)
dist_per_release$dist = geosphere::distHaversine(cbind(dist_per_release$Longitude, dist_per_release$Latitude), cbind(dist_per_release$XcoordGPS, dist_per_release$YcoordGPS))

res <- dist_per_release %>% group_by(MRR, relSite) %>% summarise(avg_dist = sum(dist*NumMoz) / sum(NumMoz))

distBoundary2Swarm = 400 # distance (in m)
xlimmy.orig <- range(MRR_swarm$X_m)
xlimmy <- xlimmy.orig + c( -1,1) * distBoundary2Swarm #diff( range( xlimmy.orig))
ylimmy.orig <- range(MRR_swarm$Y_m)
ylimmy <- ylimmy.orig + c(-1,1) * distBoundary2Swarm #diff( range( ylimmy.orig))
#define the spatial extent of the area with more dense grid cells
#NOTE THAT REAL APPLICATIONS SHOULD USE WIDER INNER GRIDS
distInnerBoundary2Swarm = 10 # distance inner (in m)
inner.xlimmy <- xlimmy.orig +c( -1,1) * distInnerBoundary2Swarm#diff( range( xlimmy.orig))
inner.ylimmy <- ylimmy.orig + c( -1,1) * distInnerBoundary2Swarm#diff( range( ylimmy.orig))
#define the grid around Bana Village (where swarms are located)
#note that there are more boundaries 
# in the Y direction as there are in the X -- an example
grid2D <- prepare.grid( covariates=list(MRR_swarm[,c("X_m","Y_m")]), 
                        N=c(100,150), xlim=xlimmy, ylim=ylimmy, 
                        inner.xlim=inner.xlimmy, inner.ylim=inner.ylimmy)

ggplot() + 
    geom_hline(aes(yintercept = grid2D$y.int), alpha = 0.25, col = "red", size = 0.5)+
    geom_vline(aes(xintercept = grid2D$x.int), alpha = 0.25, col = "red", size = 0.5)+
    geom_point(data = MRR_swarm, aes(x = X_m, y = Y_m), size = 1) + 
    geom_point(data = MRR_compounds, aes(x = X_m, y = Y_m), size = 1.75, shape = 0) + 
    theme_bw() + xlim(xlimmy.orig) + ylim(ylimmy.orig) + coord_fixed()

# III - Recode the data 

names(LocMRR)[2] <- "relSite"
tmpLoc = as.data.frame(unique(LocMRR[,c(2:4)]))
MRR_Post = merge(MRR, tmpLoc[-4,], by = "relSite", all.x = T)
MRR_Post = merge(MRR_Post, DateMRR, by = "MRR")
MRR_Post$DayDiff =  as.Date(MRR_Post$Cap_Date)- as.Date(MRR_Post$`Release date`)

MRR_Post$Col_Meth = as.factor(MRR_Post$Col_Meth)

idxCFR = which(MRR_Post$ColMeth == "Pot ")
idxPSC = which(MRR_Post$Col_Meth %in% "PSC")
idxSwarm = which(MRR_Post$Col_Meth %in% "Swarm")

MRR_Post$DayDiff[idxSwarm] = MRR_Post$DayDiff[idxSwarm] + 0.25
MRR_Post$DayDiff[idxCFR] = MRR_Post$DayDiff[idxCFR] + 0.5
MRR_Post$DayDiff[idxPSC] = MRR_Post$DayDiff[idxPSC] + 0.5

lastDay = as.numeric(max(MRR_Post$DayDiff))+1
times = seq(0,lastDay,0.25)

CP_Coord_UTM = LongLatToUTM(MRR_Post$Long_CP, MRR_Post$Lat_CP, 30)

colnames(CP_Coord_UTM) <- c("ID", "Cap_UTM_x", "Cap_UTM_y")

MRR_Post = cbind(MRR_Post, CP_Coord_UTM[,2:3])

locations = MRR_Post[,c("Cap_UTM_x", "Cap_UTM_y")]
IDs <- as.numeric( class::knn1( expand.grid( grid2D$x.mid, grid2D$y.mid), locations, 1:(grid2D$x.N*grid2D$y.N)))
MRR_Post$Loc_IDs = IDs

LocMRR = LongLatToUTM(tmpLoc$Longitude[1:3], tmpLoc$Latitude[1:3], 30)
LocMRR$ID = c("A", "B", "C")
names(LocMRR)[1] <- c("relSite")

releaseSites = LocMRR[c(rep(1:3,4)),2:3]

IDs_ReleaseSites <- as.numeric( class::knn1( expand.grid( grid2D$x.mid, grid2D$y.mid), releaseSites, 1:(grid2D$x.N*grid2D$y.N)))

releasedMoz = c(1146,1103,1158,
                1878,1734,1655,
                1665,1684,1673,
                1807,1653,1813
                #5992
)

releaseSitesCoded = NULL
for(rel in 1:12){
    mat = matrix(0, nrow = grid2D$x.N, ncol = grid2D$y.N)
    mat[IDs_ReleaseSites[rel]]<- releasedMoz[rel]
    releaseSitesCoded = rbind(releaseSitesCoded ,which(mat>0, arr.ind=T))
}
grid2D_release_sites <- releaseSitesCoded


mrr_data = plyr::ddply(MRR_Post, plyr::.(MRR, relSite,Col_Meth, DayDiff, Loc_IDs),
                 summarize, N_mos = sum(N_mos))

ggplot(mrr_data) + geom_point(aes(x = DayDiff, y = N_mos, col = Col_Meth)) + facet_grid(MRR ~ relSite)

# IV - Save the data to re-use in later scripts

saveRDS(grid2D, file = "inst/extdata/grid2D.rds")
saveRDS(grid2D_release_sites, file = "inst/extdata/grid2D_release_sites.rds")
saveRDS(mrr_data, file = "inst/extdata/mrr_data.rds")
