##--- Run this script in order to obtain prior distributions from elicited data ---##

library(rRiskDSMspread)
library(ElicitR)


data_elicited_path <- system.file("extdata", "elicitedDEID.csv", package="rRiskDSMspread")

elicitedDEID <- readr::read_csv(data_elicited_path)

# II - Get the prior distribution for the diffusion parameter of Anopheles Gambiae (Wild Type) male

eventDiffusion = c("FT1-6a")#,"FT1-6c") # Ag(WT) 6a, G3 6c (males in both cases)
elicitedDiff = elicitedDEID[intersect(
    which(elicitedDEID$basicEvent %in% eventDiffusion),
    which(!elicitedDEID$discard)),]
expertDiff = elicitedDEID$expertID[intersect(
    which(elicitedDEID$basicEvent %in% eventDiffusion),
    which(!elicitedDEID$discard))]
nExp = nrow(elicitedDiff)
par = NULL
for(i in 1:nExp){
    e.dataE <- elicit_data(
        cumul.probs = c(0.5-elicitedDiff$quant3[i]/2,0.5+elicitedDiff$quant3[i]/2), 
        quants = elicitedDiff[i,c('quant1','quant2')], family = 'lognormal')
    e.fitE <- dens_fit(e.dataE, crit = "SS")
    par = rbind(par, e.fitE$par)
}
diffusionParWT = as.data.frame(par)
diffusionParWT$expertID = as.character(expertDiff)
diffusionParWT$family = 'lognormal'
#saveRDS(diffusionPar, "../rds/diffusionParWT.rds")

dispersalVarWT <- par2sample(par = diffusionParWT, type = "diffusion", 
                             n = 1e5, weights = rep(1, nrow(diffusionParWT)))
diffusionVarWT = (dispersalVarWT*1000)^2 / (4 * 1) # m2 / day

# III - Get the prior distribution for the diffusion parameter of Anopheles Gambiae (DSM) male

eventDiffusion = c("FT1-6c") # Ag(WT) 6a, G3 6c (males in both cases)
elicitedDiff = elicitedDEID[intersect(
    which(elicitedDEID$basicEvent %in% eventDiffusion),
    which(!elicitedDEID$discard)),]
nExp = nrow(elicitedDiff)
par = NULL
for(i in 1:nExp){
    e.dataE <- elicit_data(
        cumul.probs = c(0.5-elicitedDiff$quant3[i]/2,0.5+elicitedDiff$quant3[i]/2), 
        quants = elicitedDiff[i,c('quant1','quant2')], family = 'lognormal')
    e.fitE <- dens_fit(e.dataE, crit = "SS")
    par = rbind(par, e.fitE$par)
}
diffusionParSM = as.data.frame(par)
diffusionParSM$expertID = as.character(expertDiff)
diffusionParSM$family = 'lognormal'
#saveRDS(diffusionPar, "../rds/diffusionParGE.rds")
#Generate dispersal:
dispersalVarSM <- par2sample(par = diffusionParSM, type = "diffusion", 
                             n = 1e5, weights = rep(1, nrow(diffusionParSM)))
diffusionVarSM = (dispersalVarSM*1000)^2 / (4 * 1) # m2 / day

# IV - Plot the resulting prior distributions

postPlot(diffusionVarWT,diffusionVarSM, logF=  TRUE, xlab = expression(delta), xRange = c(1e-3,1-1e-3))
axis(1, at = log(c(1 %o% 10^(-3:8))), labels = c(1 %o% 10^(-3:8)))
legend("topleft", lty=c(1, 1), lwd=c(1, 1),col=c("darkgreen", "darkblue"), 
       legend = c("Linear pool prior for Ag(WT) dispersal", "Linear pool prior for Ag(DSM) dispersal"), 
       border = NA, box.lwd = 0, cex = 0.75, bty = "n")

# V - Get the prior distribution for the survival parameter of Anopheles Gambiae (Wild Type) male

eventMortRate = c("FT1-4a-mortality", "FT1-4a-survival") # Wild type (male)
#  eventMortRate = c("FT1-4c-mortality", "FT1-4c-survival") # G3 strain (male)
elicitedMR = elicitedDEID[intersect(
    which(elicitedDEID$basicEvent %in% eventMortRate),
    which(!elicitedDEID$discard)),]
nExp = nrow(elicitedMR)
nameExpertsMort = as.character(elicitedMR$expertID)
par = NULL
for(i in 1:nrow(elicitedMR)){
    e.dataE <- elicit_data(
        cumul.probs = c(0.5-elicitedMR$quant3[i]/2,0.5+elicitedMR$quant3[i]/2), 
        quants = elicitedMR[i,c('quant1','quant2')], family = elicitedMR$family[i])
    e.fitE <- dens_fit(e.dataE, crit = "SS")
    par = rbind(par, c(e.fitE$par, elicitedMR$family[i]))
}
mortalityParWT = data.frame(alpha = as.numeric(as.character(par[,'alpha'])), beta = as.numeric(as.character(par[,'beta'])))
mortalityParWT$family = elicitedMR$family[1:nrow(elicitedMR)]
mortalityParWT$event = elicitedMR$basicEvent[1:nrow(elicitedMR)]
sampMortRateWT <- par2sample(par = mortalityParWT, type = "mortality", 
                             n = 1e5, weights = rep(1, nrow(mortalityParWT)))

# VI - Get the prior distribution for the survival parameter of Anopheles Gambiae (DSM) 

eventMortRate = c("FT1-4c-mortality", "FT1-4c-survival") # G3 strain (male)
elicitedMR = elicitedDEID[intersect(
    which(elicitedDEID$basicEvent %in% eventMortRate),
    which(!elicitedDEID$discard)),]
nExp = nrow(elicitedMR)
nameExpertsMort = as.character(elicitedMR$expertID)
par = NULL
for(i in 1:nrow(elicitedMR)){
    e.dataE <- elicit_data(
        cumul.probs = c(0.5-elicitedMR$quant3[i]/2,0.5+elicitedMR$quant3[i]/2), 
        quants = elicitedMR[i,c('quant1','quant2')], family = elicitedMR$family[i])
    e.fitE <- dens_fit(e.dataE, crit = "SS")
    par = rbind(par, c(e.fitE$par, elicitedMR$family[i]))
}
mortalityParSM = data.frame(alpha = as.numeric(as.character(par[,'alpha'])), beta = as.numeric(as.character(par[,'beta'])))
mortalityParSM$family = elicitedMR$family[1:nrow(elicitedMR)]
mortalityParSM$event = elicitedMR$basicEvent[1:nrow(elicitedMR)]
#saveRDS(mortalityPar, file = "../rds/mortalityParGE.rds")
#Generate dispersal:
sampMortRateSM <- par2sample(par = mortalityParSM, type = "mortality", 
                             n = 1e5, weights = rep(1, nrow(mortalityParSM)))

# VII - Plot the resulting prior distributions

postPlot(sampMortRateWT,sampMortRateSM, logF= TRUE, xlab = expression(mu), xRange = c(1e-3, 1-1e-15))
axis(1, at = log(c(seq(0,0.2,0.05),seq(0.2,2,0.2))), labels = c(seq(0,0.2,0.05), seq(0.2,2,0.2)))
legend("topright", lty=c(1, 1), lwd=c(1, 1),col=c("darkgreen", "darkblue"), 
       legend = c("Linear pool prior for Ag(WT) mortality", "Linear pool prior for Ag(DSM) mortality"), 
       border = NA, box.lwd = 0, cex = 0.75, bty = "n")

# VIII - Save the data to re-use in later scripts

saveRDS(mortalityParSM, file = "inst/extdata/mortalityParDSM.rds")
saveRDS(mortalityParWT, file = "inst/extdata/mortalityParWT.rds")
saveRDS(diffusionParSM, file = "inst/extdata/diffusionParDSM.rds")
saveRDS(diffusionParWT, file = "inst/extdata/diffusionParWT.rds")

# IX - Plotting all priors together

swarmAvarSM = rlnorm(1e4, -2, 1)
swarmSvarSM = rlnorm(1e4, 3.5, 4)

par(mfrow=c(2,2))
postPlot(diffusionVarWT, logF = TRUE, xlab = expression(paste("D - [",paste("m"^"2","/day"), "]")), vertL = "mean")
axis(1, at = log(c(1 %o% 10^(-3:8))), labels = c(1 %o% 10^(-3:8)))
postPlot(sampMortRateWT, logF = TRUE, xlab = expression(paste(mu," - [",paste("Day"^"-1"), "]")), vertL = "mean")
axis(1, at = log(c(seq(0,0.2,0.05),seq(0.2,2,0.2))), labels = c(seq(0,0.2,0.05), seq(0.2,2,0.2)))
postPlot(swarmAvarSM, logF = TRUE, xlab = expression(paste(alpha," - [",paste("m"^"-1"),paste("R"^"-1"), "]")), vertL = "mean")
axis(1, at = log(c(1 %o% 10^(-3:8))), labels = c(1 %o% 10^(-3:8)))
postPlot(swarmSvarSM, logF = TRUE, xlab = expression(paste(sigma," - [",paste("m"), "]")), vertL = "mean")
axis(1, at = log(c(1 %o% 10^(-3:8))),  labels = c(1 %o% 10^(-3:8)))