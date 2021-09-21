##--- Run this script in order to retrieve prior distributions from elicited data ---##

diffusion_distWT <- system.file("extdata", "diffusionParWT.rds", package="rRiskDSMspread")
diffusionParWT <- readRDS(diffusion_distWT)

dispersalVarWT <- par2sample(par = diffusionParWT, type = "diffusion", 
                             n = 1e5, weights = rep(1, nrow(diffusionParWT)))
diffusionVarWT = (dispersalVarWT*1000)^2 / (4 * 1) # m2 / day

# III - Get the prior distribution for the diffusion parameter of Anopheles Gambiae (DSM) male

diffusion_distSM <- system.file("extdata", "diffusionParSM.rds", package="rRiskDSMspread")
diffusionParSM <- readRDS(diffusion_distSM)

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

mortality_distWT <- system.file("extdata", "mortalityParWT.rds", package="rRiskDSMspread")
mortalityParWT <- readRDS(mortality_distWT)

sampMortRateWT <- par2sample(par = mortalityParWT, type = "mortality", 
                             n = 1e5, weights = rep(1, nrow(mortalityParWT)))

# VI - Get the prior distribution for the survival parameter of Anopheles Gambiae (DSM) 

mortality_distSM <- system.file("extdata", "mortalityParSM.rds", package="rRiskDSMspread")
mortalityParSM <- readRDS(mortality_distSM)

#Generate mortality rate:
sampMortRateSM <- par2sample(par = mortalityParSM, type = "mortality", 
                             n = 1e5, weights = rep(1, nrow(mortalityParSM)))

# VII - Plot the resulting prior distributions

postPlot(sampMortRateWT,sampMortRateSM, logF= TRUE, xlab = expression(mu), xRange = c(1e-3, 1-1e-15))
axis(1, at = log(c(seq(0,0.2,0.05),seq(0.2,2,0.2))), labels = c(seq(0,0.2,0.05), seq(0.2,2,0.2)))
legend("topright", lty=c(1, 1), lwd=c(1, 1),col=c("darkgreen", "darkblue"), 
       legend = c("Linear pool prior for Ag(WT) mortality", "Linear pool prior for Ag(DSM) mortality"), 
       border = NA, box.lwd = 0, cex = 0.75, bty = "n")


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
