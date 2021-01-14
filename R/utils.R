#' @export
LongLatToUTM<-function(x,y,zone){
    xy <- data.frame(ID = 1:length(x), X = x, Y = y)
    sp::coordinates(xy) <- c("X", "Y")
    proj4string(xy) <- CRS("+proj=longlat +datum=WGS84")  ## for example
    res <- spTransform(xy, CRS(paste("+proj=utm +zone=",zone," ellps=WGS84",sep='')))
    return(as.data.frame(res))
}

#' @export
par2sample <- function(par, type, n, weights){

    if(type == "mortality"){
        nExp = nrow(par)
        nGen = n
        boolVar = sample(size=nGen,x= nExp, prob = weights, replace = TRUE)
        idxMort = grepl(pattern = "mortality", par$event)
        tmp = rbeta(nGen, par$alpha[boolVar], par$beta[boolVar])
        sample_var = tmp
        sample_var[boolVar %in% which(idxMort)] <- -log(1-tmp[boolVar %in% which(idxMort)])
        sample_var[boolVar %in% which(!idxMort)] <- -log(tmp[boolVar %in% which(!idxMort)])
    }

    if(type == "diffusion"){
        nExp = nrow(par)
        nGen = n
        boolVar = sample(size=nGen,x= nExp, prob = weights, replace = TRUE)
        sample_var = rlnorm(nGen, par$mu[boolVar], par$sigma[boolVar]) # km / day
        # Get diffusion distribution from dispersal through
    }
return(sample_var)
}

myScalebar = function(units_label, yadj=1.5) {

    # Get plot coordinates
    pc = par("usr")

    # Position scale line between last two major x-axis tick marks
    # and 1/10th of the total y-range above the lower y-axis coordinate
    lines(c(floor(pc[2]-5),floor(pc[2])-4),
          rep(pc[3] + 0.1*(pc[4] - pc[3]), 2))

    # Place the units label at the midpoint of and just below the scale line
    text(x=mean(c(floor(pc[2]-5), floor(pc[2])-4)),
         y=pc[3] + 0.1*(pc[4] - pc[3]),
         label=units_label, adj=c(0.5, yadj))
}

#' @export
postPlot = function(prior,posterior = NULL, logF = FALSE, samples = TRUE, distPrior = NULL, distPosterior = NULL,
                    xlab = "", colA = "darkgreen", colB = "darkblue", shade = TRUE, vertL = NULL,
                    CI = c(0.05, 0.95), xRange = c(0.0001,0.99999),support = NULL){
    
    if(samples){
        TrFun = function(x){x}
        if(logF){TrFun = function(x){log(x)}}
        
        pPlot = TRUE
        if(is.null(posterior)){posterior = prior
        pPlot = FALSE}
        
        # if positive real support
        hdV = density(TrFun(prior), from = min(TrFun(prior)), to = max(TrFun(prior)), n = 2^11)
        hdV_Post = density(TrFun(posterior), from = min(TrFun(posterior)), to = max(TrFun(prior)), n = 2^11)
        # if bounded support, say (0, 1)
        
        if(!is.null(support)){
            hdV = density(TrFun(prior), from = support[1], to = support[2], n = 2^11)
            hdV_Post = density(TrFun(posterior), from = support[1], to = support[2], n = 2^11)
        }
        
        xPrior = hdV$x #hdV$mids
        xPosterior = hdV_Post$x #mids
        
        yPrior = hdV$y #density
        yPosterior = hdV_Post$y #density
        
        yLIM = range(c(yPrior, yPosterior), na.rm=T)
        xLIM = range(c(quantile(TrFun(prior),xRange), quantile(TrFun(posterior),xRange)))
        plot(xPrior, yPrior, col=colA, type="l",ylim = yLIM,
             xlab=xlab,ylab="",
             xlim = xLIM, axes = F, main = "")
        qPrior = quantile(TrFun(prior), CI)
        qPosterior = quantile(TrFun(posterior), CI)
        idxVar = sapply(qPrior, function(x) which.min((x - xPrior)^2))
        idxVar_Post = sapply(qPosterior, function(x) which.min((x - xPosterior)^2))
        
        
        if(shade){polygon(x = c(xPrior[idxVar[1]:idxVar[2]], rev(xPrior[idxVar[1]:idxVar[2]])), c(yPrior[idxVar[1]:idxVar[2]], rep(0, idxVar[2] - idxVar[1] + 1)), col=adjustcolor(colA, alpha.f = 0.1), border=NA)}
        
        if(pPlot){
            points(xPosterior, yPosterior, type="l", col = colB)
            if(shade){polygon(x = c(xPosterior[idxVar_Post[1]:idxVar_Post[2]], rev(xPosterior[idxVar_Post[1]:idxVar_Post[2]])), c(yPosterior[idxVar_Post[1]:idxVar_Post[2]], rep(0, idxVar_Post[2] - idxVar_Post[1] + 1)), col=adjustcolor(colB, alpha.f = 0.1), border=NA)}
        }
        axis(2)
        if(!is.null(vertL)){
            if(vertL == "mean"){
                idxMaxVar = which.min((xPrior-mean(TrFun(prior)))^2)
            }
            if(vertL == "mode"){idxMaxVar = which.max(hdV$density)}
            if(vertL == "median"){
                idxMaxVar = which.min((xPrior-quantile(TrFun(prior), 0.5))^2)
            }
            segments(xPrior[idxMaxVar], 0, xPrior[idxMaxVar], yPrior[idxMaxVar], col = colA)
            points(xPrior[idxMaxVar], yPrior[idxMaxVar],col=colA, pch=16)
        }
        if(pPlot){
            if(!is.null(vertL)){
                if(vertL == "mean"){
                    idxMaxVar = which.min((xPosterior-mean(TrFun(posterior)))^2)
                }
                if(vertL == "mode"){idxMaxVar = which.max(hdV_Post$density)}
                if(vertL == "median"){
                    idxMaxVar = which.min((xPosterior-quantile(TrFun(posterior), 0.5))^2)
                }
                segments(xPosterior[idxMaxVar], 0, xPosterior[idxMaxVar], yPosterior[idxMaxVar], col = colB)
                points(xPosterior[idxMaxVar], yPosterior[idxMaxVar],col=colB, pch=16)
            }
        }
        
    }else{
        
        
        xPrior = eval(parse(text = paste0("q",distPrior,"(seq(xRange[1],xRange[2],length.out=2000),",prior[1],",",prior[2],")")))
        xPosterior = NULL
        if(!is.null(posterior)){xPosterior = eval(parse(text = paste0("q",distPosterior,"(seq(xRange[1],xRange[2],length.out=2000),",posterior[1],",",posterior[2],")")))}
        
        yPrior = eval(parse(text = paste0("d",distPrior,"(xPrior,",prior[1],",",prior[2],")")))
        yPosterior = NULL
        if(!is.null(posterior)){yPosterior = eval(parse(text = paste0("d",distPosterior,"(xPosterior,",posterior[1],",",posterior[2],")")))}
        
        xPriorCI = eval(parse(text = paste0("q",distPrior,"(seq(CI[1],CI[2],length.out=2000),",prior[1],",",prior[2],")")))
        xPosteriorCI = NULL
        if(!is.null(posterior)){xPosteriorCI = eval(parse(text = paste0("q",distPosterior,"(seq(0.025,0.975,length.out=2000),",posterior[1],",",posterior[2],")")))}
        
        yPriorCI = eval(parse(text = paste0("d",distPrior,"(xPriorCI,",prior[1],",",prior[2],")")))
        yPosteriorCI = NULL
        if(!is.null(posterior)){yPosteriorCI = eval(parse(text = paste0("d",distPosterior,"(xPosteriorCI,",posterior[1],",",posterior[2],")")))}
        
        if(logF){
            xPrior = log(xPrior)
            if(!is.null(posterior)){xPosterior = log(xPosterior)}
            xPriorCI = log(xPriorCI)
            if(!is.null(posterior)){xPosteriorCI = log(xPosteriorCI)}
        }
        
        yLIM = range(c(yPrior, yPosterior), na.rm=T)
        xLIM = range(c(xPrior, xPosterior))
        plot(xPrior, yPrior, col=colA, type="l",ylim = yLIM,
             xlab=xlab,ylab="",
             xlim = xLIM, axes = F, main = "")
        
        
        
        polygon(x = c(xPriorCI, rev(xPriorCI)), c(yPriorCI, rep(0, length(yPriorCI))), col=adjustcolor(colA, alpha.f = 0.1), border=NA)
        
        if(!is.null(posterior)){
            points(xPosterior, yPosterior, type="l", col = colB)
            polygon(x = c(xPosteriorCI, rev(xPosteriorCI)), c(yPosteriorCI, rep(0, length(yPosteriorCI))), col=adjustcolor(colB, alpha.f = 0.1), border=NA)
        }
        axis(2)
        idxMaxVar = which.max(yPrior)
        segments(xPrior[idxMaxVar], 0, xPrior[idxMaxVar], yPrior[idxMaxVar], col = colA)
        points(xPrior[idxMaxVar], yPrior[idxMaxVar],col=colA, pch=16)
        
        if(!is.null(posterior)){
            idxMaxVar_Post = which.max(yPosterior)
            segments(xPosterior[idxMaxVar_Post], 0, xPosterior[idxMaxVar_Post], yPosterior[idxMaxVar_Post], col = colB)
            points(xPosterior[idxMaxVar_Post], yPosterior[idxMaxVar_Post],col=colB, pch=16)
        }
    }
    
    
}