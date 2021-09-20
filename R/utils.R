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


elicitedColorScale <- RColorBrewer::brewer.pal(5,"Dark2")

colPrior = "darkgreen"
colPosterior = "darkblue"
colPred = "darkorange"
colTraj = "black"

#' @export
cdfPlot = function(prior,posterior = NULL, logF = FALSE, samples = TRUE, distPrior = NULL, distPosterior = NULL,
                   xlab = "", colA = "darkgreen", colB = "darkblue",
                   CI = c(0.025, 0.975), xRange = c(0.0001,0.99999),yLIM = NULL,shade = TRUE,
                   quantBreaks = NULL, regBreaks = NULL, smoothCurve = F, spanLoess = 0.5){
    
    if(samples){
        TrFun = function(x){x}
        if(logF){TrFun = function(x){log(x)}}
        
        pPlot = TRUE
        if(is.null(posterior)){posterior = prior
        pPlot = FALSE}
        if(is.null(quantBreaks)){
            qX = unique(c(seq(0,0.01,0.005),seq(0.01,0.05,0.005),seq(0.25,0.75,0.005),seq(0.95,0.99,0.005), seq(0.99,1,0.005)))}else{
                qX = quantBreaks
            }
        
        qVar = quantile(TrFun(prior), qX)
        qVar_Post = quantile(TrFun(posterior), qX)
        if(!is.null(regBreaks)){
            if(length(regBreaks) == 1){
                qVar = seq(min(TrFun(prior)), max(TrFun(prior)), length.out = regBreaks)
                qVar_Post = seq(min(TrFun(posterior)), max(TrFun(posterior)), length.out = regBreaks)}
            if(length(regBreaks) > 1){
                qVar = regBreaks
                qVar_Post = regBreaks}
        }
        
        qVarA = sort(unique(c(qVar, qVar_Post)))
        
        yPrior = sapply(qVarA, function(x) mean(TrFun(prior) < x))
        yPosterior = sapply(qVarA, function(x) mean(TrFun(posterior) < x))
        
        #hdV = hist(TrFun(prior), breaks = qVar, plot=F)
        #hdV_Post = hist(TrFun(posterior), breaks = qVar_Post, plot=F)
        
        xPrior = qVarA #hdV$mids
        xPosterior = qVarA#hdV_Post$mids
        
        if(smoothCurve){
            yPrior = predict(loess(yPrior~xPrior, span = spanLoess))
            yPosterior = predict(loess(yPosterior~xPosterior, span = spanLoess))
        }
        
        
        
        if(is.null(yLIM)){yLIM = range(c(yPrior, yPosterior), na.rm=T)}
        xLIM = range(c(quantile(TrFun(prior),xRange), quantile(TrFun(posterior),xRange)))
        plot(xPrior, yPrior, col=colA, type="l",ylim = yLIM,
             xlab=xlab,ylab="",
             xlim = xLIM, axes = F, main = "")
        
        qPrior = quantile(TrFun(prior), CI)
        qPosterior = quantile(TrFun(posterior), CI)
        idxVar = sapply(qPrior, function(x) which.min((x - qVarA)^2))
        idxVar_Post = sapply(qPosterior, function(x) which.min((x - qVarA)^2))
        
        if(shade){polygon(x = c(xPrior[idxVar[1]:idxVar[2]], rev(xPrior[idxVar[1]:idxVar[2]])), c(yPrior[idxVar[1]:idxVar[2]], rep(0, idxVar[2] - idxVar[1] + 1)), col=adjustcolor(colA, alpha.f = 0.1), border=NA)}
        
        if(pPlot){
            points(xPosterior, yPosterior, type="l", col = colB)
            if(shade){polygon(x = c(xPosterior[idxVar_Post[1]:idxVar_Post[2]], rev(xPosterior[idxVar_Post[1]:idxVar_Post[2]])), c(yPosterior[idxVar_Post[1]:idxVar_Post[2]], rep(0, idxVar_Post[2] - idxVar_Post[1] + 1)), col=adjustcolor(colB, alpha.f = 0.1), border=NA)}
        }
        axis(2)
        idxMaxVar = which.min((yPrior - 0.5)^2)
        segments(xPrior[idxMaxVar], 0, xPrior[idxMaxVar], yPrior[idxMaxVar], col = colA)
        points(xPrior[idxMaxVar], yPrior[idxMaxVar],col=colA, pch=16)
        #abline(v = mean(prior), lty= 2, col = colA)
        
        if(pPlot){
            idxMaxVar_Post =  which.min((yPosterior - 0.5)^2)
            segments(xPosterior[idxMaxVar_Post], 0, xPosterior[idxMaxVar_Post], yPosterior[idxMaxVar_Post], col = colB)
            points(xPosterior[idxMaxVar_Post], yPosterior[idxMaxVar_Post],col=colB, pch=16)
            #abline(v = mean(posterior), lty= 2, col = colB)
        }
        
    }else{
        
        
        xPrior = eval(parse(text = paste0("q",distPrior,"(seq(xRange[1],xRange[2],length.out=2000),",prior[1],",",prior[2],")")))
        xPosterior = NULL
        if(!is.null(posterior)){xPosterior = eval(parse(text = paste0("q",distPosterior,"(seq(xRange[1],xRange[2],length.out=2000),",posterior[1],",",posterior[2],")")))}
        
        yPrior = eval(parse(text = paste0("p",distPrior,"(xPrior,",prior[1],",",prior[2],")")))
        yPosterior = NULL
        if(!is.null(posterior)){yPosterior = eval(parse(text = paste0("p",distPosterior,"(xPosterior,",posterior[1],",",posterior[2],")")))}
        
        xPriorCI = eval(parse(text = paste0("q",distPrior,"(seq(0.025,0.957,length.out=2000),",prior[1],",",prior[2],")")))
        xPosteriorCI = NULL
        if(!is.null(posterior)){xPosteriorCI = eval(parse(text = paste0("q",distPosterior,"(seq(0.025,0.975,length.out=2000),",posterior[1],",",posterior[2],")")))}
        
        yPriorCI = eval(parse(text = paste0("p",distPrior,"(xPriorCI,",prior[1],",",prior[2],")")))
        yPosteriorCI = NULL
        if(!is.null(posterior)){yPosteriorCI = eval(parse(text = paste0("p",distPosterior,"(xPosteriorCI,",posterior[1],",",posterior[2],")")))}
        
        if(logF){
            xPrior = log(xPrior)
            xPosterior = log(xPosterior)
            xPriorCI = log(xPriorCI)
            xPosteriorCI = log(xPosteriorCI)
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
        idxMaxVar = which.min((yPrior-0.5)^2)
        segments(xPrior[idxMaxVar], 0, xPrior[idxMaxVar], yPrior[idxMaxVar], col = colA)
        points(xPrior[idxMaxVar], yPrior[idxMaxVar],col=colA, pch=16)
        
        if(!is.null(posterior)){
            idxMaxVar_Post = which.min((yPosterior-0.5)^2)
            segments(xPosterior[idxMaxVar_Post], 0, xPosterior[idxMaxVar_Post], yPosterior[idxMaxVar_Post], col = colB)
            points(xPosterior[idxMaxVar_Post], yPosterior[idxMaxVar_Post],col=colB, pch=16)
        }
    }
    
    
}

#' @export
corMargPlot = function(dataF, nbins = 20, likeValue = NULL,
                       ncol = 5, marg = c(3,2,2,1),
                       nContour = 100,cexP = 0.5,
                       colPal = rev(c("#084594", "#2171B5", "#4292C6", "#6BAED6", "#9ECAE1", "#C6DBEF", "#DEEBF7", "#F7FBFF", "white")),mainNames = NULL, h = NULL) {
    
    nrowP = ncolP = ncol(dataF)
    par(mfrow = c(nrowP, ncolP), mar = marg)
    if(is.null(mainNames)){mainNames = names(dataF)}
    for(i in 1:nrowP){
        for(j in 1:ncolP){
            if(i == j){
                delta = diff(range(dataF[,i]))/10
                hist(dataF[,i],
                     breaks = seq(min(dataF[,i])-delta, max(dataF[,i])+delta, length.out = nbins),
                     main = mainNames[i])
            }
            if(i < j){
                if(!is.null(likeValue)){
                    ii <- cut(likeValue, breaks = seq(min(likeValue), max(likeValue), len = ncol),
                              include.lowest = TRUE)
                    colP = colorRampPalette(colPal)(ncol-1)[ii]
                }else{
                    colP = "gray"
                }
                plot(dataF[,j], dataF[,i], col = colP, cex = cexP, pch = 16)
            }
            if(j < i){
                if(is.null(h)){hP = c(diff(range(dataF[,j]))/10, diff(range(dataF[,i]))/10)}else{hP = c(h[j], h[i])}
                contourData= kde2d(dataF[,j], dataF[,i], n = nContour, h = hP)
                image(contourData$x, contourData$y, contourData$z, col = colorRampPalette(colPal)(ncol))
                #      contour(contourData$x, contourData$y, contourData$z, col = colorRampPalette(colPal)(ncol), nlevels = ncol-1,add=T)
            }
        }
    }
}




#' @export
elicit_data = function (e.data = NULL, expert = "expert", theta = "theta", 
                        cumul.probs = NULL, quants = NULL, method = NULL, family = NULL, 
                        support = NULL, df = NULL) 
{
    if (is.null(e.data)) {
        e.vec <- elicit_vec(cumul.probs, quants)
        e.data <- list(expert = expert, theta = theta, data = e.vec, 
                       method = method, family = family)
    }
    else {
        if (!is.null(cumul.probs)) {
            e.probs <- c(sapply(strsplit(names(e.data$data), 
                                         split = "_"), function(x) as.numeric(x[2])), 
                         cumul.probs)
            e.quants <- c(e.data$data, quants)
            e.data$data <- elicit_vec(e.probs, e.quants)
            e.data$par <- e.data$obj.value <- e.data$obj.method <- e.data$fitted.quantiles <- e.data$quant.dev <- NULL
        }
        if (!is.null(method)) 
            e.data$method <- method
        if (!is.null(e.data$family) && is.null(support) && (is.null(family) || 
                                                            identical(e.data$family, family))) {
            support <- e.data$support
        }
        if (!is.null(family) && !identical(e.data$family, family)) {
            e.data$family <- family
            e.data$par <- e.data$obj.value <- e.data$obj.method <- e.data$fitted.quantiles <- e.data$quant.dev <- NULL
        }
    }
    if (is.null(e.data$family)) {
        support <- NULL
    }
    else {
        if (is.null(support)) {
            if (e.data$family == "beta") {
                if ((any(e.data$data < 0)) || (any(e.data$data > 
                                                   1))) 
                    stop("support is [0, 1]")
                support <- c(0, 1)
            }
            if ((e.data$family == "lognormal") || (e.data$family == 
                                                   "logst") || (e.data$family == "gamma")) {
                if (any(e.data$data <= 0)) 
                    stop("support is positive reals")
                support <- c(0, Inf)
            }
            if ((e.data$family == "normal") || (e.data$family == 
                                                "st")) 
                support <- c(-Inf, Inf)
        }
        if ((any(e.data$data < support[1])) || (any(e.data$data > 
                                                    support[2]))) 
            stop("at least one elicited quantile is outside restricted support")
        if (support[1] >= support[2]) 
            stop("check support a < b")
        if (e.data$family == "beta") {
            if (!is.finite(support[1]) || !is.finite(support[2])) 
                stop("finite support required for scaled beta")
        }
        if ((e.data$family == "lognormal") || (e.data$family == 
                                               "logst") || (e.data$family == "gamma")) {
            if ((support[1] != 0) || (support[2] != Inf)) 
                warning("modifed support for density with support on positive reals")
            support <- c(0, Inf)
        }
    }
    e.data$support <- support
    if (!is.null(e.data) && !is.null(e.data$family)) {
        if (!is.null(df) && ((e.data$family != "st") && (e.data$family != 
                                                         "logst"))) 
            stop("df argument only needed for (log) Students t")
        if ((e.data$family == "st") || (e.data$family == "logst")) {
            if (is.null(e.data$df)) {
                if (is.null(df)) {
                    stop("df must be specified for (log) Students t")
                }
                else {
                    e.data$df <- df
                }
            }
        }
    }
    e.data
}


#' @export
elicit_vec = function (cumul.probs, quantiles) 
{
    e <- as.numeric(quantiles[order(cumul.probs)])
    if (any(cumul.probs < 0)) 
        stop("cumulative probabilities must be nonnegative")
    if (any(cumul.probs > 1)) 
        stop("cumulative probabilities must be <= 1")
    if (length(cumul.probs) > length(unique(cumul.probs))) 
        stop("non-unique cumulative probabilities specified")
    if (any(diff(e) < 0)) 
        stop("impossible quantile specification")
    if (length(cumul.probs) != length(quantiles)) 
        stop("unequal number of quantiles and cumulative probs")
    names(e) <- paste("q_", sort(cumul.probs), sep = "")
    e
}

#' @export
dens_fit = function(e.data, crit = "KL", pars = NULL, central.probs = NULL){
    q.elicit <- e.data$data
    dens <- e.data$family
    support <- e.data$support
    if (!is.null(e.data$df) && ((dens != "st") && (dens != "logst"))) 
        stop("df argument only needed for (log) Students t")
    if ((dens == "st") || (dens == "logst")) 
        df <- e.data$df
    cumul.probs <- sapply(strsplit(names(q.elicit), split = "_"), 
                          function(x) as.numeric(x[2]))
    cumul.probs.elicit <- c(0, cumul.probs, 1)
    prob.intervals.elicit <- diff(cumul.probs.elicit)
    if (dens == "beta") {
        if (is.null(pars)) {
            if (is.null(central.probs)) {
                pars <- c(1, 1)
            }
            else {
                if ((!(any(grepl(central.probs[1], names(q.elicit))))) || 
                    (!(any(grepl(central.probs[2], names(q.elicit)))))) 
                    stop("central.probs: no matching args")
                central.probs <- sort(central.probs)
                cen.qs <- q.elicit[paste("q_", central.probs, 
                                         sep = "")]
                z.q <- qnorm(central.probs[2]) - qnorm(central.probs[1])
                sd.n <- (cen.qs[2] - cen.qs[1])/z.q
                if (any(grepl("0\\.5$", names(q.elicit)))) {
                    mu.n <- q.elicit["q_0.5"]
                }
                else {
                    mu.n <- sum(cen.qs)/2
                }
                m <- (mu.n - support[1])/(support[2] - support[1])
                v <- sd.n^2/(support[2] - support[1])^2
                if (v < m * (1 - m)) {
                    alpha <- m * (m * (1 - m)/v - 1)
                    beta <- (1 - m) * (m * (1 - m)/v - 1)
                    pars <- c(alpha, beta)
                }
                else {
                    pars <- c(1, 1)
                }
            }
        }
        pars <- log(pars)
        names(pars) <- c("alpha", "beta")
    }
    if ((dens == "lognormal")) {
        if (is.null(pars)) {
            if (is.null(central.probs)) {
                pars <- c(0, 1)
            }
            else {
                if ((!(any(grepl(central.probs[1], names(q.elicit))))) || 
                    (!(any(grepl(central.probs[2], names(q.elicit)))))) 
                    stop("central.probs: no matching args")
                central.probs <- sort(central.probs)
                cen.qs <- log(q.elicit[paste("q_", central.probs, 
                                             sep = "")])
                z.q <- qnorm(central.probs[2]) - qnorm(central.probs[1])
                sd.n <- (cen.qs[2] - cen.qs[1])/z.q
                if (any(grepl("0\\.5$", names(q.elicit)))) {
                    mu.n <- log(q.elicit["q_0.5"])
                }
                else {
                    mu.n <- sum(cen.qs)/2
                }
                pars <- c(mu.n, sd.n)
            }
        }
        pars[2] <- log(pars[2])
        names(pars) <- c("mu", "sigma")
    }
    if ((dens == "normal")) {
        if (is.null(pars)) {
            if (is.null(central.probs)) {
                pars <- c(0, 1)
            }
            else {
                if ((!(any(grepl(central.probs[1], names(q.elicit))))) || 
                    (!(any(grepl(central.probs[2], names(q.elicit)))))) 
                    stop("central.probs: no matching args")
                central.probs <- sort(central.probs)
                cen.qs <- q.elicit[paste("q_", central.probs, 
                                         sep = "")]
                z.q <- qnorm(central.probs[2]) - qnorm(central.probs[1])
                sd.n <- (cen.qs[2] - cen.qs[1])/z.q
                if (any(grepl("0\\.5$", names(q.elicit)))) {
                    mu.n <- q.elicit["q_0.5"]
                }
                else {
                    mu.n <- sum(cen.qs)/2
                }
                pars <- c(mu.n, sd.n)
            }
        }
        pars[2] <- log(pars[2])
        names(pars) <- c("mu", "sigma")
    }
    if (dens == "st") {
        if (is.null(df)) 
            stop("must specify degrees of freedom")
        if (is.null(pars) || any(is.na(pars[1:2]))) {
            if (is.null(central.probs)) {
                pars <- c(0, 1)
            }
            else {
                if ((!(any(grepl(central.probs[1], names(q.elicit))))) || 
                    (!(any(grepl(central.probs[2], names(q.elicit)))))) 
                    stop("central.probs: no matching args")
                central.probs <- sort(central.probs)
                cen.qs <- q.elicit[paste("q_", central.probs, 
                                         sep = "")]
                if (any(grepl("0\\.5$", names(q.elicit)))) {
                    mu.n <- q.elicit["q_0.5"]
                }
                else {
                    mu.n <- sum(cen.qs)/2
                }
                tz.q <- qt(central.probs[2], df = df) - qt(central.probs[1], 
                                                           df = df)
                sd.n <- (cen.qs[2] - cen.qs[1])/tz.q
                pars <- c(mu.n, sd.n)
            }
        }
        pars[2] <- log(pars[2])
        names(pars) <- c("mu", "sigma")
    }
    if (dens == "logst") {
        if (is.null(df)) 
            stop("must specify degrees of freedom")
        if (is.null(pars) || any(is.na(pars[1:2]))) {
            if (is.null(central.probs)) {
                pars <- c(0, 1)
            }
            else {
                if ((!(any(grepl(central.probs[1], names(q.elicit))))) || 
                    (!(any(grepl(central.probs[2], names(q.elicit)))))) 
                    stop("central.probs: no matching args")
                central.probs <- sort(central.probs)
                cen.qs <- log(q.elicit[paste("q_", central.probs, 
                                             sep = "")])
                if (any(grepl("0\\.5$", names(q.elicit)))) {
                    mu.n <- log(q.elicit["q_0.5"])
                }
                else {
                    mu.n <- sum(cen.qs)/2
                }
                tz.q <- qt(central.probs[2], df = df) - qt(central.probs[1], 
                                                           df = df)
                sd.n <- (cen.qs[2] - cen.qs[1])/tz.q
                pars <- c(mu.n, sd.n)
            }
        }
        pars[2] <- log(pars[2])
        names(pars) <- c("mu", "sigma")
    }
    if (dens == "gamma") {
        if (is.null(pars)) {
            if (is.null(central.probs)) {
                pars <- c(1, 1)
            }
            else {
                if ((!(any(grepl(central.probs[1], names(q.elicit))))) || 
                    (!(any(grepl(central.probs[2], names(q.elicit)))))) 
                    stop("central.probs: no matching args")
                central.probs <- sort(central.probs)
                cen.qs <- q.elicit[paste("q_", central.probs, 
                                         sep = "")]
                z.q <- qnorm(central.probs[2]) - qnorm(central.probs[1])
                sd.n <- (cen.qs[2] - cen.qs[1])/z.q
                if (any(grepl("0\\.5$", names(q.elicit)))) {
                    mu.n <- q.elicit["q_0.5"]
                }
                else {
                    mu.n <- sum(cen.qs)/2
                }
                pars <- c((mu.n/sd.n)^2, sd.n^2/mu.n)
            }
        }
        pars <- log(pars)
        names(pars) <- c("alpha", "beta")
    }
    if (dens == "beta") {
        prob_fun <- function(x, qelicit) {
            xx <- exp(x)
            c(0, pscaled_beta(qelicit, alpha = xx[1], beta = xx[2], 
                              a = support[1], b = support[2]), 1)
        }
    }
    if ((dens == "lognormal")) {
        prob_fun <- function(x, qelicit) {
            xx <- c(x[1], exp(x[2]))
            c(0, pnorm(log(qelicit), xx[1], xx[2]), 1)
        }
    }
    if (dens == "normal") {
        prob_fun <- function(x, qelicit) {
            xx <- c(x[1], exp(x[2]))
            c(0, ptruncnorm(qelicit, a = support[1], b = support[2], 
                            mean = xx[1], sd = xx[2]), 1)
        }
    }
    if ((dens == "st")) {
        prob_fun <- function(x, qelicit) {
            xx <- c(x[1], exp(x[2]))
            c(0, pt((qelicit - xx[1])/xx[2], df = df), 1)
        }
    }
    if ((dens == "logst")) {
        prob_fun <- function(x, qelicit) {
            xx <- c(x[1], exp(x[2]))
            c(0, pt((log(qelicit) - xx[1])/xx[2], df = df), 1)
        }
    }
    if ((dens == "gamma")) {
        prob_fun <- function(x, qelicit) {
            xx <- c(exp(x[1]), exp(x[2]))
            c(0, pgamma(qelicit, shape = xx[1], scale = xx[2]), 
              1)
        }
    }
    if (crit == "KL") {
        obj.fun <- function(x, q.elicit, p.elicit) {
            cumul.prob.fit <- prob_fun(x, q.elicit)
            prob.intervals.fit <- diff(cumul.prob.fit)
            sum((log(prob.intervals.elicit) - log(prob.intervals.fit)) * 
                    prob.intervals.elicit)
        }
    }
    if (crit == "SS") {
        obj.fun <- function(x, q.elicit, p.elicit) {
            cumul.prob.fit <- prob_fun(x, q.elicit)
            prob.intervals.fit <- diff(cumul.prob.fit)
            sum((prob.intervals.elicit - prob.intervals.fit)^2)
        }
    }
    out <- optim(pars, function(x) obj.fun(x, q.elicit, p.elicit), 
                 method = "Nelder-Mead")
    if (dens == "beta") {
        par <- exp(out$par)
        quant.fits <- qscaled_beta(cumul.probs, alpha = par[1], 
                                   beta = par[2], a = support[1], b = support[2])
        names(quant.fits) <- names(q.elicit)
        d_fun <- function(x) dscaled_beta(x, par[1], par[2], 
                                          support[1], support[2])
        p_fun <- function(x) pscaled_beta(x, par[1], par[2], 
                                          support[1], support[2])
        q_fun <- function(x) qscaled_beta(x, par[1], par[2], 
                                          support[1], support[2])
    }
    if (dens == "lognormal") {
        par <- c(out$par[1], exp(out$par[2]))
        quant.fits <- qlnorm(cumul.probs, par[1], par[2])
        names(quant.fits) <- names(q.elicit)
        d_fun <- function(x) dlnorm(x, par[1], par[2])
        p_fun <- function(x) plnorm(x, par[1], par[2])
        q_fun <- function(x) qlnorm(x, par[1], par[2])
    }
    if (dens == "normal") {
        par <- c(out$par[1], exp(out$par[2]))
        quant.fits <- qtruncnorm(cumul.probs, a = support[1], 
                                 b = support[2], mean = par[1], sd = par[2])
        names(quant.fits) <- names(q.elicit)
        d_fun <- function(x) dtruncnorm(x, a = support[1], b = support[2], 
                                        mean = par[1], sd = par[2])
        p_fun <- function(x) ptruncnorm(x, a = support[1], b = support[2], 
                                        mean = par[1], sd = par[2])
        q_fun <- function(x) qtruncnorm(x, a = support[1], b = support[2], 
                                        mean = par[1], sd = par[2])
    }
    if (dens == "st") {
        par <- c(out$par[1], exp(out$par[2]))
        names(par) <- c("mu", "sigma")
        quant.fits <- qt(cumul.probs, df = df) * par["sigma"] + 
            par["mu"]
        names(quant.fits) <- names(q.elicit)
        d_fun <- function(x) dt((x - par["mu"])/par["sigma"], 
                                df = df)/par["sigma"]
        p_fun <- function(x) pt((x - par["mu"])/par["sigma"], 
                                df = df)
        q_fun <- function(x) qt(x, df = df) * par["sigma"] + 
            par["mu"]
    }
    if (dens == "logst") {
        par <- c(out$par[1], exp(out$par[2]))
        names(par) <- c("mu", "sigma")
        quant.fits <- exp(qt(cumul.probs, df = df) * par["sigma"] + 
                              par["mu"])
        names(quant.fits) <- names(q.elicit)
        d_fun <- function(x) dt((log(x) - par["mu"])/par["sigma"], 
                                df = df)/(x * par["sigma"])
        p_fun <- function(x) pt((log(x) - par["mu"])/par["sigma"], 
                                df = f)
        q_fun <- function(x) exp(qt(log(x), df = df) * par["sigma"] + 
                                     par["mu"])
    }
    if (dens == "gamma") {
        par <- exp(out$par)
        quant.fits <- qgamma(cumul.probs, shape = par[1], scale = par[2])
        names(quant.fits) <- names(q.elicit)
        d_fun <- function(x) dgamma(x, shape = par[1], scale = par[2])
        p_fun <- function(x) pgamma(x, shape = par[1], scale = par[2])
        q_fun <- function(x) qgamma(x, shape = par[1], scale = par[2])
    }
    e.data$par <- par
    e.data$obj.value <- out$value
    e.data$obj.method <- crit
    e.data$fitted.quantiles <- quant.fits
    e.data$quant.dev <- quant.fits - q.elicit
    e.data$family <- dens
    e.data$support <- support
    e.data$d_fun <- d_fun
    e.data$p_fun <- p_fun
    e.data$q_fun <- q_fun
    e.data
}

#' @export
pscaled_beta <- function (q, alpha, beta, a = 0, b = 1) 
{
    if (!is.finite(a) || !is.finite(b)) 
        stop("finite support required for scaled beta")
    y <- (q - a)/(b - a)
    pbeta(y, alpha, beta)
}

#' @export
qscaled_beta <- function (p, alpha, beta, a = 0, b = 1) 
{
    if (!is.finite(a) || !is.finite(b)) 
        stop("finite support required for scaled beta")
    q <- qbeta(p, alpha, beta)
    q * (b - a) + a
}
