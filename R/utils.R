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

x = sort(rnorm(100, 0,1))
y = sort(rnorm(100,2,1))

M_xy = sapply(y, function(t) (x-t)^2)
image(log(M_xy))

Wass= Barycenter::Greenkhorn(rep(1/100,100), rep(1/100,100), costm = M_xy)

image(log(Wass$Transportplan))


wb  = NULL
tot_cost = 10
for(i in 1:10){
    w = sort(rnorm(100, 1, 1))
    M_wx = sapply(w, function(t) (x-t)^2)
    M_wy = sapply(w, function(t) (y-t)^2)
    Wassx= Barycenter::Greenkhorn(rep(1/100,100), rep(1/100,100), costm = M_wx)
    Wassy= Barycenter::Greenkhorn(rep(1/100,100), rep(1/100,100), costm = M_wy)
    tot_cost_new = Wassx$Distance + Wassy$Distance
    if(tot_cost_new < tot_cost){
        wb = w
        tot_cost = tot_cost_new
    }

}
tot_cost
plot(x=wb, y = 1:100)
points(x = seq(-2,3,0.1), y = pnorm(seq(-2,3,0.1), 1, 1)*100, type = "l")



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

