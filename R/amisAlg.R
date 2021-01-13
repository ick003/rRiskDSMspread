#' @export
#' @example
#' y = mvtnorm::rmvnorm(50,mean = rep(2,7),sigma = diag(0.05, 7))
#' proposal = list(rp = function(N){mvtnorm::rmvnorm(N,mean = rep(0.5,7),sigma = diag(0.1, 7))},
#'                 dp = function(x){mvtnorm::dmvnorm(x,mean = rep(0.5,7),sigma = diag(0.1, 7))})
#' prior = function(x){mvtnorm::dmvnorm(x,mean = rep(1,7),sigma = diag(0.01, 7))}
#' L = function(x) {sum(apply(x,1, function(s) sum(mvtnorm::dmvnorm(y,mean = s,sigma = diag(10, 7), log = TRUE))))}
amisALG <- function(L, prior, proposal, Np, K, ...){


# Initialisation

x <- proposal$rp(Np)
delta = Np * proposal$dp(x)
w = L(x, ...) + log(prior(x)) - log(proposal$dp(x))

w = w - min(w)
w = exp(w) / sum(exp(w))

# Iterations



}
