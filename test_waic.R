require(rstan)

N=10
D=100

dat=list(N=N, D=D, K=10, y=matrix(runif(D*N),D,N), aShape=1.01, aRate=0.01, bShape=1.01, bRate=0.01)

source("waic.R")

get_waic=function(stan_fit) {
  samples=extract(stan_fit, "lp")
  waic(samples$lp)$waic
}

require(foreach)

waics=foreach(k=1:3, .combine=c) %do% {
  dat$K=k
  sf=rstan::stan("beta_cluster.stan", data=dat, cores = 4)
  get_waic(sf)    
}

waics
