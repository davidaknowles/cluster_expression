
require(rstan)

N=100
D=1000

dat=list(N=N, D=D, K=3, y=matrix(runif(D*N),D,N))

sm=stan_model("beta_cluster.stan")

model_fit=optimizing(sm, dat=dat, init=list(a=0.5+runif(dat$K),b=0.5+runif(dat$K)) as_vector=F)

model_fit$par$a
model_fit$par$b
model_fit$par$p

assignments=apply(model_fit$par$logprob,1,which.max)
