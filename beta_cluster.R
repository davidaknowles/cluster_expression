
require(rstan)

N=10
D=100

dat=list(N=N, D=D, K=10, y=matrix(runif(D*N),D,N), aShape=1.01, aRate=0.01, bShape=1.01, bRate=0.01)

sm=stan_model("beta_cluster.stan")

model_fit=optimizing(sm, dat=dat, init=list(a=0.5+runif(dat$K),b=0.5+runif(dat$K), p=1.0/dat$K+numeric(dat$K)), as_vector=F)

model_fit$par$a
model_fit$par$b
model_fit$par$p

assignments=apply(model_fit$par$logprob,1,which.max)
