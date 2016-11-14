data {
  int<lower=0> N; # samples
  int<lower=0> D; # genes
  int<lower=0> K; 
  matrix[D,N] y;
  real<lower=0> aShape; 
  real<lower=0> aRate; 
  real<lower=0> bShape; 
  real<lower=0> bRate; 
}
parameters {
  simplex[K] p; 
  real<lower=0> a[K];
  real<lower=0> b[K];
}
transformed parameters {
  vector[K] logprob[D];
  vector[K] logp; 
  vector[D] lp; 
  
  for (k in 1:K) {
    logp[k] <- log(p[k]);
  }
  for (d in 1:D) {
    for (k in 1:K) {
      real l[N]; 
      for (n in 1:N)
        l[n] <- beta_log(y[d,n], a[k], b[k]);
      logprob[d][k] <- sum(l) + logp[k];
    }
    lp[d] <- log_sum_exp(logprob[d]); 
  }
}
model {
  p ~ dirichlet( rep_vector(1.0/K,K) ); # optional, if K=3 I don't think this matters
  a ~ gamma(aShape,aRate);
  b ~ gamma(bShape,bRate);
  for (d in 1:D) {
    increment_log_prob(lp[d]); 
  }
}
