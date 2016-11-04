data {
  int<lower=0> N; # samples
  int<lower=0> D; # genes
  int<lower=0> K; 
  matrix[D,N] y;
}
parameters {
  simplex[K] p; 
  real<lower=0> a[K];
  real<lower=0> b[K];
}
transformed parameters {
  vector[K] logprob[D];
  vector[K] logp; 
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
  }
}
model {
  # p ~ dirichlet( rep_vector(1.0/K,K) ); # optional, if K=3 I don't think this matters
  a ~ gamma(1.01,0.01);
  b ~ gamma(1.01,0.01);
  for (d in 1:D) {
    increment_log_prob(log_sum_exp(logprob[d])); 
  }
}
