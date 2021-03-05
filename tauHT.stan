data {
  int<lower=1> n; //number of subjects
  int<lower=1> J; // number of timepoints
  int<lower=1, upper = n*J> N; //number of observations
  int<lower=1, upper=n> subject[N]; //subject id
  int<lower = 1> k; // scales beta variance
  real mu0; // beta0 prior mean
  real mu1; //beta1 prior mean
  real<lower=0> nu; // half-t prior df parameter
  real<lower=0> s; // half-t prior scale parameter
  real<lower=0> a_sig; // sigma inv-gamma prior
  real<lower=0> b_sig; // sigma inv-gamma prior
  real Y[N]; // height measurements
  real X[N]; // month data

}

parameters {
  real beta0; //fixed intercept
  real beta1; //fixed slope
  real<lower=0> sigma2; //residual variance
  real<lower=0> tau; // subject sd
  vector[n] alpha; //random subject intercept
}

transformed parameters{
  real<lower = 0> sigma;
  sigma = sqrt(sigma2);
}

model {
  real mu;
  alpha ~ normal(0, tau);
  beta0 ~ normal(mu0, sqrt(k)*sigma);
  beta1 ~ normal(mu1, sqrt(k)*sigma);
  sigma2 ~ inv_gamma(a_sig, b_sig);
  tau ~ student_t(nu, 0, s);
  for(i in 1:N){
    mu = beta0 + beta1*X[i] + alpha[subject[i]];
    Y[i] ~ normal(mu, sigma);
  }
}
