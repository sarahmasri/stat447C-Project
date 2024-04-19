data {
  int<lower=1> n; // Number of samples
  int<lower=1> p; // Number of predictors
  matrix[n,p] X;  // Design matrix
  real y[n];      // n-dimensional response vector
}

parameters {
  real beta0;                  // Intercept
  vector[p] beta;              // Coefficients
  vector<lower=0>[p] lambda;   // Global shrinkage parameter
  real<lower=0> tau;           // Local shrinkage parameter
}

transformed parameters {
  vector[n] theta ;
  theta = X * beta + beta0;
}

model {
  // Priors
  beta0 ~ normal(0, 10);
  beta ~ normal(0, 10);
  lambda ~ cauchy(0, 1);
  tau ~ cauchy(0, 1);
  
  // Horseshoe prior
  beta ~ normal(0, 1 * tau * lambda); 
  
  // Lieklihood
  y ~ normal(theta, 1);    // Gaussian likelihood
}

generated quantities {
  real y_pred[n];
  y_pred = normal_rng(theta, 1); // posterior predictive
}


