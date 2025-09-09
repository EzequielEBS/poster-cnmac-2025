data {
  // int<lower=1> K; // number of historical data
  int<lower=0> n0; // number of observations of each hist data
  int<lower=0> n; // number of current observations 
  int<lower=0> p; // number of current covariates
  matrix[n0,p] X0;
  vector[n0] y0;
  matrix[n,p] X;
  vector[n] y;
  real<lower=0> a; // prior shape
  real<lower=0> b; // prior scale
  matrix[p,p] V0; // prior precision
  vector[p] mu0; // prior mean
  real<lower=0> tilde_a;
  real<lower=0> tilde_b;
}

transformed data {
  vector[p] hat_beta = mdivide_left_spd(  X' * X,  X' * y);
  real<lower=0> S = (y - X * hat_beta)' * (y - X * hat_beta);
  vector[p] hat_beta0 = mdivide_left_spd( X0' * X0, X0' * y0);
  real<lower=0> S0 = (y0 - X0 * hat_beta0)' * (y0 - X0 * hat_beta0);
  vector[p] prod_tX0_y0 = X0'*y0;
}

parameters {
  real<lower=0, upper=1> delta;
  vector[p] beta;
  real<lower=0> sigma2;
}

transformed parameters {
  // Build default parameters
  real<lower=0> nu0 = (n0 * delta - p) * 0.5 + a - 1;
  matrix[p,p] Lambda0 = V0 + X0' * X0 * delta;
  vector[p] tilde_beta00 = V0 * mu0 + prod_tX0_y0 * delta;

  matrix[p,p] inv_Lambda0 = inverse_spd(Lambda0);
  vector[p] tilde_beta0 = inv_Lambda0 * tilde_beta00;
  real<lower=0> H0 = b + 0.5 * ( S0 * delta + 
                                hat_beta0'*prod_tX0_y0 * delta +
                                mu0' * V0 * (mu0 - tilde_beta0) -
                                (prod_tX0_y0 * delta)'*tilde_beta0
                                );
}

model {
  // prior
  target += beta_lpdf(delta| tilde_a, tilde_b);
  
  target += nu0 * log(H0) + 0.5 * log_determinant(Lambda0) -
            0.5 / sigma2 * ( (beta - tilde_beta0)' * 
                                    Lambda0 *
                                    (beta - tilde_beta0) +
                                    2*H0
                                  ) -
            (0.5*n0*delta + a)*log(sigma2) -
            lgamma(nu0);
}