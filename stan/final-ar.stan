functions {
#include utils/pois_rng_safe.stan
#include utils/pcp_ar.stan
#include utils/pcp_iid.stan
}
data {
  // total number of observations
  int<lower = 1> N;
  // response variable
  array[N] int y;
  // log of the offset term
  vector[N] logoff;
  // Number of fixed effects
  int<lower = 0> K;
  // data for group-level effects of ID 1
  // indicator variables to allow for increased complexity
  // * Covariates?
  int<lower = 0, upper = 1> toggle_reg;
  // Design matrix
  matrix[toggle_reg ? N : 1, K] X;
  // * indicator for unstructured random intercept
  int<lower = 0, upper = 1> toggle_un;
  // * indicator for AR1 effect
  int<lower = 0, upper = 1> toggle_ar;
  // * indicator for QR parametrized regression
  int<lower = 0, upper = 1> toggle_qr;
  // * indicator for computing marginal predictions
  int<lower = 0, upper = 1> toggle_pred;
  array[toggle_pred] int<lower = 0> N_pred;
  array[toggle_pred] real lo_pred;
  matrix[toggle_pred ? N_pred[1] : 1,
         toggle_pred ? K : 1] X_pred;
  array[toggle_pred ? N_pred[1] : 0] int<lower = 0> time_pred;
  // auxiliary data for iid re
  // * number of random intercepts
  int<lower = 0> N_i;
  // * indicator of grouping (from 1 to N_i)
  array[toggle_un ? N : 0] int<lower = 1, upper = N_i> group;
  //  auxiliary data for time re
  // * number of time points
  int<lower = 0> N_t;
  // * indicator of time (from 1 to N_t)
  array[toggle_ar ? N : 0] int<lower = 1, upper = N_t> time;
  // * hyperparameters for eta's prior
  real<lower=0, upper = 1> p_eta;
  real<lower=-1, upper = 1> eta_0;
}
transformed data {
  // QR parametrization
  matrix[toggle_qr ? N : 0, toggle_qr ? K : 0] Q_ast;
  matrix[toggle_qr ? K : 0, toggle_qr ? K : 0] R_ast;
  matrix[toggle_qr ? K : 0, toggle_qr ? K : 0] R_ast_inverse;
  if (toggle_reg) {
    if (toggle_qr) {
      // thin and scale the QR decomposition
      Q_ast = qr_thin_Q(X) * sqrt(N - 1);
      R_ast = qr_thin_R(X) / sqrt(N - 1);
      R_ast_inverse = inverse(R_ast);
    }
  }
}
parameters {
  // vector of regression coefficients
  vector[toggle_reg ? K : 0] beta0;
  // intercept
  real alpha;
  // random effects
  // * iid
  vector[toggle_un ? N_i : 0] w_i;
  array[toggle_un ? 1 : 0] real<lower = 0> sigma_i;
  // * ar
  vector[toggle_ar ? N_t : 0] w_t;
  array[toggle_ar ? 1 : 0] real<lower = 0> sigma_t;
  array[toggle_ar ? 1 : 0] real<lower = -1, upper = 1> eta;
}
transformed parameters {
  // initialize "intercept"
  vector[N] mu = rep_vector(0.0, N);
  mu += logoff;
  mu += alpha;
  // scaled random effects
  vector[toggle_un ? N_i : 0] z_i;
  if (toggle_un)
    z_i = sigma_i[1] * w_i;
  vector[toggle_ar ? N_t : 0] z_t;
  vector[toggle_ar ? N_t : 0] z_lag = rep_vector(0.0, toggle_ar ? N_t : 0);
  if (toggle_ar) {
    z_t = sigma_t[1] * w_t;
    for (k in 2:N_t) {
      z_lag[k] = z_t[k - 1];
      // z_t[k] += eta[1] * z_t[k - 1];
    }
  }
  // --- for prior sensitivity ---
  // * rep paramter
  array[toggle_reg] real lprior_beta;
  lprior_beta[1] =
    student_t_lpdf(beta0 | 3, 0, 1);
  // * iid re parameter
  array[toggle_un ? 1 : 0] real lprior_sigma_i;
  if (toggle_un) {
    lprior_sigma_i[1] =
      student_t_lpdf(sigma_i | 3, 0, 1)
      - 1 * student_t_lccdf(0 | 3, 0, 1);
    for (n in 1:N) {
      mu[n] += z_i[group[n]];
    }
  }
  // * Time parameters
  array[toggle_ar ? 1 : 0] real lprior_eta;
  array[toggle_ar ? 1 : 0] real lprior_sigma_t;
  if (toggle_ar) {
    lprior_sigma_t[1] =
      student_t_lpdf(sigma_t | 3, 0, 1)
      - 1 * student_t_lccdf(0 | 3, 0, 1);
    lprior_eta[1] = pcp_ar0_lpdf(eta[1] | p_eta, eta_0);
    for (n in 1:N) {
      mu[n] += z_t[time[n]] + eta[1] * z_lag[time[n]];
      // mu[n] += z_t[time[n]];
    }
  }
}
model {
  // Likelihood
  if (toggle_reg) {
    target += poisson_log_glm_lpmf(y | toggle_qr ? Q_ast : X, mu, beta0);
    target += lprior_beta[1];
  } else {
    for (n in 1:N)
      target += poisson_log_lpmf(y[n] | mu[n]);
  }
  // Priors and random effects
  target += student_t_lpdf(alpha | 3, 0, 1);
  if (toggle_un) {
    target += std_normal_lpdf(w_i);
    target += lprior_sigma_i;
  }
  if (toggle_ar) {
    target += std_normal_lpdf(w_t);
    target += lprior_sigma_t;
    target += lprior_eta;
  }
}
generated quantities {
  // log-likelihood for model comparisons
  vector[N] log_lik;
  // samples from posterior predictive dist
  vector[N] y_rep;
  vector[N] linpred = mu;
  vector[toggle_reg ? K : 0] beta;
  vector[toggle_pred ? N_pred[1] : 0] y_pred;
  vector[toggle_pred ? N_pred[1] : 0] mu_pred;
  if (toggle_reg) {
    linpred += toggle_qr ? Q_ast * beta0 : X * beta0;
    if (toggle_qr) {
      beta = R_ast_inverse * beta0;
    } else {
      beta = beta0;
    }
  }
  for (i in 1:N) {
    log_lik[i] = poisson_log_lpmf(y[i] | linpred[i]);
    y_rep[i] = pois_log_safe_rng(linpred[i]);
  }
    vector[toggle_ar ? N_t : 0] z_ar;
  if (toggle_ar) {
    z_ar[1] = z_t[1];
    for (k in 2:N_t) {
      z_ar[k] = z_t[k] + eta[1] * z_t[k - 1];
    }
  }
  if (toggle_pred) {
    for (n in 1:N_pred[1]) {
      if (toggle_ar) {
        mu_pred[n] = exp(lo_pred[1] +
                         alpha +
                         z_ar[time_pred[n]] +
                         X_pred[n] * beta);
      } else {
        mu_pred[n] = exp(lo_pred[1] +
                         alpha +
                         X_pred[n] * beta);
      }
      y_pred[n] = pois_log_safe_rng(log(mu_pred[n]));
      // mu[n] += z_t[time[n]];
    }
  }
}
