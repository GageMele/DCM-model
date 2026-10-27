// generated with brms 2.23.0
functions {
  /* integer sequence of values
   * Args:
   *   start: starting integer
   *   end: ending integer
   * Returns:
   *   an integer sequence from start to end
   */
  array[] int sequence(int start, int end) {
    array[end - start + 1] int seq;
    for (n in 1:num_elements(seq)) {
      seq[n] = n + start - 1;
    }
    return seq;
  }
  // compute partial sums of the log-likelihood
  real partial_log_lik_lpmf(array[] int seq, int start, int end, data vector Y, data matrix Xc, vector b, real Intercept, data matrix Xs, vector bs, data matrix Zs_1_1, vector s_1_1, data matrix Zs_2_1, vector s_2_1, data matrix Zs_3_1, vector s_3_1, data matrix Zs_4_1, vector s_4_1, data matrix Zs_5_1, vector s_5_1, data matrix Zs_6_1, vector s_6_1, data matrix Zs_7_1, vector s_7_1, data matrix Zs_8_1, vector s_8_1, data matrix Zs_9_1, vector s_9_1, data matrix Zs_10_1, vector s_10_1, data matrix Zs_11_1, vector s_11_1, data matrix Zs_12_1, vector s_12_1, data matrix Zs_13_1, vector s_13_1, data matrix Zs_14_1, vector s_14_1, data matrix Zs_15_1, vector s_15_1, data matrix Zs_16_1, vector s_16_1, data matrix X_sigma, vector b_sigma, data array[] int J_1, data vector Z_1_1, vector r_1_1, data array[] int J_2, data vector Z_2_1, vector r_2_1) {
    real ptarget = 0;
    int N = end - start + 1;
    // initialize linear predictor term
    vector[N] mu = rep_vector(0.0, N);
    // initialize linear predictor term
    vector[N] sigma = rep_vector(0.0, N);
    mu += Intercept + Xc[start:end] * b + Xs[start:end] * bs + Zs_1_1[start:end] * s_1_1 + Zs_2_1[start:end] * s_2_1 + Zs_3_1[start:end] * s_3_1 + Zs_4_1[start:end] * s_4_1 + Zs_5_1[start:end] * s_5_1 + Zs_6_1[start:end] * s_6_1 + Zs_7_1[start:end] * s_7_1 + Zs_8_1[start:end] * s_8_1 + Zs_9_1[start:end] * s_9_1 + Zs_10_1[start:end] * s_10_1 + Zs_11_1[start:end] * s_11_1 + Zs_12_1[start:end] * s_12_1 + Zs_13_1[start:end] * s_13_1 + Zs_14_1[start:end] * s_14_1 + Zs_15_1[start:end] * s_15_1 + Zs_16_1[start:end] * s_16_1;
    sigma += X_sigma[start:end] * b_sigma;
    for (n in 1:N) {
      // add more terms to the linear predictor
      int nn = n + start - 1;
      mu[n] += r_1_1[J_1[nn]] * Z_1_1[nn] + r_2_1[J_2[nn]] * Z_2_1[nn];
    }
    sigma = exp(sigma);
    ptarget += normal_lpdf(Y[start:end] | mu, sigma);
    return ptarget;
  }
}
data {
  int<lower=1> N;  // total number of observations
  vector[N] Y;  // response variable
  int<lower=1> K;  // number of population-level effects
  matrix[N, K] X;  // population-level design matrix
  int<lower=1> Kc;  // number of population-level effects after centering
  // data for splines
  int Ks;  // number of linear effects
  matrix[N, Ks] Xs;  // design matrix for the linear effects
  // data for spline 1
  int nb_1;  // number of bases
  array[nb_1] int knots_1;  // number of knots
  // basis function matrices
  matrix[N, knots_1[1]] Zs_1_1;
  // data for spline 2
  int nb_2;  // number of bases
  array[nb_2] int knots_2;  // number of knots
  // basis function matrices
  matrix[N, knots_2[1]] Zs_2_1;
  // data for spline 3
  int nb_3;  // number of bases
  array[nb_3] int knots_3;  // number of knots
  // basis function matrices
  matrix[N, knots_3[1]] Zs_3_1;
  // data for spline 4
  int nb_4;  // number of bases
  array[nb_4] int knots_4;  // number of knots
  // basis function matrices
  matrix[N, knots_4[1]] Zs_4_1;
  // data for spline 5
  int nb_5;  // number of bases
  array[nb_5] int knots_5;  // number of knots
  // basis function matrices
  matrix[N, knots_5[1]] Zs_5_1;
  // data for spline 6
  int nb_6;  // number of bases
  array[nb_6] int knots_6;  // number of knots
  // basis function matrices
  matrix[N, knots_6[1]] Zs_6_1;
  // data for spline 7
  int nb_7;  // number of bases
  array[nb_7] int knots_7;  // number of knots
  // basis function matrices
  matrix[N, knots_7[1]] Zs_7_1;
  // data for spline 8
  int nb_8;  // number of bases
  array[nb_8] int knots_8;  // number of knots
  // basis function matrices
  matrix[N, knots_8[1]] Zs_8_1;
  // data for spline 9
  int nb_9;  // number of bases
  array[nb_9] int knots_9;  // number of knots
  // basis function matrices
  matrix[N, knots_9[1]] Zs_9_1;
  // data for spline 10
  int nb_10;  // number of bases
  array[nb_10] int knots_10;  // number of knots
  // basis function matrices
  matrix[N, knots_10[1]] Zs_10_1;
  // data for spline 11
  int nb_11;  // number of bases
  array[nb_11] int knots_11;  // number of knots
  // basis function matrices
  matrix[N, knots_11[1]] Zs_11_1;
  // data for spline 12
  int nb_12;  // number of bases
  array[nb_12] int knots_12;  // number of knots
  // basis function matrices
  matrix[N, knots_12[1]] Zs_12_1;
  // data for spline 13
  int nb_13;  // number of bases
  array[nb_13] int knots_13;  // number of knots
  // basis function matrices
  matrix[N, knots_13[1]] Zs_13_1;
  // data for spline 14
  int nb_14;  // number of bases
  array[nb_14] int knots_14;  // number of knots
  // basis function matrices
  matrix[N, knots_14[1]] Zs_14_1;
  // data for spline 15
  int nb_15;  // number of bases
  array[nb_15] int knots_15;  // number of knots
  // basis function matrices
  matrix[N, knots_15[1]] Zs_15_1;
  // data for spline 16
  int nb_16;  // number of bases
  array[nb_16] int knots_16;  // number of knots
  // basis function matrices
  matrix[N, knots_16[1]] Zs_16_1;
  int<lower=1> K_sigma;  // number of population-level effects
  matrix[N, K_sigma] X_sigma;  // population-level design matrix
  int grainsize;  // grainsize for threading
  // data for group-level effects of ID 1
  int<lower=1> N_1;  // number of grouping levels
  int<lower=1> M_1;  // number of coefficients per level
  array[N] int<lower=1> J_1;  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_1_1;
  // data for group-level effects of ID 2
  int<lower=1> N_2;  // number of grouping levels
  int<lower=1> M_2;  // number of coefficients per level
  array[N] int<lower=1> J_2;  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_2_1;
  int prior_only;  // should the likelihood be ignored?
}
transformed data {
  matrix[N, Kc] Xc;  // centered version of X without an intercept
  vector[Kc] means_X;  // column means of X before centering
  array[N] int seq = sequence(1, N);
  for (i in 2:K) {
    means_X[i - 1] = mean(X[, i]);
    Xc[, i - 1] = X[, i] - means_X[i - 1];
  }
}
parameters {
  vector[Kc] b;  // regression coefficients
  real Intercept;  // temporary intercept for centered predictors
  vector[Ks] bs;  // unpenalized spline coefficients
  // parameters for spline 1
  // standardized penalized spline coefficients
  vector[knots_1[1]] zs_1_1;
  vector<lower=0>[nb_1] sds_1;  // SDs of penalized spline coefficients
  // parameters for spline 2
  // standardized penalized spline coefficients
  vector[knots_2[1]] zs_2_1;
  vector<lower=0>[nb_2] sds_2;  // SDs of penalized spline coefficients
  // parameters for spline 3
  // standardized penalized spline coefficients
  vector[knots_3[1]] zs_3_1;
  vector<lower=0>[nb_3] sds_3;  // SDs of penalized spline coefficients
  // parameters for spline 4
  // standardized penalized spline coefficients
  vector[knots_4[1]] zs_4_1;
  vector<lower=0>[nb_4] sds_4;  // SDs of penalized spline coefficients
  // parameters for spline 5
  // standardized penalized spline coefficients
  vector[knots_5[1]] zs_5_1;
  vector<lower=0>[nb_5] sds_5;  // SDs of penalized spline coefficients
  // parameters for spline 6
  // standardized penalized spline coefficients
  vector[knots_6[1]] zs_6_1;
  vector<lower=0>[nb_6] sds_6;  // SDs of penalized spline coefficients
  // parameters for spline 7
  // standardized penalized spline coefficients
  vector[knots_7[1]] zs_7_1;
  vector<lower=0>[nb_7] sds_7;  // SDs of penalized spline coefficients
  // parameters for spline 8
  // standardized penalized spline coefficients
  vector[knots_8[1]] zs_8_1;
  vector<lower=0>[nb_8] sds_8;  // SDs of penalized spline coefficients
  // parameters for spline 9
  // standardized penalized spline coefficients
  vector[knots_9[1]] zs_9_1;
  vector<lower=0>[nb_9] sds_9;  // SDs of penalized spline coefficients
  // parameters for spline 10
  // standardized penalized spline coefficients
  vector[knots_10[1]] zs_10_1;
  vector<lower=0>[nb_10] sds_10;  // SDs of penalized spline coefficients
  // parameters for spline 11
  // standardized penalized spline coefficients
  vector[knots_11[1]] zs_11_1;
  vector<lower=0>[nb_11] sds_11;  // SDs of penalized spline coefficients
  // parameters for spline 12
  // standardized penalized spline coefficients
  vector[knots_12[1]] zs_12_1;
  vector<lower=0>[nb_12] sds_12;  // SDs of penalized spline coefficients
  // parameters for spline 13
  // standardized penalized spline coefficients
  vector[knots_13[1]] zs_13_1;
  vector<lower=0>[nb_13] sds_13;  // SDs of penalized spline coefficients
  // parameters for spline 14
  // standardized penalized spline coefficients
  vector[knots_14[1]] zs_14_1;
  vector<lower=0>[nb_14] sds_14;  // SDs of penalized spline coefficients
  // parameters for spline 15
  // standardized penalized spline coefficients
  vector[knots_15[1]] zs_15_1;
  vector<lower=0>[nb_15] sds_15;  // SDs of penalized spline coefficients
  // parameters for spline 16
  // standardized penalized spline coefficients
  vector[knots_16[1]] zs_16_1;
  vector<lower=0>[nb_16] sds_16;  // SDs of penalized spline coefficients
  vector[K_sigma] b_sigma;  // regression coefficients
  vector<lower=0>[M_1] sd_1;  // group-level standard deviations
  array[M_1] vector[N_1] z_1;  // standardized group-level effects
  vector<lower=0>[M_2] sd_2;  // group-level standard deviations
  array[M_2] vector[N_2] z_2;  // standardized group-level effects
}
transformed parameters {
  // penalized spline coefficients
  vector[knots_1[1]] s_1_1;
  // penalized spline coefficients
  vector[knots_2[1]] s_2_1;
  // penalized spline coefficients
  vector[knots_3[1]] s_3_1;
  // penalized spline coefficients
  vector[knots_4[1]] s_4_1;
  // penalized spline coefficients
  vector[knots_5[1]] s_5_1;
  // penalized spline coefficients
  vector[knots_6[1]] s_6_1;
  // penalized spline coefficients
  vector[knots_7[1]] s_7_1;
  // penalized spline coefficients
  vector[knots_8[1]] s_8_1;
  // penalized spline coefficients
  vector[knots_9[1]] s_9_1;
  // penalized spline coefficients
  vector[knots_10[1]] s_10_1;
  // penalized spline coefficients
  vector[knots_11[1]] s_11_1;
  // penalized spline coefficients
  vector[knots_12[1]] s_12_1;
  // penalized spline coefficients
  vector[knots_13[1]] s_13_1;
  // penalized spline coefficients
  vector[knots_14[1]] s_14_1;
  // penalized spline coefficients
  vector[knots_15[1]] s_15_1;
  // penalized spline coefficients
  vector[knots_16[1]] s_16_1;
  vector[N_1] r_1_1;  // actual group-level effects
  vector[N_2] r_2_1;  // actual group-level effects
  // prior contributions to the log posterior
  real lprior = 0;
  // compute penalized spline coefficients
  s_1_1 = sds_1[1] * zs_1_1;
  // compute penalized spline coefficients
  s_2_1 = sds_2[1] * zs_2_1;
  // compute penalized spline coefficients
  s_3_1 = sds_3[1] * zs_3_1;
  // compute penalized spline coefficients
  s_4_1 = sds_4[1] * zs_4_1;
  // compute penalized spline coefficients
  s_5_1 = sds_5[1] * zs_5_1;
  // compute penalized spline coefficients
  s_6_1 = sds_6[1] * zs_6_1;
  // compute penalized spline coefficients
  s_7_1 = sds_7[1] * zs_7_1;
  // compute penalized spline coefficients
  s_8_1 = sds_8[1] * zs_8_1;
  // compute penalized spline coefficients
  s_9_1 = sds_9[1] * zs_9_1;
  // compute penalized spline coefficients
  s_10_1 = sds_10[1] * zs_10_1;
  // compute penalized spline coefficients
  s_11_1 = sds_11[1] * zs_11_1;
  // compute penalized spline coefficients
  s_12_1 = sds_12[1] * zs_12_1;
  // compute penalized spline coefficients
  s_13_1 = sds_13[1] * zs_13_1;
  // compute penalized spline coefficients
  s_14_1 = sds_14[1] * zs_14_1;
  // compute penalized spline coefficients
  s_15_1 = sds_15[1] * zs_15_1;
  // compute penalized spline coefficients
  s_16_1 = sds_16[1] * zs_16_1;
  r_1_1 = (sd_1[1] * (z_1[1]));
  r_2_1 = (sd_2[1] * (z_2[1]));
  lprior += normal_lpdf(b | 0,2);
  lprior += student_t_lpdf(Intercept | 3,0.73,5.448);
  lprior += normal_lpdf(bs | 0,2);
  lprior += student_t_lpdf(sds_1 | 3,0,1.5)
    - 1 * student_t_lccdf(0 | 3,0,1.5);
  lprior += student_t_lpdf(sds_2 | 3,0,1.5)
    - 1 * student_t_lccdf(0 | 3,0,1.5);
  lprior += student_t_lpdf(sds_3 | 3,0,1.5)
    - 1 * student_t_lccdf(0 | 3,0,1.5);
  lprior += student_t_lpdf(sds_4 | 3,0,1.5)
    - 1 * student_t_lccdf(0 | 3,0,1.5);
  lprior += student_t_lpdf(sds_5 | 3,0,1.5)
    - 1 * student_t_lccdf(0 | 3,0,1.5);
  lprior += student_t_lpdf(sds_6 | 3,0,1.5)
    - 1 * student_t_lccdf(0 | 3,0,1.5);
  lprior += student_t_lpdf(sds_7 | 3,0,1.5)
    - 1 * student_t_lccdf(0 | 3,0,1.5);
  lprior += student_t_lpdf(sds_8 | 3,0,1.5)
    - 1 * student_t_lccdf(0 | 3,0,1.5);
  lprior += student_t_lpdf(sds_9 | 3,0,1.5)
    - 1 * student_t_lccdf(0 | 3,0,1.5);
  lprior += student_t_lpdf(sds_10 | 3,0,1.5)
    - 1 * student_t_lccdf(0 | 3,0,1.5);
  lprior += student_t_lpdf(sds_11 | 3,0,1.5)
    - 1 * student_t_lccdf(0 | 3,0,1.5);
  lprior += student_t_lpdf(sds_12 | 3,0,1.5)
    - 1 * student_t_lccdf(0 | 3,0,1.5);
  lprior += student_t_lpdf(sds_13 | 3,0,1.5)
    - 1 * student_t_lccdf(0 | 3,0,1.5);
  lprior += student_t_lpdf(sds_14 | 3,0,1.5)
    - 1 * student_t_lccdf(0 | 3,0,1.5);
  lprior += student_t_lpdf(sds_15 | 3,0,1.5)
    - 1 * student_t_lccdf(0 | 3,0,1.5);
  lprior += student_t_lpdf(sds_16 | 3,0,1.5)
    - 1 * student_t_lccdf(0 | 3,0,1.5);
  lprior += normal_lpdf(b_sigma | 0,1);
  lprior += student_t_lpdf(sd_1 | 3,0,1.5)
    - 1 * student_t_lccdf(0 | 3,0,1.5);
  lprior += student_t_lpdf(sd_2 | 3,0,1.5)
    - 1 * student_t_lccdf(0 | 3,0,1.5);
}
model {
  // likelihood including constants
  if (!prior_only) {
    target += reduce_sum(partial_log_lik_lpmf, seq, grainsize, Y, Xc, b, Intercept, Xs, bs, Zs_1_1, s_1_1, Zs_2_1, s_2_1, Zs_3_1, s_3_1, Zs_4_1, s_4_1, Zs_5_1, s_5_1, Zs_6_1, s_6_1, Zs_7_1, s_7_1, Zs_8_1, s_8_1, Zs_9_1, s_9_1, Zs_10_1, s_10_1, Zs_11_1, s_11_1, Zs_12_1, s_12_1, Zs_13_1, s_13_1, Zs_14_1, s_14_1, Zs_15_1, s_15_1, Zs_16_1, s_16_1, X_sigma, b_sigma, J_1, Z_1_1, r_1_1, J_2, Z_2_1, r_2_1);
  }
  // priors including constants
  target += lprior;
  target += std_normal_lpdf(zs_1_1);
  target += std_normal_lpdf(zs_2_1);
  target += std_normal_lpdf(zs_3_1);
  target += std_normal_lpdf(zs_4_1);
  target += std_normal_lpdf(zs_5_1);
  target += std_normal_lpdf(zs_6_1);
  target += std_normal_lpdf(zs_7_1);
  target += std_normal_lpdf(zs_8_1);
  target += std_normal_lpdf(zs_9_1);
  target += std_normal_lpdf(zs_10_1);
  target += std_normal_lpdf(zs_11_1);
  target += std_normal_lpdf(zs_12_1);
  target += std_normal_lpdf(zs_13_1);
  target += std_normal_lpdf(zs_14_1);
  target += std_normal_lpdf(zs_15_1);
  target += std_normal_lpdf(zs_16_1);
  target += std_normal_lpdf(z_1[1]);
  target += std_normal_lpdf(z_2[1]);
}
generated quantities {
  // actual population-level intercept
  real b_Intercept = Intercept - dot_product(means_X, b);
}

