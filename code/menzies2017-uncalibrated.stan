data {
  /* Prevalence */
  int<lower=1> T;
  vector<lower=0, upper=1>[T] prev;
  
  /* Births */

  /* Annual birth rate for steady state based on 1M pop & mu_B. */
  int<lower=0> a;

  /* Mortality rates */
  real<lower=0, upper=1>
    /* Background mortality rate. */
    mu_B,
    /* Treatment uptake for early disease. */
    r_E;
}

parameters {
  real<lower=0>
    /* Fraction of births entering non-susceptible state. */
    b_alpha,
    b_beta;

  /* Mortality rates */
  real
    mu_E_mean,			/* Early disease. */
    mu_L_mean,			/* Late disease. */
    mu_T_mean;			/* Treatment. */
  real<lower=0>
    mu_E_sd,
    mu_L_sd,
    mu_T_sd;

  /* Effective rate for transmission */
  real rho_mean;
  real<lower=0> rho_sd;

  /* Disease rates */
  real
    p_mean,	 /* Progression from early to late disease. */
    r_L_mean;	 /* Treatment uptake for late disease. */
  real<lower=0>
    p_sd,
    r_L_sd;

  /* Annual cost of treatment */
  real c_T_mean;
  real<lower=0>;
}

model {
  /* Parameter distributions. */
  real b ~ beta(b_alpha, b_beta);
  real mu_E ~ log_normal(mu_E_mean, mu_E_sd);
  real mu_L ~ log_normal(mu_L_mean, mu_L_sd);
  real mu_T ~ log_normal(mu_T_mean, mu_T_sd);
  real rho ~ log_normal(rho_mean, rho_sd);
  real p ~ log_normal(p_mean, p_sd);
  real r_L ~ log_normal(r_L_mean, r_L_sd);
  real c_T ~ log_normal(c_T_mean, c_T_sd);

  /* Prevalence. */
  for () {
  prev[t] ~ (E[t] + L[t] + T[t]) / (N[t] + S[t] + E[t] + L[t] + T[t]);
  }
}

generated quantities {
  array[T] vector[3] y_sim = ode_rk45(tumor, y0, t0, ts, a, b, c, g1, g2, g3,
				      mu2, mu3, p1, p2, r2, s1, s2);
}

functions {
  /* ODEs. */
  vector hiv_pop(/* Timepoints */
		 real t,
		 /* Random variables to solve, y := (N, S, E, L, T, D) */
		 vector y,
		 /* Parameters start here */
		 ) {
    vector[6] dydt;
    real
      N = y[1],
      S = y[2],
      E = y[3],
      L = y[4],
      T = y[5],
      D = y[6];
    /* The force of infection "lambda" is an subexpression used in some the
       differential equations below.  "lambda" is defined in section "3.2 Study
       model" of the main paper and also at the end of the section "3 Model
       equations" of the PDF supplement.  Note that we call this a subexpression
       rather than use the author's term in the supplement of "endogenous time
       varying parameter" because in the context of pure ODEs, parameters are
       fixed and cannot depend on variables. */
    lambda = rho * (E + L) / (S + E + L + T);
    /* Non-susceptible (N) */
    dydt[1] = a*b - mu_B*N;
    /* Susceptible (S) */
    dydt[2] = a*(1 - b) - mu_B*S - lambda*S;
    /* Early disease (E) */
    dydt[3] = lambda*S - (mu_B + mu_E)*E - r_E*E - c_T*E;
    /* Late disease (L) */
    dydt[4] = c_T*E - (mu_B + mu_L)*L - r_L*L;
    /* Treatment (T) */
    dydt[5] = r_E*E + r_L*L - (mu_B + mu_T)*T;
    /* Dead (D) */
    dydt[6] = mu_B*(N + S + E + L + T) + mu_E*E + mu_L*L + mu_T*T;
    return dydt;
}
