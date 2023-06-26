## Perelson AS, Kirschner DE, De Boer R. Dynamics of HIV infection of CD4+ T
## cells. Math Biosci. 1993 Mar;114(1):81-125.
## doi: 10.1016/0025-5564(93)90043-a. PMID: 8096155.
##
## Note these notation changes compared to the paper:
## Paper | Here | Meaning
## ------+------+--------
## T*    | T_li | Latently infected CD4+ T cells.
## T**   | T_ai | Actively infected CD4+ T cells.

library(deSolve)
library(lhs)
library(plyr) # round_any
library(tidyverse)
library(broom) # tidy.optim
library(abind)

initial_values <- c( # Table 1, page 88.
    T = 1e3,  # mm^{-3}, Uninfected CD4+ cells.
    T_li = 0, # mm^{-3}, Latently infected CD4+ cells.
    T_ai = 0, # mm^{-3}, Actively infected CD4+ cells.
    V = 1e-3  # mm^{-3}, HIV cells.
)
VARIES <- 774
parameters <- c( # Table 1, page 88.
    s = 10,   # day^{-1}mm^{-3}, Rate of supply of CD4+ cells from precursors.
    r = 0.03, # day^{-1}, Rate of growth for the CD4+ cells.
    T_max = 1.5e3, # mm^{-3}, Maximum CD4+ cells.
    mu_T = 0.02,   # day^{-1}, Death rate of uninfected and latently CD4+ cells.
    mu_b = 0.24,   # day^{-1}, Death rate of actively infected CD4+ cells.
    mu_V = 2.4,    # day^{-1}, Death rate of free virus.
    k_1 = 2.4e-5,  # mm^{3}day^{-1}, Rate constant for CD4+ becoming infected.
    k_2 = 3e-3,    # day^{-1}, Rate latently to actively infected conversion.
    N = VARIES     # Number of free virus produced by lysing a CD4+ cell.
)
times <- seq(0, 10, by = 0.01) * 365 # years converted to days.

## Model 1: constant 's'.
hiv_s_const <- function(t, initial_values, parameters) {
    ## Conveniently access members by name by constructing an environment from
    ## the data using with() and as.list() as shown in the deSolve package
    ## vignette of the Lorenz equation.
    with(as.list(c(initial_values, parameters)), {
        ## Rates of change (equations 5a-5d, page 87).
        dT.dt <- s - mu_T*T + r*T*(1 - (T + T_li + T_ai)/T_max) - k_1*V*T
        dT_li.dt <- k_1*V*T - mu_T*T_li - k_2*T_li
        dT_ai.dt <- k_2*T_li - mu_b*T_ai
        dV.dt <- N*mu_b*T_ai - k_1*V*T - mu_V*V
        ## Return the rates of change.
        list(c(dT.dt, dT_li.dt, dT_ai.dt, dV.dt))
    })
}

## Read in the raw experimental data.
df_exp <-
    read_csv("../data/pantaleo1995-figure1.csv", col_types = "cdd")

## Set the experimental data boundaries.
cd4 <-
    df_exp %>%
    pull(cd4_cells_per_mm3) %>%
    range()
accuracy <- 100
cd4[1] <- round_any(cd4[1], accuracy = accuracy, f = floor)
cd4[2] <- round_any(cd4[2], accuracy = accuracy, f = ceiling)
df_boundaries <-
    tibble(year = 0:10,
           lower = cd4[1],
           upper = cd4[2])

## Use the default values in Table 1 as the medians, and to help fit the distributions
## add quantiles as 50% and 150% of the median.
df_prior_ranges <- tibble::tribble(
    ~param, ~x, ~distr,
    "s", c(5, 10, 15), "gamma",
    "r", c(0.02, 0.03, 0.04), "gamma",
    "mu_T", c(0.01, 0.02, 0.03), "gamma",
    "mu_b", c(0.12, 0.24, 0.36), "gamma",
    "mu_V", c(1.2, 2.4, 3.6), "gamma",
    "k_1", c(1.2e-5, 2.4e-5, 3.6e-5), "gamma",
    "k_2", c(2e-3, 3e-3, 4e-3), "gamma",
    "N", c(536, 878, 1340), "nbinom",
    ) %>%
    mutate(q = list(c(0.25, .5, 0.75)), .before = x)##  %>%
    ## unnest(cols = c(q, x))

## Objective function to minimize to fit the distributions, that returns the
## distance between the original and estimated values.
fit_error <- function(distr, params, q, x) {
    with(as.list(params, q, x), {
        x_est <-
            switch(distr,
                   gamma = qgamma(q, shape = shape, scale = scale),
                   nbinom = qnbinom(q, mu = mu, size = size))
        mat <- matrix(c(x, x_est), nrow = length(x))
        c(dist(t(mat)))
    })
}
df_prior_fit <-
    df_prior_ranges %>%
    group_by(param) %>%
    mutate(rescale = any(log10(unlist(x)) < -2) & distr == "gamma",
           par_init = case_when(distr == "gamma" ~ list(c(shape = 0.1, scale = 1)),
                                distr == "nbinom" ~ list(c(mu = 1000, size = 10))),
           ## Trying to fit a gamma distribution to values close to zero will
           ## fail, so increase the values and then divide the gamma scale
           ## parameter later.
           rescale_multiplier = ifelse(rescale,
                                       10^floor(-log10(mean(unlist(x)))),
                                       1)) %>%
    group_by(param, rescale_multiplier) %>%
    reframe(tidy(optim(par = unlist(par_init),
                       fn = function(...) fit_error(distr, ...),
                       q = unlist(q),
                       x = unlist(x) * rescale_multiplier,
                       ## All the parameters of the gamma and
                       ## negative-binomial distributions must be
                       ## greater than zero.
                       lower = setNames(rep(0 + .Machine$double.eps,
                                            length(unlist(par_init))),
                                        names(unlist(par_init))),
                       method = "L-BFGS-B"))) %>%
    mutate(value = case_when(parameter == "scale" &
                             rescale_multiplier != 1 ~
                                 value / rescale_multiplier,
                             .default = value))
## Validate the fitted distributions using both distance between the values and
## quantiles.
df_prior_validate <-
    df_prior_fit %>%
    pivot_wider(names_from = parameter, values_from = value)
args_list <-
    df_prior_validate %>%
    select(-param, -rescale_multiplier) %>%
    as.matrix() %>%
    apply(1, list) %>%
    lapply(unlist, recursive = FALSE) %>%
    lapply(na.omit) %>%
    ## We don't need to know which values were dropped.
    lapply(`attr<-`, "na.action", NULL)
df_prior_validate %>%
    select(param) %>%
    inner_join(df_prior_ranges) %>%
    group_by(param) %>%
    mutate(q_est = list(do.call(str_c("p", distr),
                                c(args_list[[cur_group_id()]], q = x))),
           x_est = list(do.call(str_c("q", distr),
                                c(args_list[[cur_group_id()]], p = q)))) %>%
    unnest(cols = c(x, x_est, q, q_est)) %>%
    print(n = Inf)

## Calibration adjustment function: Alternative Density Subtraction (ADS).
##
## Joslyn LR, Kirschner DE, Linderman JJ. CaliPro: A Calibration Protocol That
## Utilizes Parameter Density Estimation to Explore Parameter Space and
## Calibrate Complex Biological Models. Cell Mol Bioeng. 2020 Sep
## 15;14(1):31-47. doi: 10.1007/s12195-020-00650-z
ads <- function(df, n = 512) {
    range <- df %>% pull(value) %>% range()
    ## Fixed value / no passing values.
    if (all(range == range[1]) ||
        nrow(filter(df, pass)) == 0) {
        return(range)
    }
    bins_all <-
        df %>%
        group_by(pass) %>%
        reframe(density(value,
                        from = range[1],
                        to = range[2],
                        n = n) %>%
                    unclass() %>%
                    .[1:2] %>%
                as_tibble())
    bins_ads <-
        bins_all %>%
        mutate(pass = ifelse(pass, "dens_pass", "dens_fail")) %>%
        pivot_wider(id_cols = x, names_from = pass, values_from = y) %>%
        filter(dens_pass - dens_fail >= 0)
    if (nrow(bins_ads) == 0)
        return(range) ## No better values.
    bins_ads %>% pull(x) %>% range()
}

## Calibration loop.
n_iter <- 1
set.seed(123)
seeds <- sample(1:1000, size = n_iter)
iter <- 1
for (iter in 1:n_iter) {
    ## Generate LHS samples.
    ##
    ## First sample random uniform LHS values for each replicate.
    set.seed(seeds[iter])
    n_replicates <- 5
    n_samples_per_replicate <- 100
    lhs_cdf <- replicate(n_replicates,
                         randomLHS(n = n_samples_per_replicate,
                                   k = length(parameters)))
    colnames(lhs_cdf) <- names(parameters)
    ## Combine replicate dimensions.
    lhs_cdf <- array(aperm(lhs_cdf, c(1, 3, 2)),
                     dim = c(n_replicates * n_samples_per_replicate,
                             length(parameters)),
                     dimnames = list(NULL, names(parameters)))    
    ## Final LHS samples applying .
    lhs <- array(-1, dim = dim(lhs_cdf), dimnames = dimnames(lhs_cdf))
    lhs[, "s"] <- qgamma(lhs_cdf[, "s"], shape = 1.985656, scale = 5.681687)
    lhs[, "r"] <- qgamma(lhs_cdf[, "r"], shape = 4.530347876, scale = 0.006990707)
    lhs[, "T_max"] <- qnbinom(lhs_cdf[, "T_max"], mu = 562, size = 14.0126)
    lhs[, "mu_T"] <- qgamma(lhs_cdf[, "mu_T"], shape = 2.10552523, scale = 0.01068658)
    lhs[, "mu_b"] <- qgamma(lhs_cdf[, "mu_b"], shape = 1.9856561, scale = 0.1363606)
    lhs[, "mu_V"] <- qgamma(lhs_cdf[, "mu_V"], shape = 1.985657, scale = 1.363605)
    lhs[, "k_1"] <- qgamma(lhs_cdf[, "k_1"], shape = 1.985657, scale = 1.363605e-5)
    lhs[, "k_2"] <- qgamma(lhs_cdf[, "k_2"], shape = 1.594566171, scale = 0.002008847)
    ## Fix the bifurcation parameter for now instead of varying it.
    lhs[, "N"] <- 1400L

    ## Run the simulations with the LHS parameter samples.
    out <- list()
    i <- 1
    system.time({
        for (i in 1:nrow(lhs)) {
            out[[i]] <- ode(y = initial_values,
                            times = times,
                            func = hiv_s_const,
                            parms = lhs[i, ])
            i <- i + 1
        }
    })

    ## Gather the output from the LHS parameter samples.
    df_sim <-
        lapply(out, as.data.frame) %>%
        bind_rows(.id = "param_set") %>%
        as_tibble() %>%
        mutate(T_sum = T + T_li + T_ai,
               time_years = time / 365,
               param_set = as.integer(param_set),
               replicate = ceiling(param_set / n_samples_per_replicate))

    ## Classify model outputs as pass or fail using the data boundaries.
    df_class <-
        df_sim %>%
        group_by(param_set) %>%
        mutate(pass = all(T_sum >= cd4[1] & T_sum <= cd4[2]))

    message(sprintf("Overall pass rate across all replicates: %.1f%%",
                    df_class %>% pull(pass) %>% mean() * 100))

    message("Pass rate per replicate:")
    print(
        df_class %>%
        group_by(replicate) %>%
        summarize(pass = mean(pass) * 100)
    )

    ## Old parameter ranges.
    df_lhs <-
        as_tibble(lhs, rownames = "param_set") %>%
        bind_cols(as_tibble(lhs_cdf) %>%
                  rename_with(~ str_c(.x, "-cdf"))) %>%
        mutate(param_set = as.integer(param_set)) %>%
        inner_join(df_class %>% distinct(param_set, replicate, pass),
                   by = "param_set") %>%
        pivot_longer(cols = ! c(param_set, replicate, pass),
                     names_to = "param", values_to = "value") %>%
        separate_wider_delim(param, delim = "-",
                             names = c("param", "cdf"),
                             too_few = "align_start") %>%
        mutate(value_is_cdf = ! is.na(cdf)) %>%
        select(-cdf)
    ranges_prev <-
        df_lhs %>%
        select(param, value) %>%
        pivot_wider(names_from = param,
                    values_from = value,
                    values_fn = list) %>%
        reframe(across(everything(), ~ range(.x[[1]])))

    ## Calculate the new parameter ranges using ADS.
    ranges_new <-
        df_lhs %>%
        nest(.by = param) %>%
        pivot_wider(names_from = param,
                    values_from = data) %>%
        reframe(across(everything(), ~ ads(.x[[1]])))

    ## Tabulate parameter ranges before and after ADS adjutment.
    ranges_both <-
        bind_rows(before = ranges_prev, after = ranges_new, .id = "range")

    ## Find parameter ranges to change.
    ranges_changed <-
        ranges_both %>%
        mutate(bound = rep(c("lower", "upper"), 2), .before = 2) %>%
        pivot_longer(cols = ! c(range, bound), names_to = "param", values_to = "value") %>%
        pivot_wider(id_cols = c("param", "bound"), names_from = "range") %>%
        mutate(changed = before != after)

    ## If no ranges changed, end calibration.
    if (! any(pull(ranges_changed, changed))) {
        message(sprintf(c("No parameter ranges changed in iteration %d. ",
                          "Ended calibration loop."), iter))
        break
    }

    ## TODO: Fit priors to new parameter ranges.
}
