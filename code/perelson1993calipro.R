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
theme_set(theme_bw())
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

## Visualize the raw data.
gg_style <- function(gg) {
    gg +
        scale_x_continuous(breaks = 1:10) +
        labs(x = "Year",
             y = expression(CD4^"+" ~ "T-cells" ~ per ~ mm^3)) +
        guides(color = "none")
}

gg_style(
    df_exp %>%
    ggplot(aes(x = year, y = cd4_cells_per_mm3, group = group, color = group)) +
    geom_line()
)

scale <- 0.5
width <- 13.33 * scale
height <- 6.06 * scale
ggsave("../results/hiv-raw.pdf", width = width, height = height)

## Visualize the range.
gg_style(
    df_exp %>%
    group_by(year) %>%
    reframe(cd4 = range(cd4_cells_per_mm3)) %>%
    group_by(year) %>%
    mutate(boundary = c("lower", "upper")[row_number()]) %>%
    pivot_wider(id_cols = year, names_from = boundary, values_from = cd4) %>%
    ggplot(aes(x = year)) +
    geom_ribbon(alpha = 0.3, aes(ymin = lower, ymax = upper)) +
    geom_point(data = df_exp, aes(y = cd4_cells_per_mm3, color = group))
)

ggsave("../results/hiv-range.pdf", width = width, height = height)

## Visualize the smoothing function.
gg_style(
    df_exp %>%
    ggplot(aes(x = year, y = cd4_cells_per_mm3)) +
    geom_point(aes(color = group)) +
    geom_smooth()
)

ggsave("../results/hiv-smooth.pdf", width = width, height = height)

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
gg_style(
    df_boundaries %>%
    ggplot(aes(x = year)) +
    geom_ribbon(alpha = 0.2, aes(ymin = lower, ymax = upper)) +
    geom_point(data = df_exp, aes(y = cd4_cells_per_mm3, color = group)) +
    geom_label(aes(y = lower, label = lower, x = mean(year))) +
    geom_label(aes(y = upper, label = upper, x = mean(year)))
)

ggsave("../results/hiv-boundaries.pdf", width = width, height = height)

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

    diagnostics(out[[1]])
    summary(out[[1]])

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

    ## Plot the trajectories.
    gg_style(
        df_class %>%
        ggplot(aes(x = time_years, y = T_sum)) +
        geom_line(alpha = 0.1, aes(group = param_set, color = pass)) +
        geom_point(data = df_exp, aes(x = year, y = cd4_cells_per_mm3))
    )

    ggsave(sprintf("../results/hiv-calipro-traj-iter-%d.pdf", iter),
           width = width, height = height)

    ## Plot the trajectory 5%-50%-95% quantiles.
    ##
    ## See https://stackoverflow.com/a/51993331
    quantiles <- function(x, q = c(0.05, 0.95), na.rm = TRUE) {
        tibble(y = median(x),
               ymin = quantile(x, probs = q[1]),
               ymax = quantile(x, probs = q[2]))
    }
    gg_style(
        df_class %>%
        ggplot(aes(x = time_years, y = T_sum)) +
        facet_wrap(~ pass, nrow = 2) +
        geom_smooth(stat = "summary", fun.data = quantiles,
                    aes(group = pass, color = pass, fill = pass)) +
        geom_point(data = df_exp, aes(x = year, y = cd4_cells_per_mm3)) +
        guides(color = "none", fill = "none")
    )

    ggsave(sprintf("../results/hiv-calipro-quantiles-iter-%d.pdf", iter),
           width = width, height = height)

    ## Plot the parameter rug plots.
    df_lhs <-
        as_tibble(lhs, rownames = "param_set") %>%
        mutate(param_set = as.integer(param_set)) %>%
        inner_join(df_class %>% distinct(param_set, replicate, pass),
                   by = "param_set") %>%
        pivot_longer(s:N, names_to = "param", values_to = "value")

    df_lhs %>%
        mutate(param = str_replace(param, "_(.*)", "[\\1]")) %>%
        ggplot(aes(x = value, color = pass)) +
        facet_wrap(~ param, nrow = 2, scales = "free",
                   labeller = label_parsed) +
        geom_density() +
        geom_rug(sides = "b") +
        guides(color = FALSE) +
        theme(axis.text.y = element_blank(),
              axis.ticks.y = element_blank())

    ggsave(sprintf("../results/hiv-calipro-param-iter-%d.pdf", iter),
           width = width, height = height)

    ## Old ranges.
    ranges_prev <-
        df_lhs %>%
        select(param, value) %>%
        pivot_wider(names_from = param,
                    values_from = value,
                    values_fn = list) %>%
        reframe(across(everything(), ~ range(.x[[1]])))

    ## Calculate new parameter ranges.
    ranges_new <-
        df_lhs %>%
        nest(.by = param) %>%
        pivot_wider(names_from = param,
                    values_from = data) %>%
        reframe(across(everything(), ~ ads(.x[[1]])))

    ## Print ranges before and after.
    print(bind_rows(before = ranges_prev, after = ranges_new, .id = "range"))

    ## Overlay purple ranges before and after on the density plots.
    df_lhs %>%
        mutate(param = str_replace(param, "_(.*)", "[\\1]")) %>%
        ggplot(aes(x = value, color = pass)) +
        facet_wrap(~ param, nrow = 2, scales = "free",
                   labeller = label_parsed) +
        geom_density() +
        geom_rug(sides = "b") +
        geom_vline(data =
                       ranges_new %>%
                       pivot_longer(everything(),
                                    names_to = "param",
                                    values_to = "value") %>%
                       mutate(param = str_replace(param, "_(.*)", "[\\1]")),
                   color = "purple",
                   aes(xintercept = value)) +
        guides(color = FALSE) +
        theme(axis.text.y = element_blank(),
              axis.ticks.y = element_blank())

    ggsave(sprintf("../results/hiv-calipro-param-ads-iter-%d.pdf", iter),
           width = width, height = height)

    ## TODO: Fit priors to new parameter ranges.
}
