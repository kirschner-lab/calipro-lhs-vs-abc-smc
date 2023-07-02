library(RNetCDF)
library(deSolve)
library(lhs)
library(extraDistr) # {d,p,q}hnorm
library(plyr)       # round_any
library(tidyverse)
library(broom) # tidy.optim
library(DBI)
library(abind)

initial_values <- c(
    x = 10.0, # Prey population.
    y = 5.0   # Predator population.
)
parameters <- c(
    a = 1.0,  # Prey growth rate.
    b = 1.0   # Prey death rate (uncalibrated).
)
times <- seq(0, 15, length.out = 100) # Days.

## Lotka-Volterra model (predator-prey).
lv <- function(t, initial_values, parameters) {
    ## Conveniently access members by name by constructing an environment from
    ## the data using with() and as.list() as shown in the deSolve package
    ## vignette of the Lorenz equation.
    with(as.list(c(initial_values, parameters)), {
        ## Fixed parameters.
        c <- 1.50 # Predator death rate.
        d <- 0.75 # Predator growth rate.
        ## Rates of change.
        dx.dt <- +a*x - b*x*y
        dy.dt <- -c*y + d*b*x*y
        ## Return the rates of change.
        list(c(dx.dt, dy.dt))
    })
}

## Read in the raw experimental data.
df_observed_new <- function(name, group) {
    nc <- open.nc(name)
    grp <- grp.inq.nc(nc, group)$self
    mat <- var.get.nc(grp, "sim")
    dimnames(mat) <- list(c("x", "y"), times)
    t(mat) %>%
        as_tibble(rownames = "time") %>%
        mutate(time = as.numeric(time)) %>%
        gather(key = var, value = count, x, y)
}
nc_observed <- "lotkavolterra1910-01_output-idata_lv.nc"
df_observed <- df_observed_new(nc_observed, "observed_data")

## Set the experimental data boundaries.
##
## Smooth the data to apply [ymin, ymax] tolerance bars to calibrate to.
df_boundaries_set <-
    df_observed %>%
    group_by(var) %>%
    reframe(smooth = ksmooth(time, count, kernel = "normal",
                             x.points = time)) %>%
    mutate(dim = names(smooth)) %>%
    pivot_wider(names_from = dim, values_from = smooth) %>%
    unnest(c(x, y)) %>%
    dplyr::rename(time = x,
                  count = y) %>%
    bind_rows(smooth = ., obs = df_observed, .id = "type") %>%
    group_by(var) %>%
    pivot_wider(names_from = type, values_from = count) %>%
    mutate(abs_diff = abs(obs - smooth)) %>%
    mutate(err = 2 * ceiling(max(abs_diff))) %>%
    ## Count values cannot be negative.
    mutate(min = pmax(smooth - err, 0L),
           max = pmax(smooth + err, 0L))

df_boundaries_set %>%
    ggplot(aes(time, obs, fill = var, ymin = min, ymax = max,
               group = var)) +
    geom_ribbon(alpha = 0.5) +
    geom_point(aes(color = var))

df_boundaries <-
    df_boundaries_set %>%
    select(var, time, min, max) %>%
    pivot_wider(names_from = var,
                values_from = c(min, max))

## Objective function to minimize to fit the distributions, that returns the
## distance between the original and estimated values.
fit_error <- function(distr, params, q, x) {
    with(as.list(params, q, x), {
        x_est <-
            switch(distr,
                   hnorm = qhnorm(q, sigma = sigma)) # nolint
        mat <- matrix(c(x, x_est), nrow = length(x))
        c(dist(t(mat)))
    })
}
## Initial, relaxed prior distributions.
df_prior_fit <- tibble::tribble(
    ~param, ~distr, ~parameter, ~value,
    "a", "hnorm", "sigma", 1,
    "b", "hnorm", "sigma", 1,
)
## For now, sample ranges from the prior instead of fitting to ranges.
## Fitting to ranges will happen in the calibration loop.
args_list <-
    df_prior_fit %>%
    pivot_wider(names_from = parameter, values_from = value) %>%
    select(-param, -distr) %>%
    as.matrix() %>%
    apply(1, list) %>%
    lapply(unlist, recursive = FALSE) %>%
    lapply(na.omit) %>%
    ## We don't need to know which values were dropped.
    lapply(`attr<-`, "na.action", NULL)
df_prior_ranges <-
    df_prior_fit %>%
    select(-parameter, -value) %>%
    group_by(param) %>%
    mutate(q = list(c(.03, .97))) %>%
    reframe(distr = distr,
            x = do.call(str_c("q", distr), c(args_list[[cur_group_id()]],
                                             p = q)),
            q = unlist(q)) %>%
    group_by(param) %>%
    reframe(distr = head(distr, 1),
            q = list(q),
            x = list(x))

## Calibration adjustment function: Alternative Density Subtraction (ADS).
##
## Joslyn LR, Kirschner DE, Linderman JJ. CaliPro: A Calibration Protocol That
## Utilizes Parameter Density Estimation to Explore Parameter Space and
## Calibrate Complex Biological Models. Cell Mol Bioeng. 2020 Sep
## 15;14(1):31-47. doi: 10.1007/s12195-020-00650-z
ads <- function(df, n = 512) { # nolint start
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
} # nolint end

## Initialize the database to store the calibration data.
timestamp <- format(Sys.time(), "%Y%m%dT%H%M%S")
db_path <- str_c("../results/lotkavolterra1910calipro-",
                 timestamp,
                 ".sqlite")
db <- dbConnect(RSQLite::SQLite(), db_path)
header_with_iter <- function(df) {
    df %>% mutate(iter = 0L, .before = 1) %>% slice()
}
headers <- list(
    prior_ranges = header_with_iter(df_prior_ranges %>% unnest(c(q, x))),
    prior_fit = header_with_iter(df_prior_fit),
    lhs = tibble(iter = integer(),
                 param_set = integer(),
                 replicate = integer(),
                 param = character(),
                 value = numeric(),
                 q = numeric(),
                 pass = logical()),
    sim = tibble(iter = integer(),
                 param_set = integer(),
                 replicate = integer(),
                 time = numeric(),
                 time_years = numeric(),
                 x = numeric(),
                 y = numeric(),
                 pass = logical())
)
for (i in seq_along(headers)) {
    dbWriteTable(db, names(headers)[i], headers[[i]])
}

## Insert records prefixed with iteration.
db_write_iter <- function(table, df) {
    dbWriteTable(db,
                 name = table,
                 value = df %>% mutate(iter = iter, .before = 1),
                 append = TRUE)
}

## Calibration loop.
n_iter <- 50
set.seed(123)
seeds <- sample(1:1000, size = n_iter)
iter <- 1

for (iter in 1:n_iter) {
    ## Store the priors into the database.
    db_write_iter("prior_ranges", df_prior_ranges %>% unnest(c(q, x)))
    db_write_iter("prior_fit", df_prior_fit)

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
    ## Final LHS samples.
    lhs <- array(-1, dim = dim(lhs_cdf), dimnames = dimnames(lhs_cdf))
    for (param_ in df_prior_fit %>% distinct(param) %>% pull()) {
        call <-
            with(df_prior_fit %>% filter(param == param_), {
                distr <-
                    df_prior_ranges %>%
                    filter(param == param_) %>%
                    pull(distr)
                list(what = str_c("q", distr),
                     args = c(list(p = lhs_cdf[, param_]),
                              setNames(value, parameter)))
            })
        lhs[, param_] <- do.call(call[["what"]], call[["args"]])
    }

    ## Run the simulations with the LHS parameter samples.
    out <- list()
    i <- 1
    system.time({
        for (i in seq_len(nrow(lhs))) {
            out[[i]] <- ode(y = initial_values,
                            times = times,
                            func = lv,
                            parms = lhs[i, ])
            i <- i + 1
        }
    })

    ## Gather the output from the LHS parameter samples.
    df_sim <-
        lapply(out, as.data.frame) %>%
        bind_rows(.id = "param_set") %>%
        as_tibble() %>%
        mutate(param_set = as.integer(param_set),
               replicate = ceiling(param_set / n_samples_per_replicate))

    ## Classify model outputs as pass or fail using the data boundaries.
    df_class <-
        df_sim %>%
        group_by(param_set) %>%
        dplyr::mutate(id = 1:n()) %>%
        ## Time has slight differences that causes inner_join() to fail, so use
        ## the row number instead.
        left_join(df_boundaries %>%
                  select(-time) %>%
                  dplyr::mutate(id = 1:n()),
                  by = "id") %>%
        group_by(param_set) %>%
        mutate(pass = all(x >= min_x & x <= max_x &
                          y >= min_y & y <= max_y)) %>%
        select(-id, -starts_with("min_"), -starts_with("max_"))

    message(sprintf("%2d: Overall pass rate across all replicates: %.1f%%",
                    iter,
                    df_class %>% pull(pass) %>% mean() * 100))

    message(do.call(sprintf,
                    c(fmt = str_c("    Pass rate per replicate:",
                                  str_c(rep(" %.0f%%", n_replicates),
                                        collapse = "")),
                      as.list(df_class %>%
                              group_by(replicate) %>%
                              dplyr::summarize(pass = mean(pass) * 100) %>%
                              pull(pass)))))

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

    ## Store the samples and simulations into the database.
    db_write_iter("lhs",
                  df_lhs %>%
                  mutate(value_is_cdf = ifelse(value_is_cdf, "q", "value")) %>%
                  pivot_wider(names_from = value_is_cdf,
                              values_from = value))
    db_write_iter("sim", df_class)

    ranges_prev_src <-
        df_lhs %>%
        mutate(value_is_cdf = ifelse(value_is_cdf, "cdf", "before")) %>%
        pivot_wider(id_cols = c(param_set, replicate, param),
                    names_from = "value_is_cdf",
                    values_from = "value") %>%
        group_by(param) %>%
        arrange(param, before) %>%
        slice(c(1, n())) %>%
        mutate(bound = c("lower", "upper"))
    ranges_prev <-
        ranges_prev_src %>%
        select(param, before, bound) %>%
        pivot_wider(id_cols = bound, names_from = param,
                    values_from = before) %>%
        select(-bound)

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
        pivot_longer(cols = ! c(range, bound), names_to = "param",
                     values_to = "value") %>%
        pivot_wider(id_cols = c("param", "bound"), names_from = "range") %>%
        mutate(changed = before != after)

    ## If no ranges changed, end calibration.
    if (! any(pull(ranges_changed, changed))) {
        message(sprintf(c("No parameter ranges changed in iteration %d. ",
                          "Ended calibration loop."), iter))
        break
    }

    ## Subset to parameters with new ranges.
    df_prior_ranges_new <-
        ranges_changed %>%
        arrange(param, bound) %>%
        ## Add quantile from old ranges.
        left_join(ranges_prev_src, by = c("param", "bound", "before")) %>%
        ## Fit new ranges.
        rename(x = after, q = cdf) %>%
        group_by(param) %>%
        ## Only re-fit ranges that have changed.
        filter(any(changed)) %>%
        select(-bound, -before, -changed, -param_set, -replicate) %>%
        reframe(q = list(q), x = list(x)) %>%
        left_join(df_prior_ranges %>% select(param, distr),
                  by = "param")

    ## Fit priors to new parameter ranges.
    df_prior_fit_new <-
        df_prior_ranges_new %>%
        group_by(param) %>%
        mutate(par_init = case_when(distr == "hnorm" ~
                                        list(c(sigma = 1))),
               ## All the parameters of the gamma and negative-binomial
               ## distributions must be greater than zero.
               lower = list(setNames(rep(0 + .Machine$double.eps,
                                         length(unlist(par_init))),
                                     names(unlist(par_init)))),
               ## The probability parameter of the negative-binomial
               ## distribution must be less than or equal to 1.
               upper = list(setNames(rep(Inf,
                                         length(unlist(par_init))),
                                     names(unlist(par_init)))),
               upper = ifelse("prob" %in% names(unlist(par_init)), {
                   upper[[1]]["prob"] <- 1
                   upper
               }, upper)) %>%
        group_by(param, distr) %>%
        reframe(tidy(optim(par = unlist(par_init),
                           fn = function(...) fit_error(distr, ...),
                           q = unlist(q),
                           x = unlist(x),
                           lower = unlist(lower),
                           upper = unlist(upper),
                           method = "L-BFGS-B")))

    ## Carry over any missing parameter ranges from the previous run.
    df_prior_ranges_new <-
        bind_rows(df_prior_ranges_new,
                  anti_join(df_prior_ranges, df_prior_ranges_new, by = "param"))
    df_prior_fit_new <-
        bind_rows(df_prior_fit_new,
                  anti_join(df_prior_fit, df_prior_fit_new, by = "param"))

    df_prior_ranges <- df_prior_ranges_new
    df_prior_fit <- df_prior_fit_new
}

## Close the database.
dbDisconnect(db)
