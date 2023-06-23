library(RNetCDF)
library(deSolve)
library(tidyverse)
library(scales)  # trans_format, math_format
library(cowplot) # plot_grid

initial_values <- c(
    x = 10.0, # Prey population.
    y = 5.0   # Predator population.
)
parameters <- c(
    a = 1.0,  # Prey growth rate.
    b = 0.5,  # Prey death rate.
    c = 1.50, # Predator death rate.
    d = 0.75  # Predator growth rate.
)
times <- seq(0, 15, length.out = 100) # Days.
## Lotka-Volterra model (predator-prey).
lv <- function(t, initial_values, parameters) {
    ## Conveniently access members by name by constructing an environment from
    ## the data using with() and as.list() as shown in the deSolve package
    ## vignette of the Lorenz equation.
    with(as.list(c(initial_values, parameters)), {
        ## Rates of change.
        dx.dt <- +a*x - b*x*y
        dy.dt <- -c*y + d*b*x*y
        ## Return the rates of change.
        list(c(dx.dt, dy.dt))
    })
}

nc_prior <- "lotkavolterra1910-01_output-idata_untrained.nc"
nc_posterior <- "lotkavolterra1910-01_output-idata_lv.nc"

df_prior_new <- function(name, group) {
    nc <- open.nc(name)
    ## print.nc(nc)
    grp <- grp.inq.nc(nc, group)$self
    ## print.nc(grp)
    df <-
        tibble(a = var.get.nc(grp, "a"),
               b = var.get.nc(grp, "b")) %>%
        mutate(draw = row_number(),
               out = map2(a, b,
                          ~ ode(y = initial_values,
                                times = times,
                                func = lv,
                                parms = replace(parameters, 1:2, c(.x, .y))) %>%
                              `class<-`("double") %>% as_tibble())) %>%
        unnest(out) %>%
        gather(key = var, value = count, x, y)
    close.nc(nc)
    df
}
df_prior <- df_prior_new(nc_prior, "prior")

df_prior %>%
    mutate(var = case_when(var == "x" ~ "Prey",
                           var == "y" ~ "Predator")) %>%
    ggplot(aes(time, count,
               color = var,
               group = interaction(var, draw))) +
    geom_line(alpha = 0.3) +
    scale_y_log10() +
    labs(x = "Time", y = "Population", color = "") +
    theme_bw() +
    theme(legend.position = "top")

df_prior %>%
    distinct(draw, a, b) %>%
    gather("param", "value", -draw) %>%
    ## Use unicode characters for greek letters:
    ## https://stackoverflow.com/a/27741196
    mutate(param = case_when(param == "a" ~ "\u03B1",
                             param == "b" ~ "\u03B2")) %>%
    ## Reverse order.
    mutate(param = factor(param, levels = c("\u03B2", "\u03B1"))) %>%
    ggplot(aes(value, param)) +
    geom_boxplot(outlier.shape = NA) +
    geom_point(shape = 73, size = 5) +
    theme_bw() +
    labs(x = "Value", y = "Parameter")

df_data_new <- function(name, group) {
    nc <- open.nc(name)
    grp <- grp.inq.nc(nc, group)$self
    mat <- var.get.nc(grp, "sim")
    dimnames(mat) <- list(c("x", "y"), times)
    t(mat) %>%
        as_tibble(rownames = "time") %>%
        mutate(time = as.numeric(time)) %>%
        gather(key = var, value = count, x, y)
}
df_data <- df_data_new(nc_posterior, "observed_data")

df_data %>%
    mutate(var = case_when(var == "x" ~ "Prey",
                           var == "y" ~ "Predator")) %>%
    ggplot(aes(time, count,
               color = var,
               group = var)) +
    geom_point() +
    scale_y_log10() +
    labs(x = "Time", y = "Population", color = "") +
    theme_bw() +
    theme(legend.position = "top")

df_posterior_new <- function(name, group) {
    nc <- open.nc(name)
    grp <- grp.inq.nc(nc, group)$self
    df <-
        tibble(a = var.get.nc(grp, "a")[, 1],
               b = var.get.nc(grp, "b")[, 1]) %>%
        sample_n(75) %>%
        mutate(draw = row_number(),
               out = map2(a, b,
                          ~ ode(y = initial_values,
                                times = times,
                                func = lv,
                                parms = replace(parameters, 1:2, c(.x, .y))) %>%
                              `class<-`("double") %>% as_tibble())) %>%
        unnest(out) %>%
        gather(key = var, value = count, x, y)
    close.nc(nc)
    df
}
df_posterior <- df_posterior_new(nc_posterior, "posterior")

df_posterior %>%
    mutate(var = case_when(var == "x" ~ "Prey",
                           var == "y" ~ "Predator")) %>%
    ggplot(aes(time, count,
               color = var,
               group = interaction(var, draw))) +
    geom_line(alpha = 0.3) +
    scale_y_log10() +
    labs(x = "Time", y = "Population", color = "") +
    theme_bw() +
    theme(legend.position = "top")

df_posterior %>%
    distinct(draw, a, b) %>%
    gather("param", "value", -draw) %>%
    ## Use unicode characters for greek letters:
    ## https://stackoverflow.com/a/27741196
    mutate(param = case_when(param == "a" ~ "\u03B1",
                             param == "b" ~ "\u03B2")) %>%
    ## Reverse order.
    mutate(param = factor(param, levels = c("\u03B2", "\u03B1"))) %>%
    ggplot(aes(value, param)) +
    geom_boxplot(outlier.shape = NA) +
    geom_point(shape = 73, size = 5) +
    theme_bw() +
    labs(x = "Value", y = "Parameter")

label_10exp <- trans_format("log10", math_format(10^.x))
ggtraj <-
    bind_rows("Uncalibrated (n = 75)" = df_prior,
              "Noisy data (n = 1)" = df_data %>% mutate(draw = 1),
              "Calibrated (n = 75)" = df_posterior,
              .id = "sample") %>%
    mutate(var = case_when(var == "x" ~ "Prey",
                           var == "y" ~ "Predator"),
           sample = fct_inorder(sample)) %>%
    ggplot(aes(time, count,
               color = var,
               group = interaction(var, draw))) +
    facet_wrap(~ sample) +
    geom_line(alpha = 0.3) +
    scale_y_log10(labels = label_10exp) +
    labs(x = "Time", y = "Population", color = "") +
    theme_bw() +
    theme(legend.position = "top")
ggsave("../results/lv-traj.pdf", plot = ggtraj, width = 7, height = 4)

ggparam <-
    bind_rows("Uncalibrated" = df_prior,
              "Calibrated" = df_posterior,
              .id = "sample") %>%
    distinct(sample, draw, a, b) %>%
    gather("param", "value", a, b) %>%
    mutate(param = fct_rev(param),
           sample = fct_inorder(sample)) %>%
    ggplot(aes(value, param)) +
    facet_wrap(~ sample) +
    geom_boxplot(outlier.shape = NA) +
    geom_point(shape = 73, size = 5) +
    theme_bw() +
    labs(x = "Value", y = "Parameter") +
    ## ## PDF export doesn't allow using UTF-8 characters for greek
    ## ## letters, so we have to use the native R rendering method.
    scale_y_discrete(label = parse(text = c(expression(beta),
                                            expression(alpha))))
ggsave("../results/lv-param.pdf", plot = ggparam, width = 7, height = 2.5)

plot_grid(ggtraj, ggparam, align = "v", labels = "AUTO")
