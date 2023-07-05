## Plot results.
library(RNetCDF)
library(extraDistr)
library(tidyverse)
library(scales)
library(DBI)
library(Rtsne)
library(cowplot)

## Plot the original data and the calibration pass-fail criterion.
##
## Read in the raw experimental data.
times <- seq(0, 15, length.out = 100) # Days.
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
    mutate(err = 8 * ceiling(max(abs_diff))) %>%
    ## Count values cannot be negative.
    mutate(min = pmax(smooth - err, 0L),
           max = pmax(smooth + err, 0L))
plot.margin <- unit(c(30, 5.5, 5.5, 5.5), "points")
(
    plot_calipro <-
        df_boundaries_set %>%
        ggplot(aes(time, obs, fill = var, ymin = min, ymax = max,
                   group = var)) +
        theme_bw() +
        theme(plot.margin = plot.margin) +
        geom_ribbon(alpha = 0.2) +
        geom_point(aes(color = var)) +
        labs(x = "Population", y = "Time (years)") +
        guides(color = "none", fill = "none")
)

## Connect to the most recent database.
db_path <- Sys.glob("../results/lotkavolterra1910calipro-*.sqlite") %>% tail(1)
db <- dbConnect(RSQLite::SQLite(), db_path)

dbListTables(db)

df_sim <- dbReadTable(db, "sim") %>% as_tibble()

## Simulation trajectories for each iteration.
label_10exp <- trans_format("log10", math_format(10^.x))
set.seed(123) # For sample_n()
(
    plot_sim <-
        df_sim %>%
        group_by(iter, replicate) %>%
        mutate(pass_percent = mean(pass)) %>%
        mutate(pass = as.logical(pass)) %>%
        group_by(iter, param_set) %>%
        slice_sample(n = 25, weight_by = pass_percent) %>%
        pivot_longer(cols = c(x, y), names_to = "var", values_to = "values") %>%
        ggplot(aes(time, values, color = pass,
                   group = interaction(iter, param_set, var))) +
        theme_bw() +
        facet_wrap(~iter) +
        geom_line(alpha = 0.1) +
        labs(x = "Time (years)",
             y = "Counts") +
        guides(color = "none") +
        scale_y_log10(labels = label_10exp, limits = c(1, 50),
                      breaks = c(1, 10))
)

## Pass percentage line plot.
(
    plot_sim_pass <-
        df_sim %>%
        group_by(iter) %>%
        summarize(pass = sum(pass, na.rm = TRUE) * 100 / length(pass)) %>%
        ggplot(aes(x = iter, y = pass)) +
        theme_bw() +
        theme(plot.margin = plot.margin) +
        geom_line() +
        labs(x = "CaliPro iteration",
             y = "Pass %") +
        scale_y_continuous(breaks = c(0, 50))
)

(
    plot_top <-
        plot_grid(plot_calipro,
                  plot_sim_pass,
                  nrow = 1,
                  labels = c("(A)", "(B)"),
                  label_x = c(0.01, -0.05),
                  label_fontface = "plain")
)

plot_grid(plot_top,
          plot_sim,
          ncol = 1,
          labels = c("", "(C)"),
          label_fontface = "plain",
          label_x = 0.01,
          label_y = c(0.5, 1),
          rel_heights = c(1, 4))

ggsave("../results/lv-calipro-traj.pdf", width = 7, height = 9.5)

df_prior_fit <-
    dbReadTable(db, "prior_fit") %>%
    as_tibble() %>%
    arrange(iter, param)
df_sim_pass_max <-
    df_sim %>%
    group_by(iter, replicate) %>%
    mutate(pass_percent = mean(pass)) %>%
    ungroup() %>%
    slice_max(pass_percent, n = 1, with_ties = FALSE) %>%
    select(iter, replicate)
args_list <-
    df_prior_fit %>%
    select(-distr) %>%
    pivot_wider(names_from = "parameter", values_from = "value") %>%
    select(-iter, -param) %>%
    as.matrix() %>%
    apply(1, list) %>%
    lapply(unlist, recursive = FALSE) %>%
    lapply(na.omit) %>%
    ## We don't need to know which values were dropped.
    lapply(`attr<-`, "na.action", NULL)
(
    plot_param <-
        df_prior_fit %>%
        distinct(iter, param, distr) %>%
        group_by(iter, param) %>%
        reframe(x = list(do.call(what = str_c("q", distr),
                                 args = c(args_list[[cur_group_id()]],
                                          p = list(c(.03, .97)))))) %>%
        unnest(cols = x) %>%
        mutate(param = str_replace(param, "_(.+)$", "[\\1]")) %>%
        group_by(iter, param) %>%
        reframe(ymin = min(x), ymax = max(x)) %>%
        ungroup() %>%
        left_join(df_sim_pass_max %>%
                  select(-replicate) %>%
                  mutate(pass_max = TRUE)) %>%
        ggplot(aes(x = iter, ymin = ymin, ymax = ymax, color = pass_max)) +
        theme_bw() +
        facet_wrap(~ param, scales = "free_y",
                   labeller = labeller(param = label_parsed)) +
        geom_linerange() +
        labs(x = "CaliPro iteration",
             y = "Parameter value") +
        guides(color = FALSE)
)

(
    plot_sim_pass_max <-
        df_sim %>%
        inner_join(df_sim_pass_max, by = c("iter", "replicate")) %>%
        pivot_longer(cols = c(x, y), names_to = "var", values_to = "values") %>%
        mutate(pass = as.logical(pass)) %>%
        ggplot(aes(time, values, color = pass,
                   group = interaction(iter, param_set, var))) +
        theme_bw() +
        facet_wrap(~var) +
        geom_line(alpha = 0.2) +
        labs(x = "Time (years)",
             y = "Counts") +
        guides(color = "none") +
        scale_y_log10(labels = label_10exp, limits = c(0.01, 50),
                      breaks = c(0.1, 1, 10))
)

df_lhs <- dbReadTable(db, "lhs") %>% as_tibble()

df_tsne_pre <-
    df_lhs %>%
    select(-q) %>%
    pivot_wider(names_from = param, values_from = value)
mat_tsne_pre <-
    df_tsne_pre[, -c(1:4)] %>%
    as.matrix()
set.seed(123) # Rtsne is stochastic
system.time({
    tsne <- Rtsne(mat_tsne_pre)
})

(
    plot_tsne_point <-
        df_tsne_pre %>%
        mutate(tsne_x = tsne$Y[, 1],
               tsne_y = tsne$Y[, 2],
               pass = as.logical(pass)) %>%
        ggplot(aes(x = tsne_x, y = tsne_y, color = pass)) +
        theme_bw() +
        geom_point(alpha = 0.3, size = 1) +
        guides(color = "none") +
        coord_fixed() +
        theme(axis.text.x = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks.x = element_blank(),
              axis.ticks.y = element_blank(),
              axis.title.x = element_blank(),
              axis.title.y = element_blank())
)

## 2D density plots of pass and fail simulations.
(
    plot_tsne_density <-
        df_tsne_pre %>%
        mutate(tsne_x = tsne$Y[, 1],
               tsne_y = tsne$Y[, 2],
               pass = case_when(
                   pass == 1L ~ "Pass",
                   pass == 0L ~ "Fail")) %>%
        drop_na() %>%
        ggplot(aes(x = tsne_x, y = tsne_y)) +
        theme_bw() +
        facet_wrap(~pass) +
        geom_density_2d_filled() +
        guides(fill = "none") +
        coord_fixed() +
        theme(axis.text.x = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks.x = element_blank(),
              axis.ticks.y = element_blank(),
              axis.title.x = element_blank(),
              axis.title.y = element_blank())
)

(
    plot_bottom <-
        plot_grid(plot_tsne_point,
                  plot_tsne_density,
                  nrow = 1,
                  labels = c("(C)", "(D)"),
                  label_fontface = "plain",
                  label_x = c(0, -0.1),
                  rel_widths = c(1, 2),
                  align = "v",
                  axis = "b")
)
plot_grid(plot_param,
          plot_sim_pass_max,
          plot_bottom,
          nrow = 3,
          labels = c("(A)", "(B)", ""),
          label_fontface = "plain")

ggsave("../results/lv-calipro-param.pdf", width = 7, height = 9)

dbDisconnect(db)
