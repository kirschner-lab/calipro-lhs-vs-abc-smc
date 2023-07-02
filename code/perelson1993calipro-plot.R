## Plot results.
library(tidyverse)
library(DBI)
library(Rtsne)
library(cowplot)

## Connect to the most recent database.
db_path <- Sys.glob("../results/perelson1993calipro-*.sqlite") %>% tail(1)
db <- dbConnect(RSQLite::SQLite(), db_path)

dbListTables(db)

df_sim <- dbReadTable(db, "sim") %>% as_tibble()

## Simulation trajectories for each iteration.
set.seed(123) # For sample_n()
(
    plot_sim <-
        df_sim %>%
        mutate(pass = as.logical(pass)) %>%
        drop_na() %>%
        group_by(iter, param_set) %>%
        sample_n(50) %>%
        ggplot(aes(time_years, T_sum, color = pass,
                   group = interaction(iter, param_set))) +
        theme_bw() +
        facet_wrap(~iter) +
        geom_line(alpha = 0.2) +
        labs(x = "Time (years)",
             y = bquote(sum(CD4^"+" ~ "T-cells") ~ (counts))) +
        guides(color = "none") +
        scale_x_continuous(breaks = c(0, 5, 10)) +
        scale_y_continuous(breaks = c(0, 1000, 2000), limits = c(0, 2400))
)

## Pass percentage line plot.
(
    plot_sim_pass <-
        df_sim %>%
        group_by(iter) %>%
        summarize(pass = sum(pass, na.rm = TRUE) * 100 / length(pass)) %>%
        ggplot(aes(x = iter, y = pass)) +
        theme_bw() +
        geom_line() +
        labs(x = "CaliPro iteration",
             y = "Pass %") +
        scale_y_continuous(breaks = c(0, 50, 100))
)

plot_grid(plot_sim_pass,
          plot_sim,
          ncol = 1,
          labels = c("(A)", "(B)"),
          label_fontface = "plain",
          label_x = 0.03,
          label_y = c(0.5, 1),
          rel_heights = c(1, 8))

ggsave("../results/hiv-calipro-traj.pdf", width = 7, height = 9.5)

## NA values don't show up in pass-fail.
df_sim %>%
    group_by(iter, replicate) %>%
    summarize(pass_1 = sum(pass, na.rm = TRUE) * 100 / length(pass),
              pass_0 = sum(pass == 0L, na.rm = TRUE) * 100 / length(pass),
              pass_NA = sum(is.na(pass)) * 100 / length(pass)) %>%
    pivot_longer(starts_with("pass"),
                 names_to = "output",
                 values_to = "pass") %>%
    mutate(output =
               str_remove(output, "^pass_") %>%
               type.convert(as.is = TRUE) %>%
               as.logical()) %>%
    ## filter(replicate == 1) %>%
    ggplot(aes(x = iter, y = pass, fill = output,
               group = interaction(iter, replicate))) +
    theme_bw() +
    geom_bar(stat = "identity")

## Finding NA values.
df_sim %>%
    filter(if_any(T_ui:V, is.na)) %>%
    group_by(iter, replicate) %>%
    count(param_set) %>%
    print(n = Inf)

## Plot simulations with NA values.
df_sim %>%
    filter(if_any(T_ui:V, is.na)) %>%
    distinct(iter, param_set, replicate) %>%
    inner_join(df_sim, by = c("iter", "param_set", "replicate")) %>%
    mutate(pass = as.logical(pass)) %>%
    ggplot(aes(time_years, T_sum, color = pass,
               group = interaction(iter, param_set))) +
    theme_bw() +
    facet_wrap(~iter) +
    geom_line(alpha = 0.5) +
    ylim(x = c(0, 2400)) +
    labs(x = "Time (years)",
         y = bquote(sum(CD4^"+" ~ "T-cells") ~ (counts))) +
    guides(color = "none")

df_param <- dbReadTable(db, "param_fit") %>% as_tibble()

df_prior_fit <-
    dbReadTable(db, "prior_fit") %>%
    as_tibble() %>%
    arrange(iter, param)
args_list <-
    df_prior_fit %>%
    select(-distr, -rescale_multiplier) %>%
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
        ggplot(aes(x = iter, ymin = ymin, ymax = ymax)) +
        theme_bw() +
        facet_wrap(~ param, scales = "free_y",
                   labeller = labeller(param = label_parsed)) +
        geom_linerange() +
        labs(x = "CaliPro iteration",
             y = "Parameter value")
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
                  labels = c("(B)", "(C)"),
                  label_fontface = "plain",
                  label_x = c(0, -0.1),
                  rel_widths = c(1, 2),
                  align = "v",
                  axis = "b")
)
plot_grid(plot_param,
          plot_bottom,
          nrow = 2,
          labels = c("(A)", ""),
          label_fontface = "plain",
          rel_heights = c(3.5, 2))

ggsave("../results/hiv-calipro-param.pdf", width = 7, height = 8)

dbDisconnect(db)
