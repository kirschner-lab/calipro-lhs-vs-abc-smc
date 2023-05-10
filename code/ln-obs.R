## Convert flowcytometry data into aggregated counts.

library(readxl)
library(tidyverse)
library(modelr)

df <-
    read_excel("~/Downloads/ICS_Thoracic_LN_Ganchua_PLOSpath2018_pm_jlf.xlsx",
               sheet = 1L) %>%
    select(c(2, 5:7, 10:ncol(.))) %>%
    filter(LN_Gran == "GRAN",
           CFU > 0,
           Stims == "EC",
           ! is.na(Cell_Count),
           ## ! is.na(Effacement),
           ## Effacement != "none"
           ) %>%
    select(-Stims, -LN_Gran, -"CFU log 10+1", -"CFU status", -Histo) %>%
    ## mutate(
    transmute(
        nhp_id = NHP,
        days = Days_PI,
        E4 =
            Cell_Count * 2 *
            Live / 100 *
            `Live/Lymphocytes` / 100 *
            CD3 / 100 *
            CD4 / 100 *
            `CD4/IFNg` / 100,
        M = 
            Cell_Count * 2 *
            Live / 100 *
            `Live/Myeloid` / 100 *
            `Live/Myeloid/CD11b` / 100,
        B =
            CFU * 2,
        effacement = Effacement)

print(df, n = Inf)

write_csv(df, "../data/ICS_Thoracic_LN_Ganchua_PLOSpath2018_pm_jlf_aggregated.csv")

## Doesn't look like imputation would work for Macrophage (M) data.
breaks <- scales::trans_breaks("log10", function(x) 10^x, n <- 6)
labels <- scales::trans_format("log10", scales::math_format(10^.x))
limits <- c(10^1, 10^6)
df %>%
    rename(Bacteria = B,
           Macrophages = M,
           "Effector CD4+ T-cells" = E4) %>%
    gather("variable", "count", -nhp_id, -days, -effacement) %>%
    ggplot(aes(days, count, color = variable)) +
    labs(x = "Day", y = "Count", color = "Cell") +
    geom_point() +
    scale_y_log10(breaks = breaks, labels = labels, limits = limits)

ggsave("../results/ln-plot-points.pdf",
       width = width, height = height)

breaks <- scales::trans_breaks("log10", function(x) 10^x, n <- 6)
labels <- scales::trans_format("log10", scales::math_format(10^.x))
limits <- c(10^1, 10^6)
df %>%
    rename(Bacteria = B,
           Macrophages = M,
           "Effector CD4+ T-cells" = E4) %>%
    gather("variable", "count", -nhp_id, -days, -effacement) %>%
    mutate(days = factor(days)) %>%
    ggplot(aes(days, count, color = variable)) +
    labs(x = "Day", y = "Count", color = "Cell") +
    ## Preserve width of boxplots even with missing Macrophage data per
    ## https://stackoverflow.com/a/52216686
    geom_boxplot(position = position_dodge(preserve = "single")) +
    ## To overlay points on the boxplots, one has to use position_jitterdodge()
    ## per https://stackoverflow.com/a/50179098
    geom_point(position = position_jitterdodge(jitter.width = 0)) +
    scale_y_log10(breaks = breaks, labels = labels, limits = limits)

## Plot sizes in inches for PowerPoint.
width <- 13.33
height <- 6.06

ggsave("../results/ln-plot-boxplots.pdf",
       width = width, height = height)

breaks <- scales::trans_breaks("log10", function(x) 10^x, n <- 6)
labels <- scales::trans_format("log10", scales::math_format(10^.x))
limits <- c(10^1, 10^6)
df %>%
    rename(Bacteria = B,
           Macrophages = M,
           "Effector CD4+ T-cells" = E4) %>%
    mutate(gran_id = row_number()) %>%
    gather("variable", "count", -nhp_id, -gran_id, -days, -effacement) %>%
    mutate(days = factor(days)) %>%
    ggplot(aes(days, count, color = variable)) +
    facet_wrap(~gran_id) +
    theme(strip.background = element_blank(),
          strip.text.x = element_blank(),
          axis.text.x = element_text(angle = 45, vjust = 0.7)) +
    labs(x = "Day", y = "Count", color = "Cell") +
    geom_point() +
    scale_y_log10(labels = labels, limits = limits)

ggsave("../results/ln-plot-indiv.pdf",
       width = width, height = height)

breaks <- scales::trans_breaks("log10", function(x) 10^x, n <- 6)
labels <- scales::trans_format("log10", scales::math_format(10^.x))
limits <- c(10^1, 10^6)
df %>%
    rename(Bacteria = B,
           Macrophages = M,
           "Effector CD4+ T-cells" = E4) %>%
    gather("variable", "count", -nhp_id, -days, -effacement) %>%
    mutate(days = factor(days)) %>%
    ggplot(aes(days, count, color = variable)) +
    facet_wrap(~nhp_id) +
    theme(strip.background = element_blank(),
          strip.text.x = element_blank(),
          axis.text.x = element_text(angle = 45, vjust = 0.7)) +
    labs(x = "Day", y = "Count", color = "Cell") +
    geom_point() +
    scale_y_log10(labels = labels, limits = limits)

ggsave("../results/ln-plot-animal.pdf",
       width = width, height = height)

## Missing data in granulomas.
df %>%
    select(E4, M, B) %>%
    mutate(across(.fns = ~!is.na(.x))) %>%
    group_by(E4, M, B) %>%
    count()

## Missing data in animals.
df %>%
    group_by(nhp_id) %>%
    select(E4, M, B) %>%
    summarize(across(.fns = ~sum(.x, na.rm = TRUE))) %>%
    mutate(across(-1, .fns = ~. > 0)) %>%
    group_by(E4, M, B) %>%
    count()

## Ratio approach.
df %>%
    ggplot(aes(E4, M)) +
    geom_point() +
    geom_smooth(method = lm) +
    coord_fixed() +
    scale_x_log10(labels = labels) +
    scale_y_log10(labels = labels)

df %>%
    transmute(days, E4_M_ratio = E4 / M) %>%
    ggplot(aes(days, E4_M_ratio)) +
    labs(x = "Day", y = "Ratio of Effector CD4+ T-cells to Macrophages") +
    geom_point() +
    ## geom_smooth(method = lm) +
    geom_smooth() +
    scale_y_log10(labels = labels)        

ggsave("../results/ln-plot-ratio.pdf",
       width = width, height = height)

## Try predicting missing T-cell and Macrophage counts using the ratio model. 
dfr <- transmute(df, days, E4, M, ratio = E4 / M)
ratio_model <- loess(ratio ~ days - 1, data = dfr)
dfr %>%
    add_predictions(ratio_model) %>%
    mutate(E4_pred = ifelse(is.na(E4), M * pred, E4),
           M_pred = ifelse(is.na(M), E4 / pred, M)) %>%
    print(n = Inf)
