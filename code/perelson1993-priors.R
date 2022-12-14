## Perelson AS, Kirschner DE, De Boer R. Dynamics of HIV infection of CD4+ T
## cells. Math Biosci. 1993 Mar;114(1):81-125.
## doi: 10.1016/0025-5564(93)90043-a. PMID: 8096155.

library(dplyr) # coalesce
library(tidyr) # unnest
library(ggplot2)

## Can't facet stat_function, therefore need to sample from the
## distributions instead.
## https://stackoverflow.com/questions/1376967/using-stat-function-and-facet-wrap-together-in-ggplot2-in-r

df_compact <- tibble::tribble(
    ~param, ~median, ~iqr, ~distr, ~params,
    "s", 10, c(5, 15), "gamma", c(shape = 1.985656, scale = 5.681687),
    "r", 0.03, c(0.02, 0.04), "gamma", c(shape = 4.530347876, scale = 0.006990707),
    "mu[T]", 0.02, c(0.01, 0.03), "gamma", c(shape = 2.10552523, scale = 0.01068658),
    "mu[b]", 0.24, c(0.12, 0.36), "gamma", c(shape = 1.9856561, scale = 0.1363606),
    "mu[V]", 2.4, c(1.2, 3.6), "gamma", c(shape = 1.985657, scale = 1.363605),
    "k[1]", 2.4e-5, c(1.2e-5, 3.6e-5), "gamma", c(shape = 1.985657, scale = 1.363605e-5),
    "T", 562, c(485, 716), "nbinom", c(mu = 562, size = 14.0126),
    "N", 900, c(NA, NA), "nbinom", c(mu = 900, size = 13.5),
    )

## Remap values for plotmath.
distr_ <- c(gamma = "Gamma",
            nbinom = "NegBinom")
params_ <- c(shape = "k",
             scale = "theta",
             mu = "mu",
             size = "r")
df <-
    df_compact %>%
    group_by(param) %>%
    ## See ?plotmath for the label construction.
    mutate(label = paste0(param, "%~%", c(distr_[distr]),
                          " (", paste(params_[names(params[[1]])],
                                     setNames(params[[1]], NULL),
                                     sep = "==",
                                     collapse = ","), ")"),
           x = list(do.call(paste0("q", distr[[1]]),
                            c(p = list(seq(0.1, 0.9, 0.01)),
                              as.list(params[[1]]))))) %>%
    unnest(x) %>%
    mutate(p = do.call(paste0("d", distr[[1]]),
                       c(x = list(x),
                         as.list(params[[1]])))) %>%
    ungroup()
df_vlines <- 
    ## Merge median and IQR to overlay them as vertical lines.
    df %>%
    group_by(label) %>%
    distinct(median, iqr) %>%
    mutate(median = list(median)) %>%
    mutate(fit = list(sort(unlist(c_across(median:iqr))))) %>%
    unnest(fit)
label_ <- setNames(df %>% distinct(param, label) %>% pull(2),
                   df %>% distinct(param, label) %>% pull(1))
## Probability densities fitted to priors.
df %>%
    ggplot(aes(x, p, group = label)) +
    geom_line() +
    geom_vline(data = df_vlines, aes(xintercept = fit)) +
    facet_wrap(~ label, scales = "free",
               labeller = labeller(label = label_parsed)) +
    labs(x = "Parameter value",
         y = "Probability")

ggsave("../results/perelson1993-priors.pdf", width = 7, height = 4, scale = 1.2)
