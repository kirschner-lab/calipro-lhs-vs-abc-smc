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
library(dplyr)
library(tidyr)
library(ggplot2)

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
    N = VARIES,    # Number of free virus produced by lysing a CD4+ cell.
    theta = 1      # mm^{-3}, Viral concentration need to decrease s to s/2.
)
times <- seq(0, 10, by = 0.0001) * 365 # years converted to days.

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

## Model 2: decreasing 's'.
hiv_s_decr <- function(t, initial_values, parameters) {
    with(as.list(c(initial_values, parameters)), {
        s = 
    })
}

Ns <- c(774, 1500, 2000, 3000)
out <- list()
i <- 1
for (N in Ns) {
    parameters["N"] <- N
    out[[i]] <- ode(y = initial_values, times = times, func = hiv_s_const,
                    parms = parameters)
    i <- i + 1
}
names(out) <- Ns

diagnostics(out[[1]])
summary(out[[1]])

df <-
    lapply(out, as.data.frame) %>%
    setNames(Ns) %>%
    bind_rows(.id = "N") %>%
    as_tibble() %>%
    mutate(N = as.integer(N)) %>%
    group_by(N) %>%
    mutate(T_sum = T + T_li + T_ai,
           time_years = time / 365)

df_long <-
    df %>%
    gather("variable", "value", -time, -N) %>%
    mutate(time_years = time / 365) %>%
    mutate(variable = case_when(
               variable == "T_ai" ~ "CD4+ actively infected",
               variable == "T_li" ~ "CD4+ latently infected",
               variable == "T" ~ "CD4+ uninfected",
               variable == "T_sum" ~ "CD4+ uninfected & infected",
               variable == "V" ~ "Free HIV"))

## Useless plot, because T is steady state.  Need to run the model multiple
## times for different values of N, and then plot N against the final steady
## state of T.
## 
## TODO: Read my notes in Intro. to Mathematical Modeling to see what linear
## algebra operations I used to analytically find the steady state of a system
## of equations, if such a solution exists.

df %>%
    filter(time == max(time)) %>%
    ggplot(aes(N, T)) +
    geom_line()

## Plot variables of model 1.
df_long %>%
    filter(N == 774L) %>%
    ggplot(aes(time_years, value, group = variable)) +
    facet_wrap(~ variable, scales = "free_y") +
    scale_x_log10() +
    ## coord_cartesian(xlim = c(NA, 0.01)) +
    geom_line() +
    ggtitle("HIV model 1 with constant s")

## Reproduce figure 5 using model 1.
(plot_fig5 <-
     df %>%
     ggplot(aes(time_years, T_sum, group = N)) +
     geom_line() +
     geom_label(data = df %>% filter(time == max(time) / 2),
                aes(label = paste("N =", N))) +
     scale_x_continuous(breaks = 0:10,
                        sec.axis = dup_axis(labels = NULL, name = NULL)) +
     scale_y_continuous(breaks = seq(650, 1000, 25),
                        labels = c(650, NA, 700, NA, 750, NA, 800, NA, 850, NA,
                                   900, NA, 950, NA, 1000) %>%
                            as.character() %>% replace_na(""),
                        sec.axis = dup_axis(labels = NULL, name = NULL)) +
     labs(x = "Years", y = "Total T4 cells") +
     theme_classic() +
     theme(panel.border = element_rect(fill = NA, size = 1),
           axis.ticks.length.x = unit(-.25, "cm"),
           axis.ticks.length.y = unit(-.25, "cm"),
           axis.line = element_blank()
           ))

ggsave("~/Sync/lab-reports/2022-11-09/img/plot_fig5.png",
       width = 5, height = 5)

## Reproduce figure 11 using model 1 for now (instead of model 2) to be able to overlay viral data.
df %>%
    ggplot(aes(time_years, V, group = N)) +
    geom_line() +
    geom_label(data = df %>% filter(time == max(time) / 2),
               aes(label = paste("N =", N))) +
    scale_x_continuous(breaks = 0:10) +
    scale_y_continuous(n.breaks = 8) +
    labs(x = "Years", y = "Free Virus")

## Overlay data from 
