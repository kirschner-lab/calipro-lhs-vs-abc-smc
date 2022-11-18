## Menzies NA, Soeteman DI, Pandya A, Kim JJ. Bayesian Methods for Calibrating
## Health Policy Models: A Tutorial. Pharmacoeconomics. 2017
## Jun;35(6):613-624. doi: 10.1007/s40273-017-0494-4. PMID: 28247184;
## PMCID: PMC5448142.

library(deSolve)
library(stringr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(tibble)

pop_hiv_parameters <- c( # "Implied mean" values from Table 1 in the paper:
    a = 15e3,     # Annual birth rate for steady state based on 1M population
                  # size and mu_B.
    b = 0.20,     # Fraction of births entering non-susceptible state.
    mu_B = 0.015, # Background mortality rate.
    mu_E = 0.05,  # Disease-specific mortality for early disease.
    mu_L = 0.25,  # Disease-specific mortality for late disease.
    mu_T = 0.025, # Disease-specific mortality on treatment.
    rho = 0.025,  # Effective contact rate for transmission.
    p = 0.10,     # Rate of progression from early to late disease.
    #r_E = 0.0,    # Rate of treatment uptake for early disease.
    r_L = 0.50,   # Rate of treatment uptake for late disease.
    c_T = 1e3     # Annual cost of treatment.
)
b <- pop_hiv_parameters["b"]
pop_hiv_initial_values <- c( # "init_pop" defined in supplement "4 Code ..."
    N = 1 - b,
    S = b - 0.001,
    E = 0.001,
    L = 0,
    T = 0,
    D = 0
) %>%
    ## For some reason, the named vector generates N.b and S.b instead of N and
    ## S, maybe because of the dependence on b!?  Therefore we need to remove
    ## the .b unwanted suffix.
    setNames(names(.) %>% str_remove(fixed(".b")))
stopifnot(sum(pop_hiv_initial_values) == 1)
pop_hiv_initial_values <- pop_hiv_initial_values * 1e6 # Initial population.

## Population-level HIV model in a high-burden setting.
pop_hiv_model <- function(t, pop_hiv_initial_values, pop_hiv_parameters) {
    ## Conveniently access vector elements by name by constructing an
    ## environment from the elements using with() and as.list().  This was shown
    ## in the deSolve package vignette describing the Lorenz equation.
    with(as.list(c(pop_hiv_initial_values, pop_hiv_parameters)), {
        ## The force of infection "lambda" is an subexpression used in some the
        ## differential equations below.  "lambda" is defined in section "3.2
        ## Study model" of the main paper and also at the end of the section "3
        ## Model equations" of the PDF supplement.  Note that we call this a
        ## subexpression rather than use the author's term in the supplement of
        ## "endogenous time varying parameter" because in the context of pure
        ## ODEs, parameters are fixed and cannot depend on variables.
        lambda_t <- rho * (E + L) / (S + E + L + T)
        r_E_t <- ifelse(t > 30, r_L, r_L*t/30)

        ## Rates of change (Supplement PDF section "3 Model equations").
        ##
        ## Non-susceptible (N).
        dN.dt <- a*b - mu_B*N
        ## Susceptible (S).
        dS.dt <- a*(1 - b) - mu_B*S - lambda_t*S
        ## Early disease (E).
        dE.dt <- lambda_t*S - (mu_B + mu_E)*E - r_E_t*E - c_T*E
        ## Late disease (L).
        dL.dt <- c_T*E - (mu_B + mu_L)*L - r_L*L
        ## Treatment (T).
        dT.dt <- r_E_t*E + r_L*L - (mu_B + mu_T)*T
        ## Dead (D).
        dD.dt <- mu_B*(N + S + E + L + T) + mu_E*E + mu_L*L + mu_T*T

        ## Return the rates of change.
        list(c(dN.dt, dS.dt, dE.dt, dL.dt, dT.dt, dD.dt))
    })
}

years <- seq(0, 35, by = 1/12)
out <- ode(y = pop_hiv_initial_values, times = years, func = pop_hiv_model,
           parms = pop_hiv_parameters)

summary(out)

df <-
    as.data.frame(out) %>%
    as_tibble() %>%
    rename(years = time)

## The total population is increasing instead of remaining steady!?
df %>%
    mutate(total = N + S + E + L + T + D) %>%
    ggplot(aes(years, total)) +
    geom_line()

df %>%
    mutate(no_death = N + S + E + L + T) %>%
    ggplot(aes(years, no_death)) +
    geom_line()

df_long <-
    df %>%
    rename("Non-susceptible" = N,
           Susceptible = S,
           "Early disease" = E,
           "Late disease" = L,
           "Treatment" = T,
           "Dead" = D) %>%
    gather("variable", "count", -years)

ggplot(df_long, aes(years, count)) +
    facet_wrap(~ variable, scales = "free_y") +
    geom_line() +
    ggtitle("Population categories")

ggplot(df_long, aes(years, count, fill = variable)) +
    geom_area(position = "stack") +
    ## scale_y_log10() +
    ggtitle("Population categories")

data_prevalence <- tribble(
    ~years, ~prevalence, ~ymin, ~ymax,
    10,  5.0,  3.3,  7.1,
    20, 15.0, 12.0, 18.3,
    30, 10.0,  7.5, 12.8)
(grob <-
     df %>%
     ## Prevalence is defined in Table 2 of the paper.
     mutate(prevalence = (E + L + T) / (N + S + E + L + T) * 100) %>%
     ggplot(aes(years, prevalence)) +
     geom_line(color = "darkgreen") +
     theme_classic() +
     theme(plot.background = element_rect(color = "black", fill = NA,
                                          linetype = 2)) +
     labs(x = "", y = ""))
(plot_prevalence <-
     df %>%
     ## Prevalence is defined in Table 2 of the paper.
     mutate(prevalence = (E + L + T) / (N + S + E + L + T) * 100) %>%
     ggplot(aes(years, prevalence)) +
     geom_line(color = "darkgreen") +
     geom_pointrange(data = data_prevalence,
                     aes(ymin = ymin, ymax = ymax)) +
     theme_bw() +
     labs(x = "Year", y = "") +
     ggtitle("Disease prevalence (%)") +
     annotation_custom(ggplotGrob(grob),
                       xmin = 11, xmax = 29,
                       ymin = 1, ymax = 9) +
     geom_line(data = tibble(years = c(0, 11, 29, 30),
                             prevalence = c(0, 1, 1, 0),
                             group = c(1, 1, 2, 2)),
               aes(group = group),
               linetype = 2))

## "HIV survival without treatment" defined in the supplmentary code under
## "Report results".
survival_years <-
    with(as.list(pop_hiv_parameters),
         1 / (mu_B + mu_E) +
         p / (p + mu_B + mu_E) *
         1 / (mu_B + mu_L))
survival_years

data_survival <- tibble(years_min = 8,
                        years = 10,
                        years_max = 12,
                        name = 1)
(plot_survival <-
     tibble::enframe(survival_years, name = "name", value = "years") %>%
     ggplot(aes(years, name)) +
     geom_label(aes(label = format(survival_years, digits = 3)),
                color = "darkgreen") +
     geom_pointrange(data = data_survival,
                     aes(xmin = years_min, xmax = years_max)) +
     theme_bw() +
     theme(axis.ticks.y = element_blank(),
           axis.text.y = element_blank()) +
     labs(x = "Year", y = "") +
     ggtitle("Average survival (years)"))

## Treatment volume (000s).
(plot_treatment <-
     df %>%
     ggplot(aes(years, cumsum(T) / 1000)) +
     geom_line(color = "darkgreen") +
     geom_pointrange(data = tibble(x = 30, y = 75, ymin = 70, ymax = 80),
                     aes(x = x, y = y, ymin = ymin, ymax = ymax)) +
     labs(x = "Year", y = "") +
     ggtitle("Treatment volume (000s)") +
     theme_bw())
## TODO Also overlay data.

## Population hazard rates (per 100 person years).
(plot_hazard <-
     df %>%
     ## Hazard rate not defined in the paper.  Used
     ## https://en.wikipedia.org/wiki/Hazard_ratio#Definition_and_derivation
     mutate(incidence = c(0, diff(T) / S[-length(S)] / diff(years)) * 100,
            mortality = c(0, diff(D) / S[-length(S)] / diff(years)) * 100) %>%
     select(years, incidence, mortality) %>%
     gather("type", "hazard_rate", -years) %>%
     ggplot(aes(years, hazard_rate, color = type)) +
     geom_line() +
     theme_bw() +
     scale_color_manual(values = c(incidence = "blue", mortality = "red")) +
     guides(color = "none") +
     labs(x = "Year", y = "") +
     ggtitle("Population hazard rates (per 100 PY)"))

(cow <- plot_grid(plot_prevalence, plot_survival, plot_treatment, plot_hazard))

ggsave("~/Sync/lab-reports/2022-11-09/img/40273_2017_494_Fig4_uncal.png", cow,
       width = 7, height = 6)
