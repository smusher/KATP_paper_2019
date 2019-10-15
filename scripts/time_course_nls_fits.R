library(tidyverse)
library(broom)
library(nls.multstart)

time_course <-
    read_csv("/home/sam/previous_analysis/2019_manuscript/data/unroofed_timecourse_data.csv") %>%
    filter(dye < 480, test_concentration >= 1e-5) %>%
    select(unique_experiment_id, construct, repetition, time, concentration, test_concentration, response, washout, washon, steady_state)

washout_equation <- function(time, a, tau){
    a * exp(-time * tau)
}

washout <-
    time_course %>%
    filter(washout == TRUE) %>%
    group_by(construct, test_concentration, repetition) %>%
    mutate(time = time - min(time))

population_washout_fits <-
    washout %>%
    group_by(construct, test_concentration) %>%
    nest() %>%
    mutate(fit = map(data, ~ nls_multstart(
            response ~ washout_equation(time, a, tau),
            data = data.frame(.),
            iter = 500,
            lower = c(a = 0, tau = 0),
            start_lower = c(a = 0, tau = 0),
            start_upper = c(a = 2, tau = 2),
            na.action = na.omit
              )
    )
    )

population_washout_fits_glance <-
    population_washout_fits %>%
    unnest(fit %>% map(glance))

population_washout_fits_tidy <-
    population_washout_fits_glance %>%
    unnest(fit %>% map(tidy))

washon_equation <- function(time, a, tau){
    1 - ((a * exp(-time * tau)) + (1 - a))
}

washon <-
    time_course %>%
    filter(washon == TRUE) %>%
    group_by(construct, test_concentration, repetition) %>%
    mutate(time = time - min(time))

population_washon_fits <-
    washon %>%
    group_by(construct, test_concentration) %>%
    nest() %>%
    mutate(fit = map(data, ~ nls_multstart(
            response ~ washon_equation(time, a, tau),
            data = data.frame(.),
            iter = 500,
            lower = c(a = 0, tau = 0),
            start_lower = c(a = 0, tau = 0),
            start_upper = c(a = 2, tau = 2),
            na.action = na.omit
              )
    )
    )

population_washon_fits_glance <-
    population_washon_fits %>%
    unnest(fit %>% map(glance))

population_washon_fits_tidy <-
    population_washon_fits_glance %>%
    unnest(fit %>% map(tidy))

taus <-
    bind_rows(
        population_washon_fits_tidy %>%
        filter(term == "tau") %>%
        select(construct, test_concentration, estimate, std.error) %>%
        mutate(term = "Kobs"),
        population_washout_fits_tidy %>%
        filter(term == "tau") %>%
        select(construct, test_concentration, estimate, std.error) %>%
        mutate(term = "Koff")
        )

population_taus_fits <-
    taus %>%
    group_by(construct, term) %>%
    nest() %>%
    mutate(fit = map(data, ~lm(estimate ~ test_concentration, data = data.frame(.))))

saveRDS(population_washon_fits, "/home/sam/previous_analysis/2019_manuscript/data/population_washon_fits.rds")
saveRDS(population_washout_fits, "/home/sam/previous_analysis/2019_manuscript/data/population_washout_fits.rds")
saveRDS(population_taus_fits, "/home/sam/previous_analysis/2019_manuscript/data/population_taus_fits.rds")
