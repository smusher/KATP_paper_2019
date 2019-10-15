library(tidyverse)
library(broom)
library(nls.multstart)

currents <-
    read_csv("/home/sam/previous_analysis/2019_manuscript/data/electrophys_data.csv") %>%
    select(unique_experiment_id, method, measure, construct, nucleotide, concentration, response)

fluorescence <-
    read_csv("/home/sam/previous_analysis/2019_manuscript/data/unroofed_concresp_data.csv") %>%
    filter(dye < 480) %>%
    select(unique_experiment_id, method, measure, construct, concentration, response) %>%
    mutate(nucleotide = "TNP-ATP")

pcf_data <-
    read_csv("/home/sam/previous_analysis/2019_manuscript/data/pcf_data.csv") %>%
    filter(measure == "fluorescence", dye < 480) %>%
    select(unique_experiment_id, method, measure, construct, concentration, response) %>%
    mutate(nucleotide = "TNP-ATP")

concresp <-
    bind_rows(currents, fluorescence, pcf_data) %>%
    mutate(log_concentration = log10(concentration)) %>%
    filter(log_concentration > -Inf)

hill_equation <- function(log_concentration, ec50, hill, floor){
    floor + ((1 - floor) / (1 + 10^((ec50 - log_concentration) * -hill)))
}

population_hill_fits_free_floor <-
    concresp %>%
    group_by(method, measure, construct, nucleotide) %>%
    nest() %>%
    mutate(fit = map(data, ~ nls_multstart(
            response ~ hill_equation(log_concentration, ec50, hill, floor),
            data = data.frame(.),
            iter = 500,
            lower = c(-Inf, 0, 0),
            start_lower = c(ec50 = -8, hill = 0, floor= 0),
            start_upper = c(ec50 = 0, hill = 2, floor = 1),
            na.action = na.omit
              )
    )
    )

individual_hill_fits_free_floor <-
    concresp %>%
    group_by(unique_experiment_id, method, measure, construct, nucleotide) %>%
    nest() %>%
    mutate(fit = map(data, ~ nls_multstart(
            response ~ hill_equation(log_concentration, ec50, hill, floor),
            data = data.frame(.),
            iter = 500,
            lower = c(-Inf, 0, 0),
            start_lower = c(ec50 = -8, hill = 0, floor= 0),
            start_upper = c(ec50 = 0, hill = 2, floor = 1),
            na.action = na.omit
              )
    )
    )

population_hill_fits_free_floor_glance <-
    population_hill_fits_free_floor %>%
    unnest(fit %>% map(glance))

population_hill_fits_free_floor_tidy <-
    population_hill_fits_free_floor_glance %>%
    unnest(fit %>% map(tidy))

x <- seq(-8, -2, length.out = 51)
frame <- tibble(log_concentration = x)

population_hill_fits_free_floor_augment <-
    population_hill_fits_free_floor_glance %>%
    unnest(fit %>% map(augment, newdata = frame))

fixed_floor <-
    population_hill_fits_free_floor_tidy %>%
    filter(construct == "W311*-GFP+SUR", method == "unroofed", term == "floor") %>%
    pull(estimate)

hill_equation_fixed_floor <- function(log_concentration, ec50, hill){
    fixed_floor + ((1 - fixed_floor) / (1 + 10^((ec50 - log_concentration) * -hill)))
}

population_hill_fits_fixed_floor <-
    concresp %>%
    group_by(method, measure, construct, nucleotide) %>%
    nest() %>%
    mutate(fit = map(data, ~ nls_multstart(
            response ~ hill_equation_fixed_floor(log_concentration, ec50, hill),
            data = data.frame(.),
            iter = 500,
            start_lower = c(ec50 = -8, hill = 0),
            start_upper = c(ec50 = 0, hill = 2),
            na.action = na.omit
              )
    )
    )

population_hill_fits_fixed_floor_glance <-
    population_hill_fits_fixed_floor %>%
    unnest(fit %>% map(glance))

population_hill_fits_fixed_floor_tidy <-
    population_hill_fits_fixed_floor_glance %>%
    unnest(fit %>% map(tidy))

x <- seq(-8, -2, length.out = 51)
frame <- tibble(log_concentration = x)

population_hill_fits_fixed_floor_augment <-
    population_hill_fits_fixed_floor_glance %>%
    unnest(fit %>% map(augment, newdata = frame))

saveRDS(population_hill_fits_free_floor, "/home/sam/previous_analysis/2019_manuscript/data/population_hill_fits_free_floor.rds")
saveRDS(individual_hill_fits_free_floor, "/home/sam/previous_analysis/2019_manuscript/data/individual_hill_fits_free_floor.rds")
saveRDS(population_hill_fits_fixed_floor, "/home/sam/previous_analysis/2019_manuscript/data/population_hill_fits_fixed_floor.rds")
