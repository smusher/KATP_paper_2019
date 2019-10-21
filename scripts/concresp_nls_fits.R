library(tidyverse)
library(broom)
library(nls.multstart)

currents <-
    read_csv("data/electrophys_data.csv") %>%
    select(unique_experiment_id, method, measure, construct, nucleotide, concentration, response)

fluorescence <-
    read_csv("data/unroofed_concresp_data.csv") %>%
    filter(dye < 480) %>%
    select(unique_experiment_id, method, measure, construct, concentration, response) %>%
    mutate(nucleotide = "TNP-ATP")

pcf_data <-
    read_csv("data/pcf_data.csv") %>%
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

saveRDS(population_hill_fits_free_floor, "data/population_hill_fits_free_floor.rds")
