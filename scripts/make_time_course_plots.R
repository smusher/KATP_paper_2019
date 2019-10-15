library(tidyverse)
library(broom)
library(cowplot)
library(ggbeeswarm)
library(ggforce)
library(knitr)
library(kableExtra)
extrafont::loadfonts()

population_washon_fits <- readRDS("data/population_washon_fits.rds") %>%
    filter(test_concentration >= 1e-5)
population_washout_fits <- readRDS("data/population_washout_fits.rds") %>%
    filter(test_concentration >= 1e-5)
population_taus_fits <- readRDS("data/population_taus_fits.rds")

time_course <-
    read_csv("data/unroofed_timecourse_data.csv") %>%
    filter(dye < 480) %>%
    group_by(unique_experiment_id, construct, repetition) %>%
    mutate(time = time - min(time)) %>%
    filter(time <= 120, test_concentration >= 1e-5)

time_course_summary <-
    time_course %>%
    group_by(construct, unique_experiment_id, time, concentration, test_concentration) %>%
    summarise(response = mean(response)) %>%
    group_by(construct, time, concentration, test_concentration) %>%
    summarise(
        se = sd(response) / sqrt(length(response)),
        response = mean(response)
        )

washon_equation <- function(time, a, tau){
    1 - ((a * exp(-time * tau)) + (1 - a))
}

washout_equation <- function(time, a, tau){
    a * exp(-time * tau)
}

x <- seq(0, 30, length.out = 51)
frame <- tibble(time = x)
population_washout_fits_augment <-
    population_washout_fits %>%
    unnest(fit %>% map(augment, newdata = frame)) %>%
    mutate(time = time + 90)

population_washon_fits_augment <-
    population_washon_fits %>%
    unnest(fit %>% map(augment, newdata = frame)) %>%
    mutate(time = time + 30)

ggplot() +
geom_beeswarm(data = time_course, aes(x = time, y = response, fill = factor(concentration)), alpha = 0.5, shape = 21) +
geom_errorbar(data = time_course_summary, aes(x = time, ymin = response - se, ymax = response + se)) +
geom_point(data = time_course_summary, aes(x = time, y = response, fill = factor(concentration)), size = 2, shape = 21) +
geom_line(data = population_washon_fits_augment, aes(x = time, y = .fitted, colour = factor(test_concentration))) +
geom_line(data = population_washout_fits_augment, aes(x = time, y = .fitted, colour = factor(test_concentration))) +
facet_grid(rows = vars(test_concentration), cols = vars(construct)) +
theme_cowplot(font_family = "IBM Plex Sans Condensed Light", font_size = 9) +
theme(strip.background  = element_rect(fill = "white"), legend.position = "top") +
scale_colour_brewer(palette = "PuOr", aesthetics = c("colour", "fill"), direction = -1) +
coord_cartesian(ylim = c(0, 1))

population_washon_fits_tidy <-
    population_washon_fits %>%
    unnest(fit %>% map(tidy))

population_washout_fits_tidy <-
    population_washout_fits %>%
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

population_taus_fits_augment <-
    population_taus_fits %>%
    unnest(fit %>% map(augment))

ggplot() +
geom_line(data = population_taus_fits_augment, aes(x = test_concentration, y = .fitted, colour = term)) +
geom_errorbar(data = taus, aes(x = test_concentration, ymin = estimate - std.error, ymax = estimate + std.error)) +
geom_point(data = taus, aes(x = test_concentration, y = estimate, fill = term), shape = 21, size = 2) +
theme_cowplot(font_family = "IBM Plex Sans Condensed Light", font_size = 9) +
theme(strip.background  = element_rect(fill = "white"), legend.position = "top") +
facet_grid(cols = vars(construct))

population_taus_fits_tidy <-
    population_taus_fits %>%
    unnest(fit %>% map(tidy))

kon <-
    population_taus_fits_tidy %>%
    filter(term == "Kobs", term1 == "test_concentration") %>%
    rename(Kon = estimate, Kon_se = std.error) %>%
    select(construct, Kon, Kon_se)

koff <-
    population_taus_fits_tidy %>%
    filter(term == "Koff", term1 == "(Intercept)") %>%
    rename(Koff = estimate, Koff_se = std.error) %>%
    select(construct, Koff, Koff_se)

time_course_kd <-
    left_join(kon, koff) %>%
    mutate(
        Kd = Koff / Kon,
        Kd_min = (Koff - Koff_se) / (Kon + Kon_se),
        Kd_max = (Koff + Koff_se) / (Kon - Kon_se)
    )
