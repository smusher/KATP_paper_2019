library(tidyverse)
library(broom)
library(brms)
library(tidybayes)
library(cowplot)
library(ggbeeswarm)
library(ggforce)
library(knitr)
library(kableExtra)
extrafont::loadfonts()

filenames <-
    list.files(path = "data/mwc_fits", full.names = T)

models <-
    filenames %>%
    plyr::llply(readRDS)

names <-
    filenames %>%
    str_extract_all("(?<=mwc_fits\\/)[^.]*") %>%
    unlist()

model_fits <-
    tibble(name = names, fit = models) %>%
    separate(name, c("construct", "model", "hierarchy"), sep = "_")

posterior_spread <-
    model_fits %>%
    group_by(construct, model, hierarchy) %>%
    unnest(fit %>% map(spread_draws, `b_.*`, regex = T)) %>%
    rename(log_L = b_logL_Intercept, D = b_D_Intercept, log_Ka = b_logKa_Intercept, c = b_c_Intercept) %>%
    mutate(L = 10^log_L, Ka = 10^log_Ka, Po = L / (L+1)) %>%
    mutate(probability = "posterior")

prior_spread <-
    model_fits %>%
    group_by(construct, model, hierarchy) %>%
    unnest(fit %>% map(spread_draws, `prior_b_.*`, regex = T)) %>%
    rename(log_L = prior_b_logL, D = prior_b_D, log_Ka = prior_b_logKa, c = prior_b_c) %>%
    mutate(L = 10^log_L, Ka = 10^log_Ka, Po = L / (L+1)) %>%
    mutate(probability = "prior")

spreaded_draws <-
    bind_rows(posterior_spread, prior_spread) %>%
    select(construct, model, hierarchy, probability, c, D, log_Ka, log_L)

gathered_draws <-
    spread_draws %>%
    gather(parameter, estimate, c, D, log_Ka, log_L)

ggplot(spreaded_draws %>% filter(hierarchy == "population", probability == "posterior", construct == "W311*-GFP+SUR")) +
# geom_autopoint(shape = ".", aes(colour = construct)) +
geom_autodensity(colour = "black", position = "identity", alpha = 0.5, aes(fill = model)) +
stat_density_2d(geom = "contour", aes(x = .panel_x, y = .panel_y, colour = model)) +
facet_matrix(vars(c:log_L), layer.lower = 2, layer.diag = 1, layer.upper = 3) +
theme_cowplot(font_family = "IBM Plex Sans Condensed", font_size = 9) +
theme(strip.background  = element_rect(fill = "white"), legend.position = "top")

toplot <-
    gathered_draws %>%
    filter(probability == "posterior", hierarchy == "population")

toplot_2 <-
    gathered_draws %>%
    filter(probability == "prior", hierarchy == "population")

ggplot() +
ggtitle("mwc") +
geom_density(data = toplot %>% filter(construct == "W311*-GFP+SUR"), aes(estimate, fill = model), colour = "black", alpha = 0.5) +
geom_density(data = toplot_2 %>% filter(construct == "W311*-GFP+SUR"), aes(estimate), colour = "black", fill = "grey50", alpha = 0.5) +
scale_fill_brewer(palette = "Set2") +
facet_wrap(vars(parameter), scales = "free", ncol = 4) +
theme_cowplot(font_family = "IBM Plex Sans Condensed Light", font_size = 9) +
theme(strip.background  = element_rect(fill = "white"), legend.position = "top") -> dens_1

ggplot(toplot %>% filter(model == "single"), aes(estimate, fill = construct)) +
ggtitle("single") +
geom_density(colour = "black", alpha = 0.5) +
scale_fill_brewer(palette = "Set2") +
facet_wrap(vars(parameter), scales = "free", ncol = 4) +
theme_cowplot(font_family = "IBM Plex Sans Condensed Light", font_size = 9) +
theme(strip.background  = element_rect(fill = "white"), legend.position = "top") -> dens_2

ggplot() +
ggtitle("coop") +
geom_density(data = toplot %>% filter(model == "coop"), aes(estimate, fill = construct), colour = "black", alpha = 0.5) +
geom_density(data = toplot_2 %>% filter(model == "coop"), aes(estimate), colour = "black", linetype = 2, alpha = 0.5) +
scale_fill_brewer(palette = "Set2") +
facet_wrap(vars(parameter), scales = "free", ncol = 4) +
theme_cowplot(font_family = "IBM Plex Sans Condensed Light", font_size = 9) +
theme(strip.background  = element_rect(fill = "white"), legend.position = "top") -> dens_3

pcf_data <-
    read_csv("/home/sam/previous_analysis/2019_manuscript/data/pcf_data.csv") %>%
    filter(dye < 480 | is.na(dye), concentration > 0) %>%
    select(unique_experiment_id, construct, measure, method, concentration, response)

pcf_data_summary <-
    pcf_data %>%
    group_by(construct, measure, concentration) %>%
    summarise(
        se = sd(response)/sqrt(length(response)),
        response = mean(response)
    )

x <- c(10^seq(-7, -2, length.out = 51))
z <- c(0, 1)

frame <-
    expand.grid(
        concentration = x,
        binding_mask = z,
        unique_experiment_id = 1
    ) %>%
    as_tibble()

fits <-
    model_fits %>%
    unnest(fit %>% map(fitted_draws, newdata = frame, re_formula = NA, allow_new_levels = TRUE)) %>%
    group_by(construct, model, hierarchy, binding_mask, concentration) %>%
    point_interval(.value, .point = median, .interval = qi) %>%
    mutate(
        measure = case_when(binding_mask == 1 ~ "fluorescence", TRUE ~ "current"),
        response = case_when(measure == "fluorescence" ~ 1 - .value, TRUE ~ .value),
        .lower = case_when(measure == "fluorescence" ~ 1 - .lower, TRUE ~ .lower),
        .upper = case_when(measure == "fluorescence" ~ 1 - .upper, TRUE ~ .upper)
    )

fancy_scientific <- function(l) {
    l <- format(l, scientific = TRUE)
    l <- gsub("^(.*)e", "10^", l)
    parse(text=l)
}

split_construct <- function(l) {
    gsub("[+]", "\n +", l)
}

fits <-
    fits %>%
    ungroup %>%
    mutate(model = case_when(model == "coop" ~ "Negative Cooperativity", model == "mwc" ~ "Full MWC", model == "single" ~ "Single binding"))

ggplot() +
geom_ribbon(data = fits %>% filter(hierarchy == "population"), aes(x = concentration, ymin = .lower, ymax = .upper, fill = measure), alpha = 0.25) +
geom_quasirandom(data = pcf_data, aes(x = concentration, y = response, fill = measure), shape = 21, alpha = 0.25, size = 1.5, width = 0.2) +
geom_line(data = fits %>% filter(hierarchy == "population"), aes(x = concentration, y = response, colour = measure), size = 1) +
geom_line(data = fits %>% filter(hierarchy == "fixed"), aes(x = concentration, y = response, colour = measure), size = 1, linetype = 2) +
geom_errorbar(data = pcf_data_summary, aes(x = concentration, ymin = response - se, ymax = response + se), width = 0.2) +
geom_point(data = pcf_data_summary, aes(x = concentration, y = response, fill = measure), shape = 21, size = 3) +
facet_grid(rows = vars(construct), cols = vars(model), labeller = labeller(construct = split_construct)) +
theme_cowplot(font_family = "IBM Plex Sans Condensed", font_size = 12) +
theme(strip.background  = element_rect(fill = "white"), legend.position = "top") +
scale_colour_manual(values = c("#1E88E5", "#FB8C00"), aesthetics = c("colour", "fill")) +
scale_x_log10(labels = fancy_scientific) +
scale_y_continuous(breaks = c(0.0, 0.5, 1.0)) +
labs(
    x = "[TNP-ATP] (M)",
    y = expression(F/F[max] ~ or ~ I/I[max])
    ) -> fit_plot

legend <- dens_1 %>% get_legend()
dens_1 <- dens_1 + theme(legend.position = "none")
dens_2 <- dens_2 + theme(legend.position = "none")
dens_3 <- dens_3 + theme(legend.position = "none")
top_row <- plot_grid(legend, dens_1, dens_2, dens_3, nrow = 4, align = "v", axis = "l", rel_heights = c(0.5, 1, 1, 1))
plots <- plot_grid(top_row, fit_plot, ncol = 1, align = "hv", axis = "lb")

ggsave("mike_update/population_bayes.pdf", plots, width = 200, height = 290, units = "mm")

toplot <-
    gather_draws %>%
    filter(probability == "posterior", hierarchy == "individual", parameter != "L", parameter != "Ka", parameter != "Po")

ggplot(toplot %>% filter(model == "mwc"), aes(estimate, fill = construct)) +
ggtitle("mwc") +
geom_density(colour = "black", alpha = 0.5) +
scale_fill_brewer(palette = "Set2") +
facet_wrap(vars(parameter), scales = "free", ncol = 4) +
theme_cowplot(font_family = "IBM Plex Sans Condensed Light", font_size = 9) +
theme(strip.background  = element_rect(fill = "white"), legend.position = "top") -> dens_1

ggplot(toplot %>% filter(model == "single"), aes(estimate, fill = construct)) +
ggtitle("single") +
geom_density(colour = "black", alpha = 0.5) +
scale_fill_brewer(palette = "Set2") +
facet_wrap(vars(parameter), scales = "free", ncol = 4) +
theme_cowplot(font_family = "IBM Plex Sans Condensed Light", font_size = 9) +
theme(strip.background  = element_rect(fill = "white"), legend.position = "top") -> dens_2

ggplot(toplot %>% filter(model == "coop"), aes(estimate, fill = construct)) +
ggtitle("coop") +
geom_density(colour = "black", alpha = 0.5) +
scale_fill_brewer(palette = "Set2") +
facet_wrap(vars(parameter), scales = "free", ncol = 4) +
theme_cowplot(font_family = "IBM Plex Sans Condensed Light", font_size = 9) +
theme(strip.background  = element_rect(fill = "white"), legend.position = "top") -> dens_3

ggplot() +
geom_ribbon(data = fits %>% filter(hierarchy == "individual"), aes(x = concentration, ymin = .lower, ymax = .upper, fill = measure), alpha = 0.25) +
geom_quasirandom(data = pcf_data, aes(x = concentration, y = response, fill = measure), shape = 21, alpha = 0.25, size = 1.5, width = 0.2) +
geom_line(data = fits %>% filter(hierarchy == "individual"), aes(x = concentration, y = response, colour = measure), size = 1) +
geom_errorbar(data = pcf_data_summary, aes(x = concentration, ymin = response - se, ymax = response + se), width = 0.2) +
geom_point(data = pcf_data_summary, aes(x = concentration, y = response, fill = measure), shape = 21, size = 2.5) +
facet_grid(rows = vars(construct), cols = vars(model), labeller = labeller(construct = split_construct)) +
theme_cowplot(font_family = "IBM Plex Sans Condensed Light", font_size = 9) +
theme(strip.background  = element_rect(fill = "white"), legend.position = "top") +
scale_colour_brewer(palette = "Set2", aesthetics = c("colour", "fill")) +
scale_x_log10(labels = fancy_scientific) +
scale_y_continuous(breaks = c(0.0, 0.5, 1.0)) +
labs(
    x = "[TNP-ATP] (M)",
    y = expression(F/F[max] ~ or ~ I/I[max])
    ) -> fit_plot

legend <- dens_1 %>% get_legend()
dens_1 <- dens_1 + theme(legend.position = "none")
dens_2 <- dens_2 + theme(legend.position = "none")
dens_3 <- dens_3 + theme(legend.position = "none")
top_row <- plot_grid(legend, dens_1, dens_2, dens_3, nrow = 4, align = "v", axis = "l", rel_heights = c(0.5, 1, 1, 1))
plots <- plot_grid(top_row, fit_plot, ncol = 1, align = "hv", axis = "lb")

ggsave("mike_update/individual_bayes.pdf", plots, width = 200, height = 290, units = "mm")

gather_draws %>%
group_by(construct, model, hierarchy, parameter) %>%
point_interval(estimate, .point = median, .interval = qi) %>%
select(construct, model, hierarchy, parameter, estimate, .lower, .upper) %>%
filter(parameter %in% c("L", "Po", "D", "Ka", "c")) %>%
kable("latex", booktabs = T, digits = 2) %>%
landscape() %>%
kable_styling(latex_options = "striped", font_size = 12) %>%
add_header_above(c(" " = 4, "median estimates" = 3)) -> fit_table

save_kable(fit_table, "mike_update/bayes_table.pdf", latex_header_includes = "\\setmainfont[UprightFont = *-Light]{IBM Plex Serif}")

 wt <-
    model_fits %>%
    filter(construct == "W311*-GFP+SUR", hierarchy != "individual") %>%
    mutate(loo_score = fit %>% map(loo, reloo = TRUE))
