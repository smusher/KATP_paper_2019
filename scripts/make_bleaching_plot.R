library(tidyverse)
library(broom)
library(nls.multstart)
library(patchwork)
library(ggbeeswarm)

pcf_data <-
    read_csv("data/pcf_data.csv") %>%
    filter(measure == "fluorescence", dye < 480, construct == "W311*-GFP+SUR") %>%
    select(unique_experiment_id, construct, concentration, cumulative_exposure, raw_intensity) %>%
    group_by(unique_experiment_id) %>%
    mutate(
    	normalised_intensity = raw_intensity / max(raw_intensity),
    	cumulative_exposure = cumulative_exposure / 1000
    	)

bleaching_equation <- function(time, a, k){
    a * exp(time * -k) + (1 - a)
}

pcf_bleaching_fits <-
	pcf_data %>%
	filter(concentration == 0) %>%
	group_by(construct, unique_experiment_id) %>%
	nest() %>%
	mutate(fit = map(data, ~ nls_multstart(
            normalised_intensity ~ bleaching_equation(cumulative_exposure, a, k),
            data = data.frame(.),
            iter = 500,
            lower = c(a = 0, k = 0),
            start_lower = c(a = 0, k = 0.01),
            start_upper = c(a = 2, k = 1),
            na.action = na.omit
              )
    )
    )

frame <-
	tibble(
		cumulative_exposure = seq(min(pcf_data$cumulative_exposure), max(pcf_data$cumulative_exposure), length.out = 51)
	)

pcf_bleaching_plots <-
	pcf_bleaching_fits %>%
	mutate(fit = map(fit, augment, newdata = frame)) %>%
	unnest(fit) %>%
	left_join(pcf_data, by = c("unique_experiment_id", "cumulative_exposure", "construct")) %>%
	group_by(unique_experiment_id) %>%
	mutate(
        bleaching_curve = .fitted * max(raw_intensity, na.rm = TRUE),
        corrected_intensity = raw_intensity / .fitted
        )

fancy_scientific <- function(l) {
    l <- format(l, scientific = TRUE)
    l <- gsub("^(.*)e", "10^", l)
    l <- gsub("0.001", "10^-3", l)
    parse(text=l)
}

ggplot() +
geom_point(data = pcf_bleaching_plots %>% filter(unique_experiment_id > 7), aes(x = cumulative_exposure, y = corrected_intensity, colour = factor(concentration)), fill = NA, shape = 21, size = 3) +
geom_point(data = pcf_data %>% filter(unique_experiment_id > 7), aes(x = cumulative_exposure, y = raw_intensity, fill = factor(concentration)), shape = 21, size = 3) +
geom_line(data = pcf_bleaching_plots %>% filter(unique_experiment_id > 7), aes(x = cumulative_exposure, y = bleaching_curve), linetype = 2, size = 1) +
scale_fill_brewer(
	palette = "PuOr",
	direction = -1,
	labels = fancy_scientific,
    aesthetics = c("fill", "colour")
	) +
facet_wrap(vars(unique_experiment_id)) +
labs(
	x = "Cumulative exposure (seconds)",
	y = "Raw fluorescence intensity (A.U.)",
	fill = "[TNPATP] (M)"
	) +
hrbrthemes::theme_ipsum() +
theme(legend.position = "bottom") -> bleaching_plot

pcf_bleaching_tidy <-
    pcf_bleaching_fits %>%
    mutate(fit = map(fit, tidy)) %>%
    unnest(fit)

ggplot() +
geom_quasirandom(data = pcf_bleaching_tidy %>% filter(unique_experiment_id !=6), aes(x = term, y = estimate, fill = term), shape = 21, size = 3, width = 0.2) +
scale_fill_brewer(palette = "Set2") +
hrbrthemes::theme_ipsum() +
theme(legend.position = "bottom") -> terms_plot

pcf_data_min <-
    pcf_data %>%
    group_by(unique_experiment_id) %>%
    filter(concentration == 0) %>%
    filter(cumulative_exposure == max(cumulative_exposure))

ggplot() +
geom_quasirandom(data = pcf_data_min, aes(x = construct, y = normalised_intensity), shape = 21, size = 3, width = 0.2) +
scale_fill_brewer(palette = "Set2") +
hrbrthemes::theme_ipsum() +
theme(legend.position = "bottom") +
coord_cartesian(ylim = c(0, 1)) +
geom_hline(data = pcf_data_min, aes(yintercept = mean(normalised_intensity))) -> max_bleaching_plot

mean(pcf_data_min$normalised_intensity)