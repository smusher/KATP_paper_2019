library(tidyverse)
library(broom)
library(cowplot)
library(ggbeeswarm)
library(knitr)
library(kableExtra)
extrafont::loadfonts()

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

concresp_summary <-
    concresp %>%
    group_by(construct, method, measure, concentration, nucleotide) %>%
    summarise(
        se = sd(response)/sqrt(length(response)),
        response = mean(response)
    )

x <- seq(-7, -2, length.out = 51)
frame <- tibble(log_concentration = x)

hill_equation <- function(log_concentration, ec50, hill, floor){
    floor + ((1 - floor) / (1 + 10^((ec50 - log_concentration) * -hill)))
}

free_floor <-
    readRDS("data/population_hill_fits_free_floor.rds") %>%
    unnest(fit %>% map(augment, newdata = frame)) %>%
    mutate(concentration = 10^log_concentration)

free_floor_tidy <-
    readRDS("data/population_hill_fits_free_floor.rds") %>%
    unnest(fit %>% map(tidy))

fancy_scientific <- function(l) {
    l <- format(l, scientific = TRUE)
    l <- gsub("^(.*)e", "10^", l)
    parse(text=l)
}

split_construct <- function(l) {
    gsub("[+]", "\n +", l)
}

ggplot() +
geom_quasirandom(data = concresp %>% filter(method == "pcf"), aes(x = concentration, y = response, fill = measure), shape = 21, alpha = 0.5, size = 1.5, width = 0.2) +
geom_line(data = free_floor %>% filter(method == "pcf"), aes(x = concentration, y = .fitted, colour = measure), size = 1) +
geom_errorbar(data = concresp_summary %>% filter(method == "pcf"), aes(x = concentration, ymin = response - se, ymax = response + se), width = 0.2) +
geom_point(data = concresp_summary %>% filter(method == "pcf"), aes(x = concentration, y = response, fill = measure), shape = 21, size = 2.5) +
facet_wrap(vars(construct), ncol = 2, labeller = labeller(construct = split_construct)) +
theme_cowplot(font_family = "IBM Plex Sans Condensed Light", font_size = 9) +
theme(strip.background  = element_rect(fill = "white")) +
scale_colour_brewer(palette = "Set2", aesthetics = c("colour", "fill")) +
scale_x_log10(labels = fancy_scientific) +
scale_y_continuous(breaks = c(0.0, 0.5, 1.0)) +
labs(
    x = "[TNP-ATP] (M)",
    y = expression(F/F[max] ~ or ~ I/I[max])
    ) -> pcf_plot

ggplot() +
geom_quasirandom(data = free_floor_tidy %>% filter(method == "pcf"), aes(x = construct, y = estimate, fill = measure), shape = 21, size = 2, stroke = 0.5, dodge.width = 0.8) +
facet_grid(cols = vars(term), scales = "free") +
scale_fill_brewer(palette = "Set2", guide = guide_legend(title.position = "top")) +
coord_flip() +
scale_x_discrete(labels = split_construct) +
theme_cowplot(font_family = "IBM Plex Sans Condensed Light", font_size = 9) +
theme(axis.title.y = element_blank(), legend.position = "top", strip.background  = element_rect(fill = "white")) -> pcf_swarm

ggplot() +
geom_line(data = free_floor %>% filter(method == "unroofed", construct == "W311*-GFP+SUR") %>% select(-construct), aes(x = concentration, y = .fitted), size = 1, colour = "grey50") +
geom_quasirandom(data = concresp %>% filter(method == "unroofed"), aes(x = concentration, y = response, fill = construct), shape = 21, alpha = 0.5, size = 1.5, width = 0.2) +
geom_line(data = free_floor %>% filter(method == "unroofed"), aes(x = concentration, y = .fitted, colour = construct), size = 1) +
geom_errorbar(data = concresp_summary %>% filter(method == "unroofed"), aes(x = concentration, ymin = response - se, ymax = response + se), width = 0.2) +
geom_point(data = concresp_summary %>% filter(method == "unroofed"), aes(x = concentration, y = response, fill = construct), shape = 21, size = 2.5) +
facet_wrap(vars(construct), ncol = 3, labeller = labeller(construct = split_construct)) +
theme_cowplot(font_family = "IBM Plex Sans Condensed Light", font_size = 9) +
theme(strip.background  = element_rect(fill = "white"), legend.position = "none") +
scale_colour_brewer(palette = "Set2", aesthetics = c("colour", "fill")) +
scale_x_log10(labels = fancy_scientific) +
scale_y_continuous(breaks = c(0.0, 0.5, 1.0)) +
labs(
    x = "[TNP-ATP] (M)",
    y = expression(F/F[max])
    ) -> unroofed_plot

ggplot() +
geom_quasirandom(data = free_floor_tidy %>% filter(method == "unroofed", construct != "W311*,G334D-GFP+SUR", construct != "W311*-GFP+SUR-K205E"), aes(x = construct, y = estimate, fill = construct), shape = 21, size = 2, stroke = 0.5, dodge.width = 0.8) +
facet_grid(cols = vars(term), scales = "free") +
scale_fill_brewer(palette = "Set2", guide = guide_legend(title.position = "top")) +
coord_flip() +
scale_x_discrete(labels = split_construct) +
theme_cowplot(font_family = "IBM Plex Sans Condensed Light", font_size = 9) +
theme(axis.title.y = element_blank(), legend.position = "none", strip.background  = element_rect(fill = "white")) -> unroofed_swarm

ggplot() +
geom_line(data = free_floor %>% filter(measure == "current", construct == "W311*-GFP+SUR") %>% select(-construct), aes(x = concentration, y = .fitted), size = 1, colour = "grey50") +
geom_quasirandom(data = concresp %>% filter(measure == "current"), aes(x = concentration, y = response, fill = construct), shape = 21, alpha = 0.5, size = 1.5, width = 0.2) +
geom_line(data = free_floor %>% filter(measure == "current"), aes(x = concentration, y = .fitted, colour = construct), size = 1) +
geom_errorbar(data = concresp_summary %>% filter(measure == "current"), aes(x = concentration, ymin = response - se, ymax = response + se), width = 0.2) +
geom_point(data = concresp_summary %>% filter(measure == "current"), aes(x = concentration, y = response, fill = construct), shape = 21, size = 2.5) +
facet_grid(rows = vars(construct), cols = vars(nucleotide), labeller = labeller(construct = split_construct)) +
theme_cowplot(font_family = "IBM Plex Sans Condensed Light", font_size = 9) +
theme(strip.background  = element_rect(fill = "white"), legend.position = "none", strip.text.y = element_text(size = 4)) +
scale_colour_brewer(palette = "Set2", aesthetics = c("colour", "fill")) +
scale_x_log10(labels = fancy_scientific) +
scale_y_continuous(breaks = c(0.0, 0.5, 1.0)) +
labs(
    x = "[TNP-ATP] or [ATP] (M)",
    y = expression(I/I[max])
    ) -> patch_plot

ggplot() +
geom_quasirandom(data = free_floor_tidy %>% filter(measure == "current"), aes(x = construct, y = estimate, fill = nucleotide), shape = 21, size = 2, stroke = 0.5, dodge.width = 0.8) +
facet_grid(cols = vars(term), scales = "free") +
scale_fill_brewer(palette = "Set2", guide = guide_legend(title.position = "top")) +
coord_flip() +
scale_x_discrete(labels = split_construct) +
theme_cowplot(font_family = "IBM Plex Sans Condensed Light", font_size = 9) +
theme(axis.title.y = element_blank(), legend.position = "top", strip.background  = element_rect(fill = "white")) -> patch_swarm

pcf_plot_no_legend <- pcf_plot + theme(legend.position = "none")

plots <-
    plot_grid(
        pcf_plot_no_legend, pcf_swarm,
        unroofed_plot, unroofed_swarm,
        patch_plot, patch_swarm,
        ncol = 2, align = "h", axis = "b",
        rel_heights = c(2, 2.5, 3)
    )

ggsave("mike_update/concresp.pdf", plots, width = 200, height = 290, units = "mm")

free_floor_tidy %>%
select(-statistic, -p.value) %>%
rename(Method = method, Measure = measure, Nucleotide = nucleotide, Construct = construct, Term = term, Estimate = estimate, `Standard Error` = std.error) %>%
mutate(
    Method = case_when(Method == "pcf" ~ "PCF", Method == "unroofed" ~ "Unroofed", Method == "patch_clamp" ~ "Patch Clamp"),
    Measure = case_when(Measure == "current" ~ "Current", Measure == "fluorescence" ~ "Fluorescence"),
    Term = case_when(Term == "ec50" ~ "EC50 (M-1)", Term == "hill" ~ "h", Term == "floor" ~ "Emax")
    ) %>%
kable("latex", booktabs = T, digits = 2) %>%
kable_styling(latex_options = "striped", font_size = 12) %>%
collapse_rows(1:4, row_group_label_position = "stack") %>%
footnote(general = "EC50 values are expressed in log10 values.") %>%
landscape() -> concresp_table

save_kable(concresp_table, "mike_update/concresp_table.pdf", latex_header_includes = "\\setmainfont[UprightFont = *-Light]{IBM Plex Serif}")

ggplot() +
geom_quasirandom(data = concresp %>% filter(measure == "current", construct %in% c("W311*-GFP+SUR", "WT-GFP+SUR")), aes(x = concentration, y = response, fill = nucleotide), shape = 21, alpha = 0.5, size = 1.5, width = 0.2) +
geom_line(data = free_floor %>% filter(measure == "current", construct %in% c("W311*-GFP+SUR", "WT-GFP+SUR")), aes(x = concentration, y = .fitted, colour = nucleotide), size = 1) +
geom_errorbar(data = concresp_summary %>% filter(measure == "current", construct %in% c("W311*-GFP+SUR", "WT-GFP+SUR")), aes(x = concentration, ymin = response - se, ymax = response + se), width = 0.2) +
geom_point(data = concresp_summary %>% filter(measure == "current", construct %in% c("W311*-GFP+SUR", "WT-GFP+SUR")), aes(x = concentration, y = response, fill = nucleotide), shape = 21, size = 2.5) +
facet_grid(rows = vars(construct), labeller = labeller(construct = split_construct)) +
theme_cowplot(font_family = "IBM Plex Sans Condensed Light", font_size = 9) +
theme(strip.background  = element_rect(fill = "white"), legend.position = "none", strip.text.y = element_text(size = 4)) +
scale_colour_brewer(palette = "Set2", aesthetics = c("colour", "fill")) +
scale_x_log10(labels = fancy_scientific) +
scale_y_continuous(breaks = c(0.0, 0.5, 1.0)) +
labs(
    x = "[TNP-ATP] or [ATP] (M)",
    y = expression(I/I[max])
    )
