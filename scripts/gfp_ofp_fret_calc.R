library(tidyverse)

import <- function(directory, pattern) {

    oldwd <- getwd()

    setwd(directory)

    filenames <-
        list.files(pattern = pattern)

    tempdf <-
        filenames %>%
        plyr::llply(
            read_csv, col_types = cols(
                construct = col_character(),
                date = col_date(),
                cumulative_exposure = col_integer(),
                excitation_wavelength = col_double(),
                wavelength = col_double(),
                raw_spectra = col_double(),
                background_corrected_spectra = col_double()
            )
        ) %>%
        bind_rows(.id = "n") %>%
        type_convert()

    setwd(oldwd)

    return(tempdf)
}

spectra_1 <-
    import("/home/sam/current_analysis/kir_sur_fret/whole_cell", "spectra") %>%
    filter(wavelength > 450, wavelength < 700) %>%
    select(n, construct, wavelength, background_corrected_spectra, excitation_wavelength) %>%
    mutate(method = "whole_cell")

spectra_2 <-
    import("/home/sam/current_analysis/kir_sur_fret/unroofed", "spectra") %>%
    filter(wavelength > 450, wavelength < 700) %>%
    select(n, construct, wavelength, background_corrected_spectra, excitation_wavelength) %>%
    mutate(method = "unroofed")

spectra <- bind_rows(spectra_1, spectra_2)

ggplot(spectra, aes(x = wavelength, y = background_corrected_spectra, colour = factor(excitation_wavelength), group = interaction(excitation_wavelength, n))) +
geom_line() +
facet_grid(rows = vars(construct), cols = vars(method)) +
geom_vline(xintercept = 516) +
geom_vline(xintercept = 518)

gfp_peak <-
    spectra %>%
    filter(wavelength > 516, wavelength < 518) %>%
    spread(excitation_wavelength, background_corrected_spectra) %>%
    group_by(n, construct, method) %>%
    summarise(gfp_intensity_max = mean(`490`))

spectra <-
    spectra %>%
    spread(excitation_wavelength, background_corrected_spectra) %>%
    left_join(gfp_peak) %>%
    group_by(n, construct, method) %>%
    mutate(
        norm_490 = `490` / gfp_intensity_max,
        norm_565 = `565` / gfp_intensity_max
        )

gfp_idealised <-
    spectra %>%
    filter(construct == "WT-GFP+SUR") %>%
    group_by(wavelength) %>%
    summarise(gfp_490 = mean(norm_490))

spectra <-
    left_join(spectra, gfp_idealised) %>%
    group_by(n, construct, method) %>%
    mutate(
        norm_490_subtracted = norm_490 - gfp_490,
        `490_subtracted` = norm_490_subtracted * gfp_intensity_max,
        ratio_490_565 = `490_subtracted` / `565`
        ) %>%
    drop_na() %>%
    filter(ratio_490_565 < Inf) %>%
    ungroup() %>%
    mutate(construct = case_when(
        construct == "W311ANAP-GFP+SUR_mOrange" ~ "W311*-GFP+SUR_mOrange",
        construct == "W311ANAP-GFP+TMD0_232_mOrange" ~ "W311*-GFP+TMD0_232_mOrange",
        TRUE ~ construct
        )
    )

spectra_summary <-
    spectra %>%
    group_by(construct, method, wavelength) %>%
    summarise(
        se = sd(ratio_490_565) / sqrt(length(ratio_490_565)),
        mean = mean(ratio_490_565)
        ) %>%
    ungroup()

write_csv(spectra, "/home/sam/previous_analysis/2019_manuscript/data/gfp_ofp_fret_spectra.csv")

ratio <-
    spectra %>%
    filter(wavelength > 575, wavelength < 615) %>%
    group_by(n, construct, method) %>%
    summarise(ratio_490_565 = mean(ratio_490_565))

morange_ratio <-
    ratio %>%
    filter(construct == "SUR_mOrange") %>%
    ungroup() %>%
    summarise(ratio = mean(ratio_490_565)) %>%
    pull(ratio)

ratio <-
    ratio %>%
    mutate(ratio = ratio_490_565 / morange_ratio )

ratio_summary <-
    ratio %>%
    group_by(construct, method) %>%
    summarise(
        se = sd(ratio) / sqrt(length(ratio)),
        ratio = mean(ratio)
        )

ggplot(spectra) +
geom_line(aes(x = wavelength, y = `490`, group = interaction(construct, n)), colour = "black") +
geom_line(aes(x = wavelength, y = `490_subtracted`, group = interaction(construct, n)), colour = "red") +
geom_line(aes(x = wavelength, y = `565`, group = interaction(construct, n)), colour = "green") +
facet_grid(rows = vars(method), cols = vars(construct), scales = "free")

ggplot(spectra) +
geom_line(aes(x = wavelength, y = `490`, group = interaction(construct, n, method), colour = construct), linetype = 1) +
geom_line(aes(x = wavelength, y = `565`, group = interaction(construct, n, method), colour = construct), linetype = 2)

ggplot() +
geom_ribbon(
    data = spectra_summary %>% filter(construct == "SUR_mOrange") %>% select(-construct),
    aes(
        x = wavelength,
        ymin = mean - se,
        ymax = mean + se
        ),
    fill = "grey80",
    alpha = 0.25
    ) +
geom_line(
    data = spectra_summary %>% filter(construct == "SUR_mOrange"),
    aes(
        x = wavelength,
        y = mean
        ),
    colour = "grey80"
    ) +
geom_smooth(
    data = spectra %>% filter(construct != "WT-GFP+SUR", construct != "SUR_mOrange"),
    aes(
        x = wavelength,
        y = ratio_490_565,
        colour = construct,
        fill = construct,
        group = interaction(construct, method, n)
        )
    ) +
xlim(575, 615) +
ylim(0, 2) +
facet_grid(construct~method)

ggplot(spectra_summary %>% filter(construct != "WT-GFP+SUR")) +
geom_ribbon(aes(x = wavelength, ymin = mean-se, ymax = mean+se, fill = construct, linetype = method), alpha = 0.25) +
geom_line(aes(x = wavelength, y = mean, colour = construct, linetype = method)) +
xlim(575, 615) +
ylim(0, 1.6) +
cowplot::theme_cowplot()

ggplot() +
geom_point(data = ratio %>% filter(construct != "WT-GFP+SUR"), aes(x = method, y = ratio, fill = construct), size = 2, shape = 21, position = position_dodge(width = 1)) +
geom_pointrange(data = ratio_summary %>% filter(construct != "WT-GFP+SUR"), aes(x = method, ymin = ratio-se, y = ratio, ymax = ratio+se, fill = construct), size = 2.5, shape = 21, position = position_dodge(width = 1)) +
cowplot::theme_cowplot()
