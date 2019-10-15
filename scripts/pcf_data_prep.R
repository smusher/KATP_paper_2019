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
                concentration = col_double(),
                dye = col_double(),
                response = col_double(),
                raw_current = col_double(),
                wash_current = col_double(),
                barium_current = col_double(),
                raw_intensity = col_double(),
                corrected_intensity = col_double()
            )
        ) %>%
        bind_rows(.id = "n") %>%
        type_convert()

    setwd(oldwd)

    return(tempdf)
}

w311_current <-
    import("/home/sam/current_analysis/paper/pcf/W311-GFP+SUR/reanalyse", "currents") %>%
    mutate(
        construct = "W311*-GFP+SUR",
        measure = "current",
        )

w311c166s_current <-
    import("/home/sam/current_analysis/paper/pcf/W311-C166S-GFP+SUR/reanalyse", "currents") %>%
    mutate(
        construct = "W311*,C166S-GFP+SUR",
        measure = "current"
        )

w311_surk205e_current <-
    import("/home/sam/current_analysis/paper/pcf/W311-GFP+SUR-K205E", "currents") %>%
    mutate(
        construct = "W311*-GFP+SUR-K205E",
        measure = "current"
        )

w311_surk205a_current <-
    import("/home/sam/current_analysis/paper/pcf/W311-GFP+SUR-K205A", "currents") %>%
    mutate(
        construct = "W311*-GFP+SUR-K205A",
        measure = "current"
        )

w311_fluorescence <-
    import("/home/sam/current_analysis/paper/pcf/W311-GFP+SUR/reanalyse", "intensities") %>%
    mutate(
        construct = "W311*-GFP+SUR",
        measure = "fluorescence",
        response = 1 - response
        )

w311c166s_fluorescence <-
    import("/home/sam/current_analysis/paper/pcf/W311-C166S-GFP+SUR/reanalyse", "intensities") %>%
    mutate(
        construct = "W311*,C166S-GFP+SUR",
        measure = "fluorescence",
        response = 1 - response
        )

w311_surk205e_fluorescence <-
    import("/home/sam/current_analysis/paper/pcf/W311-GFP+SUR-K205E", "intensities") %>%
    mutate(
        construct = "W311*-GFP+SUR-K205E",
        measure = "fluorescence",
        response = 1 - response
        )

w311_surk205a_fluorescence <-
    import("/home/sam/current_analysis/paper/pcf/W311-GFP+SUR-K205A", "intensities") %>%
    mutate(
        construct = "W311*-GFP+SUR-K205A",
        measure = "fluorescence",
        response = 1 - response
        )

combined <-
    bind_rows(w311_current, w311c166s_current, w311_surk205e_current, w311_surk205a_current, w311_fluorescence, w311c166s_fluorescence, w311_surk205e_fluorescence, w311_surk205a_fluorescence) %>%
    mutate(
        construct = factor(construct),
        method = "pcf"
        ) %>%
    group_by(construct, n) %>%
    mutate(unique_experiment_id = group_indices())

write_csv(combined, "/home/sam/previous_analysis/2019_manuscript/data/pcf_data.csv", col_names = TRUE)

combined_summary <-
    combined %>%
    group_by(construct, measure, concentration, dye) %>%
    summarise(
        se = sd(response)/sqrt(length(response)),
        response = mean(response)
        )

ggplot() +
geom_point(data = combined %>% filter(dye < 480 | is.na(dye), concentration > 0), aes(x = concentration, y = response, fill = measure), shape = 21, alpha = 0.25, size = 2) +
geom_errorbar(data = combined_summary %>% filter(dye < 480 | is.na(dye), concentration > 0), aes(x = concentration, ymin = response - se, ymax = response + se), width = 0.2) +
geom_point(data = combined_summary %>% filter(dye < 480 | is.na(dye), concentration > 0), aes(x = concentration, y = response, fill = measure), shape = 21, size = 3) +
facet_grid(rows = vars(construct)) +
cowplot::theme_cowplot() +
scale_colour_brewer(palette = "Set1") +
scale_x_log10()
