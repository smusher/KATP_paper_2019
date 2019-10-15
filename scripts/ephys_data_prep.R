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
                concentration = col_double(),
                response = col_double(),
                raw_current = col_double(),
                wash_current = col_double(),
                barium_current = col_double()
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
        nucleotide = "TNP-ATP",
        method = "pcf"
        )

w311c166s_current <-
    import("/home/sam/current_analysis/paper/pcf/W311-C166S-GFP+SUR/reanalyse", "currents") %>%
    mutate(
        construct = "W311*,C166S-GFP+SUR",
        nucleotide = "TNP-ATP",
        method = "pcf"
        )

w311_surk205e_current <-
    import("/home/sam/current_analysis/paper/pcf/W311-GFP+SUR-K205E", "currents") %>%
    mutate(
        construct = "W311*-GFP+SUR-K205E",
        nucleotide = "TNP-ATP",
        method = "pcf"
        )

w311_surk205a_current <-
    import("/home/sam/current_analysis/paper/pcf/W311-GFP+SUR-K205A", "currents") %>%
    mutate(
        construct = "W311*-GFP+SUR-K205A",
        nucleotide = "TNP-ATP",
        method = "pcf"
        )

wt_current <-
    read_csv("/home/sam/current_analysis/paper/pcf/wt_currents.csv", col_types = cols(
        n = col_integer(),
        construct = col_character(),
        nucleotide = col_character(),
        concentration = col_double(),
        response = col_double(),
        tolbutamide_response = col_double()
        )
    ) %>%
    mutate(method = "patch_clamp")

currents <-
    bind_rows(w311_current, w311c166s_current, w311_surk205e_current, w311_surk205a_current, wt_current) %>%
    mutate(
        construct = factor(construct),
        measure = "current"
        ) %>%
    group_by(construct, n, nucleotide) %>%
    mutate(unique_experiment_id = group_indices())

write_csv(currents, "/home/sam/previous_analysis/2019_manuscript/data/electrophys_data.csv", col_names = TRUE)
