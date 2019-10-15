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
                raw_intensity = col_double(),
                corrected_intensity = col_double()
            )
        ) %>%
        bind_rows(.id = "n") %>%
        type_convert()

    setwd(oldwd)

    return(tempdf)
}

data_1 <-
  import("/home/sam/current_analysis/paper/tnp-atp_doseresponse_unroofed/W311-GFP+SUR", "intensities") %>%
  mutate(construct = "W311*-GFP+SUR")
data_2 <-
  import("/home/sam/current_analysis/paper/tnp-atp_doseresponse_unroofed/W311-GFP", "intensities") %>%
  mutate(construct = "W311*-GFP")
data_3 <-
  import("/home/sam/current_analysis/paper/tnp-atp_doseresponse_unroofed/W311-G334D-GFP+SUR", "intensities") %>%
  mutate(construct = "W311*,G334D-GFP+SUR")
data_5 <-
  import("/home/sam/current_analysis/paper/tnp-atp_doseresponse_unroofed/W311-C166S-GFP+SUR", "intensities") %>%
  mutate(construct = "W311*,C166S-GFP+SUR")
data_6 <-
  import("/home/sam/current_analysis/paper/tnp-atp_doseresponse_unroofed/W311-GFP+TMD0_232", "intensities") %>%
  mutate(construct = "W311*-GFP+TMD0_232")
data_7 <-
  import("/home/sam/current_analysis/paper/tnp-atp_doseresponse_unroofed/W311-GFP+TMD0_195", "intensities") %>%
  mutate(construct = "W311*-GFP+TMD0_195")
data_8 <-
  import("/home/sam/current_analysis/paper/tnp-atp_doseresponse_unroofed/W311-GFP+SUR-K205E", "intensities") %>%
  mutate(construct = "W311*-GFP+SUR-K205E")

unroofed <-
    bind_rows(data_1, data_2, data_3, data_5, data_6, data_7, data_8) %>%
    mutate(
      construct = factor(construct),
      measure = "fluorescence",
      method = "unroofed",
      response = 1 - response
      ) %>%
    group_by(construct, n) %>%
    mutate(unique_experiment_id = group_indices())

write_csv(unroofed, "/home/sam/previous_analysis/2019_manuscript/data/unroofed_concresp_data.csv", col_names = TRUE)
