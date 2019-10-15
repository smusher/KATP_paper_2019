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
  import("/home/sam/current_analysis/paper/tnp-atp_kinetics_unroofed/W311-GFP+SUR/100ms_exposure_3s_interval_60s_wash", "intensities") %>%
  mutate(construct = "W311*-GFP+SUR")
data_2 <-
  import("/home/sam/current_analysis/paper/tnp-atp_kinetics_unroofed/W311-GFP", "intensities") %>%
  mutate(construct = "W311*-GFP")

time_course <-
  bind_rows(data_1, data_2) %>%
  group_by(construct, n, dye) %>%
  mutate(
    response = 1 - (corrected_intensity / corrected_intensity[1L]),
    concentration = ifelse(is.na(concentration), 0, concentration),
    test_concentration = max(concentration),
    time = cumulative_exposure * 0.03,
    repetition = case_when(
      time < 120 ~ 1,
      time >= 120 & time < 240 ~ 2,
      TRUE ~ 3
      ),
    washout = case_when(
      time >= 90 & time < 120 ~ TRUE,
      time >= 210 & time < 240 ~ TRUE,
      time >= 330 & time < 360 ~ TRUE,
      time >= TRUE ~ FALSE),
    washon = case_when(
      time >= 30 & time < 60 ~ TRUE,
      time >= 150 & time < 180 ~ TRUE,
      time >= 270 & time < 300 ~ TRUE,
      TRUE ~ FALSE),
    steady_state = case_when(
      washout == TRUE | washon == TRUE ~ FALSE,
      TRUE ~ TRUE
      ),
    unique_experiment_id = group_indices()
    ) %>%
  filter(test_concentration > 1e-6, test_concentration < 3e-4)

write_csv(time_course, "/home/sam/previous_analysis/2019_manuscript/data/unroofed_timecourse_data.csv", col_names = TRUE)
