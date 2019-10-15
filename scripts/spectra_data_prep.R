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
                wavelength = col_double(),
                dye = col_double(),
                raw_spectra = col_double(),
                background_corrected_spectra = col_double(),
                bleaching_corrected_spectra = col_double()
            )
        ) %>%
        bind_rows(.id = "n") %>%
        type_convert()

    setwd(oldwd)

    return(tempdf)
}

data_3 <-
  import("/home/sam/current_analysis/paper/tnp-atp_doseresponse_unroofed/W311-GFP+SUR", "spectra") %>%
  mutate(construct = "W311*-GFP+SUR", method = "unroofed")
data_4 <-
  import("/home/sam/current_analysis/paper/tnp-atp_doseresponse_unroofed/W311-GFP", "spectra") %>%
  mutate(construct = "W311*-GFP", method = "unroofed")
data_5 <-
  import("/home/sam/current_analysis/paper/tnp-atp_doseresponse_unroofed/W311-G334D-GFP+SUR", "spectra") %>%
  mutate(construct = "W311*,G334D-GFP+SUR", method = "unroofed")
data_6 <-
  import("/home/sam/current_analysis/paper/tnp-atp_doseresponse_unroofed/W311-C166S-GFP+SUR", "spectra") %>%
  mutate(construct = "W311*,C166S-GFP+SUR", method = "unroofed")
data_7 <-
  import("/home/sam/current_analysis/paper/tnp-atp_doseresponse_unroofed/W311-GFP+TMD0_232", "spectra") %>%
  mutate(construct = "W311*-GFP+TMD0_232", method = "unroofed")
data_8 <-
  import("/home/sam/current_analysis/paper/tnp-atp_doseresponse_unroofed/W311-GFP+TMD0_195", "spectra") %>%
  mutate(construct = "W311*-GFP+TMD0_195", method = "unroofed")
data_9 <-
  import("/home/sam/current_analysis/paper/tnp-atp_doseresponse_unroofed/W311+SUR", "spectra") %>%
  mutate(construct = "W311*+SUR", method = "unroofed")
data_10 <-
  import("/home/sam/current_analysis/paper/tnp-atp_doseresponse_unroofed/WT-GFP+SUR", "spectra") %>%
  mutate(construct = "WT-GFP+SUR", method = "unroofed")
data_11 <-
  import("/home/sam/current_analysis/paper/pcf/W311-GFP+SUR/reanalyse", "spectra") %>%
  mutate(construct = "W311*-GFP+SUR", method = "combined")
data_12 <-
  import("/home/sam/current_analysis/paper/pcf/W311-C166S-GFP+SUR/reanalyse", "spectra") %>%
  mutate(construct = "W311*,C166S-GFP+SUR", method = "combined")
data_13 <-
  import("/home/sam/current_analysis/paper/pcf/WT-GFP+SUR", "spectra") %>%
  mutate(construct = "WT-GFP+SUR", method = "combined")
data_14 <-
  import("/home/sam/current_analysis/paper/pcf/WT-GFP", "spectra") %>%
  mutate(construct = "WT-GFP", method = "combined")

spectra <-
  bind_rows(data_3, data_4, data_5, data_6, data_7, data_8, data_9, data_10, data_11, data_12, data_13, data_14) %>%
  filter(dye < 480 | is.na(dye), wavelength > 405) %>%
  mutate(concentration = case_when(is.na(concentration) ~ 0, TRUE ~ concentration))

summed_spectra <-
  spectra %>%
  filter(cumulative_exposure == 0) %>%
  group_by(construct, method, wavelength) %>%
  summarise(
    raw_spectra = sum(raw_spectra),
    background_corrected_spectra = sum(background_corrected_spectra),
    bleaching_corrected_spectra = sum(bleaching_corrected_spectra)
    ) %>%
  group_by(construct, method) %>%
  mutate(
    normalised_raw_spectra = raw_spectra / max(raw_spectra),
    normalised_background_corrected_spectra = background_corrected_spectra / max(background_corrected_spectra),
    normalised_bleaching_corrected_spectra = bleaching_corrected_spectra / max(bleaching_corrected_spectra)
    )

unroofed <-
  summed_spectra %>%
  filter(method == "unroofed")

unroofed$gfp_spectra <-
  unroofed %>%
  filter(construct == "WT-GFP+SUR") %>%
  pull(normalised_background_corrected_spectra)

unroofed <-
  unroofed %>%
  mutate(subtracted_spectra = normalised_background_corrected_spectra - gfp_spectra)

combined <-
  summed_spectra %>%
  filter(method == "combined")

combined$gfp_spectra <-
  combined %>%
  filter(construct == "WT-GFP+SUR") %>%
  pull(normalised_background_corrected_spectra)

combined <-
  combined %>%
  mutate(subtracted_spectra = normalised_background_corrected_spectra - gfp_spectra)

subtracted_spectra <-
  bind_rows(unroofed, combined) %>%
  mutate(normalised_subtracted_spectra = subtracted_spectra / max(subtracted_spectra)) %>%
  right_join(summed_spectra)

average_spectra <-
  spectra %>%
  group_by(construct, method, concentration, wavelength) %>%
  summarise(
    se_bleaching_corrected_spectra = sd(bleaching_corrected_spectra) / sqrt(length(bleaching_corrected_spectra)),
    avg_bleaching_corrected_spectra = mean(bleaching_corrected_spectra)
    ) %>%
  group_by(construct, method) %>%
  mutate(
    norm_bleaching = avg_bleaching_corrected_spectra / max(avg_bleaching_corrected_spectra)
    )

write_csv(subtracted_spectra, "/home/sam/previous_analysis/2019_manuscript/data/summed_spectra.csv", col_names = TRUE)
write_csv(average_spectra, "/home/sam/previous_analysis/2019_manuscript/data/average_spectra.csv", col_names = TRUE)
