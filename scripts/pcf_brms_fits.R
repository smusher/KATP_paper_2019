library(tidyverse)
library(broom)
library(brms)
library(rstan)
library(tidybayes)
library(future)
source("scripts/specify_model_formulae.R")
brms_iter <- 30000
brms_warmup <- 20000
brms_chains <- 4
brms_thin <- 10
brms_seed <- 2019

pcf_data <-
    read_csv("data/pcf_data.csv") %>%
    filter(dye < 480 | is.na(dye), concentration > 0) %>%
    select(unique_experiment_id, construct, measure, method, concentration, response) %>%
    mutate(
        response = case_when(measure == "fluorescence" ~ 1 - response, TRUE ~ response),
        binding_mask = case_when(measure == "fluorescence" ~ 1, TRUE ~ 0)
    ) %>%
    group_by(unique_experiment_id, construct, measure) %>%
    mutate(response_normalised = case_when(
        measure == "fluorescence" ~ (response-min(response))/(max(response)-min(response)),
        measure == "current" ~ response
    )) %>%
    ungroup()

unroofed_data <-
    read_csv("data/unroofed_concresp_data.csv") %>%
    filter(dye < 480, concentration > 0) %>%
    select(unique_experiment_id, construct, measure, method, concentration, response) %>%
    filter(construct %in% c("W311*-GFP+SUR", "W311*,C166S-GFP+SUR", "W311*-GFP+SUR-K205A", "W311*-GFP+SUR-K205E"), method == "unroofed") %>%
    mutate(
        response = 1 - response,
        binding_mask = 1
    )

all_data <-
    bind_rows(pcf_data, unroofed_data)

mwc_formula_transformed <-
    mwc_formula %>%
    str_replace_all("Ka", "(10^logKa)") %>%
    str_replace_all("L", "(10^logL)")

single_formula_transformed <-
    single_formula %>%
    str_replace_all("Ka", "(10^logKa)") %>%
    str_replace_all("L", "(10^logL)")

coop_formula_transformed <-
    mwc_formula_with_cooperativity %>%
    str_replace_all("Ka", "(10^logKa)") %>%
    str_replace_all("L", "(10^logL)")

mwc_formula_fixed <-
    mwc_formula %>%
    str_replace_all("Ka", "(10^logKa)") %>%
    str_replace_all("L", "0.8")

single_formula_fixed <-
    single_formula %>%
    str_replace_all("Ka", "(10^logKa)") %>%
    str_replace_all("L", "0.8")

coop_formula_fixed <-
    mwc_formula_with_cooperativity %>%
    str_replace_all("Ka", "(10^logKa)") %>%
    str_replace_all("L", "0.8")

mwc_formula_fixed_cs <-
    mwc_formula %>%
    str_replace_all("Ka", "(10^logKa)") %>%
    str_replace_all("L", "6.01")

single_formula_fixed_cs <-
    single_formula %>%
    str_replace_all("Ka", "(10^logKa)") %>%
    str_replace_all("L", "6.01")

coop_formula_fixed_cs <-
    mwc_formula_with_cooperativity %>%
    str_replace_all("Ka", "(10^logKa)") %>%
    str_replace_all("L", "6.01")

mwc_priors <-
    c(
        prior(uniform(0, 1), nlpar = "D", lb = 0, ub = 1),
        prior(uniform(2, 6), nlpar = "logKa", lb = 2, ub = 6),
        prior(normal(0, 0.7), nlpar = "logL")
    )

coop_priors <-
    c(
        prior(uniform(0, 1), nlpar = "D", lb = 0, ub = 1),
        prior(uniform(2, 6), nlpar = "logKa", lb = 2, ub = 6),
        prior(normal(0, 0.7), nlpar = "logL"),
        prior(uniform(0, 1), nlpar = "c", lb = 0, ub = 1)
    )

mwc_fixed_priors <-
    c(
        prior(uniform(0, 1), nlpar = "D", lb = 0, ub = 1),
        prior(uniform(2, 6), nlpar = "logKa", lb = 2, ub = 6)
    )

coop_fixed_priors <-
    c(
        prior(uniform(0, 1), nlpar = "D", lb = 0, ub = 1),
        prior(uniform(2, 6), nlpar = "logKa", lb = 2, ub = 6),
        prior(uniform(0, 1), nlpar = "c", lb = 0, ub = 1)
    )

mwc_population <-
    bf(
        as.formula(mwc_formula_transformed),
        logL ~ 1,
        D ~ 1,
        logKa ~ 1,
        nl = TRUE
    )

coop_population <-
    bf(
        as.formula(coop_formula_transformed),
        logL ~ 1,
        D ~ 1,
        logKa ~ 1,
        c ~ 1,
        nl = TRUE
    )

single_population <-
    bf(
        as.formula(single_formula_transformed),
        logL ~ 1,
        D ~ 1,
        logKa ~ 1,
        nl = TRUE
    )

# mwc_individual <-
#     bf(
#         as.formula(mwc_formula_transformed),
#         logL ~ 1 + (1|unique_experiment_id),
#         D ~ 1,
#         logKa ~ 1,
#         nl = TRUE
#     )

# coop_individual <-
#     bf(
#         as.formula(coop_formula_transformed),
#         logL ~ 1 + (1|unique_experiment_id),
#         D ~ 1,
#         logKa ~ 1,
#         c ~ 1,
#         nl = TRUE
#     )

# single_individual <-
#     bf(
#         as.formula(single_formula_transformed),
#         logL ~ 1 + (1|unique_experiment_id),
#         D ~ 1,
#         logKa ~ 1,
#         nl = TRUE
#     )

mwc_fixed <-
    bf(
        as.formula(mwc_formula_fixed),
        D ~ 1,
        logKa ~ 1,
        nl = TRUE
    )

coop_fixed <-
    bf(
        as.formula(coop_formula_fixed),
        D ~ 1,
        logKa ~ 1,
        c ~ 1,
        nl = TRUE
    )

single_fixed <-
    bf(
        as.formula(single_formula_fixed),
        D ~ 1,
        logKa ~ 1,
        nl = TRUE
    )

mwc_fixed_cs <-
    bf(
        as.formula(mwc_formula_fixed_cs),
        D ~ 1,
        logKa ~ 1,
        nl = TRUE
    )

coop_fixed_cs <-
    bf(
        as.formula(coop_formula_fixed_cs),
        D ~ 1,
        logKa ~ 1,
        c ~ 1,
        nl = TRUE
    )

single_fixed_cs <-
    bf(
        as.formula(single_formula_fixed_cs),
        D ~ 1,
        logKa ~ 1,
        nl = TRUE
    )

do_fit <- function(construct_name, formula, priors){
    model_name <- deparse(substitute(formula))
    save_name <- paste(construct_name, model_name, sep = "_")
    print(save_name)
    data <- pcf_data %>% filter(construct == construct_name)
    fit <-
        brm(
            formula = formula,
            data = data,
            prior = priors,
            family = gaussian(),
            iter = brms_iter,
            warmup = brms_warmup,
            chains = brms_chains,
            thin = brms_thin,
            seed = brms_seed,
            control = list(adapt_delta = 0.99, max_treedepth = 15),
            file = paste("data/mwc_fits/", save_name, sep = ""),
            sample_prior = "yes",
            future = TRUE
        )
    return(fit)
}

mwc_wt_pop <- do_fit("W311*-GFP+SUR", mwc_population, mwc_priors)
# mwc_wt_ind <- do_fit("W311*-GFP+SUR", mwc_individual, mwc_priors)
mwc_wt_fix <- do_fit("W311*-GFP+SUR", mwc_fixed, mwc_fixed_priors)
mwc_cs_pop <- do_fit("W311*,C166S-GFP+SUR", mwc_population, mwc_priors)
# mwc_cs_ind <- do_fit("W311*,C166S-GFP+SUR", mwc_individual, mwc_priors)
mwc_cs_fix <- do_fit("W311*,C166S-GFP+SUR", mwc_fixed_cs, mwc_fixed_priors)
mwc_ke_pop <- do_fit("W311*-GFP+SUR-K205E", mwc_population, mwc_priors)
# mwc_ke_ind <- do_fit("W311*-GFP+SUR-K205E", mwc_individual, mwc_priors)
mwc_ke_fix <- do_fit("W311*-GFP+SUR-K205E", mwc_fixed, mwc_fixed_priors)
mwc_ka_pop <- do_fit("W311*-GFP+SUR-K205A", mwc_population, mwc_priors)
# mwc_ka_ind <- do_fit("W311*-GFP+SUR-K205A", mwc_individual, mwc_priors)
mwc_ka_fix <- do_fit("W311*-GFP+SUR-K205A", mwc_fixed, mwc_fixed_priors)

single_wt_pop <- do_fit("W311*-GFP+SUR", single_population, mwc_priors)
# single_wt_ind <- do_fit("W311*-GFP+SUR", single_individual, mwc_priors)
single_wt_fix <- do_fit("W311*-GFP+SUR", single_fixed, mwc_fixed_priors)
single_cs_pop <- do_fit("W311*,C166S-GFP+SUR", single_population, mwc_priors)
# single_cs_ind <- do_fit("W311*,C166S-GFP+SUR", single_individual, mwc_priors)
single_cs_fix <- do_fit("W311*,C166S-GFP+SUR", single_fixed_cs, mwc_fixed_priors)
single_ke_pop <- do_fit("W311*-GFP+SUR-K205E", single_population, mwc_priors)
# single_ke_ind <- do_fit("W311*-GFP+SUR-K205E", single_individual, mwc_priors)
single_ke_fix <- do_fit("W311*-GFP+SUR-K205E", single_fixed, mwc_fixed_priors)
single_ka_pop <- do_fit("W311*-GFP+SUR-K205A", single_population, mwc_priors)
# single_ka_ind <- do_fit("W311*-GFP+SUR-K205A", single_individual, mwc_priors)
single_ka_fix <- do_fit("W311*-GFP+SUR-K205A", single_fixed, mwc_fixed_priors)

coop_wt_pop <- do_fit("W311*-GFP+SUR", coop_population, coop_priors)
# coop_wt_ind <- do_fit("W311*-GFP+SUR", coop_individual, coop_priors)
coop_wt_fix <- do_fit("W311*-GFP+SUR", coop_fixed, coop_fixed_priors)
coop_cs_pop <- do_fit("W311*,C166S-GFP+SUR", coop_population, coop_priors)
# coop_cs_ind <- do_fit("W311*,C166S-GFP+SUR", coop_individual, coop_priors)
coop_cs_fix <- do_fit("W311*,C166S-GFP+SUR", coop_fixed_cs, coop_fixed_priors)
coop_ke_pop <- do_fit("W311*-GFP+SUR-K205E", coop_population, coop_priors)
# coop_ke_ind <- do_fit("W311*-GFP+SUR-K205E", coop_individual, coop_priors)
coop_ke_fix <- do_fit("W311*-GFP+SUR-K205E", coop_fixed, coop_fixed_priors)
coop_ka_pop <- do_fit("W311*-GFP+SUR-K205A", coop_population, coop_priors)
# coop_ka_ind <- do_fit("W311*-GFP+SUR-K205A", coop_individual, coop_priors)
coop_ka_fix <- do_fit("W311*-GFP+SUR-K205A", coop_fixed, coop_fixed_priors)
