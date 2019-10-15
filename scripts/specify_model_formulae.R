library(tidyverse)

all_states <-
    list(
        "(1)",
        "(4 * Ka * concentration)",
        "(6 * Ka^2 * concentration^2)",
        "(4 * Ka^3 * concentration^3)",
        "(Ka^4 * concentration^4)",
        "(L)",
        "(L * D * 4 * Ka * concentration)",
        "(L * D^2 * 6 * Ka^2 * concentration^2)",
        "(L * D^3 * 4 * Ka^3 * concentration^3)",
        "(L * D^4 * Ka^4 * concentration^4)"
    )

open_states <-
    list(
        "(L)",
        "(L * D * 4 * Ka * concentration)",
        "(L * D^2 * 6 * Ka^2 * concentration^2)",
        "(L * D^3 * 4 * Ka^3 * concentration^3)",
        "(L * D^4 * Ka^4 * concentration^4)"
    )

bound_states <-
    list(
        "(1 * 4 * Ka * concentration)",
        "(2 * 6 * Ka^2 * concentration^2)",
        "(3 * 4 * Ka^3 * concentration^3)",
        "(4 * Ka^4 * concentration^4)",
        "(1 * L * D * 4 * Ka * concentration)",
        "(2 * L * D^2 * 6 * Ka^2 * concentration^2)",
        "(3 * L * D^3 * 4 * Ka^3 * concentration^3)",
        "(4 * L * D^4 * Ka^4 * concentration^4)"
    )

binding_eq <-
    paste(
        "(",
        paste(bound_states, collapse = " + "),
        ") / (4 * (",
        paste(all_states, collapse = " + "),
        "))",
        sep = ""
    )

gating_eq <-
    paste(
        "((",
        paste(open_states, collapse = " + "),
        ") / (",
        paste(all_states, collapse = " + "),
        ") / (L / (L + 1)))",
        sep = ""
    )

mwc_formula <-
    paste(
        "response ~ (binding_mask * (",
        binding_eq,
        ")) + ((1 - binding_mask) * (",
        gating_eq,
        "))",
        sep = ""
    )

all_states <-
    list(
        "(1)",
        "(4 * Ka * concentration)",
        "(6 * Ka^2 * concentration^2)",
        "(4 * Ka^3 * concentration^3)",
        "(Ka^4 * concentration^4)",
        "(L)",
        "(L * D * 4 * Ka * concentration)",
        "(L * D * 6 * Ka^2 * concentration^2)",
        "(L * D * 4 * Ka^3 * concentration^3)",
        "(L * D * Ka^4 * concentration^4)"
    )

open_states <-
    list(
        "(L)",
        "(L * D * 4 * Ka * concentration)",
        "(L * D * 6 * Ka^2 * concentration^2)",
        "(L * D * 4 * Ka^3 * concentration^3)",
        "(L * D * Ka^4 * concentration^4)"
    )

bound_states <-
    list(
        "(1 * 4 * Ka * concentration)",
        "(2 * 6 * Ka^2 * concentration^2)",
        "(3 * 4 * Ka^3 * concentration^3)",
        "(4 * Ka^4 * concentration^4)",
        "(1 * L * D * 4 * Ka * concentration)",
        "(2 * L * D * 6 * Ka^2 * concentration^2)",
        "(3 * L * D * 4 * Ka^3 * concentration^3)",
        "(4 * L * D * Ka^4 * concentration^4)"
    )

binding_eq <-
    paste(
        "(",
        paste(bound_states, collapse = " + "),
        ") / (4 * (",
        paste(all_states, collapse = " + "),
        "))",
        sep = ""
    )

gating_eq <-
    paste(
        "((",
        paste(open_states, collapse = " + "),
        ") / (",
        paste(all_states, collapse = " + "),
        ") / (L / (L + 1)))",
        sep = ""
    )

single_formula <-
    paste(
        "response ~ (binding_mask * (",
        binding_eq,
        ")) + ((1 - binding_mask) * (",
        gating_eq,
        "))",
        sep = ""
    )

all_states <-
    list(
        "(1)",
        "(4 * Ka * concentration)",
        "(6 * Ka^2 * concentration^2 * c)",
        "(4 * Ka^3 * concentration^3 * c^2)",
        "(Ka^4 * concentration^4 * c^3)",
        "(L)",
        "(L * D * 4 * Ka * concentration)",
        "(L * D^2 * 6 * Ka^2 * concentration^2 * c)",
        "(L * D^3 * 4 * Ka^3 * concentration^3 * c^2)",
        "(L * D^4 * Ka^4 * concentration^4 * c^3)"
    )

open_states <-
    list(
        "(L)",
        "(L * D * 4 * Ka * concentration)",
        "(L * D^2 * 6 * Ka^2 * concentration^2 * c)",
        "(L * D^3 * 4 * Ka^3 * concentration^3 * c^2)",
        "(L * D^4 * Ka^4 * concentration^4 * c^3)"
    )

bound_states <-
    list(
        "(1 * 4 * Ka * concentration)",
        "(2 * 6 * Ka^2 * concentration^2 * c)",
        "(3 * 4 * Ka^3 * concentration^3 * c^2)",
        "(4 * Ka^4 * concentration^4 * c^3)",
        "(1 * L * D * 4 * Ka * concentration)",
        "(2 * L * D^2 * 6 * Ka^2 * concentration^2 * c)",
        "(3 * L * D^3 * 4 * Ka^3 * concentration^3 * c^2)",
        "(4 * L * D^4 * Ka^4 * concentration^4 * c^3)"
    )

binding_eq <-
    paste(
        "(",
        paste(bound_states, collapse = " + "),
        ") / (4 * (",
        paste(all_states, collapse = " + "),
        "))",
        sep = ""
    )

gating_eq <-
    paste(
        "((",
        paste(open_states, collapse = " + "),
        ") / (",
        paste(all_states, collapse = " + "),
        ") / (L / (L + 1)))",
        sep = ""
    )

mwc_formula_with_cooperativity <-
    paste(
        "response ~ (binding_mask * (",
        binding_eq,
        ")) + ((1 - binding_mask) * (",
        gating_eq,
        "))",
        sep = ""
    )
