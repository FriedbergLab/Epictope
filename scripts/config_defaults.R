# define species to run MSA against
species <- c("bos_taurus", "canis_lupus_familiaris", "gallus_gallus", "homo_sapiens", "mus_musculus", "takifugu_rubripes", "xenopus_tropicalis")
assign("species", species, envir = .GlobalEnv)

# weights for tagging features
h_weight = 1.5 # shannon entropy
rsa_weight = 1 # solvent accessible surface area
ss_weight = 1 # secondary structure
br_weight = 1 # disordered binding region
assign("h_weight", h_weight, envir = .GlobalEnv)
assign("rsa_weight", rsa_weight, envir = .GlobalEnv)
assign("ss_weight", ss_weight, envir = .GlobalEnv)
assign("br_weight", br_weight, envir = .GlobalEnv)

# value for secondary structures, must be 0-1.
# each letter refers to a type of secondary structure
# the number indicates the value or "suitability" for tagging.
# values should be from 0-1, with higher values indicating greater
# suitability for tagging.
ss_key <- list(
    "G" = 0,
    "H" = 0,
    "I" = 0,
    "E" = 0,
    "C" = 1,
    "T" = 0.5,
    "B" = 0.5,
    "S" = 0.5,
    "P" = 0)
assign("ss_key", ss_key, envir = .GlobalEnv)

# reference values for maximum solvent accessibility of amino acids.
# default values estimate from the following study;
# https://doi.org/10.1371/journal.pone.0080635
aa <- c("A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")
empirical <- c(121, 265, 187, 187, 148, 214, 214, 97, 216, 195, 191, 230, 203, 228, 154, 143, 163, 264, 255, 165)
max_sasa <- setNames(empirical, aa)
assign("max_sasa", max_sasa, envir = .GlobalEnv)
