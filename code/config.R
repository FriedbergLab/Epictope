# folder structure
dataFolder <- "data"
outputFolder <- "outputs"
model_folder <- paste0(dataFolder, "/models")
cds_folder <- paste0(dataFolder, "/CDS")

# inputs
# ex: "O57472"
query <- "O57472"

# define species to run MSA against
species <- c("bos_taurus", "canis_lupus_familiaris", "gallus_gallus", "homo_sapiens", "mus_musculus", "takifugu_rubripes", "xenopus_tropicalis")

# weights for tagging features
h_weight = 1.5
rsa_weight = 1
ss_weight = 1
br_weight = 1

# value for secondary structures, must be 0-1.
ss_key <- list(
    "G" = 0, 
    "H" = 0, 
    "I" = 0,
    "E" = 0, 
    "C" = 1,
    "T" = 0.5, 
    "B" = 0.5, 
    "S" = 0.5)

# reference values for maximum solvent accesibility.
# https://doi.org/10.1371/journal.pone.0080635
aa <- c("A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")
empirical <- c(121, 265, 187, 187, 148, 214, 214, 97, 216, 195, 191, 230, 203, 228, 154, 143, 163, 264, 255, 165)
max_sasa <- setNames(empirical, aa)

