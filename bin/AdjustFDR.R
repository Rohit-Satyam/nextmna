suppressPackageStartupMessages(library("argparser"))
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("plyranges"))
suppressPackageStartupMessages(library("fs"))



# PARSE ARGUMENTS -------------------------------------------------------------

## INTERACTIVE ARGUMENT PARSING -----------------------------------------------
# Just change the arguments here and source the script.
# project_path <-     # Input path. Must be a path to DANPOS output.
# output_dir <- paste0(project_path, "./adjustedFDR/") # Output directory
# chromosome_sizes_file <-   # Chrom.sizes file for wig to bigwig conversion.
n_simulations <- NA # Number of simulated p-values for FDR estimation. [default: 100x fdr-values to estimate]
sampling_type_occupancy <- "near" # sampling type for fdr estimation of occupancy change fdr. Either 'same', 'random' or 'near'. [default: "near"]
sampling_type_fuzziness <- "near" # sampling type for fdr estimation of fuzziness change fdr. Either 'same', 'random' or 'near'. [default: "near"]
# sample_names <-     # names of samples seperated by ','. example: 'sampleA,sampleB,sampleC'. If not provided, names are searched automaticly.
# control_names <-    # names of controls seperated by ','. example: 'controlA,controlB'. If not provided, names are searched automaticly.

## COMMAND LINE PARSING -------------------------------------------------------
p <- arg_parser("calculate correct FDR values from DANPOS output")
p <- add_argument(p, "input_path", short = "-i", help = "Input path. Must be a path to DANPOS output.", default = ".")
p <- add_argument(p, "--output_dir", short = "-o", help = "Output directory.", default = "./adjustedFDR/", type = "character")
p <- add_argument(p, "--chr_sizes", short = "-c", help = "Chrom.sizes file for wig to bigwig conversion.", type = "character")
p <- add_argument(p, "--n_simulations", short = "-n", help = "Number of simulated p-values for FDR estimation. [default: 100x fdr-values to estimate]", default = NA, type = "numeric")
p <- add_argument(p, "--sampling_type_occupancy", short = "-sto", help = "sampling type for fdr estimation of occupancy change fdr. Either 'same', 'random' or 'near'.", default = "near")
p <- add_argument(p, "--sampling_type_fuzziness", short = "-stf", help = "sampling type for fdr estimation of fuzziness change fdr. Either 'same', 'random' or 'near'.", default = "near")
p <- add_argument(p, "--sample_names", short = "-s", help = "names of samples seperated by ','. example: 'sampleA,sampleB,sampleC'. If not provided, names are searched automaticly.", type = "character", default = NA)
p <- add_argument(p, "--control_names", short = "-l", help = "names of controls seperated by ','. example: 'controlA,controlB'. If not provided, names are searched automaticly.", type = "character", default = NA)


if (!interactive()) {
  this_file <- commandArgs() %>%
    tibble::enframe(name = NULL) %>%
    tidyr::separate(
      col = value,
      into = c("key", "value"),
      sep = "=",
      fill = "right"
    ) %>%
    dplyr::filter(key == "--file") %>%
    dplyr::pull(value)
  source(paste0(dirname(this_file), "/functions.R"))
  argv <- parse_args(p)
} else {
  source("functions.R")
  argv <- parse_args(p, argv = project_path)
  argv$output_dir <- output_dir
  argv$chr_sizes <- chromosome_sizes_file
  argv$n_simulations <- n_simulations
  argv$sampling_type_occupancy <- sampling_type_occupancy
  argv$sampling_type_fuzziness <- sampling_type_fuzziness
  argv$sample_names <- sample_names
  argv$control_names <- control_names
}

# MAIN -------------------------------------------------------------------------
dir.create(paste0(argv$output_dir), showWarnings = FALSE, recursive = TRUE)

print(argv$input_path)
print(argv$sample_names)
print(argv$control_names)

names_list <- get_names(
  project_path = argv$input_path,
  sample_names = argv$sample_names,
  control_names = argv$control_names
)

for (i in seq_along(names_list$control_names)) {
  walk(names_list$sample_names,
    calc_FDR,
    control_name = names_list$control_names[i],
    outPath = argv$output_dir,
    inPath = argv$input_path,
    sampling_n = argv$n_simulations,
    chromosome_sizes_file = argv$chr_sizes,
    sampling_type_occ = argv$sampling_type_occupancy,
    sampling_type_fuzz = argv$sampling_type_fuzziness
  )
}
