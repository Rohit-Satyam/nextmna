# FUNCTIONS -----------------------------------------

#' get sample and control names
#'
#' This function tries to prepare sample and control names for later usage.
#' It either splits the supplied string or tries to automatically detect sample and control names if none are supplied.
#' Further it checks for plausibility of sample/control names
#'
#' @param project_path path to automatically look for sample/control names and plausibility.
#' @param sample_names string of comma separated sample names. Overwrites any automatically detected names.
#' @param control_names string of comma separated control names. Overwrites any automatically detected names.
#' @return Returns a List with a sample_name vector and a control_name vector.
get_names <- function(project_path, sample_names, control_names) {
  autodetect <- 0
  potential_names <- list.files(paste0(project_path, "/pooled"), pattern = ".positions.xls") %>% sub(x = ., "\\..*", "")
  all_names <- list.files(project_path, pattern = ".positions.integrative.xls") %>%
    sub(x = ., "\\..*", "") %>%
    str_split(pattern = "-", simplify = TRUE, n = 2)
  s_names <- all_names[, 1]
  c_names <- unique(all_names[, 2])
  if (is.na(sample_names)) {
    autodetect <- 1
    sample_names <- s_names
  } else {
    sample_names <- sample_names %>%
      str_split(pattern = ",", simplify = TRUE) %>%
      c()
  }
  if (is.na(control_names)) {
    autodetect <- 1
    control_names <- c_names
  } else {
    control_names <- control_names %>%
      str_split(pattern = ",", simplify = TRUE) %>%
      c()
  }
  print(potential_names)
  if (!all(sample_names %in% potential_names) | !all(control_names %in% potential_names)) {
    if (autodetect == 1) {
      stop("There is some issue with sample and contrast names.\n There should be no '-' or '.' in the sample/contrast name when automaticly detecting names automaticly. \n Either rename or provide names manually.")
    } else {
      stop("There is some issue with sample and contrast names.\n Make sure all names provided are correct and present in the DANPOS results.")
    }
  }

  list("sample_names" = sample_names, "control_names" = control_names)
}

#' Load Danpos Bigwig
#'
#' This function loads a pooled bigwig file ("XYZ.Fnor.smooth.bw") and
#' converts it to a GRanges(default) or tibble. If there is only a wig file, it converts it to a bigwig before loading.
#'
#' @param path base path to search for danpos output bigwigs
#' @param name name of sample to import
#' @param chromosome_sizes_file path to chromosome.sizes file.
#' This is a tab separated file with chromosome name and chromosome length in each line
#' @param as_gr `TRUE` (default) function returns a GRanges-Object,
#'   else returns a tibble.
#' @return returns a GRanges object (default) or tibble.
load_danpos_bigwig <- function(path,
                               name,
                               chromosome_sizes_file = NULL,
                               as_gr = TRUE) {
  if (!file.exists(paste0(path, "/pooled/", name, ".Fnor.smooth.bw"))) {
    if (file.exists(paste0(path, "/pooled/", name, ".Fnor.smooth.wig"))) {
      convert_wig_to_bigwig(
        wig_file = paste0(path, "/pooled/", name, ".Fnor.smooth.wig"),
        chr_sizes_file = chromosome_sizes_file
      )
    } else {
      stop("There is no pooled wig or bigwig output file in the provided path. Please make sure to provide a path to DANPOS output!")
    }
  }

  coverage_tibble <- read_bigwig(paste0(path, "/pooled/", name, ".Fnor.smooth.bw")) %>%
    as_tibble() %>%
    select(seqnames, start, end, score)


  # adjust for different wiggle step sizes
  maxSize <- coverage_tibble %>%
    group_by(seqnames) %>%
    summarise(max = max(end))
  starts <- maxSize$max %>%
    map(function(x) {
      1:x
    }) %>%
    unlist()
  coverage_tibble <- GRanges(
    seqnames = rep(coverage_tibble$seqnames, times = width(coverage_tibble %>%
      as_granges())),
    ranges = IRanges(start = starts, width = 1),
    score = rep(coverage_tibble$score, times = width(coverage_tibble %>%
      as_granges()))
  ) %>%
    as_tibble()
  ########

  if (as_gr == TRUE) {
    coverage_gr <- coverage_tibble %>% as_granges()
  } else {
    coverage_tibble
  }
}

#' Convert wig to bigwig
#'
#' This functions converts a wig file to a bigwig file
#'
#' @param wig_file Path to wig-File.
#' @param chr_sizes_file path to tab delimited file indicating chromosome sizes.
#' @return Writes a bigwig file named as the input file into the directory of the input.
convert_wig_to_bigwig <- function(wig_file, chr_sizes_file) {
  if (is.null(chr_sizes_file) | is.na(chr_sizes_file)) {
    stop("You have to provide a chromosome.sizes file for conversion of wig to bigwigs!")
  }
  chr_sizes <- read_delim(chr_sizes_file, delim = "\t", col_names = FALSE)
  chr_sizes <- Seqinfo(seqnames = chr_sizes$X1, seqlengths = chr_sizes$X2) %>% sortSeqlevels()
  new_name <- wig_file %>%
    path_ext_remove()

  wig_gr <- read_wig(wig_file)
  wig_gr <- sortSeqlevels(wig_gr)
  seqinfo(wig_gr) <- chr_sizes
  write_bigwig(wig_gr, file = paste0(new_name, ".bw"))
}

#' Load Danpos Positions
#'
#' This function loads a danpos nucleosome positions file ("XYZ.positions.integrative.xls") and
#' converts it to a GRanges(default) or tibble.
#'
#' @param path base path to search for danpos output bigwigs
#' @param name name of sample to import
#' @param contrast name of contrast/control sample
#' @param as_gr `TRUE` (default) function returns a GRanges-Object,
#'   else returns a tibble.
#' @return returns a GRanges object (default) or tibble.
load_danpos_positions <- function(path,
                                  name,
                                  contrast,
                                  as_gr = TRUE) {
  positions_path <- paste0(path, "/", name, "-", contrast, ".positions.integrative.xls")

  # Load positions as granges
  # Add +1 to all positions for true positions.
  # Not sure if its only required because of 0-based and 1-based genome issue or because of original wig step-size.
  #
  # POTENTIAL ISSUE:
  # -[] Check if positions are dependent on wig-step-size!
  positions_tibble <- read_delim(positions_path, delim = "\t") %>%
    mutate(
      start = start + 1,
      end = end + 1,
      diff_smt_loca = diff_smt_loca + 1,
      center = center + 1,
      control_smt_loca = as.numeric(control_smt_loca) + 1,
      treat_smt_loca = treat_smt_loca + 1,
      ID = paste0("Nuc", 1:length(start))
    ) %>%
    rename(chr = "seqnames")

  if (as_gr == TRUE) {
    positions_gr <- positions_tibble %>% as_granges()
  } else {
    positions_tibble
  }
}

#' calculate variance of nucleosome
#'
#' This function calculates the variance of a supposed nucleosome defined by its arguments.
#' It is based on the implementation of calculating fuzziness in danpos (uses the standard deviation of read positions in each peak as an estimate of nucleosome fuzziness),
#' but ultimately origins in Jiang and Pugh (2009).
#'
#' @details
#' `start_pos`, `end_pos` and `dyad_position` define a nucleosome. The `scores` dataframe is then used to calculate the variance in this boundary.
#' Therefore it is IMPORTANT to make sure the positions and scores match in chromosome!
#'
#' @param start_pos start position for calculating variance.
#' @param end_pos end position for calculating variance.
#' @param dyad_position dyad position of nucleosome.
#' @param scores data.frame or tibble with scores for each bp of the chromosome.
#'   Score columns must be named "score_treat" and/or "score_control".
#' @param type either "treat" or "control".
#' @return Returns a data.frame with variance and count value.
#' @references
#' Jiang, C., & Pugh, B. F. (2009). A compiled and systematic reference map of nucleosome positions across the Saccharomyces cerevisiae genome. Genome biology, 10(10), R109. https://doi.org/10.1186/gb-2009-10-10-r109
#' Chen, K., Xi, Y., Pan, X., Li, Z., Kaestner, K., Tyler, J., Dent, S., He, X., & Li, W. (2013). DANPOS: dynamic analysis of nucleosome position and occupancy by sequencing. Genome research, 23(2), 341–351. https://doi.org/10.1101/gr.142067.112

calculate_variance <- function(start_pos,
                               end_pos,
                               dyad_position,
                               scores,
                               type) {
  start <- start_pos - dyad_position
  end <- end_pos - dyad_position
  peak <- start:end
  if (type == "control") {
    dev_squared_base <- sum(peak^2)
    n <- length(peak)
    dev_squared_sum <- sum(peak^2 * scores$score_control[dyad_position + peak])
    counts <- sum(scores$score_control[dyad_position + peak])
  } else if (type == "treat") {
    dev_squared_base <- sum(peak^2)
    n <- length(peak)
    dev_squared_sum <- sum(peak^2 * scores$score_treat[dyad_position + peak])
    counts <- sum(scores$score_treat[dyad_position + peak])
  } else {
    stop("Please provide 'type' as either 'treat' or 'control'.")
  }

  counts <- counts + n
  dev_squared_sum <- dev_squared_sum + dev_squared_base
  variance <- sqrt(dev_squared_sum / counts)
  setNames(c(variance, counts), c("variance", "counts")) %>%
    as.list() %>%
    as.data.frame()
}

#' calculate point_diff_fdr values
#'
#' calculate FDR values for difference between nucleosomes.
#'
#' @param treat_scores coverage scores for treatment
#' @param control_scores coverage scores for control
#' @param positions table with nucleosome positions
#' @param sampling_n number of pvalues which will be simulated additionally to the ones for nucleosomes.
#' @param random_seed random seed for sampling. default: `NULL`
adjust_point_diff_fdr <- function(treat_scores, control_scores, positions, sampling_n, random_seed = NULL) {
  if (sampling_n > nrow(treat_scores)) {
    warning("n higher than possible scores. Using all scores instead")
    sampling_n <- nrow(treat_scores)
  }


  ####### DEBUG
  # treat_scores = treat_tib
  # control_scores = control_tib
  # positions = positions_gr
  # sampling_n = sampling_n
  ######################################


  joined_tib <- full_join(
    x = control_scores,
    y = treat_scores,
    by = c("seqnames", "start", "end"),
    suffix = c("_control", "_treat")
  )

  if (!is.null(random_seed)) {
    set.seed(random_seed)
  }

  sampled_tib <- joined_tib %>%
    select(score_control, score_treat) %>%
    sample_n(size = sampling_n, replace = FALSE) %>%
    mutate(ID = "simulated")

  positions_diff_scores <- positions %>%
    as_tibble() %>%
    select(control_point_val, treat_point_val, ID) %>%
    rename(control_point_val = "score_control", treat_point_val = "score_treat")

  pval_tib <- bind_rows(positions_diff_scores, sampled_tib) %>%
    mutate(
      score_control = score_control,
      score_treat = score_treat
    ) %>%
    rowwise() %>%
    mutate(
      pval_lambda_control = (ppois(score_treat,
        score_control,
        lower.tail = FALSE,
        log.p = FALSE
      )),
      pval_lambda_treat = (ppois(score_control,
        score_treat,
        lower.tail = FALSE,
        log.p = FALSE
      ))
    )



  # simDat <- padj_tib %>% dplyr::filter(ID != "simulated")

  # ggplot(simDat[1:50000,]) +
  #  geom_point(aes(x = score_treat, y = score_control), alpha = 0.05)

  # ppois(q = 0, lambda = 601, lower.tail = FALSE)

  # dpois(x = , lambda = 60.1)

  # hist(simDat$pval_lambda_treat )

  # sum(simDat$padj_lambda_treat < 0.05 | simDat$padj_lambda_control < 0.05)
  # sum(simDat$pval_lambda_treat < 0.05 | simDat$pval_lambda_control < 0.05)

  padj_tib <- pval_tib %>%
    ungroup() %>%
    mutate(
      padj_lambda_treat = p.adjust(pval_lambda_treat, method = "BH"),
      padj_lambda_control = p.adjust(pval_lambda_control, method = "BH")
    )

  point_diff_fdr_adjusted <- padj_tib %>%
    rowwise() %>%
    mutate(point_diff_fdr_adj = min(padj_lambda_treat, padj_lambda_control)) %>%
    select(ID, point_diff_fdr_adj) %>%
    filter(ID != "simulated")
}


#' calculate smt_diff_fdr values
#'
#' calculate FDR values for difference between nucleosome summits.
#'
#' @param treat_scores coverage scores for treatment
#' @param control_scores coverage scores for control
#' @param positions table with nucleosome positions
#' @param sampling_n number of pvalues which will be simulated additionally to the ones for nucleosomes.
#' @param sampling_type can be one of "standard"(default), "random" or "near".
#'   "same" samples a random position and calculates pvalue for score between treatment and control on this position. Used in original DANPOS.
#'   "random" samples two random positions and calculates pvalue for treatment score on one position and control score on the other position.
#'   "near" samples a random positions and a second position near (+-150bp) that position. It calculates pvalue for treatment score on one position and control score on the other position. (default)
#' @param sampling_near_range range in which the second dyad will be simulated around the first dyad when simulating positions with the `sampling_type` = "near". default: 150
#' @param random_seed random seed for sampling. default: `NULL`
adjust_smt_diff_fdr <- function(treat_scores,
                                control_scores,
                                positions,
                                sampling_n,
                                sampling_type = "near",
                                sampling_near_range = 150,
                                random_seed = NULL) {
  # Check Input Arguments
  if (sampling_type != "same" & sampling_type != "near" & sampling_type != "random") {
    stop("sampling_type has to be 'standard', 'near' or 'random'!")
  }

  if (sampling_n > nrow(treat_scores)) {
    warning("Sampling n higher than possible scores. Using maximum amount of scores instead.")
    sampling_n <- nrow(treat_scores)
  }

  joined_tib <- full_join(
    x = control_scores,
    y = treat_scores,
    by = c("seqnames", "start", "end"),
    suffix = c("_control", "_treat")
  )

  if (!is.null(random_seed)) {
    set.seed(random_seed)
  }

  if (sampling_type == "same") {
    sampled_tib <- joined_tib %>%
      select(score_control, score_treat) %>%
      sample_n(size = sampling_n, replace = FALSE) %>%
      mutate(ID = "simulated")
  } else if (sampling_type == "random") {
    sampled_tib <- bind_cols(
      joined_tib$score_control %>%
        sample(size = sampling_n, replace = TRUE),
      joined_tib$score_treat %>%
        sample(size = sampling_n, replace = TRUE)
    ) %>%
      rename(...1 = "score_control", ...2 = "score_treat") %>%
      mutate(ID = "simulated")
  } else if (sampling_type == "near") {
    control_score_positions <- sample(1:nrow(joined_tib), size = sampling_n, replace = TRUE)
    treat_score_positions <- control_score_positions + sample(-sampling_near_range:sampling_near_range, size = sampling_n, replace = TRUE)

    # Make sure none of the positions is outside of the possible range
    while ((!all(treat_score_positions > 0)) || (!all(treat_score_positions <= nrow(joined_tib)))) {
      change_list <- c(
        which(treat_score_positions <= 0),
        which(treat_score_positions > nrow(joined_tib))
      )
      treat_score_positions[change_list] <- control_score_positions[change_list] + sample(-sampling_near_range:sampling_near_range, size = length(change_list), replace = TRUE)
    }

    # Make sure none of the scores is on another chromosome and the position stays in the possible range
    while ((!all(joined_tib[control_score_positions, ]$seqnames == joined_tib[treat_score_positions, ]$seqnames))) {
      change_list <- which(joined_tib[control_score_positions, ]$seqnames != joined_tib[treat_score_positions, ]$seqnames)
      treat_score_positions[change_list] <- control_score_positions[change_list] + sample(-sampling_near_range:sampling_near_range, size = length(change_list), replace = TRUE)

      while ((!all(treat_score_positions > 0)) || (!all(treat_score_positions <= nrow(joined_tib)))) {
        change_list <- c(
          which(treat_score_positions <= 0),
          which(treat_score_positions > nrow(joined_tib))
        )
        treat_score_positions[change_list] <- control_score_positions[change_list] + sample(-sampling_near_range:sampling_near_range, size = length(change_list), replace = TRUE)
      }
    }

    sampled_tib <- bind_cols(
      joined_tib$score_control[control_score_positions],
      joined_tib$score_treat[treat_score_positions]
    ) %>%
      rename(...1 = "score_control", ...2 = "score_treat") %>%
      mutate(ID = "simulated")
  }

  positions_smt_scores <- positions %>%
    as_tibble() %>%
    select(control_smt_val, treat_smt_val, ID) %>%
    rename(control_smt_val = "score_control", treat_smt_val = "score_treat")

  pval_tib <- bind_rows(positions_smt_scores, sampled_tib) %>%
    rowwise() %>%
    mutate(
      pval_lambda_control = (ppois(score_treat,
        score_control,
        lower.tail = FALSE,
        log.p = FALSE
      )),
      pval_lambda_treat = (ppois(score_control,
        score_treat,
        lower.tail = FALSE,
        log.p = FALSE
      ))
    )

  padj_tib <- pval_tib %>%
    ungroup() %>%
    mutate(
      padj_lambda_treat = p.adjust(pval_lambda_treat, method = "BH"),
      padj_lambda_control = p.adjust(pval_lambda_control, method = "BH")
    )

  point_smt_fdr_adjusted <- padj_tib %>%
    rowwise() %>%
    mutate(smt_diff_fdr_adj = min(padj_lambda_treat, padj_lambda_control)) %>%
    select(ID, smt_diff_fdr_adj) %>%
    filter(ID != "simulated")
}

#' calculate fuzziness_diff_fdr values
#'
#' calculate FDR values for difference between nucleosome fuzziness.
#'
#' @param treat_scores coverage scores for treatment
#' @param control_scores coverage scores for control
#' @param positions table with nucleosome positions
#' @param sampling_n number of pvalues which will be simulated additionally to the ones for nucleosomes.
#' @param sampling_type can be one of "standard"(default), "random" or "near".
#'   "same" samples a random position and calculates pvalue for fuzziness between treatment and control nucleosomes with dyad on this position. Used in original DANPOS.
#'   "random" samples two random positions and calculates pvalue between treatment fuzziness of a nucleosome with one position as dyad and control fuzziness of a nucleosome with dyad on the other position.
#'   "near" samples a random positions and a second position near that position (maximum distance can be supplied via `sampling_near_distance`, default: +-150bp). It calculates pvalue between treatment fuzziness with one position as dyad and control fuzzinesswith dyad on the other position. (default)
#' @param sampling_near_range range in which the second dyad will be simulated around the first dyad when simulating positions with the `sampling_type` = "near". default: 150
#' @param random_seed random seed for sampling. default: `NULL`
adjust_fuzziness_diff_fdr <- function(treat_scores,
                                      control_scores,
                                      positions,
                                      sampling_n,
                                      sampling_type = "near",
                                      sampling_near_range = 150,
                                      random_seed = NULL) {
  if (sampling_type != "same" & sampling_type != "near" & sampling_type != "random") {
    stop("sampling_type has to be 'standard', 'near' or 'random'!")
  }

  if (sampling_n > nrow(treat_scores)) {
    warning("Sampling n higher than possible scores. Using maximum amount of scores instead.")
    sampling_n <- nrow(treat_scores)
  }

  # join tables with coverage data
  joined_tib <- full_join(
    x = control_scores,
    y = treat_scores,
    by = c("seqnames", "start", "end"),
    suffix = c("_control", "_treat")
  )

  # Get maximum sizes of all chromosomes for trimming later
  seqinfo_pre <- joined_tib %>%
    group_by(seqnames) %>%
    summarise(n = n())
  sequence_sizes <- Seqinfo(
    seqnames = as.character(seqinfo_pre$seqnames),
    seqlengths = seqinfo_pre$n
  ) %>% sortSeqlevels()

  if (!is.null(random_seed)) {
    set.seed(random_seed)
  }

  if (sampling_type == "same") {
    sampled_tib <- joined_tib %>%
      mutate(
        dyadpos_control = start,
        dyadpos_treat = start,
        seqnames_control = seqnames,
        seqnames_treat = seqnames
      ) %>%
      select(seqnames_control, seqnames_treat, dyadpos_control, dyadpos_treat) %>%
      dplyr::slice_sample(n = sampling_n, replace = FALSE) %>%
      mutate(ID = "simulated")
  } else if (sampling_type == "random") {
    sampled_tib <- bind_cols(
      joined_tib %>%
        select(seqnames, start) %>%
        dplyr::slice_sample(n = sampling_n, replace = TRUE),
      joined_tib %>%
        select(seqnames, start) %>%
        dplyr::slice_sample(n = sampling_n, replace = TRUE)
    ) %>%
      rename(
        seqnames...1 = "seqnames_control",
        start...2 = "dyadpos_control",
        seqnames...3 = "seqnames_treat",
        start...4 = "dyadpos_treat"
      ) %>%
      mutate(ID = "simulated")
  } else if (sampling_type == "near") {
    control_dyad_positions <- sample(1:nrow(joined_tib), size = sampling_n, replace = TRUE)
    treat_dyad_positions <- control_dyad_positions + sample(-sampling_near_range:sampling_near_range, size = sampling_n, replace = TRUE)

    # Make sure none of the positions is outside of the possible range
    while ((!all(treat_dyad_positions > 0)) || (!all(treat_dyad_positions <= nrow(joined_tib)))) {
      change_list <- c(
        which(treat_dyad_positions <= 0),
        which(treat_dyad_positions > nrow(joined_tib))
      )
      treat_dyad_positions[change_list] <- control_dyad_positions[change_list] + sample(-sampling_near_range:sampling_near_range, size = length(change_list), replace = TRUE)
    }

    # Make sure none of the dyads is on another chromosome and the position stays in the possible range
    while ((!all(joined_tib[control_dyad_positions, ]$seqnames == joined_tib[treat_dyad_positions, ]$seqnames))) {
      change_list <- which(joined_tib[control_dyad_positions, ]$seqnames != joined_tib[treat_dyad_positions, ]$seqnames)
      treat_dyad_positions[change_list] <- control_dyad_positions[change_list] + sample(-sampling_near_range:sampling_near_range, size = length(change_list), replace = TRUE)

      while ((!all(treat_dyad_positions > 0)) || (!all(treat_dyad_positions <= nrow(joined_tib)))) {
        change_list <- c(
          which(treat_dyad_positions <= 0),
          which(treat_dyad_positions > nrow(joined_tib))
        )
        treat_dyad_positions[change_list] <- control_dyad_positions[change_list] + sample(-sampling_near_range:sampling_near_range, size = length(change_list), replace = TRUE)
      }
    }

    sampled_tib <- bind_cols(
      joined_tib %>%
        select(seqnames, start) %>%
        dplyr::slice(control_dyad_positions),
      joined_tib %>%
        select(seqnames, start) %>%
        dplyr::slice(treat_dyad_positions)
    ) %>%
      rename(
        seqnames...1 = "seqnames_control",
        start...2 = "dyadpos_control",
        seqnames...3 = "seqnames_treat",
        start...4 = "dyadpos_treat"
      ) %>%
      mutate(ID = "simulated")
  }

  # smt_loca are required for fuzziness calculation.
  # If there are any nucleosomes without a smt_loca in either treat or control the corresponding smt_loca of the other will be used.
  positions <- positions %>%
    as_tibble() %>%
    mutate(
      control_smt_loca = ifelse(is.na(control_smt_loca), treat_smt_loca, control_smt_loca),
      treat_smt_loca = ifelse(is.na(treat_smt_loca), control_smt_loca, treat_smt_loca)
    )

  # Conversion of position data to bind simulated data later.
  positions_smt_pos <- positions %>%
    as_tibble() %>%
    select(control_smt_loca, treat_smt_loca, ID, seqnames) %>%
    rename(control_smt_loca = "dyadpos_control", treat_smt_loca = "dyadpos_treat", seqnames = "seqnames_control") %>%
    mutate(seqnames_treat = seqnames_control)

  # bind simulated and position data. Convert to long format for GRanges operations
  bind_rows(sampled_tib, positions_smt_pos) %>%
    mutate(data_couple = paste0("id_", 1:length(ID))) %>%
    pivot_longer(!c(ID, data_couple),
      names_to = c(".value", "type"),
      names_sep = "_",
      values_to = "dyad_pos",
      values_drop_na = TRUE
    ) -> longer

  # convert dyad position to nucleosome and make sure nucleosome ranges are in bounds of data
  longer_gr <- longer %>%
    filter(!is.na(dyadpos)) %>%
    mutate(start = dyadpos - 75, end = dyadpos + 75) %>%
    as_granges() %>%
    sortSeqlevels()

  seqinfo(longer_gr) <- sequence_sizes
  longer_gr <- trim(longer_gr)

  # split data into chromosomes for better matching with nucleosome positioning
  joined_list <- joined_tib %>% base::split(.$seqnames)

  # calculate variance of nucleosome positions.
  variance_longer <- longer_gr %>%
    as_tibble() %>%
    rowwise() %>%
    mutate(test_vec = list(calculate_variance(
      start_pos = start,
      end_pos = end,
      dyad_position = dyadpos,
      scores = joined_list[[seqnames]],
      type = type
    ))) %>%
    unnest(cols = test_vec)

  # Calculate p-values with ftest
  pval_tib <- variance_longer %>%
    pivot_wider(id_cols = c(ID, data_couple), names_from = type, values_from = c(seqnames, dyadpos, variance, counts)) %>%
    rowwise() %>%
    mutate(
      pval_a = pf(variance_treat / variance_control,
        counts_treat,
        counts_control,
        log.p = FALSE
      ),
      pval_b = pf(variance_control / variance_treat,
        counts_control,
        counts_treat,
        log.p = FALSE
      )
    )

  # Correct pvalues
  padj_tib <- pval_tib %>%
    ungroup() %>%
    mutate(
      padj_a = p.adjust(pval_a, method = "BH"),
      padj_b = p.adjust(pval_b, method = "BH")
    )

  # select minimal FDR for nucleosome and return nucleosome ID + FDR
  fuzziness_diff_fdr_adjusted <- padj_tib %>%
    rowwise() %>%
    mutate(fuzziness_diff_fdr_adj = min(padj_a, padj_b)) %>%
    select(ID, fuzziness_diff_fdr_adj) %>%
    filter(ID != "simulated")
}



#' get smoothed maxima
#'
#' finds highest maximum in the vector after loess smoothing. If no maxima is found, NA is returned
#'
#' @param x vector to find maximum in.
#' @param s loess smooth span. default: 0.6
get_smooth_max <- function(x, s = 0.6) {
  if (all(is.na(x))) {
    return(NA)
  }
  fit <- loess(x ~ c(1:length(x)), span = s) %>%
    predict()
  maxima <- which(diff(sign(diff(fit))) == -2)
  if (length(maxima) > 1) {
    maxima <- maxima[which.max(fit[which(diff(sign(diff(fit))) == -2) + 1])]
  } else if (length(maxima) < 1) {
    maxima <- NA
  }
  return(maxima)
}


#' find treat to control dyad distances
#'
#' find treat to control distances for smoothed nucleosome tracks
#'
#' @param treat_scores coverage scores for treatment
#' @param control_scores coverage scores for control
#' @param positions table with nucleosome positions
adjust_treat2control_dis <- function(treat_scores,
                                     control_scores,
                                     positions) {
  treat_cov <- treat_scores %>%
    as_granges() %>%
    coverage(., weight = .$score) %>%
    as.list(.) %>%
    map(~ as.vector(.))
  control_cov <- control_scores %>%
    as_granges() %>%
    coverage(., weight = .$score) %>%
    as.list(.) %>%
    map(~ as.vector(.))

  treat2control_dis_adjusted <- positions %>%
    plyranges::select(ID) %>%
    as_tibble() %>%
    mutate(start = ifelse(start < 0, 1, start)) %>%
    group_by(ID) %>%
    nest() %>%
    mutate(
      nucVec_treat = map(data, ~ treat_cov[[.x$seqnames]] %>% .[.x$start:.x$end]),
      nucVec_control = map(data, ~ control_cov[[.x$seqnames]] %>% .[.x$start:.x$end])
    ) %>%
    mutate(
      treat_smt_loca_smooth = map(nucVec_treat, ~ get_smooth_max(.)),
      control_smt_loca_smooth = map(nucVec_control, ~ get_smooth_max(.))
    ) %>%
    dplyr::select(-c("nucVec_treat", "nucVec_control")) %>%
    unnest(cols = c(data, treat_smt_loca_smooth, control_smt_loca_smooth)) %>%
    mutate(
      treat2control_dis_adj = abs(treat_smt_loca_smooth - control_smt_loca_smooth),
      treat_smt_loca_smooth = start + treat_smt_loca_smooth - 1,
      control_smt_loca_smooth = start + control_smt_loca_smooth - 1
    ) %>%
    dplyr::select(ID, treat_smt_loca_smooth, control_smt_loca_smooth, treat2control_dis_adj)
}


#'  calculate FDR values
#'
#' Function that takes a sample name, a contrast/control name and a path to danpos differential nucleosome calling output
#' and calculates true FDR values for those nucleosome positions. Results are saved in a new csv-file.
#' Additionally a number for simulated pvalues can be provided (n=5000000 takes approximately 40-50min and it is scaling nearly linearly).
#' When pooled danpos results are wig files a chromosome.sizes file must be supplied.
#'
#' @param sample_name name of the sample as provided to danpos
#' @param control_name name of the control as provided to danpos
#' @param outPath output Path
#' @param inPath input Path, must be a DANPOS output directory
#' @param sampling_n number of p-values simulated for FDR estimation.
#'  If not provided sampling_n will be 100x the number of nucleosomes in sample.
#' @param chromosome_sizes_file path to a chromosome.sizes file. required for conversion from wig to bigwig.
calc_FDR <- function(sample_name,
                     control_name,
                     outPath,
                     inPath,
                     sampling_n,
                     chromosome_sizes_file = NULL,
                     sampling_type_occ = "near",
                     sampling_type_fuzz = "near") {
  message("####### ", sample_name, " #######")
  message(Sys.time(), ": load danpos output")

  ###### DEBUG
  # sample_name = "EPI"
  # control_name = "TCT"
  # outPath <- argv$output_dir
  # inPath <- argv$input_path
  # sampling_n <- argv$n_simulations
  # chromosome_sizes_file = argv$chr_sizes
  # sampling_type_occ = "near"
  # sampling_type_fuzz = "near"
  ###############################################

  positions_gr <- load_danpos_positions(inPath,
    name = sample_name,
    contrast = control_name
  )

  if (is.null(sampling_n) | is.na(sampling_n)) {
    sampling_n <- length(positions_gr) * 100
  }
  message("calculating FDR values with ", sampling_n, " simulated p-values.")

  treat_tib <- load_danpos_bigwig(inPath, name = sample_name, as_gr = FALSE, chromosome_sizes_file)
  control_tib <- load_danpos_bigwig(inPath, name = control_name, as_gr = FALSE, chromosome_sizes_file)

  message(Sys.time(), ": adjust point diff FDR")
  point_diff_fdr_adjusted <- adjust_point_diff_fdr(
    treat_scores = treat_tib,
    control_scores = control_tib,
    positions = positions_gr,
    sampling_n = sampling_n
  )

  message(Sys.time(), ": adjust smt diff FDR")
  smt_diff_fdr_adjusted <- adjust_smt_diff_fdr(
    treat_scores = treat_tib,
    control_scores = control_tib,
    positions = positions_gr,
    sampling_n = sampling_n,
    sampling_type = sampling_type_occ
  )

  message(Sys.time(), ": adjust fuzziness diff FDR")
  fuzziness_diff_fdr_adjusted <- adjust_fuzziness_diff_fdr(
    treat_scores = treat_tib,
    control_scores = control_tib,
    positions = positions_gr,
    sampling_n = sampling_n,
    sampling_type = sampling_type_fuzz
  )

  message(Sys.time(), ": adjust treat2control distance")
  treat2control_dis_adjusted <- adjust_treat2control_dis(
    treat_scores = treat_tib,
    control_scores = control_tib,
    positions = positions_gr
  )

  positions_gr_out <- positions_gr %>%
    as_tibble() %>%
    left_join(point_diff_fdr_adjusted, by = "ID") %>%
    left_join(smt_diff_fdr_adjusted, by = "ID") %>%
    left_join(fuzziness_diff_fdr_adjusted, by = "ID") %>%
    left_join(treat2control_dis_adjusted, by = "ID") %>%
    as_granges()

  message(Sys.time(), ": write output")
  positions_gr_out %>%
    as_tibble() %>%
    write_delim(delim = "\t", file = paste0(outPath, "/", sample_name, "-", control_name, ".positions.FDR_corrected.xls"))
}
