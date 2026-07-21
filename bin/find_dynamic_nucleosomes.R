# Find nucleosomes with changes in occupancy, fuzziness and shifts

## Dependancies 

library(tidyverse)
library(plyranges)
library(ggpubr)

## Functions 

get_high_var_nucs <- function(ranked_dat, slope_at_cutoff = 1, plot_cutoff = FALSE) {
  # Fit curve using smooth.spline
  fitx <- smooth.spline(x = ranked_dat$fdr_rank, y = ranked_dat$norm_fdr)
  f <- predict(fitx)
  f1 <- predict(fitx, deriv = 1)
  
  # Find the cutoff position where the slope is > slope_at_cutoff
  slope_above_cutoff <- which(f1$y > slope_at_cutoff)
  
  # positions are only valid after the first 1000 lowest variances. 
  # This is necessary as sometimes the first ones have a quite high predicted slope. 
  valid_cutoff_positions <- slope_above_cutoff[slope_above_cutoff > 1000]
  cutOff <- min(valid_cutoff_positions) / nrow(ranked_dat)
  
  
  if (plot_cutoff) {
    # Plot 
    ranked_dat %>% ggplot(aes(x = fdr_rank, y = norm_fdr)) +
      geom_point() +
      geom_point(data = ranked_dat %>%
                   filter(fdr_rank > cutOff),
                 aes(x = fdr_rank, y = norm_fdr),
                 color = "red") +
      # geom_smooth(aes(x = fdr_rank, y = norm_fdr), method = "loess", span = 0.01, color = "black", size = 0.5) +
      geom_vline(xintercept = cutOff, linetype = "dashed") +
      theme_classic() +
      labs(x = "normalized rank", 
           y = "norm. -log10(FDR)") +
      labs_pubr() +
      border()
  }

  # Extract CutOff FDR
  
  cutOff_FDR <- ranked_dat %>%
    filter(fdr_rank == cutOff) %>% 
    pull(fdr)
  }

# Load Data ------------------------------------------------------------------------

# nucleosome positions - DANPOS results adjusted
# can be downloaded via GEO or found in the pipeline results.
nucs_i24 <- read_delim("../2_pipeline/04_differential_nucleosome_positioning/results/adjustedFDR/iKO_24h-noKO_24h.positions.FDR_corrected.xls",
                       delim = "\t",
                       escape_double = FALSE,
                       trim_ws = TRUE
)

nucs_i0 <- read_delim("~/Desktop/GSE228947_iKO_0h-noKO_0h.positions.FDR_corrected.tsv",
                      delim = "\t",
                      escape_double = FALSE,
                      trim_ws = TRUE
)

out_dir <- "../03_analysis/"
dir.create(out_dir)

# ## OCCUPANCY----------------------------------------------------------------------
# transform fdr and calculate rank and normalize both 
ranked_dat <- nucs_i24 %>%
  mutate(fdr = -log10(smt_diff_fdr_adj)) %>% 
  mutate(fdr = ifelse(fdr == Inf,310, fdr)) %>% #Values are becoming Inf because of floating point precision (max/min values in R are around 10^-308).  
  arrange(fdr) %>% 
  mutate(fdr_rank = 1:dplyr::n(),
         fdr_rank = fdr_rank/max(fdr_rank),
         norm_fdr = fdr/max(fdr))

occCutOff <- get_high_var_nucs(ranked_dat = ranked_dat, slope_at_cutoff = 1, plot_cutoff = FALSE)

# Join Data into one table and clean
outTable <- ranked_dat %>% 
  mutate(highOccFDR = if_else(fdr > occCutOff), TRUE, FALSE) %>% 
  dplyr::select(c("seqnames", "start", "end", "width", "strand", "center",
                  "control_smt_loca", "treat_smt_loca", "control_smt_val",
                  "treat_smt_val", "smt_log2FC", "smt_diff_log10pval", 
                  "smt_diff_FDR", "control_fuzziness_score", "ID", # control_fuzziness_score is included to filter for well positioned nucleosomes
                  "smt_diff_fdr_adj", "fdr_rank", "norm_fdr", "highOccFDR"))

# Save Data table with all relevant datacolumns
outTable %>% 
  write_delim(file = paste0(out_dir, "fdr_cutoff_Occupancy_table_i24.tsv", delim = "\t"))

# Save Bed file with high variance in occupancy nucleosome positions
outTable %>%
  mutate(name = ID, score = smt_diff_fdr_adj) %>% 
  filter(highOccFDR) %>% 
  as_granges() %>% 
  write_bed(file = paste0(out_dir,"/HighOccNuc.bed"))

# ## FUZZINESS----------------------------------------------------------------------
# transform fdr and calculate rank and normalize both 
ranked_dat <- nucs_i24 %>%
  filter(!is.na(fuzziness_diff_fdr_adj)) %>%
  mutate(fdr = -log10(fuzziness_diff_fdr_adj)) %>% 
  mutate(fdr = ifelse(fdr == Inf,310, fdr)) %>% #Values are becoming Inf because of floating point precision (max/min values in R are around 10^-308).  
  arrange(fdr) %>% 
  mutate(fdr_rank = 1:dplyr::n(),
         fdr_rank = fdr_rank/max(fdr_rank),
         norm_fdr = fdr/max(fdr))

fuzzCutOff <- get_high_var_nucs(ranked_dat = ranked_dat, slope_at_cutoff = 1, plot_cutoff = FALSE)

# Join Data into one table and clean
outTable <- ranked_dat %>% 
  mutate(highFuzzFDR = if_else((fdr > fuzzCutOff), TRUE, FALSE)) %>% 
  dplyr::select(c("seqnames", "start", "end", "width", "strand", "center",
                  "control_fuzziness_score", "treat_fuzziness_score",
                  "fuzziness_log2FC", "fuzziness_diff_log10pval",
                  "fuzziness_diff_FDR", "ID", "fuzziness_diff_fdr_adj", 
                  "fdr_rank", "norm_fdr", "highFuzzFDR"))
 
# Save Data table with all relevant datacolumns
outTable %>% 
  write_delim(file = paste0(out_dir, "fdr_cutoff_Fuzziness_table_i24.tsv", delim = "\t"))

# Save Bed file with high variance in fuzziness nucleosome positions
outTable %>% mutate(name = ID,score = fuzziness_diff_fdr_adj) %>% 
  filter(highFuzzFDR) %>% as_granges() %>% 
  write_bed(file = paste0(out_dir,"/HighFuzzNuc.bed"))

# ## SHIFT--------------------------------------------------------------------------
# transform fdr and calculate  rank and normalize both 

ranked_dat <- nucs_i24 %>% 
  filter(treat2control_dis_adj != "NA") %>%
  mutate(fdr = (treat2control_dis_adj)) %>% # using nucleosome distance as substitute for significance value
  mutate(fdr = ifelse(fdr == Inf,310, fdr)) %>% #Values are becoming Inf because of floating point precision (max/min values in R are around 10^-308).  
  arrange(fdr) %>% 
  mutate(fdr_rank = 1:dplyr::n(),
         fdr_rank = fdr_rank/max(fdr_rank),
         norm_fdr = fdr/max(fdr))

shiftCutOff <- get_high_var_nucs(ranked_dat = ranked_dat, slope_at_cutoff = 1, plot_cutoff = FALSE)

# Join Data into one table and clean
outTable <- ranked_dat %>% 
  mutate(highShiftFDR = if_else((fdr > shiftCutOff), TRUE, FALSE)) %>% 
  dplyr::select(c("seqnames", "start", "end", "width", "strand", "center",
                  "control_smt_loca", "treat_smt_loca", "treat2control_dis",
                  "treat_smt_loca_smooth", "control_smt_loca_smooth", 
                  "treat2control_dis_adj", "control_fuzziness_score", "ID", # control_fuzziness_score is included to filter for well positioned nucleosomes
                  "fdr_rank", "norm_fdr", "highShiftFDR"))
  

# Save Data table with all relevant datacolumns
outTable %>% 
  write_delim(file = paste0(out_dir, "fdr_cutoff_shift_table_i24.tsv", delim = "\t"))

# Save Bed file with high variance in fuzziness nucleosome positions
outTable %>%
  mutate(name = ID, score = treat2control_dis_adj) %>% 
  filter(highShiftFDR) %>% 
  as_granges() %>% 
  write_bed(file = paste0(out_dir,"/HighShiftNuc.bed"))


# transfer cutoff to schizont stage ------------------------------------------------

#occcutoff: 
outTable <- nucs_i0 %>%
  mutate(fdr = -log10(smt_diff_fdr_adj)) %>% 
  mutate(fdr = ifelse(fdr == Inf,310, fdr)) %>% #Values are becoming Inf because of floating point precision (max/min values in R are around 10^-308).  
  arrange(fdr) %>% 
  mutate(fdr_rank = 1:dplyr::n(),
         fdr_rank = fdr_rank/max(fdr_rank),
         norm_fdr = fdr/max(fdr)) %>% 
  mutate(highOccFDR = if_else(fdr > occCutOff), TRUE, FALSE) %>% 
  dplyr::select(c("seqnames", "start", "end", "width", "strand", "center",
                  "control_smt_loca", "treat_smt_loca", "control_smt_val",
                  "treat_smt_val", "smt_log2FC", "smt_diff_log10pval", 
                  "smt_diff_FDR", "control_fuzziness_score", "ID", # control_fuzziness_score is included to filter for well positioned nucleosomes
                  "smt_diff_fdr_adj", "fdr_rank", "norm_fdr", "highOccFDR"))
# Save Data table with all relevant datacolumns
outTable %>% 
  write_delim(file = paste0(out_dir, "fdr_cutoff_Occupancy_table_i0.tsv", delim = "\t"))

#fuzzcutoff: 
outTable <- nucs_i0 %>%
  filter(!is.na(fuzziness_diff_fdr_adj)) %>%
  mutate(fdr = -log10(fuzziness_diff_fdr_adj)) %>% 
  mutate(fdr = ifelse(fdr == Inf,310, fdr)) %>% #Values are becoming Inf because of floating point precision (max/min values in R are around 10^-308).  
  arrange(fdr) %>% 
  mutate(fdr_rank = 1:dplyr::n(),
         fdr_rank = fdr_rank/max(fdr_rank),
         norm_fdr = fdr/max(fdr)) %>% 
  mutate(highFuzzFDR = if_else((fdr > fuzzCutOff), TRUE, FALSE)) %>% 
  dplyr::select(c("seqnames", "start", "end", "width", "strand", "center",
                  "control_fuzziness_score", "treat_fuzziness_score",
                  "fuzziness_log2FC", "fuzziness_diff_log10pval",
                  "fuzziness_diff_FDR", "ID", "fuzziness_diff_fdr_adj", 
                  "fdr_rank", "norm_fdr", "highFuzzFDR"))
# Save Data table with all relevant datacolumns
outTable %>% 
  write_delim(file = paste0(out_dir, "fdr_cutoff_Fuzziness_table_i0.tsv", delim = "\t"))

#shiftcutoff: 
outTable <- nucs_i0 %>% 
  filter(treat2control_dis_adj != "NA") %>%
  mutate(fdr = (treat2control_dis_adj)) %>% 
  mutate(fdr = ifelse(fdr == Inf,310, fdr)) %>% #Values are becoming Inf because of floating point precision (max/min values in R are around 10^-308).  
  arrange(fdr) %>% 
  mutate(fdr_rank = 1:dplyr::n(),
         fdr_rank = fdr_rank/max(fdr_rank),
         norm_fdr = fdr/max(fdr)) %>% 
  mutate(highShiftFDR = if_else((fdr > shiftCutOff), TRUE, FALSE)) %>% 
  dplyr::select(c("seqnames", "start", "end", "width", "strand", "center",
                  "control_smt_loca", "treat_smt_loca", "treat2control_dis",
                  "treat_smt_loca_smooth", "control_smt_loca_smooth", 
                  "treat2control_dis_adj", "control_fuzziness_score", "ID", # control_fuzziness_score is included to filter for well positioned nucleosomes
                  "fdr_rank", "norm_fdr", "highShiftFDR"))


# Save Data table with all relevant datacolumns
outTable %>% 
  write_delim(file = paste0(out_dir, "fdr_cutoff_shift_table_i0.tsv", delim = "\t"))

