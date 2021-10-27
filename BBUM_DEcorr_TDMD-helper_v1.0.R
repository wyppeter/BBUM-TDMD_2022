# Defines wrapper for BBUM correction of sRNA seq data of ZSWIM8 KO/KD.
# If running this as part of a larger single R code file, please insert this
#   code at the head.
# Peter Y. Wang 2021
# Bartel Lab, Whitehead Institute/MIT

# Import tidyverse libraries ----
library(tibble)
library(tidyr)
library(dplyr)
library(purrr)
library(ggplot2)
# ## Install (if necessary) and import bbum ----
if(!require("bbum")) {
  if(!require("devtools")) { install.packages("devtools") }
  devtools::install_github("wyppeter/bbum")
}
library(bbum)

# BBUM_DEcorr_TDMD function ----
#
# Wrapper for the main BBUM_DEcorr function for correction and multiple testing
#   adjustment by BBUM modeling; specifically for TDMD analysis
#
# See documentation of BBUM_DEcorr for more details: ?BBUM_DEcorr
# Argument names were largely preserved from older version of BBUM_DEcorr_TDMD
#   for back-compatibility.
# Args share mostly the same names as in BBUM_DEcorr():
# - deseq_results --> df.deseq
# - miRNA.col     --> classCol
# - fails.cutoff  --> excluded
# - other args carry the exact same names as in BBUM_DEcorr()
#
# For guide/passenger strand identification, please name star strands in any of
#   the following ways:
# - let-7a <--> let-7a*
# - let-7a <--> let-7a_star
# - let-7a-5p <--> let-7a-3p
#
# Outputs a data.frame of the results table
#   with added columns from BBUM correction and significance calling.
# Cols are added as described in BBUM_DEcorr().
# Additional cols specific to TDMD analysis are:
# - miRNA: extracted name of miRNAs (str)
# - miRNA.pre: miRNA name without 3p/5p arms/star annotations (for calling
#     guide-passsenger strand pairs) (str)
# - FCdir: direction of fold change (factor, -1 or +1)
# - FCdir.up: whether fold change went up i.e. sample class (bool)
# - tdmd.miR, tdmd.star: TDMD-sensitive miR and its star (bool)
# - tdmd.fct: summary of tdmd.miR and tdmd.star as factor strs (factor)
#
# Example usage:
# ### DESeq2 calls
# dds = DESeqDataSetFromMatrix(countData = cts, ...)
# dds = DESeq(dds)
# res = results(dds, tidy=T) %>%
#   as.data.frame() %>%
#   rename(miRNA = row)
# ### Then call BBUM_DEcorr_TDMD
# res.BBUMcorr = BBUM_DEcorr_TDMD(
#   deseq_results = res,
#   miRNA.col = "miRNA"
# )
#
BBUM_DEcorr_TDMD = function(
  deseq_results,
  miRNA.col = NULL,
  pBBUM.alpha = 0.05,
  fails.cutoff = c(),
  outliers = c(),
  add_starts = list(), only_start = F,
  limits = list(),
  auto_outliers = T, rthres = 1,
  rtrimmax = 0.05, atrimmax = 10,
  quiet = F
) {

  ## Process input and call bbum::BBUM_DEcorr() ----
  df.bbum.out = deseq_results %>%
    tibble::as_tibble(rownames = NA, .name_repair = "minimal") %>%
    dplyr::mutate(
      FCdir = factor(sign(log2FoldChange)),
      FCdir.up = FCdir == +1
    ) %>%
    bbum::BBUM_DEcorr(
      df.deseq = .,
      classCol = "FCdir.up",
      geneName = miRNA.col,
      pBBUM.alpha = pBBUM.alpha,
      excluded = fails.cutoff,
      outliers = outliers,
      add_starts = add_starts,
      only_start = only_start,
      limits = limits,
      auto_outliers = auto_outliers,
      rthres = rthres,
      rtrimmax = rtrimmax,
      atrimmax = atrimmax,
      quiet = quiet
    ) %>%
    dplyr::mutate(
      miRNA = as.character(geneName),
      miRNA.pre = sub("(\\-(3p|5p))$|(_star)$|(\\*)$", "", miRNA)  # guide/star pair detector
    )

  ## Call miRNAs, prepare output, etc. ----
  df.bbum = df.bbum.out %>%
    dplyr::group_by(miRNA.pre) %>%
    dplyr::mutate(
      tdmd.miR = BBUM.hits) %>% # Call sensitive miRNAs (primary effect)
    dplyr::mutate(
      tdmd.star = any(tdmd.miR) & !tdmd.miR,  # Call stars of sensitive miRNAs
      tdmd.fct = factor(
        dplyr::if_else(tdmd.miR,
                       "guide",
                       dplyr::if_else(tdmd.star,
                                      "star",
                                      dplyr::if_else(outlier,
                                                     "outlier",
                                                     "none")
                       )), levels = c("none","star","guide","outlier"))  # summaritive factor for categories
    ) %>%
    dplyr::ungroup()

  ## Done ----
  return(df.bbum)
}
