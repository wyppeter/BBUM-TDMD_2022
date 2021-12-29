# DESeq2 run with BBUM correction and calling for Shi et al. data
# Peter Y. Wang 2021
# Bartel Lab, Whitehead Institute/MIT

library(tibble)
library(dplyr)
library(tidyr)
library(forcats)
library(ggplot2)

shortlist = c(
  # "MEF",
  # "iMN",
  # "HeLa"

  "MEF",
  "iMN",
  "S2",
  "K562-KO",
  "K562-KDa",
  "K562-KDb",
  "A549",
  "HeLa",
  "MCF7"

)

# Setting data up
cellsv = c(
  "MEF",
  "iMN",
  "S2",

  "K562-KO",
  "K562-KDa",
  "K562-KDb",

  "A549",
  "HeLa",
  "MCF7"
)
# p_adj thresholds defined in Shi et al. 2020
# FDRp_v = c(1E-7, 1E-4,  1E-8,
#            0.01, 1E-20, 1E-20,
#            0.05,  0.05, 0.05
# )  # original
FDRp_v = c(1E-7, 1E-4, 1E-5,
           0.01, 1E-7, 1E-4,
           0.05,  0.05, 0.05
)  # current run
names(FDRp_v) = cellsv

##########################

collatedata = function(cellshere){
  # Values table from my own rerun of DESeq2 (for raw pvals and other values)
  df.raw = read.csv(
    paste0("./dfout-", cellshere, ".csv"),
    stringsAsFactors = T)

  # Repair names
  # Calculate relevant value fields
  df.raw %>%
    mutate(cells = cellshere)
}
df = lapply(cellsv, collatedata) %>%
  bind_rows() %>%

  filter(cells %in% shortlist) %>%

  mutate(FCdir = factor(FCdir))

df.corr = df %>%
  filter(!is.na(padj), !excluded) %>%

  mutate(tdmd_new = tdmd.miR,
         tdmd_old = padj < FDRp_v[cells]) %>%

  mutate(cells = factor(cells, levels = cellsv)) %>%

  # Identify status
  mutate(status = factor(
    if_else(
      tdmd_old,
      if_else(
        tdmd_new,
        "agreed",
        "removed"
      ),
      if_else(
        tdmd_new,
        "added",
        "none"
      )
    ), levels = c("agreed","removed","added","none"))
  ) %>%

  mutate(tdmd.fct = factor(tdmd.fct, levels = c("none","star","guide","outlier"))) %>%

  arrange(cells,
          tdmd.fct,
          -abs(log2FoldChange),
          miRNA)

# Log output ----
# write.csv(df.corr,
#           file = paste0(
#             "./dfout_compare.csv"
#           ),
#           row.names = F, quote = F)


####################

# # Verify cutoff
# df.corr %>%
#   filter(!FCdir.up, !outlier, !excluded) %>%
#   group_by(cells) %>%
#   summarize_at("padj", min) %>%
#   mutate(padj = signif(padj, 3)) %>%
#   mutate(cutoff.used = FDRp_v[cells]) %>%
#   View()

############PLOT###############

# Fig. 1A
df.corr %>%
  mutate(oldsecon = !tdmd_old & padj < 0.05,
         fig1a.fct = factor(if_else(
           tdmd_old,
           "prim",
           if_else(
             oldsecon,
             "seco",
             "null"
           )
         ), levels = c("prim","seco","null"))) %>%
  arrange(fct_rev(fig1a.fct)) %>%
  ggplot2::ggplot(ggplot2::aes(x = WTmean,
                               y = log2FoldChange,
                               color = fig1a.fct)) +
  ggplot2::geom_hline(yintercept = 0, color = "black", alpha = 0.5,
                      size = 0.5, linetype = "dashed") +
  ggplot2::geom_point(alpha = 0.6, size = 1.5, shape = 16) +
  facet_wrap("cells", scales = "free", ncol = 3) +
  ggplot2::scale_color_manual(
    breaks = c("prim","seco","null"),
    values = c(
      "red3",
      "darkblue",
      "gray80"
    )) +
  ggplot2::scale_x_continuous(trans = "log10", breaks = 10^seq(0,100,2)) +
  ggplot2::scale_y_continuous(breaks = seq(-100,100,2)) +
  ggplot2::coord_cartesian(ylim = c(-5, 5),
                           xlim = c(1, NA)) +
  ggplot2::labs(y = "Fold change (log2)", x = "Normalized WT expression",
                title = "MA plot", color = "Gene category") +
  ggplot2::theme_classic(base_size = 12)

# Fig. 1B
## Histogram ----
binN = 20
binwidth = 1/binN
df.corr_plot_bin = df.corr %>%
  group_by(cells) %>%
  dplyr::mutate(pvalue.binned.raw = ggplot2::cut_interval(pvalue, n = binN),
                pvalue.binned = as.numeric(
                  gsub("(^[\\(\\[])|(,.*)", "", pvalue.binned.raw)
                  # extract number from range text
                )
  ) %>%
  dplyr::group_by(cells, BBUM.th, BBUM.class, pvalue.binned) %>%
  dplyr::tally(.) %>%
  dplyr::group_by(cells, BBUM.th, BBUM.class) %>%
  dplyr::mutate(freq = n/sum(n)/binwidth) %>%
  dplyr::mutate(freq = dplyr::if_else(BBUM.class, freq, freq * (1-BBUM.th)))
df.corr_plot_bin %>%
  ggplot2::ggplot(ggplot2::aes(
    x = pvalue.binned,
    y = freq,
    color = factor(BBUM.class, levels = c(TRUE,FALSE)))) +
  facet_wrap("cells", scales = "free", ncol = 3) +
  ggplot2::geom_hline(yintercept = 0, color = "gray60", alpha = 0.75,
                     size = 0.5) +
  ggplot2::geom_vline(xintercept = c(0,1), color = "gray60", alpha = 0.75,
                     size = 0.5) +
  ggplot2::geom_step(stat = "identity", alpha = 0.75, size = 1) +
  ggplot2::scale_color_manual(breaks = c(FALSE, TRUE),
                             values = c("goldenrod3","turquoise4"),
                             labels = c("Background", "Signal")
  ) +
  ggplot2::labs(x = "p value", y = "Density of probability",
               title = "Histogram of p values", color = "Data set") +
  ggplot2::coord_cartesian(ylim = c(0, 4)) +
  ggplot2::theme_classic(base_size = 12)

# Fig. 2
df.coefs = df.corr %>%
  ungroup() %>%
  select(cells, BBUM.l, BBUM.a, BBUM.th, BBUM.r) %>%
  distinct()
bbum.model.graph = tibble::tibble(
  p = sort(c(10^seq(-300,-3,0.05), seq(1E-3,1,1E-3)))
) %>%
  crossing(tibble(
    cells = cellsv[cellsv %in% shortlist]
  )) %>%
  left_join(df.coefs, by = "cells") %>%
  dplyr::mutate(dbum.model  = dbum(p, BBUM.l, BBUM.a),
                pbum.model  = pbum(p, BBUM.l, BBUM.a),
                dbbum.model = dbbum(p, BBUM.l, BBUM.a, BBUM.th, BBUM.r),
                pbbum.model = pbbum(p, BBUM.l, BBUM.a, BBUM.th, BBUM.r)
  ) %>%
  mutate(cells = factor(cells, levels = cellsv))

df.corr %>%
  ggplot2::ggplot(ggplot2::aes(
    x = pvalue,
    color = factor(BBUM.class, levels = c(TRUE,FALSE)))) +
  facet_wrap("cells", scales = "free", ncol = 3) +
  ggplot2::geom_hline(yintercept = 0, color = "gray60",
                      alpha = 0.75, size = 0.5) +
  ggplot2::geom_vline(xintercept = c(1), color = "gray60",
                      alpha = 0.75, size = 0.5) +
  ggplot2::geom_line(data = tibble::tibble(
    pvalue     = 10^(seq(-10,0,0.01)),
    pvalue.lin = 10^(seq(-10,0,0.01))
  ), ggplot2::aes(y = pvalue.lin), color = "gray60", alpha = 0.75,
  size = 0.7, linetype = "dashed") +
  ggplot2::geom_step(stat = "ecdf", alpha = 0.75, size = 0.5) +
  ggplot2::geom_line(data = bbum.model.graph, inherit.aes = F,
                     ggplot2::aes(x = p, y = pbum.model),
                     color = "goldenrod3",
                     alpha = 0.5, size = 0.7) +
  ggplot2::geom_line(data = bbum.model.graph, inherit.aes = FALSE,
                     ggplot2::aes(x = p, y = pbbum.model),
                     color = "turquoise4",
                     alpha = 0.5, size = 0.7) +
  ggplot2::scale_color_manual(breaks = c(FALSE, TRUE),
                              values = c("goldenrod3","turquoise4"),
                              labels = c("Background", "Signal")
  ) +
  ggplot2::scale_x_continuous(trans = "log10") +
  ggplot2::coord_cartesian(xlim = c(1E-6,1)) +
  ggplot2::labs(x = "p value", y = "Cumulative frequency",
                title = "ECDF of p values, in log", color = "Data set") +
  ggplot2::theme_classic(base_size = 12)

# Fig 3A
pBBUM.alpha = 0.05
df.pdir = df.corr %>%
  dplyr::select(cells, geneName, tdmd.fct, status, pvalue, pBBUM, BBUM.class) %>%
  tidyr::pivot_longer(cols = tidyr::starts_with("p"),
                      names_to = "p", values_to = "val") %>%
  dplyr::mutate(p = factor(p, levels = c("pvalue","pBBUM"))) %>%
  dplyr::mutate(p.dir = -log10(val) * dplyr::if_else(BBUM.class,+1,-1)) %>%
  dplyr::arrange(tdmd.fct) %>%
  mutate(cells = factor(cells, levels = cellsv))
df.pdir %>%
  ggplot2::ggplot(ggplot2::aes(y = p, x = p.dir,
                               color = tdmd.fct,
                               group = geneName
  )) +
  facet_grid(cells~., scales = "free") +
  ggplot2::geom_point(alpha = 0.6, shape = 16) +
  geom_point(data = df.pdir %>%
               filter(status %in% c("added", "removed")),
             color = "black", alpha = 0.6, shape = 1) +
  ggplot2::geom_path(alpha = 0.05, size = 0.25) +
  scale_color_manual(
    breaks = c("none","star","guide","outlier"),
    values = c(
      "gray80",
      "cyan3",
      "red3",
      "gray25"
    )) +
  ggplot2::geom_vline(xintercept = 0, color = "black", size = 0.5,
                      linetype = "dashed") +
  ggplot2::geom_vline(xintercept = c(log10(pBBUM.alpha),
                                     -log10(pBBUM.alpha)),
                      color = "salmon1", alpha = 0.5, size = 0.5,
                      linetype = "dashed") +
  ggplot2::scale_x_continuous(breaks = seq(-100,100,5),
                              labels = 10^-abs(seq(-100,100,5))) +
  ggplot2::coord_cartesian(xlim = c(-7.5, 20)) +
  ggplot2::labs(x = "Value", y = "Statistic",
                title = "Correction of p values",
                color = "Gene category") +
  ggplot2::theme_classic(base_size = 12)

# Fig. 3B
df.corr %>%
  ggplot2::ggplot(ggplot2::aes(x = WTmean,
                               y = log2FoldChange,
                               color = tdmd.fct)) +
  ggplot2::geom_hline(yintercept = 0, color = "black", alpha = 0.5,
                      size = 0.5, linetype = "dashed") +
  ggplot2::geom_point(alpha = 0.6, size = 1.5, shape = 16) +
  facet_wrap("cells", scales = "free", ncol = 3) +
  ggplot2::scale_color_manual(
    breaks = c("none","star","guide","outlier"),
    values = c(
      "gray80",
      "cyan3",
      "red3",
      "gray25"
    )) +
  geom_point(data = df.corr %>%
               filter(status %in% c("added", "removed")),
             color = "black", size = 1.5, stroke = 1, shape = 1, alpha = 0.6) +
  ggplot2::scale_x_continuous(trans = "log10", breaks = 10^seq(0,100,2)) +
  ggplot2::scale_y_continuous(breaks = seq(-100,100,2)) +
  ggplot2::coord_cartesian(ylim = c(-5, 5),
                           xlim = c(1, NA)) +
  ggplot2::labs(y = "Fold change (log2)", x = "Normalized WT expression",
                title = "MA plot", color = "Gene category") +
  ggplot2::theme_classic(base_size = 12)
