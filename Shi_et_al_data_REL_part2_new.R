# DESeq2 run with BBUM correction and calling for Shi et al. data
# Plots
# Peter Y. Wang 2023
# Bartel Lab, Whitehead Institute/MIT

library(tibble)
library(dplyr)
library(tidyr)
library(forcats)
library(ggplot2)
library(ggrepel)
library(ggnewscale)
library(bbum)

shortlist = c(
  "MEF",
  "iMN",
  "HeLa"

  # "MEF",
  # "iMN",
  # "S2",
  # "A549",
  # "HeLa",
  # "MCF7",
  # "K562-KO",
  # "K562-KDa",
  # "K562-KDb"

)

# Setting data up
cellsv = c(
  "MEF",
  "iMN",
  "S2",

  "A549",
  "HeLa",
  "MCF7",

  "K562-KO",
  "K562-KDa",
  "K562-KDb"
)
# p_adj thresholds as defined in Shi et al. 2020
# >> original
# FDRp_v = c(1E-7, 1E-4,  1E-8,
#            0.05,  0.05, 0.05,
#            0.01, 1E-20, 1E-20
# )
# >> current run
FDRp_v = c(1E-7, 1E-4, 1E-5,
           0.05,  0.05, 0.05,
           0.01, 1E-7, 1E-4
)
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

  # filter(cells %in% shortlist) %>%

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

  mutate(tdmd.fct = factor(tdmd.fct, levels = c("none","passenger","guide","outlier"))) %>%

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

############PLOTS###############

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
  ggplot(aes(x = WTmean,
                               y = log2FoldChange,
                               color = fig1a.fct)) +
  geom_hline(yintercept = 0, color = "black", alpha = 0.5,
                      size = 0.5, linetype = "dashed") +
  geom_point(alpha = 0.6, size = 1.5, shape = 16) +
  facet_wrap("cells", scales = "free", ncol = 3) +
  scale_color_manual(
    breaks = c("prim","seco","null"),
    values = c(
      "red3",
      "darkblue",
      "gray80"
    )) +
  scale_x_continuous(trans = "log10", breaks = 10^seq(0,100,2)) +
  scale_y_continuous(breaks = seq(-100,100,2)) +
  coord_cartesian(ylim = c(-5, 5),
                           xlim = c(1, NA)) +
  labs(y = "Fold change (log2)", x = "Normalized WT expression",
                title = "MA plot", color = "Gene category") +
  theme_classic(base_size = 12)

# Fig. 1B
## Histogram ----
binN = 20
binwidth = 1/binN
df.corr_plot_bin = df.corr %>%
  group_by(cells) %>%
  mutate(pvalue.binned.raw = cut_interval(pvalue, n = binN),
                pvalue.binned = as.numeric(
                  gsub("(^[\\(\\[])|(,.*)", "", pvalue.binned.raw)
                  # extract number from range text
                )
  ) %>%
  group_by(cells, BBUM.th, BBUM.class, pvalue.binned) %>%
  tally(.) %>%
  group_by(cells, BBUM.th, BBUM.class) %>%
  mutate(freq = n/sum(n)/binwidth) %>%
  mutate(freq = if_else(BBUM.class, freq, freq * (1-BBUM.th)))
df.corr_plot_bin %>%
  filter(cells %in% shortlist) %>%
  ggplot(aes(
    x = pvalue.binned,
    y = freq,
    color = factor(BBUM.class, levels = c(TRUE,FALSE)))) +
  facet_wrap("cells", scales = "free", ncol = 3) +
  geom_hline(yintercept = 0, color = "gray60", alpha = 0.75,
                     size = 0.5) +
  geom_vline(xintercept = c(0,1), color = "gray60", alpha = 0.75,
                     size = 0.5) +
  geom_step(stat = "identity", alpha = 0.75, size = 1) +
  scale_color_manual(breaks = c(FALSE, TRUE),
                             values = c("goldenrod3","darkslategray3"),
                             labels = c("Background", "Signal")
  ) +
  labs(x = "p value", y = "Density of probability",
               title = "Histogram of p values", color = "Data set") +
  coord_cartesian(ylim = c(0, 4)) +
  scale_x_continuous(breaks = c(0,0.5,1)) +
  theme_classic(base_size = 12)

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
  mutate(dbum.model  = dbum(p, BBUM.l, BBUM.a),
                pbum.model  = pbum(p, BBUM.l, BBUM.a),
                dbbum.model = dbbum(p, BBUM.l, BBUM.a, BBUM.th, BBUM.r),
                pbbum.model = pbbum(p, BBUM.l, BBUM.a, BBUM.th, BBUM.r)
  ) %>%
  mutate(cells = factor(cells, levels = cellsv))

df.corr %>%
  filter(cells %in% shortlist) %>%
  ggplot(aes(
    x = pvalue,
    color = factor(BBUM.class, levels = c(TRUE,FALSE)))) +
  facet_wrap("cells", scales = "free", ncol = 3) +
  geom_hline(yintercept = 0, color = "gray60",
             alpha = 0.75, size = 0.5) +
  geom_vline(xintercept = c(1), color = "gray60",
             alpha = 0.75, size = 0.5) +
  geom_ribbon(data = bbum.model.graph, inherit.aes = F,
              aes(x = p, ymin = p, ymax = pbum.model),
              fill = "palegoldenrod") +
  geom_line(data = tibble::tibble(
    pvalue     = 10^(seq(-10,0,0.01)),
    pvalue.lin = 10^(seq(-10,0,0.01))
  ), aes(y = pvalue.lin), color = "gray60", alpha = 0.75,
  size = 0.7, linetype = "dashed") +
  geom_step(stat = "ecdf", alpha = 0.75, size = 0.5) +
  geom_line(data = bbum.model.graph, inherit.aes = F,
            aes(x = p, y = pbum.model),
            color = "goldenrod3",
            alpha = 0.5, size = 0.7) +
  geom_line(data = bbum.model.graph, inherit.aes = FALSE,
            aes(x = p, y = pbbum.model),
            color = "turquoise4",
            alpha = 0.5, size = 0.7) +
  scale_color_manual(breaks = c(FALSE, TRUE),
                     values = c("goldenrod3","turquoise4"),
                     labels = c("Background", "Signal")
  ) +
  scale_x_continuous(trans = "log10") +
  coord_cartesian(xlim = c(1E-6,1)) +
  labs(x = "p-value", y = "Cumulative frequency",
       title = "ECDF of p values, in log", color = "Data set") +
  theme_classic(base_size = 12)

# Fig. 3A
df.corr %>%
  ggplot(aes(x = WTmean,
             y = log2FoldChange,
             color = tdmd.fct)) +
  geom_hline(yintercept = 0, color = "black", alpha = 0.5,
                      size = 0.5, linetype = "dashed") +
  geom_point(alpha = 0.6, size = 1.5, shape = 16) +
  facet_wrap("cells", scales = "free", ncol = 3) +
  scale_color_manual(
    breaks = c("none","passenger","guide","outlier"),
    values = c(
      "gray80",
      "cyan3",
      "red3",
      "gray25"
    )) +
  geom_point(data = df.corr %>%
               filter(status %in% c("added", "removed")),
             color = "black", size = 1.5, stroke = 1, shape = 1, alpha = 0.6) +
  geom_text_repel(data = df.corr %>%
               filter(miRNA %in% c(
                 "miR-7-5p",
                 "miR-7a-5p",
                 "miR-7b-5p")),
             aes(label = miRNA), size = 2, min.segment.length = 0,
             nudge_x = 8, nudge_y = 1, alpha = 0.25, color = "black") +
  scale_x_continuous(trans = "log10", breaks = 10^seq(0,100,2)) +
  scale_y_continuous(breaks = seq(-100,100,2)) +
  coord_cartesian(ylim = c(-5, 5),
                           xlim = c(1, NA)) +
  labs(y = "Fold change (log2)", x = "Normalized WT expression",
                title = "MA plot", color = "Gene category") +
  theme_classic(base_size = 12)

# Fig 3B
df.strands = df.corr %>%
  group_by(cells, miRNA.pre) %>%
  filter(any(tdmd.miR)) %>%
  filter(n() > 1) %>%
  mutate(changed = any(status %in% c("added","removed"))) %>%
  group_by(cells) %>%
  filter(any(changed)) %>%
  ungroup() %>%
  arrange(changed) %>%
  mutate(strandhere = factor(
    if_else(tdmd.miR, "Guide", "Passenger"),
    levels = c("Passenger","Guide")
  )) %>%
  select(cells, miRNA.pre, miRNA, strandhere, changed, log2FoldChange)

df.strands.pair = full_join(
  df.strands %>% filter(strandhere == "Guide") %>%
    select(-strandhere) %>%
    rename(log2FC.guide = log2FoldChange),
  df.strands %>% filter(strandhere == "Passenger")%>%
    select(-strandhere) %>%
    rename(log2FC.pass = log2FoldChange,
           miRNA.star = miRNA),
  by = c(
    "cells", "miRNA.pre", "changed"
  )
) %>%
  pivot_longer(cols = c(log2FC.guide, log2FC.pass),
               names_prefix = "log2FC.",
               names_to = "strandhere",
               values_to = "log2FoldChange") %>%
  mutate(strandhere = factor(strandhere, levels = c("pass","guide")))

df.strands.pair %>%
  ggplot(aes(x = strandhere,
             y = log2FoldChange,
             group = miRNA.pre)) +
  geom_hline(yintercept = 0, color = "black", alpha = 0.15,
             size = 0.5, linetype = "dashed") +
  geom_path(
    aes(color = changed), alpha = 0.8, size = 0.5) +
  scale_color_manual(
    values = c(
      "gray90",
              "black"
    )) +
  new_scale_color() +
  geom_point(aes(color = strandhere),
             alpha = 0.6, size = 1.5, shape = 16) +
  facet_wrap("cells", scales = "free", ncol = 1) +
  scale_color_manual(
    breaks = c("pass","guide"),
    values = c(
      "cyan3",
             "red3"
    )) +
  geom_point(data = df.strands.pair %>%
               filter(changed),
             color = "black", size = 1.5, stroke = 1, shape = 1, alpha = 0.6) +
  geom_text(data = df.strands.pair %>%
              filter(changed, strandhere == "guide"),
            aes(label = miRNA), size = 2,
            nudge_x = 0.2, hjust = 0, vjust = 0.5, color = "black") +
  scale_y_continuous(breaks = seq(-100,100,1)) +
  coord_cartesian(ylim = c(-0.8, 4), xlim = c(NA, 3)) +
  labs(y = "Fold change (log2)", x = "Strand",
       title = "Strand plot", color = "Gene category") +
  theme_classic(base_size = 12)

# OLD Fig 3D
pBBUM.alpha = 0.05
df.pdir = df.corr %>%
  select(cells, geneName, tdmd.fct, status, pvalue, pBBUM, BBUM.class) %>%
  pivot_longer(cols = starts_with("p"),
                      names_to = "p", values_to = "val") %>%
  mutate(p = factor(p, levels = c("pvalue","pBBUM"))) %>%
  mutate(p.dir = -log10(val) * if_else(BBUM.class,+1,-1)) %>%
  arrange(tdmd.fct) %>%
  mutate(cells = factor(cells, levels = cellsv))
df.pdir %>%
  filter(cells %in% shortlist) %>%
  filter(is.finite(p.dir)) %>%
  ggplot(aes(y = p, x = p.dir,
                               color = tdmd.fct,
                               group = geneName
  )) +
  facet_grid(cells~., scales = "free") +
  geom_path(alpha = 0.15, size = 0.25) +
  geom_point(alpha = 0.6, shape = 16) +
  geom_point(data = df.pdir %>%
               filter(status %in% c("added", "removed")) %>%
               filter(cells %in% c(
                 "MEF",
                 "iMN",
                 "HeLa"
               )),
             color = "black", alpha = 0.6, shape = 1, stroke = 0.8) +
  scale_color_manual(
    breaks = c("none","passenger","guide","outlier"),
    values = c(
      "gray80",
      "cyan3",
      "red3",
      "gray25"
    )) +
  geom_vline(xintercept = 0, color = "black", size = 0.5,
                      linetype = "dashed") +
  geom_vline(xintercept = c(log10(pBBUM.alpha),
                                     -log10(pBBUM.alpha)),
                      color = "salmon1", alpha = 0.5, size = 0.5,
                      linetype = "dashed") +
  scale_x_continuous(breaks = seq(-100,100,5),
                              labels = 10^-abs(seq(-100,100,5))) +
  scale_y_discrete(limits = rev) +
  coord_cartesian(xlim = c(-7.5, 20)) +
  labs(x = "Value", y = "Statistic",
                title = "Correction of p values",
                color = "Gene category") +
  theme_classic(base_size = 12)

# NEW Fig 3D
# pBBUM.alpha = NA_real_  # toggle for different rows
pBBUM.alpha = 0.05

# VOLC.YLIM = 30  # toggle for different rows
VOLC.YLIM = 15

df.corr.3D = df.corr %>%
  filter(cells %in% c(
    "MEF",
    "iMN",
    "HeLa"
  )) %>%
  select(cells, geneName, tdmd.fct, status, log2FoldChange, pvalue, pBBUM, BBUM.class) %>%
  pivot_longer(cols = starts_with("p"),
               names_to = "p", values_to = "val") %>%
  mutate(p = factor(p, levels = c("pvalue","pBBUM"))) %>%
  arrange(tdmd.fct) %>%
  mutate(cells = factor(cells, levels = cellsv)) %>%
  mutate(val = if_else(val <= 10^-VOLC.YLIM, 0, val))
df.corr.3D %>%
  ggplot(aes(x = log2FoldChange,
             y = -log10(val),
             color = tdmd.fct,
             group = geneName)) +
  geom_vline(xintercept = 0, color = "black", alpha = 0.5,
                      size = 0.5) +
  geom_hline(yintercept = 0, color = "black", alpha = 0.5,
                      size = 0.5) +
  geom_hline(yintercept = -log10(pBBUM.alpha),
             color = "salmon1", alpha = 0.5, size = 0.5,
             linetype = "dashed") +
  geom_point(alpha = 0.75, size = 1.5, shape = 16) +
  geom_point(data = df.corr.3D %>%
               filter(status %in% c("added", "removed")),
             color = "black", alpha = 0.6, shape = 1, stroke = 0.8) +

  facet_grid(p~cells, scales = "free") +

  scale_color_manual(
    breaks = c("none","passenger","guide","outlier"),
    values = c(
      "gray80",
      "cyan3",
      "red3",
      "gray25"
    )) +
  scale_y_continuous(breaks = seq(-100,100,5),
                     labels = 10^-abs(seq(-100,100,5))) +
  coord_cartesian(
    ylim = c(0, VOLC.YLIM),
    xlim = c(-2, 3.2)) +
  labs(x = "Fold change (log2)", y = "p",
                title = "Volcano plot", color = "Gene category") +
  theme_classic(base_size = 12)
