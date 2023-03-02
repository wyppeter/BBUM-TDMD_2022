# DESeq2 run with BBUM correction and calling for Han et al. data
# + Plots
# Peter Y. Wang 2022
# Bartel Lab, Whitehead Institute/MIT

library(tidyverse)
library(DESeq2)
library(gridExtra)
library(bbum)

# Collate data ----
expconfig = read.csv("./Han-exp-config_TS7.csv")

df = expconfig %>%
  rowwise() %>%
  do(., {
    SRR = .$SRR_ID
    cells = .$cells
    cond = .$condition

    cts = read.csv(paste0(SRR, ".csv")) %>%
      mutate(cells = cells, cond = cond) %>%
      dplyr::select(cells, cond, miR, counts)

    cts
  }) %>%
  group_by(cells) %>%
  pivot_wider(names_from = cond, values_from = counts)
df[is.na(df)] = 0
df.all = df %>%
  mutate(miR = sub("(hsa|mmu|dme)\\-", "", miR)) %>%
  rowwise() %>%
  filter(min(N1,N2,N3,Z1,Z2,Z3) > 5)

coldata = data.frame(
  condition = c(rep("untreated",3),rep("treated",3)),
  cond = c("N1","N2","N3","Z1","Z2","Z3")
) %>%
  mutate(cond =      factor(cond,      levels = c("N1","N2","N3","Z1","Z2","Z3")),
         condition = factor(condition, levels = c("untreated","treated"))) %>%
  column_to_rownames(var = "cond")

# DE analysis ----
df.DE = df.all %>%
  group_by(cells) %>%
  do(., {

    DE.in = .

    print("==============================")
    cellshere = DE.in$cells %>% unique()
    print(cellshere)

    DE.in = DE.in %>%
      select(-cells) %>%
      column_to_rownames("miR")

    # DESeq2 run ----
    dds = DESeqDataSetFromMatrix(countData = DE.in,
                                 colData = coldata,
                                 design = ~ condition)
    dds = DESeq(dds,
                test = "Wald",
                minReplicatesForReplace = Inf,
                betaPrior = F,
                fitType = "parametric",
                quiet = T)

    # DESeq2 Results ----
    res = results(dds,
                  independentFiltering = F,
                  alpha = 0.05,
                  pAdjustMethod = "BH"
    )

    df.res = as.data.frame(res) %>%
      rownames_to_column("miRNA")

    # Carry out BBUM correction and signif calling
    df.res.corr = df.res %>%
      BBUM_DEcorr_TDMD(
        auto_outliers = T,  # toggle to test without outlier detection.
        miRNA.col = "miRNA",
        pBBUM.alpha = 0.05
        )

    # Get WTmean ----
    cts.n = dds %>%
      counts(normalized = T) %>%
      as.data.frame() %>% rownames_to_column("miRNA")
    df.wt = cts.n %>%
      select(miRNA, starts_with("N")) %>%
      pivot_longer(starts_with("N"), names_to = "sample", values_to = "ctN") %>%
      select(-sample) %>%
      group_by(miRNA) %>%
      summarize_at("ctN", mean) %>%
      transmute(miRNA = miRNA, WTmean = ctN)

    # Output ----
    df.out = df.res.corr %>%
      left_join(DE.in %>%
                  as.data.frame() %>%
                  rownames_to_column("miRNA"),
                by = "miRNA") %>%
      left_join(df.wt,
                by = "miRNA")

    df.out
  })

#-------------------------------------------

# Fig. 3C
df.DE %>%
  arrange(tdmd.fct, abs(log2FoldChange)) %>%
  ggplot(aes(x = WTmean, y = log2FoldChange, color = tdmd.fct)) +
  geom_hline(yintercept = 0, color = "black", alpha = 0.5,
             size = 0.5, linetype = "dashed") +
  facet_wrap("cells", scales = "free") +
  geom_point(alpha = 0.6, size = 1.5, shape = 16) +
  scale_color_manual(
    breaks = c("none","passenger","guide","outlier"),
    values = c(
      "gray80",
      "cyan3",
      "red3",
      "magenta2"
    )) +
  scale_x_continuous(trans = "log10", breaks = 10^seq(0,100,2)) +
  scale_y_continuous(breaks = seq(-100,100,2)) +
  coord_cartesian(xlim = c(1, NA), ylim = c(-6,4)) +
  labs(y = "Fold change (log2)", x = "Normalized WT expression",
       title = "MA plot", color = "Gene category") +
  theme_classic(base_size = 12)

df.DE.hits = df.DE %>%
  filter(tdmd.miR)

# bbum::BBUM_plot(df.DE %>% filter(cells == "K562"),
#           option = "pcorr")

df.DE %>%
  write.csv("./Han_DE_results.csv", quote = F, row.names = F)

#---------------------------
