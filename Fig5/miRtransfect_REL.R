# DESeq2 run with BBUM correction and calling for miRNA transfection data
# from McGeary & Lin et al. 2019
# Peter Y. Wang 2023
# Bartel Lab, Whitehead Institute/MIT

library(DESeq2)
library(bbum)
library(tidyverse)

# Parameters -----

ALPHA = 0.05

theme_plot = theme(
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  axis.line = element_line(colour = "black"),
  axis.ticks = element_line(colour = "black"),
  axis.text = element_text(color = "black", size = 7),
  axis.title = element_text(color = "black", size = 7),
  legend.background = element_blank(),
  legend.key = element_blank(),
  strip.background = element_blank(),
  strip.text.x = element_text(size = 7, colour = "black"),
  strip.text.y = element_text(size = 7, colour = "black", angle = 0)
)

# Get data -----
HeLa.dir = "hela_counts/"
HEK.dir = "hek_counts/"

filename.info = read.csv("files-info.csv")

df = filename.info %>%
  # Pull in df for each file (each line of info file)
  rowwise() %>%
  do(., {
    expe.info = .
    cells.this = expe.info$cells
    miR.this = sub("-", ".", expe.info$miR)
    repl.this = expe.info$repl
    fname.this = expe.info$filename
    df.this = read.delim(paste0(
      if_else(
        cells.this == "HeLa",
        HeLa.dir,
        HEK.dir
      ),
      fname.this
    ), header = F) %>%
      dplyr::rename(geneID = "V1", ct = "V2") %>%
      mutate(cells = cells.this, miR = miR.this, repl = repl.this)
    df.this
  }) %>%
  filter(!grepl("^__", geneID)) %>%
  mutate(expeID = paste(cells,
                        miR,
                        repl,
                        sep = "_"))

# DESeq2 + BBUM correction -----
df.de = df %>%
  # Do analysis set per cell
  group_by(cells) %>%
  do(., {
    df.cell = .

    cell.here = df.cell$cells %>% unique()
    print(paste0("Now processing data from ", cell.here))

    V.miRs = df.cell$miR %>% unique()
    N.miRs = length(V.miRs)
    V.treats = df.cell$expeID %>% unique()
    N.treats = length(V.treats)

    # Count matrix per cell
    cts = df.cell %>%
      dplyr::select(geneID, ct, expeID) %>%
      pivot_wider(names_from = "expeID", values_from = "ct") %>%
      column_to_rownames("geneID") %>%
      as.matrix()

    cts.notenoughreads = df.cell %>%
      group_by(geneID) %>%
      summarize_at("ct", min) %>%
      filter(ct < 5) %>%  # at least 5 in all
      ungroup() %>%
      pull(geneID)

    # Analysis here
    df.de.cell = V.miRs %>%
      # Do this per miRNA
      lapply(function(miR.here){

        print(paste(cell.here, miR.here, "in progress..."))

        # Setting up condition table
        miR.pos = grepl(paste0("_", miR.here, "_"), V.treats)

        this.treat.v = rep("ref", N.treats)
        this.treat.v[miR.pos] = miR.here
        this.treat.v = factor(this.treat.v, levels = c("ref", miR.here))

        conds = data.frame(
          condition = this.treat.v,
          cond = factor(V.treats)) %>%
          column_to_rownames(var = "cond")

        # print(conds)

        # Ready. DESeq2 time
        dds = DESeqDataSetFromMatrix(countData = cts,
                                     colData = conds,
                                     design = ~ condition)
        dds = DESeq(dds,
                    test = "Wald",
                    minReplicatesForReplace = Inf,
                    betaPrior = F,
                    fitType = "parametric",
                    quiet = T)
        res = results(dds,
                      independentFiltering = F,
                      alpha = ALPHA,
                      pAdjustMethod = "BH"
        )

        # Set up for BBUM
        df.res = as.data.frame(res) %>%
          rownames_to_column("geneID") %>%
          mutate(miR = miR.here) %>%
          mutate(FCdown = log2FoldChange < 0)

        # BBUM time!
        df.bbum = df.res %>%
          BBUM_DEcorr(classCol = "FCdown",
                      geneName = "geneID",
                      pBBUM.alpha = ALPHA,
                      excluded = cts.notenoughreads
          )

        df.bbum

      }) %>%
      bind_rows()

    df.de.cell %>%
      mutate(cells = cell.here)
  })


# Output polishing-----
point.classes = c(
  "FALSE.none",
  "FALSE.outlier",
  "TRUE.outlier",
  "TRUE.none",
  "FALSE.hit",
  "TRUE.hit"
)
point.classes.n = c(
  "TRUE.hit",
  "FALSE.hit",
  "TRUE.none",
  "FALSE.outlier",
  "total"
)
df.processed = df.de %>%
  mutate(expeID = paste(cells,
                        sub("\\.", "\uad", miR),
                        sep = "\n"),
         padj.safe = if_else(is.na(padj), 1.0, padj),
         DESeq.hits = padj.safe < ALPHA,
         point.class = interaction(DESeq.hits, BBUM.fct))
expeID.order = df.processed %>%
  transmute(expeID = expeID,
            FCup.BH.i = as.numeric(!FCdown & DESeq.hits),
            BBUM.hits.i = as.numeric(BBUM.hits)
            ) %>%
  group_by(expeID) %>%
  summarise_at(c("FCup.BH.i", "BBUM.hits.i"), sum) %>%
  arrange(-FCup.BH.i, -BBUM.hits.i) %>%
  pull(expeID)
df.processed = df.processed %>%
  mutate(expeID = factor(expeID, levels = expeID.order)) %>%
  mutate(point.class = factor(point.class, levels = point.classes)) %>%
  arrange(expeID, point.class, -pvalue)

# Count table for indicating gene category counts
df.n = df.processed %>%
  ungroup() %>%
  filter(!excluded) %>%
  count(expeID, point.class) %>%
  pivot_wider(names_from = point.class, values_from = n, values_fill = 0) %>%
  mutate(total = FALSE.none+TRUE.none+FALSE.outlier+TRUE.outlier+FALSE.hit+TRUE.hit,
         FALSE.outlier = TRUE.outlier + FALSE.outlier) %>%
  select(-TRUE.outlier, -FALSE.none) %>%
  pivot_longer(cols = c(-expeID), names_to = "point.class", values_to = "counts") %>%
  mutate(counts = paste0(if_else(point.class == "total", "/", ""), counts)) %>%
  mutate(point.class = factor(point.class, levels = point.classes.n),
         class.i = as.numeric(point.class)) %>%
  arrange(expeID, point.class)

# Number of omitted points in Fig. 5B
df.omitted = df.processed %>%
  ungroup() %>%
  filter(!excluded) %>%
  filter(pvalue < 10^-20) %>%
  group_by(expeID) %>%
  count()

# PLOTS -----

## MA -----
df.processed %>%
  filter(!excluded) %>%
  ggplot(aes(x = baseMean, y = log2FoldChange,

             color = point.class

             )) +

  geom_hline(yintercept = 0, color = "black", alpha = 0.5,
             size = 0.5, linetype = "dashed") +
  geom_point(alpha = 0.5, size = 0.5, shape = 16) +

  geom_text(data = df.n,
            size = 1.7,
            hjust = 1, vjust = 1,
            aes(

              label = counts,
              x = 5 * 10^5,
              y = -0.25 - (class.i-1) * 1
                ),
            show.legend = F
            ) +

  scale_color_manual(
    breaks = c(point.classes, "total"),
    values = c(
      "gray80",
      "magenta2",
      "magenta2",
      "darkblue",
      "goldenrod4",
      "red3",
      "black"
    ),
    labels = c(
      "Other points",
      "Outliers",
      "Outliers",
      "p_adj < 0.05, p_BBUM > 0.05",
      "p_adj > 0.05, p_BBUM < 0.05",
      "p_adj < 0.05, p_BBUM < 0.05",
      "Total"
    )
    ) +

  scale_x_continuous(trans = "log10", breaks = 10^seq(0,100,2)) +
  scale_y_continuous(breaks = seq(-100,100,2)) +
  coord_cartesian(ylim = c(-5, 5),
                  xlim = c(1, 5 * 10^5)) +

  facet_wrap("expeID", scales = "free", ncol = 6) +

  labs(y = "Fold change (log2)", x = "Mean normalized expression",
       color = "Gene category") +
  theme_plot

## Volcano -----
df.processed %>%
  filter(!excluded) %>%
  ggplot(aes(x = log2FoldChange, y = -log10(pvalue),

             color = point.class

  )) +

  geom_vline(xintercept = 0, color = "black", alpha = 0.5,
             size = 0.5, linetype = "dashed") +
  geom_point(alpha = 0.5, size = 0.5, shape = 16) +

  scale_color_manual(
    breaks = c(point.classes, "total"),
    values = c(
      "gray80",
      "magenta2",
      "magenta2",
      "darkblue",
      "goldenrod4",
      "red3",
      "black"
    ),
    labels = c(
      "Other points",
      "Outliers",
      "Outliers",
      "p_adj < 0.05, p_BBUM > 0.05",
      "p_adj > 0.05, p_BBUM < 0.05",
      "p_adj < 0.05, p_BBUM < 0.05",
      "Total"
    )
  ) +

  scale_x_continuous(breaks = seq(-100,100,2)) +
  coord_cartesian(xlim = c(-5, 5),
                  ylim = c(1, 20)) +

  facet_wrap("expeID", scales = "free", ncol = 6) +

  labs(x = "Fold change (log2)", y = "-log10(p)",
       color = "Gene category") +
  theme_plot
