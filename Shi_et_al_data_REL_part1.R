# DESeq2 run with BBUM correction and calling for Shi et al. data
# Peter Y. Wang 2022
# Bartel Lab, Whitehead Institute/MIT

library(tibble)
library(dplyr)
library(tidyr)
library(forcats)
library(DESeq2)
library(bbum)

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

dme_outliers = c(
  "dme-miR-137-3p",
  "dme-miR-307a-3p",
  "dme-miR-980-3p"  # CV outliers
)

sRNAseq_run = function(cellshere){

  # Get counts ----
  # cols are MIR, N_A, N_B, N_C, Z_A, Z_B, Z_C
  cts.raw = read.csv(paste0(
    "./ct-", cellshere, ".csv"
    ), stringsAsFactors = F
  ) %>%
    filter(!MIR %in% dme_outliers) %>%
    mutate(MIR = sub("(hsa|mmu|dme)\\-", "",     MIR),
           MIR = sub("mir\\-",           "miR-", MIR)
           )

  # Read cutoff ----
  cts.notenoughreads = cts.raw %>%
    pivot_longer(cols = -MIR, names_to = "expe", values_to = "ct") %>%
    group_by(MIR) %>%
    summarize_at("ct", min) %>%
    filter(ct < 5) %>%  # at least 5 in all
    ungroup() %>%
    pull(MIR)

  # Get input for DESeq2 ----
  cts = cts.raw %>%
    column_to_rownames("MIR") %>%
    as.matrix()
  condhere = colnames(cts)
  # print(head(cts))
  coldata = data.frame(
      condition = c(rep("untreated",3),rep("treated",3)),
      cond = c("N_A","N_B","N_C","Z_A","Z_B","Z_C")
    ) %>%
    filter(cond %in% condhere) %>%
    mutate(cond =      factor(cond,      levels = condhere),
           condition = factor(condition, levels = c("untreated","treated"))) %>%
    column_to_rownames(var = "cond")
  # print(coldata)

  # DESeq2 call ----
  print("==============================")
  print(cellshere)

  dds = DESeqDataSetFromMatrix(countData = cts,
                               colData = coldata,
                               design = ~ condition)
  dds = DESeq(dds,
              test = "Wald",
              minReplicatesForReplace = Inf,
              betaPrior = F,
              fitType = "parametric",
              quiet = T)
  # resultsNames(dds)

  # DESeq2 Results ----
  res = results(dds,
                independentFiltering = F,
                alpha = 0.05,
                pAdjustMethod = "BH"
                )

  # plotMA(res)  # check DESeq2 plot

  df.res = as.data.frame(res) %>%
    rownames_to_column("miRNA")

  # Carry out BBUM correction and signif calling ----
  df.res.corr = df.res %>%
    # bind_rows(., df.res%>%dplyr::slice(1:5)%>%
    #             mutate(miRNA = "ARTIFICIALOUTLIER",
    #                    log2FoldChange = -6,
    #                    pvalue = 1E-300)
    # ) %>%
    BBUM_DEcorr_TDMD(miRNA.col = "miRNA",
                     fails.cutoff = cts.notenoughreads,
                     auto_outliers = T,
                     pBBUM.alpha = 0.05)

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
    left_join(cts %>%
                as.data.frame() %>%
                rownames_to_column("miRNA"),
              by = "miRNA") %>%
    left_join(df.wt,
              by = "miRNA")

  # Log output ----
  write.csv(df.out,
            file = paste0(
              "./dfout-", cellshere, ".csv"
            ),
            row.names = F, quote = F)

  return(df.out)
}
run.res.Shi = sapply(cellsv,
                     sRNAseq_run,
                     USE.NAMES = T, simplify = F)

BBUM_plot(
  run.res.Shi$"MEF",
  # option = "symm",
  option = "MA",
  # option = "pcorr",
  # option = "ecdf_log",
  # option = "pp",
  # option = "volcano",
  expressionCol = "WTmean"
)

