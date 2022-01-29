# Simulation script for BBUM benchmarking: evaluating other cutoffs
# Peter Y. Wang 2022
# Bartel Lab, Whitehead Institute/MIT

# Import libraries ---------------
library(tidyverse)
library(bbum)
library(ggsci)
library(EnvStats)
library(boot)

######################

pBBUM.alpha.set = c(0.005, 0.01, 0.02, 0.05, 0.1, 0.2)

getsimul = function(datatype){
  bind_rows(
    read.csv(
      paste0("./sim_data_0_", datatype, "_20211019.csv"),
      stringsAsFactors = F
    ) %>%
      mutate(repN = 0),
    read.csv(
      paste0("./sim_data_1_", datatype, "_20211019.csv"),
      stringsAsFactors = F
    ) %>%
      mutate(repN = 1),
    read.csv(
      paste0("./sim_data_2_", datatype, "_20211019.csv"),
      stringsAsFactors = F
    ) %>%
      mutate(repN = 2),
    read.csv(
      paste0("./sim_data_3_", datatype, "_20211019.csv"),
      stringsAsFactors = F
    ) %>%
      mutate(repN = 3),
    read.csv(
      paste0("./sim_data_4_", datatype, "_20211019.csv"),
      stringsAsFactors = F
    ) %>%
      mutate(repN = 4),
    read.csv(
      paste0("./sim_data_5_", datatype, "_20211019.csv"),
      stringsAsFactors = F
    ) %>%
      mutate(repN = 5)
  ) %>%
    select(
      # -lambda,-a,-theta,-r,-outth,-outr,
      -N,-N.up,-N.down, -N.prim.emp,-N.outl.emp,
      -padj, -cutoff,-SHI.outl.hits,
      # -pBBUM.auto, -pBBUM.naiv,
      -BBUM.auto.hits, -BBUM.naiv.hits,
      -BBUM.auto.fct, -BBUM.naiv.fct
    ) %>%
    mutate(simultype = if_else(datatype == "+out", "with", "without"))
}

df.simul.call = bind_rows(
  getsimul("-out"),
  getsimul("+out")
) %>%
  crossing(data.frame(alphahere = pBBUM.alpha.set)) %>%
  mutate(BBUM.auto.hits = FC.up & pBBUM.auto < alphahere,
         BBUM.naiv.hits = FC.up & pBBUM.naiv < alphahere)

df.simul.eval = df.simul.call %>%
  filter(FC.up) %>%
  pivot_longer(
    cols = ends_with(".hits",ignore.case = F),
    names_to = "callmethod",
    values_to = "called"
  ) %>%
  mutate(callmethod = sub("\\.hits", "", callmethod)) %>%
  mutate(
    status = case_when(
      cate == 1 &  called ~ "TP",
      cate != 1 &  called ~ "FP",
      cate == 1 & !called ~ "FN",
      cate != 1 & !called ~ "TN",
      TRUE ~ NA_character_
    )
  )

df.simul.eval.summ = df.simul.eval %>%
  group_by(simultype, repN, sampleID, callmethod, alphahere,
           status) %>%
  tally() %>%
  pivot_wider(names_from = "status", values_from = "n") %>%
  ungroup() %>%
  bind_rows(data.frame(simultype = "null",
                       sampleID = "null",
                       callmethod  = "null",
                       alphahere = NA_real_,
                       TP = NA_real_,
                       FP = NA_real_,
                       TN = NA_real_,
                       FN = NA_real_)) %>%  # prevent missing cols
  replace(is.na(.), 0) %>%
  filter(sampleID != "null") %>%
  mutate(
    sensitivity = TP/(TP+FN),
    specificity = TN/(TN+FP),
    FDR         = FP/(FP+TP)
  ) %>%
  mutate(
    FDR         = if_else(is.nan(FDR), 0 , FDR)
  )

eval_plot = df.simul.eval.summ %>%
  pivot_longer(cols = c(sensitivity, specificity, FDR), names_to = "statistic")


mean_forboot = function(x, i){
  mean(x[i], na.rm = T)
}
cv_forboot = function(x, i){
  cv(x[i], na.rm = T)
}

eval_plot_summstats = eval_plot %>%
  group_by(simultype, alphahere, callmethod, statistic) %>%
  do(., {
    dfhere = .
    bootmean = boot(
      dfhere$value, statistic = mean_forboot,
      R = 3000, sim = "ordinary"
    ) %>%
      boot.ci(type = "basic")
    bootcv = boot(
      dfhere$value, statistic = cv_forboot,
      R = 3000, sim = "ordinary"
    ) %>%
      boot.ci(type = "basic")
    ci_mean = bootmean$basic[c(4,5)]
    ci_cv = bootcv$basic[c(4,5)]
    mean_val = mean(dfhere$value, na.rm = T)
    cv_val = sqrt(var(dfhere$value, na.rm = T))/mean_val
    dfhere %>%
      select(simultype, alphahere, callmethod, statistic) %>%
      distinct() %>%
      mutate(v_mu = mean_val,
             v_mu.lo = ci_mean[1],
             v_mu.hi = ci_mean[2],
             v_cv = cv_val,
             v_cv.lo = ci_cv[1],
             v_cv.hi = ci_cv[2]
      )
  })

# write.csv(eval_plot_summstats, "./eval_stats_OTHERALPHAS.csv", row.names = F, quote = F)

######################

eval_plot_summstats = read.csv("./eval_stats_OTHERALPHAS.csv")

# Fig. 4C
eval_plot_summstats %>%
  filter(statistic == "FDR") %>%
  # Plot
  ggplot(aes(x = alphahere, y = v_mu,
             color = callmethod,
             group = interaction(callmethod,simultype))) +
  # Diagonal
  geom_abline(color = "gray90", slope = 1, intercept = 0, size = 0.5) +
  # Shi lines
  geom_hline(yintercept = 0.035554256,
             color = "olivedrab2",
             size = 0.5, alpha = 0.75,
             linetype = "solid") +
  geom_hline(yintercept = 0.005934703,
             color = "olivedrab2",
             size = 0.5, alpha = 0.75,
             linetype = "twodash") +
  # BBUM points
  geom_path(aes(linetype = simultype),
    size = 0.5, alpha = 0.75) +
  geom_errorbar(aes(ymin = v_mu.lo, ymax = v_mu.hi),
                width = 0.1, alpha = 0.75) +
  geom_point(aes(shape = callmethod), size = 2, alpha = 0.75) +
  # Scales
  scale_shape_manual(values = c(17, 16)) +
  scale_color_manual(values = c("hotpink3", "royalblue3")) +
  scale_linetype_manual(values = c("twodash","solid")) +
  # Axes
  scale_x_continuous(trans = "log10", breaks = c(0.003,0.01,0.03,0.1,0.3)) +
  scale_y_continuous(trans = "log10", breaks = c(0.003,0.01,0.03,0.1,0.3)) +
  coord_fixed(xlim = c(0.001,0.3), ylim = c(0.001,0.3)) +
  labs(x = "pBBUM FDR threshold (aBBUM)", y = "Mean empirical FDR") +
  # Done
  theme_classic(base_size = 12)
