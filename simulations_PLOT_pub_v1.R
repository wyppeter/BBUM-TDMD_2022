# Simulation script for BBUM benchmarking: plotting results

# Import libraries ---------------
library(tidyverse)
library(ggsci)
library(EnvStats)

# Random param ranges ----

pBBUM.alpha = 0.05

r.lambda = c( 0.1,  0.9)  # runif
r.a      = c( 0.1,  0.9)  # runif
r.theta  = c(-1.5, -0.5)  # 10^runif, 0.03  ~ 0.3
r.outth  = c(-2.5, -1.5)  # 10^runif, 0.003 ~ 0.03
r.r      = c(-1.5, -0.5)  # 10^runif, 0.03  ~ 0.3
r.outr   = c(-2.0, -1.0)  # 10^runif, 0.01  ~ 0.1
r.N      = c( 200, 1000)   # uniform sample
MIN.HITS = 3
MIN.OUTS = 1


################################### htcomp done

datatype = "-out"
df.simul.fit = bind_rows(
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
)

################################### htcomp done


df.simul.eval = df.simul.fit %>%
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
  group_by(repN, sampleID, callmethod,
           status,
           lambda, a, theta, r, outth, outr,
           N, N.prim.emp, N.outl.emp) %>%
  tally() %>%
  pivot_wider(names_from = "status", values_from = "n") %>%
  ungroup() %>%
  bind_rows(data.frame(sampleID = "null",
                       callmethod  = "null",
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


eval_plot_summstats = eval_plot %>%
  group_by(repN, callmethod, statistic) %>%
  summarize_at("value", list(mean_v = mean, var_v = var)) %>%
  mutate(CV = sqrt(var_v)/mean_v)

eval_stats = eval_plot_summstats %>%
  group_by(callmethod, statistic) %>%
  mutate(mean.mean = mean(mean_v), mean.se = sqrt(var(mean_v)/length(mean_v)),
         mean.CI.lo = mean.mean - mean.se*1.96,
         mean.CI.hi = mean.mean + mean.se*1.96,
         CV.mean = mean(CV), CV.se = sqrt(var(CV)/length(CV)),
         CV.CI.lo = CV.mean - CV.se*1.96,
         CV.CI.hi = CV.mean + CV.se*1.96

  ) %>%
  ungroup() %>%
  select(callmethod, statistic, starts_with(c("mean.", "CV."))) %>%
  distinct() %>%
  arrange(statistic, callmethod)

eval_zeroes = eval_plot %>%
  filter(value %in% c(0,1)) %>%
  group_by(repN, callmethod, statistic, value) %>%
  tally() %>%
  select(repN, callmethod, statistic, value, n) %>%
  distinct() %>%
  arrange(repN, callmethod, statistic, value)

###################

repHere = 2
eval_plot %>%
  mutate(callmethod = factor(
    as.character(callmethod),
    levels = c("BBUM.auto","BBUM.naiv","SHI.outl")
    )) %>%
  filter(repN == repHere) %>%
  ggplot(aes(
    x = value
  )) +
  facet_grid(callmethod~statistic) +
  geom_histogram(stat = "bin", fill = "mediumpurple4", binwidth = 0.02, color = NA) +
  geom_text(data = eval_zeroes %>%
              filter(repN == repHere),
            aes(label = n, x = value*0.8+0.1), y = Inf, vjust = 3,
            color = "mediumpurple4", size = 2
  ) +
  coord_cartesian(ylim = c(0,100)) +
  scale_x_continuous(breaks = c(0.05,0.5,0.95)) +
  labs(title = paste(datatype, "rep", repHere)) +
  theme_classic()

View(eval_stats)
# write.csv(eval_stats,paste0("eval_stats_",datatype,"_normalCV.csv"),row.names = F,quote = F)
