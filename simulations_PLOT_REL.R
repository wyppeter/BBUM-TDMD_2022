# Simulation script for BBUM benchmarking: plotting results
# Peter Y. Wang 2022
# Bartel Lab, Whitehead Institute/MIT

# Import libraries ---------------
library(tidyverse)
library(ggsci)
library(EnvStats)
library(boot)

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

# datatype = "-out"
datatype = "+out"
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


mean_forboot = function(x, i){
  mean(x[i], na.rm = T)
}
cv_forboot = function(x, i){
  cv(x[i], na.rm = T)
}

eval_plot_summstats = eval_plot %>%
  group_by(callmethod, statistic) %>%
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
      select(callmethod, statistic) %>%
      distinct() %>%
      mutate(v_mu = mean_val,
             v_mu.lo = ci_mean[1],
             v_mu.hi = ci_mean[2],
             v_cv = cv_val,
             v_cv.lo = ci_cv[1],
             v_cv.hi = ci_cv[2]
             )
  })

eval_zeroes = eval_plot %>%
  filter(value %in% c(0,1)) %>%
  group_by(callmethod, statistic, value) %>%
  tally() %>%
  select(callmethod, statistic, value, n) %>%
  distinct() %>%
  mutate(callmethod = factor(
    as.character(callmethod),
    levels = c("SHI.outl","BBUM.naiv","BBUM.auto")
  )) %>%
  arrange(callmethod, statistic, value)

###################

# Fig 4A/B
eval_plot %>%
  mutate(callmethod = factor(
    as.character(callmethod),
    levels = c("SHI.outl","BBUM.naiv","BBUM.auto")
    )) %>%
  ggplot(aes(
    x = value
  )) +
  facet_grid(callmethod~statistic) +
  geom_histogram(stat = "bin", fill = "mediumpurple4",
                 # binwidth = 0.02,
                 binwidth = 0.01,
                 color = NA) +
  geom_text(data = eval_zeroes,
            aes(label = n, x = value*0.8+0.1), y = Inf, vjust = 2,
            color = "mediumpurple4", size = 2
  ) +
  coord_cartesian(ylim = c(0,400)) +
  scale_x_continuous(breaks = c(0,0.05,0.5,1)) +
  labs(title = paste(datatype)) +
  theme_classic(base_size = 12)

# View(eval_plot_summstats)

# write.csv(eval_plot_summstats,
#           paste0("eval_stats_",datatype,"_boot.csv"),row.names = F,quote = F)
