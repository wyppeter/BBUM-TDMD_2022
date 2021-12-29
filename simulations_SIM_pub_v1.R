# Simulation script for BBUM benchmarking

# Import libraries ---------------
library(tidyverse)
# Import bbum, pre package-publishing, on computer cluster... ----
sapply( paste0(
  "./bbum_R_static/",
  list.files(path = "./bbum_R_static")), source )

# Is this a run with outliers?
OUTLIERSETUP = T
OUTFLAG = ifelse(OUTLIERSETUP, "+", "-")
SIMID = 0

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

# Generate data ----
gen_param = function(n, outliersgen = T, hitsgen = T){
  list(
    lambda =     runif(n, min = r.lambda[1], max = r.lambda[2]),
    a      =     runif(n, min = r.a[1],      max = r.a[2]),
    theta  = if_else(!rep(hitsgen,     n), rep(0, n),
              10^runif(n, min = r.theta[1],  max = r.theta[2])),
    r      =  10^runif(n, min = r.r[1],      max = r.r[2]),
    outth  = if_else(!rep(outliersgen, n), rep(0, n),
              10^runif(n, min = r.outth[1],  max = r.outth[2])),
    outr   =  10^runif(n, min = r.r[1],      max = r.r[2]),
    N      = sample(x = r.N[1]:r.N[2], size = n, replace = T)
  )
}

gen_data = function(sampleN.targ, oversample = 20,
                    outliersgen, hitsgen){

  minHITS = if_else(hitsgen,     MIN.HITS, 0)
  minOUTS = if_else(outliersgen, MIN.OUTS, 0)

  sampleN = sampleN.targ*oversample
  sampleID_v = seq(1, sampleN)

  dt.setup = data.frame(
    sampleID = paste0(
      "s", str_pad(
        sampleID_v,
        ceiling(log10(sampleN)),
        "left", "0"
        )
      )
    ) %>%
    cbind(data.frame(gen_param(sampleN,
                               outliersgen = outliersgen,
                               hitsgen = hitsgen
    ))) %>%
    mutate(N.up = rbinom(sampleN, size = N, prob = 0.5),
           N.down = N - N.up)

  dt = dt.setup %>%
    group_by(sampleID) %>%
    do(., {
      params = .
      dt.sample = data.frame(
        geneID = paste0(
          "g", str_pad(1:params$N,
                       ceiling(log10(max(dt.setup$N))),
                       "left", "0")
          ),
        FCdir = factor(c(
          rep(+1, params$N.up),
          rep(-1, params$N.down)
          ))
        ) %>%
        cbind(data.frame(mapply(c,
                                rbbum.ID(params$N.up,
                                         params$lambda, params$a,
                                         params$theta, params$r),
                                rbbum.ID(params$N.down,
                                         params$lambda, params$a,
                                         params$outth, params$outr),
                                SIMPLIFY = F)
        ))
      dt.sample %>%
        crossing(params)  # attach params as new cols
    }) %>%

    mutate(FC.up = FCdir == +1) %>%

    # Filter down to the wanted ones
    group_by(sampleID) %>%
    filter(sum(cate == 1 &  FC.up) >= minHITS,
           sum(cate == 1 & !FC.up) >= minOUTS
    ) %>%
    ungroup() %>%
    mutate(sampleID = factor(as.character(sampleID))) %>%
    filter(sampleID %in% sample(
      levels(sampleID),
      sampleN.targ,
      replace = F)
      ) %>%

    group_by(sampleID) %>%
    mutate(padj = p.adjust(pvalue, method = "BH")) %>%
    mutate(N.prim.emp = sum(cate == 1 &  FC.up),
           N.outl.emp = sum(cate == 1 & !FC.up)
    ) %>%
    ungroup()

  dt
}

# gen_data(1) %>% view

sampleN = 1000
df.simul = gen_data(sampleN,
                    hitsgen = T,
                    outliersgen = OUTLIERSETUP
                    )

# Calc cutoffs using Shi et al. 2020-esque methods ----
pick_cutoff = function(dt){

  downpvals = dt %>%
    filter(!FC.up) %>%
    pull(padj)

  logp = sort(-log10(downpvals), decreasing = T)
  topmost = ceiling(logp[1])
  topmost = if_else(topmost < -log10(0.05), -log10(0.05), topmost)
  cutoff = 10^-topmost

  dt %>%
    mutate(
      cutoff = cutoff,
      SHI.outl.hits = padj < cutoff & FC.up
    )
}

# BBUM fitters; adds coef cols to df
BBUM_fit_wrapper = function(dt){
  fitout.auto = BBUM_DEcorr(
    df.deseq = dt,
    classCol = "FC.up",
    pBBUM.alpha = pBBUM.alpha,
    auto_outliers = T,
    quiet = T
  )
  fitout.naiv = BBUM_DEcorr(
    df.deseq = dt,
    classCol = "FC.up",
    pBBUM.alpha = pBBUM.alpha,
    auto_outliers = F,
    quiet = T
  )

  bind_cols(dt,
            fitout.auto %>%
              transmute(pBBUM.auto = pBBUM,
                        BBUM.auto.hits = BBUM.hits,
                        BBUM.auto.fct = BBUM.fct),
            fitout.naiv %>%
              transmute(pBBUM.naiv = pBBUM,
                        BBUM.naiv.hits = BBUM.hits,
                        BBUM.naiv.fct = BBUM.fct),
  )
}

df.simul.fit = df.simul %>%
  group_by(sampleID) %>%
  do(BBUM_fit_wrapper(.)) %>%
  group_by(sampleID) %>%
  do(pick_cutoff(.)) %>%
  ungroup()

write.csv(df.simul.fit,
          file = paste(
            "./sim_data",
            SIMID,
            OUTFLAG,
            "out",
            "20211019.csv",
            sep = "_"
          ),
          quote = F,
          row.names = F
          )
