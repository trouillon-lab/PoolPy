
library(binGroup2)


res1 <- OTC1(
  algorithm = "A2",      # three-stage hierarchical
  p         = 0.05,      # prevalence
  Se        = 1,      # sensitivity
  Sp        = 1,      # specificity
  group.sz  = 3:15,     # range of initial group sizes to consider
  obj.fn    ='ET'# c("ET", "MAR", "GR1")
  # ET  = expected number of tests per individual
  # MAR = misclassification-adjusted risk
  # GR1 = generalized risk function
)

# Extract optimal configuration (design) for ET, etc.
res1$opt.ET$OTC 
res1$opt.ET$value


summary(res1)
res1$opt.ET$OTC  
res1$opt.ET
res1$opt.ET$ET 


res2 <- OTC1(
  algorithm = "A2",
  p         = 0.1,
  Se        = 1,
  Sp        = 1,
  group.sz  = 5:20,
  obj.fn    = c("ET", "MAR", "GR1"),
  weights   = matrix(data = c(1, 1), nrow = 1, ncol = 2)
)

Config(res2)
summary(res2)

# ── Prevalence sweep ──────────────────────────────────────────────────────────

# User-defined parameters
algorithm   <- "A2"
Se_sweep    <- 1
Sp_sweep    <- 1
output_file <- "timing_comparison_R.csv"

# Prevalence sequence: 0.5, 0.2, 0.1, 0.05, 0.02, ..., 0.0001
prevalences <- c(0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001, 0.0005, 0.0002, 0.0001)

sweep_df <- data.frame(
  prevalence = numeric(),
  Method     = character(),
  ET         = numeric(),
  time       = numeric(),
  stringsAsFactors = FALSE
)

for (p in prevalences) {
  # group.sz lower: max(2, round(-log(p)))   upper: min(200, round(1/p))
  gz_lower <- max(2L, round(-log(p)))
  gz_upper <- min(200L, round(1 / p))
  if (algorithm=="A2") {
    gz_upper <- min(200L, round(2*sqrt(1 / p)))
    }
  if (gz_lower >= gz_upper) gz_upper <- gz_lower + 10L

  cat(sprintf("Running p = %.4f  |  group.sz = %d:%d\n", p, gz_lower, gz_upper))

  timing <- system.time({
    res <- OTC1(
      algorithm  = algorithm,
      p          = p,
      Se         = Se_sweep,
      Sp         = Sp_sweep,
      group.sz   = gz_lower:gz_upper,
      obj.fn     = "ET",
      trace      = FALSE,
      print.time = FALSE
    )
  })

  et_value <- res$opt.ET$value   # expected tests per individual (already per-sample)

  sweep_df <- rbind(sweep_df, data.frame(
    prevalence = p,
    Method     = algorithm,
    ET         = et_value,
    time       = timing["elapsed"],
    stringsAsFactors = FALSE
  ))

  write.csv(sweep_df, output_file, row.names = FALSE)
  cat(sprintf("  ET = %.4f  |  time = %.2f s\n", et_value, timing["elapsed"]))
}

sweep_df

prevalences <- c(0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001, 0.0005, 0.0002, 0.0001)

sweep_df <- data.frame(
  prevalence = numeric(),
  Method     = character(),
  ET         = numeric(),
  time       = numeric(),
  stringsAsFactors = FALSE
)

for (p in prevalences) {
  # group.sz lower: max(2, round(-log(p)))   upper: min(200, round(1/p))
  gz_lower <- max(2L, round(-log(p)))
  gz_upper <- min(200L, round(1 / p))
  if (gz_lower >= gz_upper) gz_upper <- gz_lower + 10L

  cat(sprintf("Running p = %.4f  |  group.sz = %d:%d\n", p, gz_lower, gz_upper))

  timing <- system.time({
    res <- OTC1(
      algorithm  = algorithm,
      p          = p,
      Se         = Se_sweep,
      Sp         = Sp_sweep,
      group.sz   = gz_lower:gz_upper,
      obj.fn     = "ET",
      trace      = FALSE,
      print.time = FALSE
    )
  })

  et_value <- res$opt.ET$value   # expected tests per individual (already per-sample)

  sweep_df <- rbind(sweep_df, data.frame(
    prevalence = p,
    Method     = algorithm,
    ET         = et_value,
    time       = timing["elapsed"],
    stringsAsFactors = FALSE
  ))

  write.csv(sweep_df, output_file, row.names = FALSE)
  cat(sprintf("  ET = %.4f  |  time = %.2f s\n", et_value, timing["elapsed"]))
}

sweep_df
