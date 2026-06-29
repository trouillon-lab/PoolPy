
library(binGroup2)

Se_sweep    <- 1
Sp_sweep    <- 1
output_file <- "timing_comparison_BG2_R.csv"



sweep_df
prevalences <- c(0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001, 0.0005, 0.0002, 0.0001)

sweep_df <- data.frame(
  prevalence = numeric(),
  Method     = character(),
  ET         = numeric(),
  time       = numeric(),
  stringsAsFactors = FALSE
)

for (algorithm in c("A2M", "D2", "A2", "D3")) {
  for (p in prevalences){
    # group.sz lower: max(2, round(-log(p)))   upper: min(200, round(1/p))
    gz_lower <- max(3L, round(-log(p)))
    gz_upper <- min(200L, round(2*sqrt(1 / p)))
    if (algorithm=="A2") {
      gz_upper <- min(200L, round(5*sqrt(sqrt(1 / p))))
      }
    if (gz_lower >= gz_upper -10) gz_upper <- round(gz_lower + 10L)

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
}
sweep_df
