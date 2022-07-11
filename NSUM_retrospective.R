# libraries ----
library(parallel)
library(ggplot2)
library(ggrepel)

# parameters ----
retro_data <- read.csv("./NSUM_retrospective.csv")
n_studies <- nrow(retro_data)
rel_error.avg <- rel_error.025 <- rel_error.975 <- rep(NA, n_studies)

results <- cbind(retro_data, rel_error.avg, rel_error.025, rel_error.975)
reps <- 10000

for (s in seq_len(n_studies)) {
  n <- results[s, "n"]
  M <- results[s, "M"]
  N <- results[s, "N"]
  di <- results[s, "di"]
  diu <- results[s, "diu"]
  
  start <- Sys.time()
  vals <- simplify2array(
    mclapply(seq_len(reps), mc.cores = detectCores(), FUN = function (r) {
      set.seed(r)
      
      # sample di
      di <- rbinom(n = n, size = M, prob = di/M)
      
      # sample diu
      diu <- rbinom(n = n, size = N, prob = diu/N)
      
      # NSUM
      N_hat <- M * sum(diu) / sum(di)
      
      # relative error
      return(abs(N_hat - N) / N)
    }))
  
  message(s)
  
  results[s, "rel_error.avg"] <- mean(vals)
  results[s, "rel_error.025"] <- quantile(vals, .025)
  results[s, "rel_error.975"] <- quantile(vals, .975)
}

df <- as.data.frame(results)
df$n_min <- ceiling(400 * (df$M / df$N) * (1 / df$di))
df$study <- retro_data$Study

# prevalence is consistent with degrees
plot(df$N / df$M, df$diu / df$di)
abline(coef = c(0, 1))

ggplot(df, aes(x = n, y = n_min, label = study)) +
  geom_point(aes(color = di), size = 3.5) +
  geom_abline(slope = 1, intercept = 0) +
  geom_label_repel(size = 1.5, force = 2.5
                   , box.padding = .75, min.segment.length = 0
                   , aes(segment.size = .35)) +
  scale_color_viridis_c(option = "C", direction = -1) +
  scale_x_continuous(limits = c(0, 8000)) +
  scale_y_continuous(limits = c(0, 8000)) +
  labs(x = "Study sample size"
       , y = "Minimum sample size from (6)"
       , color = "Average\ndegree") +
  guides(size = "none") +
  coord_fixed(ratio = 1) +
  theme_bw(base_size = 12)
