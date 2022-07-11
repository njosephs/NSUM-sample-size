# libraries ----
library(parallel)
library(igraph)

# simulation params ----
M <- 1000                   # |V|, size of population
N <- .15*M                  # prevalence
deviation <- seq(0, 1, .01) # deviation away from vertex exchangeable
p <- .1                     # underlying population graph density
model <- c("PA", "SBM", "ERGM (+)", "ERGM (-)")
alpha <- .05
eps <- .1
reps <- 500

grid <- expand.grid(model, deviation, stringsAsFactors = FALSE)
colnames(grid) <- c("model", "deviation")

# simulation ----
res <- matrix(nrow = nrow(grid), ncol = 6)
colnames(res) <- c("model", "deviation"
                   , "coverage", "coverage.sd"
                   , "error", "error.sd")

for (g in seq_len(nrow(grid))) {
  model <- grid[g, "model"]
  deviation <- grid[g, "deviation"]
  
  z_alpha <- qnorm(1 - alpha/2)
  
  start <- Sys.time()
  vals <- simplify2array(
    mclapply(seq_len(reps), mc.cores = detectCores(), FUN = function (r) {
      set.seed(r)
      
      if (model == "PA") {
        
        # sample true population network G
        G0 <- erdos.renyi.game(M*(1-deviation), p, type = "gnp")
        
        # change proportion of G from ER --> PA
        if (deviation == 0) {
          G <- G0
        } else if (deviation == 1) {
          G <- sample_pa(n = M, power = 1.4, m = M/20
                         , out.pref = TRUE, directed = FALSE)
        } else {
          G <- sample_pa(n = M, power = 1.4, m = M/20
                         , start.graph = G0
                         , out.pref = TRUE, directed = FALSE)
        }
        
        A <- get.adjacency(G)
        
      } else if (model == "SBM") {
        
        pm <- cbind(c(p * (1+deviation), p * (1-deviation))
                    , c(p * (1-deviation), p * (1+deviation)))
        
        G <- sample_sbm(M, pm, c(.5*M, .5*M))
        
        A <- get.adjacency(G)
        
      } else {
        
        tmp <- network(M, directed = FALSE, density = p)
        
        beta <- ifelse(model == "ERGM (+)", deviation, -deviation)
        
        G <- simulate(tmp ~ triangle
                      , coef = beta)
        
        A <- as.matrix.network.adjacency(G)
        
      }
      
      # use true density as estimate for p
      p <- sum(A) / (M*(M-1))
      
      # sample size and variance
      n <- ceiling(z_alpha^2/eps^2 * (1-p)/(N*p))
      v_n <- N/n  * (1-p)/p
      
      # compute d_i and d_i^u (sample N as hidden population)
      d_i <- Matrix::rowSums(A)
      d_i_u <- Matrix::rowSums(A[, sample.int(M, N)])
      
      # sample n from population
      s <- sample.int(M, n)
      
      # estimate N
      N_hat <- M * sum(d_i_u[s]) / sum(d_i[s])
      
      # check coverage
      coverage <- N >= N_hat - z_alpha*sqrt(v_n) && N <=  N_hat + z_alpha*sqrt(v_n)
      
      # check relative error
      error <- abs(N_hat - N) / N
      
      return(c(coverage, error))
    }))
  
  message(g, " (", round(difftime(Sys.time(), start, units = "secs")), " secs) : "
          , " model = ", model
          , " deviation = ", deviation
          , " , coverage = ", mean(vals[1, ])
          , " , error = ", round(mean(vals[2, ]), 3))
  
  res[g, "model"] <- model
  res[g, "deviation"] <- deviation
  res[g, "coverage"] <- mean(vals[1, ])
  res[g, "coverage.sd"] <- sd(vals[1, ])
  res[g, "error"] <- mean(vals[2, ])
  res[g, "error.sd"] <- sd(vals[2, ])
}

# plot ----
library(ggplot2)

df <- as.data.frame(res)
df[c(2:6)] <- lapply(df[c(2:6)], function(x) as.numeric(x))

colors <- setNames(viridis::viridis(length(unique(df$model))+1, option = "C")
                   , levels(df$model))[-(length(unique(df$model))+1)] # skip yellow

ggplot(df, aes(x = deviation, y = error, group = model, color = model)) +
  geom_line() +
  geom_point(position = position_dodge(0.005), size = .5) +
  geom_errorbar(aes(ymin = error - error.sd
                    , ymax = error + error.sd)
                , width = .005, position = position_dodge(0.005), alpha = .25) +
  geom_hline(yintercept = eps, linetype = "dashed", color = "red") +
  ylim(c(0, .25)) +
  labs(y = expression("Relative error")
       , x = "Deviation from vertex-exchangeable"
       , color = "Model"
       , shape = "") +
  scale_color_manual(values = colors) +
  theme_bw(base_size = 12)

ggplot(df, aes(x = deviation, y = coverage, group = model, color = model)) +
  geom_line() +
  geom_point(position = position_dodge(0.005), size = .5) +
  geom_errorbar(aes(ymin = pmax(0, coverage - coverage.sd)
                    , ymax = pmin(1, coverage + coverage.sd))
                , width = .005, position = position_dodge(0.005), alpha = .25) +
  geom_hline(yintercept = 1-alpha, linetype = "dashed", color = "red") +
  ylim(c(.25, 1)) +
  labs(y = "Coverage"
       , x = "Deviation from vertex-exchangeable"
       , color = "Model") +
  scale_color_manual(values = colors) +
  theme_bw(base_size = 12)