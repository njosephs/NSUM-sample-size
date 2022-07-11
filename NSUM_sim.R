# libraries ----
library(parallel)
library(igraph)
library(ergm)

# parameters ----
M <- c(1000, 5000, 10000)      # |V|, size of population
N <- seq(.01, .51, .02)        # prevalence
eps <- .1                      # relative error tolerance
alpha <- c(.01, .05, .1, .2)   # probability threshold
reps <- 500                    # number of replicates
G_model <- c("ER", "SBM", "small-world", "PA", "ERGM")

grid <- expand.grid(M, N, alpha, G_model, stringsAsFactors = FALSE)
colnames(grid) <- c("M", "N", "alpha", "G_model")

# ERGM too slow when M > 1000...
grid <- grid[!(grid$G_model == "ERGM" & grid$M > 1000), ]

res <- matrix(nrow = nrow(grid), ncol = 8)
colnames(res) <- c("M", "N", "alpha", "G_model"
                   , "coverage", "coverage.sd"
                   , "error", "error.sd")

for (g in seq_len(nrow(grid))) {
  M <- grid[g, "M"]
  N <- grid[g, "N"] * M
  G_model <- grid[g, "G_model"]
  alpha <- grid[g, "alpha"]
  
  z_alpha <- qnorm(1 - alpha/2)
  
  start <- Sys.time()
  vals <- simplify2array(
    mclapply(seq_len(reps), mc.cores = detectCores(), FUN = function (r) {
      set.seed(r)
      
      # sample true population network G
      if (G_model == "ER") {
        p <- .1
        G <- erdos.renyi.game(M, p, type = "gnp")
        
      } else if (G_model == "SBM") {
        pm <- cbind(c(.9, .75, .5)
                    , c(.75, .6, .25)
                    , c(.5, .25, .1)) / 5
        G <- sample_sbm(M, pm, c(.4*M, .35*M, .25*M))
        
      } else if (G_model == "small-world") {
        p <- .1
        G <- sample_smallworld(1, M, 50, p)
        
      } else if (G_model == "PA") {
        p <- 1.4
        G <- sample_pa(M, power = p, m = M/20, directed = FALSE)
        
      } else if (G_model == "ERGM") {
        tmp <- network(M, directed = FALSE, density = .1)
        
        beta <- -1
        G <- simulate(tmp ~ edges + kstar(2) + triangle
                      , coef = c(0, -1, beta))
      }
      
      if (G_model == "ERGM") {
        A <- as.matrix.network.adjacency(G)
      } else {
        A <- get.adjacency(G)
      }
      
      # use true density as estimate for p
      p <- sum(A) / (M*(M-1))
      
      # compute d_i
      d_i <- Matrix::rowSums(A)
      
      # compute d_i^u by assigning N to hidden population
      d_i_u <- Matrix::rowSums(A[, sample.int(M, N)])
      
      # use minimum sample size
      n <- ceiling(z_alpha^2/eps^2 * (1-p)/(N*p))
      v_n <- N/n  * (1-p)/p
      
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
          , " M = ", M
          , " , N = ", N
          , " , alpha = ", alpha
          , " , G_model = ", G_model
          , " , coverage = ", mean(vals[1, ])
          , " , error = ", round(mean(vals[2, ]), 3))
  
  res[g, "M"] <- M
  res[g, "N"] <- N
  res[g, "alpha"] <- alpha
  res[g, "G_model"] <- G_model
  res[g, "coverage"] <- mean(vals[1, ])
  res[g, "coverage.sd"] <- sd(vals[1, ])
  res[g, "error"] <- mean(vals[2, ])
  res[g, "error.sd"] <- sd(vals[2, ])
}

# plot ----
library(ggplot2)

df <- as.data.frame(res)
df[c(1:2, 5:8)] <- lapply(df[c(1:2, 5:8)], function(x) as.numeric(as.character(x)))

empty_grid <- expand.grid(M, N, alpha, "ERGM", stringsAsFactors = FALSE)
colnames(empty_grid) <- c("M", "N", "alpha", "G_model")
empty_grid$coverage <- empty_grid$coverage.sd <- empty_grid$error <- empty_grid$error.sd <- NA

df <- rbind(df, empty_grid)
df$hidden_prop <- df$N / df$M
df$col <- ifelse(df$G_model == "ERGM" & df$M > 1000, "black", NA)

colors <- setNames(viridis::viridis(length(unique(df$alpha))+1, option = "C")
                   , levels(df$alpha))[-(length(unique(df$alpha))+1)] # skip yellow

ggplot(df, aes(x = hidden_prop, y = error, group = alpha, color = alpha)) +
  geom_line() +
  geom_point(position = position_dodge(0.01), size = .5) +
  geom_errorbar(aes(ymin = error - error.sd
                    , ymax = error + error.sd)
                , width = .01, position = position_dodge(0.01), alpha = .75) +
  geom_hline(yintercept = eps, linetype = "dashed", color = "black") +
  geom_rect(fill = df$col, inherit.aes = FALSE
            , xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
  ylim(c(0, .25)) +
  labs(y = "Relative error"
       , x = "Prevalence"
       , color = expression(alpha)) +
  facet_grid(~G_model~M) +
  scale_color_manual(values = colors) +
  theme_set(theme_bw() + theme(legend.key=element_blank())) +
  theme_bw(base_size = 12)

ggplot(df, aes(x = hidden_prop, y = coverage, group = alpha, color = alpha)) +
  geom_line() +
  geom_point(position = position_dodge(0.01), size = .5) +
  geom_errorbar(aes(ymin = pmax(0, coverage - coverage.sd)
                    , ymax = pmin(1, coverage + coverage.sd))
                , width = .01, position = position_dodge(0.01), alpha = .75) +
  geom_hline(yintercept = 1-as.numeric(as.character(unique(df$alpha)))
             , linetype = "dashed", color = "black") +
  geom_rect(fill = df$col, inherit.aes = FALSE
            , xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
  ylim(c(.25, 1)) +
  labs(y = "Coverage"
       , x = "Prevalence"
       , color = expression(alpha)
       , fill = NULL) +
  facet_grid(~G_model~M) +
  scale_color_manual(values = colors) +
  theme_bw(base_size = 12)