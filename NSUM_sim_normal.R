# libraries ----
library(parallel)
library(igraph)
library(ergm)

# parameters ----
M <- 1000                      # |V|, size of population
N <- .025*M                    # prevalence
eps <- .1                      # relative error tolerance
alpha <- .05                   # probability threshold
reps <- 500                    # number of replicates
G_model <- c("ER", "SBM", "small-world", "PA", "ERGM")

z_alpha <- qnorm(1 - alpha/2)
p <- .1
n <- seq(10, .75*M, length.out = 50)

grid <- expand.grid(n, G_model, stringsAsFactors = FALSE)
colnames(grid) <- c("n", "G_model")

res <- matrix(nrow = nrow(grid), ncol = 5)
colnames(res) <- c("n", "G_model", "Shapiro", "KS", "KL")

for (g in seq_len(nrow(grid))) {
  n <- ceiling(grid[g, "n"])
  G_model <- grid[g, "G_model"]
  
  start <- Sys.time()
  vals <- mclapply(seq_len(reps), mc.cores = detectCores(), FUN = function (r) {
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
      pp <- 1.4
      G <- sample_pa(M, power = pp, m = M/20, directed = FALSE)
      
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
    
    # compute d_i
    d_i <- Matrix::rowSums(A)
    
    # compute d_i^u by assigning N to hidden population
    d_i_u <- Matrix::rowSums(A[, sample.int(M, N)])
    
    # sample n from population
    s <- replicate(reps, sample.int(M, n))
    
    # estimate N
    N_hat <- M * apply(s, 2, function(rep) sum(d_i_u[rep]) / sum(d_i[rep]))
    N_hat_sd <- sqrt(N / n * (1-p) / p)
    
    return(list(shapiro = shapiro.test(N_hat)$p.value
                , KS = ks.test(N_hat, "pnorm", N, N_hat_sd)$p.value
                , KL = FNN::KL.divergence((N_hat - N) / N_hat_sd
                                          , dnorm(seq(-3, 3, length.out = reps)))[10]))
  })
  
  message(g, " (", round(difftime(Sys.time(), start, units = "secs")), " secs)")
  
  res[g, "n"] <- n
  res[g, "G_model"] <- G_model
  res[g, "Shapiro"] <- mean(sapply(vals, "[[", "shapiro") < alpha)
  res[g, "KS"] <- mean(sapply(vals, "[[", "KS") < alpha)
  res[g, "KL"] <- mean(sapply(vals, "[[", "KL"))
}

# plot ----
library(dplyr)
library(ggplot2)

df <- as.data.frame(res)
df[c(1, 3:5)] <- lapply(df[c(1, 3:5)], function(x) as.numeric(as.character(x)))

ggplot(df, aes(n, Shapiro)) +
  geom_point() +
  geom_line() +
  geom_hline(yintercept = alpha, linetype = "dashed", color = "black") +
  geom_vline(xintercept = ceiling(z_alpha^2/eps^2 * (1-p)/(N*p)), linetype = "dashed", color = "black") +
  labs(y = "% of significant Shapiro-Wilk tests"
       , x = "n") +
  scale_x_continuous(breaks = c(0, ceiling(z_alpha^2/eps^2 * (1-p)/(N*p))
                                , 200, 400, 600)
                     , labels = c("0", "min n", "200", "400", "600")) +
  ylim(c(0, .4)) +
  facet_grid(rows = vars(G_model)) +
  theme_set(theme_bw() + theme(legend.key = element_blank())) +
  theme_bw(base_size = 12)
