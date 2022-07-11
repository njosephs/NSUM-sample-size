# libraries ----
library(ggplot2)

# q <- seq(.001, .1, .0001) # heat map
q <- c(.005, .01, .025, .05, .1, .2)
# d_bar <- seq(1, 25, .025)
d_bar <- seq(1, 100, .025)
M <- 10^(3:5)

samp_size <- function(vals, alpha = .05, eps = .1) {
  q <- vals[[1]]
  d_bar <- vals[[2]]
  M <- vals[[3]]
  
  z_alpha <- qnorm(1 - alpha/2)
  
  n <- z_alpha^2 / eps^2 * (1/q) * (1/d_bar - 1/M)
  
  return(n)
}

grid <- expand.grid(q, d_bar, M)
names(grid) <- c("q", "d_bar", "M")

df <- cbind(grid, apply(as.matrix(grid), 1, samp_size))
names(df) <- c("q", "d_bar", "M", "n")
df$M <- factor(df$M, labels = c("1,000", "10,000", "100,000"))

# line plot ----
df$q <- factor(df$q, labels = c("0.5%", "1%", "2.5%", "5%", "10%", "20%"))

colors <- setNames(viridis::viridis(length(unique(df$q))+1, option = "C")
                   , levels(df$q))[-(length(unique(df$q))+1)] # skip yellow

ggplot(df[df$M == "10,000", ], aes(x = d_bar, y = n, group = q, color = q)) +
  geom_line() +
  ylim(c(0, 10000)) +
  labs(x = "Average degree", y = "Minimum sample size", color = "Prevalence") +
  scale_color_manual(values = colors) +
  theme_bw(base_size = 12)