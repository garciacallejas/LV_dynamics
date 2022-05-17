
library(tidyverse)

# -------------------------------------------------------------------------
# functions for solving Lotka-Volterra systems
source("R/GLV_functions.R")

# -------------------------------------------------------------------------
# exponential growth
r <- c(-0.03, -0.02, 0, 0.02, 0.03)
N0 <- 2
time <- seq(1,100,1)
cont.mat <- sapply(r, function(r) N0 * exp(r * time))
cont.df <- as.data.frame(cont.mat)
names(cont.df) <- as.character(r)
cont.df$timestep <- time
cont.long <- cont.df %>% pivot_longer(cols = 1:5, names_to = "growth_rate",
                                      values_to = "density") 

ggplot(cont.long,aes(x = timestep, y = density, group = growth_rate)) + 
  geom_line(aes(color = growth_rate)) + 
  theme_bw() +
  NULL

# -------------------------------------------------------------------------
logistic_growth <- function(times, y, parms) {
  n <- y[1]
  r <- parms[1]
  K <- parms[2]
  dN.dt <- r * n * (1 - n/K)
  return(list(c(dN.dt)))
}

library(deSolve)
prms <- c(r = 1, K = 100)
init.N <- c(1)
t.s <- seq(0.1, 50, by = 0.01)
out <- ode(y = init.N, times = t.s, logistic_growth, parms = prms)
out.df <- as.data.frame(out)
names(out.df)[2] <- "density"

ggplot(out.df, aes(x = time, y = density)) + 
  geom_line() +
  theme_bw() +
  NULL

# different growth rates and initial conditions
out.list <- list()
for (i in 1:20) {
  y <- runif(n = 1, min = 0, max = 150)
  prms <- c(r = runif(1, 0.01, 2), K = 100)
  i.out <- as.data.frame(ode(y, times = t.s, logistic_growth, prms))
  names(i.out)[2] <- "density"
  i.out$replicate <- i
  out.list[[i]] <- i.out
}
out.df <- bind_rows(out.list)

ggplot(out.df, aes(x = time, y = density, group = replicate)) + 
  geom_line() +
  theme_bw() +
  NULL

# -------------------------------------------------------------------------
# Lotka-Volterra systems

# unstable system

set.seed(1) # for reproducibility
r_1 <- rep(1, 3)
A_1 <- -matrix(c(10, 9, 5, 
                 9, 10, 9, 
                 5, 9, 10), 3, 3, byrow = TRUE)
# check the existence of feasible equilibrium
# all eigenvalues <0?
print(solve(A_1, -r_1)) # not feasible
x0_1 <- runif(3)
res_1 <- integrate_GLV(r_1, A_1, x0_1)
lv1 <- ggplot(data = res_1) +
  aes(x = time, y = density, colour = species) + #ylim(0,1e5) +
  geom_line()
plot_eigenvalues(A_1)

# -------------------------------------------------------------------------
# stable

set.seed(2) # for reproducibility
r_2 <- rep(10, 3)
A_2 <- -matrix(c(10, 7, 12, 
                 15, 10, 8, 
                 7, 11, 10), 3, 3, byrow = TRUE)
# check the existence of feasible equilibrium
print(solve(A_2, -r_2)) # feasible
x0_2 <- runif(3)
res_2 <- integrate_GLV(r_2, A_2, x0_2)
lv2 <- ggplot(data = res_2) +
  aes(x = time, y = density, colour = species) + #ylim(0,1e5) +
  geom_line()
plot_eigenvalues(A_2)

# -------------------------------------------------------------------------
# stable limit cycles

set.seed(3) # for reproducibility
r_3 <- rep(1, 3)
A_3 <- -matrix(c(10, 6, 12, 
                 14, 10, 2, 
                 8, 18, 10), 3, 3, byrow = TRUE)
# check the existence of feasible equilibrium
print(solve(A_3, -r_3)) # feasible
x0_3 <- 0.1 * runif(3)
res_3 <- integrate_GLV(r_3, A_3, x0_3, maxtime = 250)
lv3 <- ggplot(data = res_3) +
  aes(x = time, y = density, colour = species) + #ylim(0,1e5) +
  geom_line()
plot_eigenvalues(A_3)

# -------------------------------------------------------------------------
# chaos

set.seed(4) # for reproducibility
r_4 <- c(1, 0.72, 1.53, 1.27)
A_4 <- -matrix(c(1, 1.09, 1.52, 0, 
                 0, 0.72, 0.3168, 0.9792, 
                 3.5649, 0, 1.53, 0.7191,
                 1.5367, 0.6477, 0.4445, 1.27), 4, 4, byrow = TRUE)
# check the existence of feasible equilibrium
print(solve(A_4, -r_4)) # feasible
x0_4 <- 0.1 * runif(4)
res_4 <- integrate_GLV(r_4, A_4, x0_4, maxtime = 500)
lv4 <- ggplot(data = res_4) +
  aes(x = time, y = density, colour = species) + #ylim(0,1e5) +
  geom_line()
plot_eigenvalues(A_4)

# -------------------------------------------------------------------------
# a slightly more ecologically realistic interaction matrix

# richness
S <- 25
# connectance
c <- 0.2
# generate interaction matrix
A_5 <- horizontal_community_matrix(S,c)
# growth rates
r <- rep(1,S)
# initial abundances
x0 <- rep(.1,S)

# solve the system
lv.dynamics <- integrate_GLV(r = r,
                             A = A_5,
                             x0 = x0,
                             positive.competition = TRUE)

print(solve(A_5, -r)) # quite a few species go extinct

# check the resulting dynamics
head(lv.dynamics)

# plot
pl <- ggplot(data = lv.dynamics) +
    aes(x = time, y = density, colour = species) + #ylim(0,1e5) +
    geom_line()

#pl
plot_eigenvalues(A_5)

# -------------------------------------------------------------------------

# a binary food web
S <- 25
c <- 0.2
A_6 <- niche(S,c)

# set interaction strengths
A_6[A_6 == 1] <- runif(sum(A_6 == 1))
# intraspecific strength. Should this be zero/one, or something else?
diag(A_6) <- 1
plot_eigenvalues(A_6)


