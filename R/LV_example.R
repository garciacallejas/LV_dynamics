
library(ggplot2)

# -------------------------------------------------------------------------
# functions for solving Lotka-Volterra systems
source("R/GLV_functions.R")

# -------------------------------------------------------------------------

# richness
S <- 25

# connectance
c <- 0.2

# generate interaction matrix
A <- horizontal_community_matrix(S,c)

# growth rates
r <- rep(1,S)

# initial abundances
x0 <- rep(.1,S)

lv.dynamics <- integrate_GLV(r = r,
                             A = A,
                             x0 = x0,
                             positive.competition = TRUE)

# check the resulting dynamics
head(lv.dynamics)

# plot
pl <- ggplot(data = lv.dynamics) +
    aes(x = time, y = density, colour = species) + #ylim(0,1e5) +
    geom_line()

#pl
