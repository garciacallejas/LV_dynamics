# adapted from https://stefanoallesina.github.io/Sao_Paulo_School/multi.html#how-many-species-will-coexist
# and my own SME repo

library(tidyverse)
library(deSolve)

# Generalized Lotka-Volterra model
GLV <- function(t, x, parameters){
    with(as.list(c(x, parameters)), {
        x[x < 10^-8] <- 0 # prevent numerical problems
        dxdt <- x * (r + A %*% x)
        list(dxdt)
    })
}

GLV.positive.competition <- function(t, x, parameters){
    with(as.list(c(x, parameters)), {
        x[x < 10^-8] <- 0 # prevent numerical problems
        dxdt <- x * (r - A %*% x)
        list(dxdt)
    })
}

# function to clean output into a long dataframe
tidy_ODE_output <- function(out){
    out <- as.data.frame(out)
    # colnames(out) <- c("time", paste("sp", 1:(ncol(out) -1), sep = "_"))
    out <- as_tibble(out) %>% gather(species, density, -time)
    # pl <- ggplot(data = out) + 
    #     aes(x = time, y = density, colour = species) + 
    #     geom_line()
    # show(pl)
    return(out)
}

# general function to integrate GLV
integrate_GLV <- function(r, A, x0, maxtime = 100, steptime = 0.5,positive.competition = FALSE){
    times <- seq(0, maxtime, by = steptime)
    parameters <- list(r = r, A = A)
    # solve numerically
    if(positive.competition){
        out <- ode(y = x0, times = times, 
                   func = GLV.positive.competition, 
                   parms = parameters, 
                   method = "ode45")
    }else{
        out <- ode(y = x0, times = times, 
                   func = GLV, parms = parameters, 
                   method = "ode45")
    }
    # make into tidy form
    out2 <- tidy_ODE_output(out)
    return(out2)
}

# from Allesina, builds a globally stable matrix
build_LDstable <- function(n){
    A <- matrix(0, n, n)
    A[upper.tri(A)] <- rnorm(n * (n - 1) / 2)
    # make symmetric
    A <- A + t(A)
    # now find the largest eigenvalue
    l1A <- max(eigen(A, only.values = TRUE, symmetric = TRUE)$values)
    if (l1A > 0){
        # set the diagonal to make it stable
        diag(A) <- diag(A) - l1A - 0.01
    }
    return(A)
}

# builds an interaction matrix with different parameters
# S = richness
# c = connectance
# tau = inverse of kurtosis (tau ~1.5 gives kurtosis ~3 for S=25)
# see https://statisticaloddsandends.wordpress.com/2019/04/15/the-sinh-arcsinh-normal-distribution/
# min.diag.dom = the minimum diagonal dominance wanted in each row
# restricted.positive = whether interactions are only >0 or not
# int.mean = mean of interaction strength (taken from a normal distribution)
# ind.sd = standard deviation of interaction strength
horizontal_community_matrix <- function(S = 5,
                                        c = 0.5,
                                        tau = 1.5,
                                        min.diag.dom = 0,
                                        restricted.positive = TRUE,
                                        int.mean = 0,
                                        int.sd = 1){
    
    a.rows <- S
    a.cols <- S
    l <- round(c * (a.rows*a.cols))
    
    A <- matrix(0,nrow = a.rows,ncol = a.cols)
    if(restricted.positive){
        ints <- abs(gamlss.dist::rSHASHo(l, mu = int.mean, 
                                         sigma = int.sd, nu = 0, tau = tau))
    }else{
        ints <- gamlss.dist::rSHASHo(l, mu = int.mean, 
                                     sigma = int.sd, nu = 0, tau = tau)
    }
    
    # randomly assign interaction strengths outside the diagonal
    for(i in 1:l){
        my.sample.row <- sample(1:a.rows,1,replace = T)
        my.sample.col <- sample(1:a.cols,1,replace = T)
        
        while(A[my.sample.row,my.sample.col] != 0 & 
              my.sample.row == my.sample.col){
            my.sample.row <- sample(1:a.rows,1,replace = T)
            my.sample.col <- sample(1:a.cols,1,replace = T)
        }
        A[my.sample.row,my.sample.col] <- ints[i]
    }# for i
    
    # diag values
    for(i.row in 1:a.rows){
        non.diag <- abs(sum(A[i.row,]))
        if(min.diag.dom > 0){
            # values around that needed to achieve dominance in this row
            A[i.row,i.row] <- abs(rnorm(1,mean = (non.diag + min.diag.dom),sd = .1))
        }else{
            # values from the same distribution as the rest
            A[i.row,i.row] <- abs(gamlss.dist::rSHASHo(1, mu = int.mean, 
                                                       sigma = int.sd, nu = 0, tau = tau))
        }
        # cat(i.row,"diag:",A[i.row,i.row],"-non diag:",non.diag,"-dominance:",A[i.row,i.row] - non.diag,"\n")
    }
    return(A)
}



