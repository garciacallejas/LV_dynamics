# adapted from https://stefanoallesina.github.io/Sao_Paulo_School/multi.html#how-many-species-will-coexist


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

jacobian_matrix_GLV <- function(r, A, x0, maxtime = 100, steptime = 0.5,positive.competition = FALSE){
  times <- seq(0, maxtime, by = steptime)
  parameters <- list(r = r, A = A)
  # solve numerically
  
  model_stode_pos_comp = function(t,y,parms=NULL,A,growth.rate) {
    dy = y*(growth.rate-A%*%y)
    return(list(dy,1))
  }
  model_stode_neg_comp = function(t,y,parms=NULL,A,growth.rate) {
    dy = y*(growth.rate+A%*%y)
    return(list(dy,1))
  }
  

  
  if(positive.competition){
    out <- rootSolve::stode(y = x0,
                            time=0,
                            func=model_stode_pos_comp,
                            parms=NULL,
                            A=A, 
                            growth.rate=r,
                            positive = TRUE)[[1]]
    # out <- ode(y = x0, times = times, 
    #            func = GLV.positive.competition, 
    #            parms = parameters, 
    #            method = "ode45")
  }else{
    out <- rootSolve::stode(y = x0,
                            time=0,
                            func=model_stode_neg_comp,
                            parms=NULL,
                            A=A, 
                            growth.rate=r,
                            positive = TRUE)[[1]]
  }
  model_J = function(t,y,parms=NULL,A,r) {
    dy = y*(r+A%*%y) 
    return(as.list(dy))
  }
  jac.matrix <- jacobian.full(y = out,func = model_J,A = A, r = r)
  return(jac.matrix)
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

plot_eigenvalues <- function(M, prediction = NULL){
  eig <- eigen(M, only.values = TRUE)$values
  dt <- tibble(Real = Re(eig), Imaginary = Im(eig))
  pl <- ggplot(dt) + aes(x = Real, y = Imaginary) + 
    geom_point() + 
    coord_equal() + 
    geom_vline(xintercept = 0, colour = "red", linetype = 2)
  if (is.null(prediction) == FALSE) {
    pl <- pl + geom_vline(xintercept = prediction, colour = "black", linetype = 2)
  }
  show(pl)
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

#' Niche Model Food Web
#'
#' @param S Number of species in the community.
#' @param C The connectance, or fraction of realized links in the food web.
#'
#' @return An adjacency matrix for a niche model food web.
#' @export
#'
#' @section Reference:
#' Williams, R. J., and N. D. Martinez. 2000. Simple rules yield complex food webs. Nature 404:180–183.
#'
#' @examples
#' niche(20, .1)
niche <- function(S, C){
  cond <- FALSE
  while(!cond){
    n.i <- sort(runif(S), decreasing = F)
    r.i <- rbeta(S,1,((1/(2*C))-1))*n.i
    c.i <- runif(S, r.i/2, n.i)
    
    a <- matrix(0, nrow = S, ncol = S)
    
    for(i in 2:S){
      for(j in 1:S){
        if(n.i[j] > (c.i[i] - (.5 * r.i[i])) & n.i[j] < (c.i[i] + .5 * r.i[i])){
          a[j, i] <- 1
        }
      }
    }
    
    cond <- igraph::is.connected(igraph::graph.adjacency(a))
  }
  
  return(a)
}

