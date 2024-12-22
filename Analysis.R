library(matrixcalc)
library(expm)

u = c(10,7,0) # Payout

# Lottery probabilities
lottery <- data.frame(c(0, 0.2,  0,    0.2),
                      c(1, 0.75, 0.25, 0),
                      c(0, 0.05, 0.75, 0.8))
lottery$exp <- c(as.matrix(lottery) %*% u) # expected payout
colnames(lottery) <- c('10','7','0','exp')

# Set the Games' Hamiltonian
H_A <- diag(lottery$exp[1:2])
H_B <- diag(lottery$exp[3:4])
H_AB <- (H_A %x% diag(2) + diag(2) %x% H_B) / 2

# Observed Game Outcomes
# Game Occupancy       Basis                   Cohort
#  A   15, 17          |1>, |2>                Student
#  A    8, 19          |1>, |2>                Trader
#  B    8, 24          |3>, |4>                Student
#  B    3, 24          |3>, |4>                Trader
#  AB   4, 13, 3, 10   |13>, |14>, |23>, |24>  Student
#  AB   8,  7, 9, 30   |13>, |14>, |23>, |24>  Trader
n_A <- data.frame(c(15, 17), c(8, 19))
n_B <- data.frame(c(8, 24), c(3, 24))
n_AB <- data.frame(c(4, 13, 3, 10), c(8,  7, 9, 30))
colnames(n_A) <- c('S','T')
colnames(n_B) <- c('S','T')
colnames(n_AB) <- c('S','T')

# Empirical Density Matrices
rho_A <- list(diag(n_A$S)/sum(n_A$S),diag(n_A$T)/sum(n_A$T))
rho_B <- list(diag(n_B$S)/sum(n_B$S),diag(n_B$T)/sum(n_B$T))
rho_AB <- list(diag(n_AB$S)/sum(n_AB$S),diag(n_AB$T)/sum(n_AB$T))

# The Canonical Economic Distribution
canon_dist <- function(beta, H){
  if (!is.numeric(beta)) stop('beta must be a number.')
  if (beta < 0) stop('beta must be greater than zero.')
  if (!is.square.matrix(H)) stop('The H matrix is not square.')
  return(expm(beta * H) / sum(diag(expm(beta * H))))
}

# Kullback-Leibler Divergence
kl_div <- function(model, data) {
  if (!is.square.matrix(data)) stop('The emperical matrix is not square.')
  if (!is.square.matrix(model)) stop('The modeled matrix is not square.')
  if (is.diagonal.matrix(data) && is.diagonal.matrix(model)) {
    return(sum(diag(data) * log(diag(data) / diag(model))))
  } else {
    return(sum(diag(data * log(data / model))))
  }
}

# R^2
mat_r2 <- function(model, data) {
  if (!is.square.matrix(model)) stop('The model matrix is not square.')
  if (!is.square.matrix(data)) stop('The data matrix is not square.')
  return(cor(diag(model), diag(data))^2)
}

mat_chisq <- function(model, data) {
  if (!is.square.matrix(model)) stop('The model matrix is not square.')
  if (!is.square.matrix(data)) stop('The data matrix is not square.')
  model <- diag(model)
  data <- diag(data)
  return(suppressWarnings(chisq.test(x = data,
                                     p = model,
                                     simulate.p.value = TRUE)))
}

# Looking at the data we estimate two different strategies for each of the
# cohorts.
#
# The students take a naive approach where they treat the combined
# game as an independent composition of the two sub games.
#
# The traders look at the combined game as a unitary whole.
#
# Thus it is clear that the Allais paradox isn't a paradox but is instead an
# issue of experience at evaluating and analyzing risk. In both cases, utility
# was always canonically distributed. What differed was how experience led to
# different evaluation of the risk in the system.

# Objective function for estimating beta
obj_func <- function(beta, H, rho) {
  if (!is.square.matrix(rho)) stop('The rho matrix is not square.')
  return(sum(diag((rho - canon_dist(beta, H))^2)))
}
min_obj <- list(optimize(obj_func, c(0, 10), rho = rho_A[[1]], H = H_A),
                optimize(obj_func, c(0, 10), rho = rho_B[[1]], H = H_B),
                optimize(obj_func, c(0, 10), rho = rho_A[[2]], H = H_A),
                optimize(obj_func, c(0, 10), rho = rho_B[[2]], H = H_B),
                optimize(obj_func, c(0, 10), rho = rho_AB[[2]], H = H_AB))
names(min_obj) <- c('SA', 'SB', 'TA', 'TB', 'TAB')

kl <- list(kl_div(rho_A[[1]] %x% rho_B[[1]], rho_AB[[1]]),
           kl_div(canon_dist(min_obj$SA$minimum, H_A) %x%
                          canon_dist(min_obj$SB$minimum, H_B), rho_AB[[1]]),
           kl_div(canon_dist(min_obj$TAB$minimum, H_AB), rho_AB[[2]]))
r2 <- list(mat_r2(rho_A[[1]] %x% rho_B[[1]], rho_AB[[1]]),
           mat_r2(canon_dist(min_obj$SA$minimum, H_A) %x%
                    canon_dist(min_obj$SB$minimum, H_B), rho_AB[[1]]),
           mat_r2(canon_dist(min_obj$TAB$minimum, H_AB), rho_AB[[2]]))
# chisq <- list(mat_chisq(rho_A[[1]] %x% rho_B[[1]], rho_AB[[1]]),
#               mat_chisq(canon_dist(min_obj$SA$minimum, H_A) %x%
#                        canon_dist(min_obj$SB$minimum, H_B), rho_AB[[1]]),
#               mat_chisq(canon_dist(min_obj$TB$minimum, H_AB), diag(4)*n_AB[[2]]))
chisq <- list(mat_chisq(rho_A[[1]] %x% rho_B[[1]], diag(4)*n_AB[[1]]),
              mat_chisq(canon_dist(min_obj$SA$minimum, H_A) %x%
                          canon_dist(min_obj$SB$minimum, H_B), diag(4)*n_AB[[1]]),
              mat_chisq(canon_dist(min_obj$TAB$minimum, H_AB), diag(4)*n_AB[[2]]))
# chisq <- list(mat_chisq(rho_A[[1]] %x% rho_B[[1]], diag(4)*n_AB[[1]]),
#               mat_chisq(canon_dist(min_obj$SA$minimum, H_A) %x%
#                           canon_dist(min_obj$SB$minimum, H_B), diag(4)*n_AB[[1]]),
              # mat_chisq(canon_dist(min_obj$TAB$minimum, H_AB), diag(4)*c(4,10,10,30)))
names(kl) <- c('SE', 'S', 'T')
names(r2) <- c('SE', 'S', 'T')
names(chisq) <- c('SE', 'S', 'T')
