### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
### ~~~~~~~~~~~~~~~~~~~~ Variance estimations ATE ATT ATC ~~~~~~~~~~~~~ ###
### ~~~~~~~~~~~~~~~~~~~~ Simulation Study          ~~~~~~~~~~~~~~~~~~~~ ###
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###

####### Summary statistics functions

### by Yi Liu
### Nov 7, 2021

library(dplyr)

PE <- function(data) {
  cols <- c("ATC", "ATT")
  data <- data[, cols] %>% as.matrix()
  apply(data, 2, mean)
}

bias <- function(truth, data) {
  cols <- c("ATC", "ATT")
  data <- data[, cols] %>% as.matrix()
  truth <- matrix(rep(as.matrix(truth[, cols]), nrow(data)), nrow=2) %>% t()
  100*apply(data/truth-1, 2, mean) %>% abs()
}

bias.med <- function(truth, data) {
  cols <- c("ATC", "ATT")
  data <- data[, cols] %>% as.matrix()
  truth <- matrix(rep(as.matrix(truth[, cols]), nrow(data)), nrow=2) %>% t()
  100*apply(data/truth-1, 2, median) %>% abs()
}

RMSE <- function(truth, data) {
  sqmean <- function(n) sqrt(mean(n))
  cols <- c("ATC", "ATT")
  data <- data[, cols] %>% as.matrix()
  truth <- matrix(rep(as.matrix(truth[, cols]), nrow(data)), nrow=2) %>% t()
  apply((data-truth)^2, 2, sqmean)
}

RE <- function(data) {
  cols <- c("ATC", "ATT")
  PE <- data[, cols] %>% as.matrix()
  var.cols <- c("ATC.wbvar", "ATT.wbvar")
  var <- data[, var.cols] %>% as.matrix()
  emp <- apply(PE, 2, var)
  emp/apply(var, 2, mean)
}

RE.med <- function(data) {
  cols <- c("ATC", "ATT")
  PE <- data[, cols] %>% as.matrix()
  var.cols <- c("ATC.wbvar", "ATT.wbvar")
  var <- data[, var.cols] %>% as.matrix()
  emp <- apply(PE, 2, var)
  emp/apply(var, 2, median)
}

CP <- function(truth, data) {
  cols <- c("ATC.ifci", "ATT.ifci")
  data <- data[, cols] %>% as.matrix()
  apply(data, 2, mean)
}
