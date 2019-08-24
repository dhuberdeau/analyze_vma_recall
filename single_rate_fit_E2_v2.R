library(aod)
library(R.matlab)
library(rlist)


setwd("~/OneDrive/Documents/Yale/NTB_lab/batters_problem/vma_recall_repo/")

### DEFINE NONLINEAR FUNCTION ###
f_exp2 <- function(beta, tr){
  ## 1-rate exponent:
  y_hat <- beta[1]*(1 - exp(beta[2]*tr))
  
  ## 2-rate exponent:
  # y_hat <- beta[1]*(1 - exp(beta[2]*tr) - exp(beta[3]*tr))
  
  ## 2-rate / different asymptotes:
  # y_hat <- (1 - beta[1]*exp(beta[2]*tr) - beta[3]*exp(beta[4]*tr))
  
  return(y_hat)
  }

### DEFINE OBJECTIVE FUNCTION ###
f <- function(beta, y, tr){
  # y_hat <- beta[1]*(1 - exp(beta[2]*tr) - exp(beta[3]*tr))
  y_hat <- f_exp2(beta, tr)
  er__ <- (y - y_hat)^2
  er_ <- sum(er__)
  # er <- inf.omit(er_)
  return(er_)
} 

### DEFINE FIT FUNCTION ###
fit_exp <- function(dat, typ){
  dat <- na.exclude(dat)
  y <- dat$pc
  tr <- dat$trials
  ty <- dat$type
  
  ty <- as.factor(ty)
  y <- as.numeric(y)
  tr <- as.numeric(tr)
  
  ## for 3-rate
  # fit <- nlm(f, c(1,-.5, 1, -.8), y = y[ty == typ], tr = tr[ty == typ])
  # fit <- nls(y ~ 1 + beta1*exp(beta2*tr) + beta3*exp(beta4*tr),
  #            data = data.frame(y = y[ty == typ], tr = tr[ty == typ]),
  #            start = list(beta1 = .5, beta2 = -0.5, beta3 = .5, beta4 = -0.8),
  #            nls.control(maxiter = 5000, minFactor = 0.000001))
  
  ## baseline:
  # fit <- nls(y ~ beta1*(1 + exp(beta2*tr)),
  #            data = data.frame(y = y[ty == typ], tr = tr[ty == typ]),
  #            start = list(beta1 = .5, beta2 = -0.5),
  #            nls.control(maxiter = 5000, minFactor = 0.000001))
  
  ## plinear
  fit <- nls(y ~ beta1*(1 - exp(beta2*tr)),
             data = data.frame(y = y[ty == typ], tr = tr[ty == typ]),
             start = list(beta1 = .9, beta2 = -0.5),
             nls.control(warnOnly = TRUE))#, alg = "plinear")
  
  ## for 2-rate
  # fit <- nlm(f, c(1,-.5,-.8), y = y[ty == typ], tr = tr[ty == typ])
  
  ## for 1-rate
  # fit <- nlm(f, c(1,-.5), y = y[ty == typ], tr = tr[ty == typ])
  
  return(fit)
}


### READ IN AND PROCESS THE DATA ###
fit_struct_single = list()
file_nums <- c(1,2,3,4) #exclude 4 b/c of failure to converge in nls
for (i_group in file_nums){
  pc_1_ <- read.csv(paste('trial_occurance_data_', i_group, sep=''),
                               header = FALSE, na.strings = 'NaN')
  # V1 = trial number, V2 = Prob. Correct, V3 = type (should be a factor)
  pc_1_ <- data.frame(pc = pc_1_$V2,
                     trials = pc_1_$V1,
                     type = pc_1_$V3)
  pc_1 <- na.exclude(pc_1_)
  
  ### Fit exponent ###
  fit_2 <- fit_exp(pc_1, 2)
  # fit_1 <- fit_exp(pc_1, 1)
  
  
  ### plot the average raw data ###
  t <- seq(1, 18 , 1)
  pc_mean_1 <- seq(1,18,1)
  pc_mean_2 <- seq(1,18,1)
  for (i_tr in t){
    pc_mean_1[i_tr] <- mean(pc_1$pc[pc_1$trials==i_tr & pc_1$type==1])
    # if (i_tr == 1){
    #   plot(temp_tr, temp_pc, col="green")
    # }else{
    #   par(new = TRUE)
    #   plot(temp_tr, temp_pc, col="green")
    # }
  
    pc_mean_2[i_tr] <- mean(pc_1$pc[pc_1$trials==i_tr & pc_1$type==2])
    # temp_pc <- mean(pc_1$pc[pc_1$trials==i_tr & pc_1$type==2])
    # temp_tr <- seq(i_tr, i_tr, length.out = length(temp_pc))
    # par(new = TRUE)
    # plot(temp_tr, temp_pc, col="blue")
  }
  plot(pc_mean_1~t, col="green", xlim = c(0, 18), ylim = c(0,1))
  par(new = TRUE)
  points(pc_mean_2~t, col="blue")
  
  # with(pc_1[pc_1$type == 1,], plot(pc~trials,col="green"))
  # par(new = TRUE)
  # with(pc_1[pc_1$type == 2,], plot(pc~trials,col="blue"))
  
  
  ### Plot fitted exponents ###
  
  # y_hat_1 <- f_exp2(fit_1$estimate, t)
  s_2 <- summary(fit_2)
  beta_fit <- c(s_2$coefficients[1,1], s_2$coefficients[2,1])
  # y_hat_2 <- f_exp2(fit_2$estimate, t)
  # y_hat_2 <- f_exp2(beta_fit, t)
  y_hat_2 <- beta_fit[1]*(1 - exp(beta_fit[2]*t)) 
  
  mat_out = data.frame(t, y_hat_2)
  writeMat(paste("exp_fit_single_E2_", i_group, ".mat",sep=''), mat_out = mat_out)
  
  fit_struct_single <- list.append(fit_struct_single, fit_2)
  
  # lines(t, y_hat_1, col="green")
  lines(t, y_hat_2, col="blue")
}