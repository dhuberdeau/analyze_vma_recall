library(aod)
library(R.matlab)

setwd(getwd())

### DEFINE FIT FUNCTION ###
fit_sigmoid <- function(dat){
  dat <- dat[which(dat$quarter > 0),]
  pt <- dat$pt
  y <- dat$pc
  ty <- dat$type
  qt <- dat$quarter
  tr <- dat$trials
  ty[ty == 4] <- 3
  ty[ty == 3] <- 3
  ty <- as.factor(ty)
  qt <- as.factor(qt)
  tr <- as.numeric(tr)

  fit <- glm(y ~ pt*ty*tr, data = data.frame(pt,ty,y,tr), family = binomial(link = logit))
  # fit <- glm(y ~ pt + ty + tr + tr:ty, data = data.frame(pt,ty,y,tr), family = binomial(link = logit))

  # for (dqrt in c(1,2,3,4)){sig_out <- plot_sigmoid_by_quarter(y, pt, ty, tr, qt, dqrt, fit)}
  # writeMat(
  #   paste(c("/Users/david/OneDrive/Documents/Yale/NTB_lab/batters_problem/vma_recall_repo/sigmoid_fit_E2_",
  #           toString(dqrt), ".mat"), sep="", collapse=""),
  #   sig_out = sig_out)

  return(list(fit, data.frame(y, pt, ty, tr, qt)))
}

fit_sigmoid_0 <- function(dat){

  dat <- dat[which(dat$quarter == 0),]
  pt <- dat$pt
  y <- dat$pc

  fit <- glm(y ~ pt, data = data.frame(pt,y), family = binomial(link = logit))
  return(list(fit, data.frame(y, pt)))
}

### DEFINE PLOT FUNCTION ###
plot_sigmoid_by_quarter <- function(pc, pt, ty, tr, qt, dqrt, fit){
  PT_MIN_MAX <- .6
  TRIAL_QUARTERS <- 60

  pc <- pc[qt == dqrt]
  pt <- pt[qt == dqrt]
  ty <- ty[qt == dqrt]
  tr <- tr[qt == dqrt]

  dat <- data.frame(pc, pt, ty)

  t <- seq(-PT_MIN_MAX, PT_MIN_MAX , .05)
  tr_ <- seq(1 + TRIAL_QUARTERS*(dqrt-1),1 + TRIAL_QUARTERS*(dqrt-1),length.out = length(t))
  tr_ <- as.numeric(tr_)

  unq_types <- unique(ty)

  with(dat[ty == 1,], plot(pc ~ pt, xlim = c(-PT_MIN_MAX, PT_MIN_MAX), col = "green"))
  p <- seq(1,1,length.out = length(t))
  p <- as.factor(p)
  y_1 <- predict(fit, newdata = data.frame(pt=t,ty=p,tr=tr_), type="response")
  lines(t, y_1, col = "green")

  par(new = TRUE)
  with(dat[ty == 2,], plot(pc ~ pt, xlim = c(-PT_MIN_MAX, PT_MIN_MAX), axes=FALSE, col = "blue"))
  p <- seq(2,2,length.out = length(t))
  p <- as.factor(p)
  y_2 <- predict(fit, newdata = data.frame(pt=t,ty=p,tr=tr_), type="response")
  lines(t, y_2, col = "blue")

  par(new = TRUE)
  with(dat[ty == 3,], plot(pc ~ pt, xlim = c(-PT_MIN_MAX, PT_MIN_MAX), axes=FALSE, col= "black"))
  p <- seq(3,3,length.out = length(t))
  p <- as.factor(p)
  y_3 <- predict(fit, newdata = data.frame(pt=t,ty=p,tr=tr_), type="response")
  lines(t, y_3, col = "black")

  return(data.frame(t, y_1, y_2, y_3))
}

### DEFINE PLOT FUNCTION ###
plot_sigmoid_0 <- function(pc, pt, fit){
  PT_MIN_MAX <- .6

  dat <- data.frame(pc, pt)

  t <- seq(-PT_MIN_MAX, PT_MIN_MAX , .05)

    with(dat, plot(pc ~ pt, xlim = c(-PT_MIN_MAX, PT_MIN_MAX), col = "red"))
    y_0 <- predict(fit, newdata = data.frame(pt=t), type="response")
    lines(t, y_0, col = "red")

  return(data.frame(t, y_0))
}

### LOAD AND ORGANIZE DATA ###

### 3-target variant ###
pc_data_quarters <- read.csv('raw_data_mat_E2_pc_quarters1',
                    header = FALSE, na.strings = 'NaN')
pc_data_trials <- read.csv('raw_data_mat_E2_pc_trials1',
                             header = FALSE, na.strings = 'NaN')

# head(pc_data_trials)
# head(pc_data_quarters)
# V1=correct or not, V2=trial type, V3=PT, V4=subject, V5=quarters/trials

dat_q3 <- data.frame(
    pc = pc_data_trials$V1,
    type = pc_data_trials$V2,
    pt = pc_data_trials$V3,
    subject = pc_data_trials$V4,
    quarter = pc_data_quarters$V5,
    trials = pc_data_trials$V5)

dat_q_3 <- na.exclude(dat_q3)

### 4-target variant ###
pc_data_quarters <- read.csv('raw_data_mat_E2_pc_quarters2',
                             header = FALSE, na.strings = 'NaN')
pc_data_trials <- read.csv('raw_data_mat_E2_pc_trials2',
                           header = FALSE, na.strings = 'NaN')

# head(pc_data_trials)
# head(pc_data_quarters)
# V1=correct or not, V2=trial type, V3=PT, V4=subject, V5=quarters/trials

dat_q4 <- data.frame(
  pc = pc_data_trials$V1,
  type = pc_data_trials$V2,
  pt = pc_data_trials$V3,
  subject = pc_data_trials$V4,
  quarter = pc_data_quarters$V5,
  trials = pc_data_trials$V5)

dat_q_4 <- na.exclude(dat_q4)

### 6-target variant ###
pc_data_quarters <- read.csv('raw_data_mat_E2_pc_quarters3',
                             header = FALSE, na.strings = 'NaN')
pc_data_trials <- read.csv('raw_data_mat_E2_pc_trials3',
                           header = FALSE, na.strings = 'NaN')

# head(pc_data_trials)
# head(pc_data_quarters)
# V1=correct or not, V2=trial type, V3=PT, V4=subject, V5=quarters/trials

dat_q6 <- data.frame(
  pc = pc_data_trials$V1,
  type = pc_data_trials$V2,
  pt = pc_data_trials$V3,
  subject = pc_data_trials$V4,
  quarter = pc_data_quarters$V5,
  trials = pc_data_trials$V5)

dat_q_6 <- na.exclude(dat_q6)

### 3-target-replication variant ###
pc_data_quarters <- read.csv('raw_data_mat_E2_pc_quarters4',
                             header = FALSE, na.strings = 'NaN')
pc_data_trials <- read.csv('raw_data_mat_E2_pc_trials4',
                           header = FALSE, na.strings = 'NaN')

# head(pc_data_trials)
# head(pc_data_quarters)
# V1=correct or not, V2=trial type, V3=PT, V4=subject, V5=quarters/trials

dat_q3b <- data.frame(
  pc = pc_data_trials$V1,
  type = pc_data_trials$V2,
  pt = pc_data_trials$V3,
  subject = pc_data_trials$V4,
  quarter = pc_data_quarters$V5,
  trials = pc_data_trials$V5)

dat_q_3b <- na.exclude(dat_q3b)

### --- Fit No Pre-Cue trials --- ###

fit_3_out <- fit_sigmoid(dat_q_3)
fit_3 <- fit_3_out[[1]]
fit_data_3 <- fit_3_out[[2]]
for (dqrt in c(1,2,3,4)){
  sig_out <- plot_sigmoid_by_quarter(fit_data_3$y,
                                      fit_data_3$pt,
                                      fit_data_3$ty,
                                      fit_data_3$tr,
                                      fit_data_3$qt,
                                      dqrt, fit_3)
  writeMat(
      paste(c("sigmoid_fit_E2_",
              "3T_", toString(dqrt), ".mat"), sep="", collapse=""), sig_out = sig_out)
  }
fit_3_out_0 <- fit_sigmoid_0(dat_q_3)
fit_3_0 <- fit_3_out_0[[1]]
fit_data_3_0 <- fit_3_out_0[[2]]
sig_out <- plot_sigmoid_0(fit_data_3_0$y, fit_data_3_0$pt, fit_3_0)
writeMat(
    paste(c("sigmoid_fit_E2_",
            "3T_0.mat"), sep="", collapse=""), sig_out = sig_out)
summary(fit_3)

fit_4_out <- fit_sigmoid(dat_q_4)
fit_4 <- fit_4_out[[1]]
fit_data_4 <- fit_4_out[[2]]
for (dqrt in c(1,2,3,4)){sig_out <- plot_sigmoid_by_quarter(fit_data_4$y,
                                                            fit_data_4$pt,
                                                            fit_data_4$ty,
                                                            fit_data_4$tr,
                                                            fit_data_4$qt,
                                                            dqrt, fit_4)
    writeMat(
      paste(c("sigmoid_fit_E2_",
              "4T_", toString(dqrt), ".mat"), sep="", collapse=""), sig_out = sig_out)
  }
fit_4_out_0 <- fit_sigmoid_0(dat_q_4)
fit_4_0 <- fit_4_out_0[[1]]
fit_data_4_0 <- fit_4_out_0[[2]]
sig_out <- plot_sigmoid_0(fit_data_4_0$y, fit_data_4_0$pt, fit_4_0)
writeMat(
    paste(c("sigmoid_fit_E2_",
            "4T_0.mat"), sep="", collapse=""), sig_out = sig_out)
summary(fit_4)

fit_6_out <- fit_sigmoid(dat_q_6)
fit_6 <- fit_6_out[[1]]
fit_data_6 <- fit_6_out[[2]]
for (dqrt in c(1,2,3,4)){sig_out <- plot_sigmoid_by_quarter(fit_data_6$y,
                                                            fit_data_6$pt,
                                                            fit_data_6$ty,
                                                            fit_data_6$tr,
                                                            fit_data_6$qt,
                                                            dqrt, fit_6)

      writeMat(
        paste(c("sigmoid_fit_E2_",
                "6T_", toString(dqrt), ".mat"), sep="", collapse=""), sig_out = sig_out)
    }
fit_6_out_0 <- fit_sigmoid_0(dat_q_6)
fit_6_0 <- fit_6_out_0[[1]]
fit_data_6_0 <- fit_6_out_0[[2]]
sig_out <- plot_sigmoid_0(fit_data_6_0$y, fit_data_6_0$pt, fit_6_0)
writeMat(
    paste(c("sigmoid_fit_E2_",
            "6T_0.mat"), sep="", collapse=""), sig_out = sig_out)
summary(fit_6)

fit_3b_out <- fit_sigmoid(dat_q_3b)
fit_3b <- fit_3b_out[[1]]
fit_data_3b <- fit_3b_out[[2]]
for (dqrt in c(1,2,3,4)){sig_out <- plot_sigmoid_by_quarter(fit_data_3b$y,
                                                            fit_data_3b$pt,
                                                            fit_data_3b$ty,
                                                            fit_data_3b$tr,
                                                            fit_data_3b$qt,
                                                            dqrt, fit_3b)
  writeMat(
    paste(c("sigmoid_fit_E2_",
            "3Tb_", toString(dqrt), ".mat"), sep="", collapse=""), sig_out = sig_out)
}
fit_3b_out_0 <- fit_sigmoid_0(dat_q_3b)
fit_3b_0 <- fit_3b_out_0[[1]]
fit_data_3b_0 <- fit_3b_out_0[[2]]
sig_out <- plot_sigmoid_0(fit_data_3b_0$y, fit_data_3b_0$pt, fit_3b_0)
writeMat(
    paste(c("sigmoid_fit_E2_",
            "3Tb_0.mat"), sep="", collapse=""), sig_out = sig_out)
summary(fit_3b)
