#library(boot)
library(effects)
library(stratification)
library(lme4)
library(car)
library(lmerTest)

setwd(getwd())

##### --- Functions --- #####
lme_func <- function(dat){
  # Use mixed effects to account for subject repeated measures:
  #
  # INPUT: dat - data frame with fields pt, type, response, and subject.
  # OUTPUT: list of anova results for interaction, pt, and type tests, and full model

  lm1_ <- lmer(response ~ type + pt + (1|subject), data = dat, REML = FALSE)
  lm2_ <- lmer(response ~ type*pt + (1|subject), data = dat, REML = FALSE) #test the interaction
  lm3_ <- lmer(response ~ type + (1|subject), data = dat, REML = FALSE) # test pt
  lm4_ <- lmer(response ~ pt + (1|subject), data = dat, REML = FALSE) # test type

  av_intx <- anova(lm1_, lm2_) # LRT (liklihood ratio test)
  av_pt <- anova(lm1_, lm3_)
  av_type <- anova(lm1_, lm4_)

  plot(allEffects(lm2_))

  # difflsmeans(lm2_, test.effs = "type")

  return(list(av_intx, av_pt, av_type, lm2_))
}

plot_func <- function(dat){
  # plot types:
  xmin <- -.1
  xmax <- .7
  ymin <- min(dat$response, na.rm=TRUE)
  ymax <- max(dat$response, na.rm=TRUE)

  with(dat[dat$type == 0,], plot(response ~ pt, main = "No Precue", xlim=c(xmin, xmax), ylim=c(ymin, ymax)))
  with(dat[dat$type == 1,], plot(response ~ pt, main = "Direct", xlim=c(xmin, xmax), ylim=c(ymin, ymax)))
  with(dat[dat$type == 2,], plot(response ~ pt, main = "Symbolic", xlim=c(xmin, xmax), ylim=c(ymin, ymax)))
}

posthoc_pt_func <- function(dat, type_){
  dat_ <- dat[dat$type == type_,]
  lm0 <- lmer(response ~ pt + (1|subject), data = dat_, REML = FALSE)
  lm_ <- lmer(response ~ 1 + (1|subject), data = dat_, REML = FALSE)
  pt_0 <- anova(lm0, lm_)

  return(list(pt_0, lm0))
}

##### --- probability correct --- ######
pc_data <- read.csv('table_E1_pc',
                    header = FALSE, na.strings = 'NaN')
# V1: PT, V2: PC, V3: trial type, V4: subject

# rename variables and put into single data frame
dat <- data.frame(subject = pc_data$V4,
                  response = pc_data$V2,
                  type = pc_data$V3,
                  pt = pc_data$V1)

# make sure variables are interpreted properly as factors or continuous
dat$pt <- as.numeric(dat$pt)
dat$type <- factor(dat$type)
dat$subject <- factor(dat$subject)
dat$response <- as.numeric(dat$response)

# plot the data:
plot_func(dat)

# Test interaction and main effects:
lme_pc = lme_func(dat)
mod_pc = lme_pc[[4]]
#mod_pc <- lmer(response ~ type*pt + (1|subject), data = dat, REML = FALSE) #test the interaction
mcomp_pc <- difflsmeans(mod_pc, test.effs = "type")
post_pc_0 <- posthoc_pt_func(dat,0)
post_pc_1 <- posthoc_pt_func(dat,1)
post_pc_2 <- posthoc_pt_func(dat,2)

##### --- Peak Velocity --- #####
pv_data <- read.csv('table_E1_pv',
                    header = FALSE, na.strings = 'NaN')
# V1: PT, V2: PC, V3: trial type, V4: subject

# rename variables and put into single data frame
dat <- data.frame(subject = pv_data$V4,
                  response = pv_data$V2,
                  type = pv_data$V3,
                  pt = pv_data$V1)

# make sure variables are interpreted properly as factors or continuous
dat$pt <- as.numeric(dat$pt)
dat$type <- factor(dat$type)
dat$subject <- factor(dat$subject)
dat$response <- as.numeric(dat$response)

plot_func(dat)
lme_pv = lme_func(dat)
mod_pv = lme_pv[[4]]
#mod_pv <- lmer(response ~ type*pt + (1|subject), data = dat, REML = FALSE) #test the interaction
mcomp_pv <- difflsmeans(mod_pv, test.effs = "type")
post_pv_0 <- posthoc_pt_func(dat,0)
post_pv_1 <- posthoc_pt_func(dat,1)
post_pv_2 <- posthoc_pt_func(dat,2)

##### --- Variability --- #####
var_data <- read.csv('table_E1_var',
                    header = FALSE, na.strings = 'NaN')
# V1: PT, V2: PC, V3: trial type, V4: subject

# rename variables and put into single data frame
dat <- data.frame(subject = var_data$V4,
                  response = var_data$V2,
                  type = var_data$V3,
                  pt = var_data$V1)

# make sure variables are interpreted properly as factors or continuous
dat$pt <- as.numeric(dat$pt)
dat$type <- factor(dat$type)
dat$subject <- factor(dat$subject)
dat$response <- as.numeric(dat$response)

plot_func(dat)
lme_var = lme_func(dat)
mod_var = lme_var[[4]]
#mod_var <- lmer(response ~ type*pt + (1|subject), data = dat, REML = FALSE) #test the interaction
mcomp_var <- difflsmeans(mod_var, test.effs = "type")
post_var_0 <- posthoc_pt_func(dat,0)
post_var_1 <- posthoc_pt_func(dat,1)
post_var_2 <- posthoc_pt_func(dat,2)


##### --- simulation test --- ######
# Does a simulation of the data look sensible? #

PT = c(-.05, 0, .05, .1, .15, .2, .25, .3, .35, .4, .45, .5, .55, .6, .65, .7)
T1 = c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
betas = fixef(mod_pc)
PC0 = betas[1] + betas[4]*PT
PC1 = betas[1] + betas[2]*T1 + betas[4]*PT + betas[5]*T1*PT
PC2 = betas[1] + betas[3]*T1 + betas[4]*PT + betas[6]*T1*PT
plot(PT, PC0, ylim = c(0, 1))
par(new = TRUE)
plot(PT, PC1, ylim = c(0, 1), axes = FALSE)
par(new = TRUE)
plot(PT, PC2, ylim = c(0, 1), axes = FALSE)
