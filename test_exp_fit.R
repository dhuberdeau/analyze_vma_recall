setwd("~/OneDrive/Documents/Yale/NTB_lab/batters_problem/vma_recall_repo/")

## run both models: single rate exponent and double rate exponent
source('single_rate_fit_E2_v2.R')
source('double_rate_fit_E2_v2.R')

## test the model types
E1_compare <- anova(fit_struct_single[[1]], fit_struct_double[[1]])
E2_compare <- anova(fit_struct_single[[2]], fit_struct_double[[2]])
E3_compare <- anova(fit_struct_single[[3]], fit_struct_double[[3]])

S1 <- summary(fit_struct_single[[1]])
S2 <- summary(fit_struct_single[[2]])
S3 <- summary(fit_struct_single[[3]])
S4 <- summary(fit_struct_single[[4]])