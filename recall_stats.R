# perform repeated measures anova on variability data from each group:

library(lme4)

# define mixed-effects anova function:
anova_custom <- function(var1) {
  #var1$V1 <- as.factor(var1$V1)
  #var1$V3 <- as.factor(var1$V3)
  
  baseline <- lmer(V4 ~ 1 + (1|V2), data = var1, REML=FALSE)
  blockModel <- lmer(V4 ~ V1 + (1|V2), data = var1, REML=FALSE)
  typeModel <- lmer(V4 ~ V1 + V3 + (1|V2), data = var1, REML=FALSE)
  fullModel <- lmer(V4 ~ V1 * V3 + (1|V2), data = var1, REML=FALSE)
  
  # fullModel <- lmer(V4 ~ V1 * V3 + (1|V2), data = var1, REML=FALSE) 
  #baseline <- lme(V4 ~ 1, random = ~1 | V2/V1/V3, data = var1, method = "ML")
  #blockModel <- lme(V4 ~ V1, random = ~1 | V2/V1/V3, data = var1, method = "ML")
  #typeModel <- lme(V4 ~ V1 + V3, random = ~1 | V2/V1/V3, data = var1, method = "ML")
  #fullModel <- lme(V4 ~ V1*V3, random = ~1 | V2/V1/V3, data = var1, method = "ML")
  
  a_result <- anova(baseline, blockModel, typeModel, fullModel)
  return(a_result) 
}

# define repeated measure anova function:
rm_anova_custom <- function(var1) {
  var1$V1 <- as.factor(var1$V1)
  #var1$V3 <- as.factor(var1$V3)
  
  baseline <- lmer(V4 ~ 1 + (1|V2), data = var1, REML=FALSE)
  blockModel <- lmer(V4 ~ V1 + (1|V2), data = var1, REML=FALSE)
  typeModel <- lmer(V4 ~ V1 + V3 + (1|V2), data = var1, REML=FALSE)
  fullModel <- lmer(V4 ~ V1 * V3 + (1|V2), data = var1, REML=FALSE)
  
  # fullModel <- lmer(V4 ~ V1 * V3 + (1|V2), data = var1, REML=FALSE) 
  #baseline <- lme(V4 ~ 1, random = ~1 | V2/V1/V3, data = var1, method = "ML")
  #blockModel <- lme(V4 ~ V1, random = ~1 | V2/V1/V3, data = var1, method = "ML")
  #typeModel <- lme(V4 ~ V1 + V3, random = ~1 | V2/V1/V3, data = var1, method = "ML")
  #fullModel <- lme(V4 ~ V1*V3, random = ~1 | V2/V1/V3, data = var1, method = "ML")
  
  a_result <- anova(baseline, blockModel, typeModel, fullModel)
  return(a_result) 
}

rm_anova_model <- function(var1) {
  
  var1$V1 <- as.factor(var1$V1)
  #var1$V3 <- as.factor(var1$V3)
  
  baseline <- lmer(V4 ~ 1 + (1|V2), data = var1, REML=FALSE)
  blockModel <- lmer(V4 ~ V1 + (1|V2), data = var1, REML=FALSE)
  typeModel <- lmer(V4 ~ V1 + V3 + (1|V2), data = var1, REML=FALSE)
  fullModel <- lmer(V4 ~ V1 * V3 + (1|V2), data = var1, REML=FALSE)
  
  # fullModel <- lmer(V4 ~ V1 * V3 + (1|V2), data = var1, REML=FALSE) 
  #baseline <- lme(V4 ~ 1, random = ~1 | V2/V1/V3, data = var1, method = "ML")
  #blockModel <- lme(V4 ~ V1, random = ~1 | V2/V1/V3, data = var1, method = "ML")
  #typeModel <- lme(V4 ~ V1 + V3, random = ~1 | V2/V1/V3, data = var1, method = "ML")
  #fullModel <- lme(V4 ~ V1*V3, random = ~1 | V2/V1/V3, data = var1, method = "ML")
  
  #a_result <- anova(baseline, blockModel, typeModel, fullModel)
  return(fullModel) 
}

# group 1:
var1 <- read.csv("~/OneDrive/Documents/Yale/NTB_lab/batters_problem/vma_recall_repo/variability_1", 
                 header = FALSE, na.strings = 'NaN')

result_1 <- rm_anova_custom(var1)
model_1 <- rm_anova_model(var1)

# group 2:
var2 <- read.csv("~/OneDrive/Documents/Yale/NTB_lab/batters_problem/vma_recall_repo/variability_2", 
                 header = FALSE, na.strings = 'NaN')
result_2 <- rm_anova_custom(var2)
model_2 <- rm_anova_model(var2)

# group 3:
var3 <- read.csv("~/OneDrive/Documents/Yale/NTB_lab/batters_problem/vma_recall_repo/variability_3", 
                 header = FALSE, na.strings = 'NaN')
result_3 <- rm_anova_custom(var3)
model_3 <- rm_anova_model(var3)

# group 4:
var4 <- read.csv("~/OneDrive/Documents/Yale/NTB_lab/batters_problem/vma_recall_repo/variability_4", 
                 header = FALSE, na.strings = 'NaN')
result_4 <- rm_anova_custom(var4)
model_4 <- rm_anova_model(var4)