# perform repeated measures anova on variability data from each group:

library(lme4)

#### DEFINE FUNCTIONS FOR DOING REPEATED MEASURES ANOVAs #######

# define mixed-effects anova for probability of success measure: 
anova_custom_PrSucc <- function(var1) {
  # V1 - preparation time
  # V2 - subject
  # V3 - trial type
  # V4 - probability of making a correct response
  
  baseline <- lmer(V4 ~ 1 + (1|V2), data = var1, REML=FALSE)
  blockModel <- lmer(V4 ~ V1 + (1|V2), data = var1, REML=FALSE)
  typeModel <- lmer(V4 ~ V1 + V3 + (1|V2), data = var1, REML=FALSE)
  fullModel <- lmer(V4 ~ V1 * V3 + (1|V2), data = var1, REML=FALSE)
  
  a_result <- anova(baseline, blockModel, typeModel, fullModel)
  return(a_result) 
}

# define mixed-effects anova for probability of succcess collapsed across individuals:
anova_allSubs_PrSucc <- function(var1) {
  # V1 - preparation time
  # V2 - subject
  # V3 - trial type
  # V4 - was the response correct or not?
  
  baseline <- lmer(V4 ~ 1 + (1|V2), data = var1, REML=FALSE)
  blockModel <- lmer(V4 ~ V1 + (1|V2), data = var1, REML=FALSE)
  typeModel <- lmer(V4 ~ V1 + V3 + (1|V2), data = var1, REML=FALSE)
  fullModel <- lmer(V4 ~ V1 * V3 + (1|V2), data = var1, REML=FALSE)
  
  a_result <- anova(baseline, blockModel, typeModel, fullModel)
  return(a_result) 
  
}

# define mixed-effects anova for probability of correct recall acorss blocks:
anova_Exp2_learning <- function(var1) {
  # V1 - Blocks
  # V2 - subjects
  # V3 - Trial Type
  # V4 - probability of success
  
  baseline <- lmer(V4 ~ 1 + (1|V2), data = var1, REML=FALSE)
  blockModel <- lmer(V4 ~ V1 + (1|V2), data = var1, REML=FALSE)
  typeModel <- lmer(V4 ~ V1 + V3 + (1|V2), data = var1, REML=FALSE)
  fullModel <- lmer(V4 ~ V1 * V3 + (1|V2), data = var1, REML=FALSE)
  
  a_result <- anova(baseline, blockModel, typeModel, fullModel)
  return(a_result) 
  
}

# define mixed effects anova for comparing variability 
# (exploratory analysis not reported in paper)
me_anova_custom <- function(var1){
  # mixed effects model, comparing variability across trial types, without a second factor
  
  baseline <- lmer(V3 ~ 1 + (1|V2), data = var1, REML=FALSE)
  fullModel <- lmer(V3 ~ V1 + (1|V2), data = var1, REML=FALSE)
  
  a_result <- anova(baseline, fullModel)
  
  return(a_result) 
}

# define repeated measure anova function for variability
# (exploratory analysis not reported in paper)
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



## probability of success measure:
ps1 <- read.csv("~/OneDrive/Documents/Yale/NTB_lab/batters_problem/vma_recall_repo/prob_succ_E1", 
                header = FALSE, na.strings = 'NaN')

ps_result_1 <- anova_custom_PrSucc(ps1)

ps1_D_S <- subset(ps1, V3 == 1 | V3 == 2, select = c(V1, V2, V3, V4))

ps_D_S_result_1 <- anova_custom_PrSucc(ps1_D_S)

## analysis collapsing across individuals:
ps_allSubs <- read.csv("~/OneDrive/Documents/Yale/NTB_lab/batters_problem/vma_recall_repo/prob_succ_E1", 
                       header = FALSE, na.strings = 'NaN')

ps_allSubs_results <- anova_allSubs_PrSucc(ps_allSubs)

## analyze probability of success across blocks (Exp2):
ps_e2_1 <- read.csv("~/OneDrive/Documents/Yale/NTB_lab/batters_problem/vma_recall_repo/prob_succ_E2_3T", 
  header = FALSE, na.strings = 'NaN')
ps_e2_2 <- read.csv("~/OneDrive/Documents/Yale/NTB_lab/batters_problem/vma_recall_repo/prob_succ_E2_4T", 
                    header = FALSE, na.strings = 'NaN')
ps_e2_3 <- read.csv("~/OneDrive/Documents/Yale/NTB_lab/batters_problem/vma_recall_repo/prob_succ_E2_6T", 
                    header = FALSE, na.strings = 'NaN')
ps_e2_4 <- read.csv("~/OneDrive/Documents/Yale/NTB_lab/batters_problem/vma_recall_repo/prob_succ_E2_3Tb", 
                    header = FALSE, na.strings = 'NaN')

e2_results_1 <- anova_Exp2_learning(ps_e2_1)
e2_results_2 <- anova_Exp2_learning(ps_e2_2)
e2_results_3 <- anova_Exp2_learning(ps_e2_3)
e2_results_4 <- anova_Exp2_learning(ps_e2_4)

## test for differences between vma recall and memory test recall
mem_e2_1 <- read.csv("~/OneDrive/Documents/Yale/NTB_lab/batters_problem/vma_recall_repo/memTest_score_E2_3T", 
                    header = FALSE, na.strings = 'NaN')
mem_e2_1$V1 <- as.factor(mem_e2_1$V1)

mem_e2_2 <- read.csv("~/OneDrive/Documents/Yale/NTB_lab/batters_problem/vma_recall_repo/memTest_score_E2_4T", 
                     header = FALSE, na.strings = 'NaN')
mem_e2_2$V1 <- as.factor(mem_e2_2$V1)

mem_e2_3 <- read.csv("~/OneDrive/Documents/Yale/NTB_lab/batters_problem/vma_recall_repo/memTest_score_E2_6T", 
                     header = FALSE, na.strings = 'NaN')
mem_e2_3$V1 <- as.factor(mem_e2_3$V1)

mem_e2_4 <- read.csv("~/OneDrive/Documents/Yale/NTB_lab/batters_problem/vma_recall_repo/memTest_score_E2_3Tb_1", 
                     header = FALSE, na.strings = 'NaN')
mem_e2_4$V1 <- as.factor(mem_e2_4$V1)

mem_e2_5 <- read.csv("~/OneDrive/Documents/Yale/NTB_lab/batters_problem/vma_recall_repo/memTest_score_E2_3Tb_2", 
                     header = FALSE, na.strings = 'NaN')
mem_e2_5$V1 <- as.factor(mem_e2_5$V1)

mem_t1 <- t.test(mem_e2_1$V2[mem_e2_1$V1==0], mem_e2_1$V2[mem_e2_1$V1==1], paired = TRUE)
mem_t2 <- t.test(mem_e2_2$V2[mem_e2_2$V1==0], mem_e2_2$V2[mem_e2_2$V1==1], paired = TRUE)
mem_t3 <- t.test(mem_e2_3$V2[mem_e2_3$V1==0], mem_e2_3$V2[mem_e2_3$V1==1], paired = TRUE)
mem_t4 <- t.test(mem_e2_4$V2[mem_e2_4$V1==0], mem_e2_4$V2[mem_e2_4$V1==1], paired = TRUE)
mem_t5 <- t.test(mem_e2_5$V2[mem_e2_5$V1==0], mem_e2_5$V2[mem_e2_5$V1==1], paired = TRUE)


######################################
######################################

######################################
## Variability measure (Experiment 1):
## exploratory analysis (not included in paper)
######################################

var0 <- read.csv("~/OneDrive/Documents/Yale/NTB_lab/batters_problem/vma_recall_repo/variability_E1", 
                 header = FALSE, na.strings = 'NaN')

result_0 <- me_anova_custom(var0)


######################################
## Variability measure (Experiment 2):
## exploratory analysis (not included in paper)
######################################

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