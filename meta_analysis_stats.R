library(R.matlab)

# Load in the meta analysis data:
pt_data_ <- read.csv('min_pt_data',
                             header = FALSE, na.strings = 'NaN')

pt_data <- na.exclude(pt_data_)

lmfit <- lm(V1 ~ V2, data = pt_data)

summary(lmfit)
# lm_out <- data.frame(lmfit)
# writeMat("/Users/david/OneDrive/Documents/Yale/NTB_lab/batters_problem/vma_recall_repo/meta_analysis_pt_fit",
#           lm_out = lm_out)
