setwd("/Users/paulzuo/Documents/Penn2018-2019/STAT 590/")
nhanesi_df <- read.csv("smoking_nhanesi_data.csv")
summary(nhanesi_df)
nhanesi_df <- subset(nhanesi_df, select = -c(diabetes.ever, heart.disease, malignant.tumor.present, malignant.tumor.past, polio.paralysis, fracture.of.hip, fracture.of.spine))
nhanesi_df <- subset(nhanesi_df, select = -c(X.1, X))
nhanesi_df <- subset(nhanesi_df, select = -c(SEQN))
nhanesi_df$age_lived_since_1971 <- nhanesi_df$yr.death.expanded - 71
nhanesi_df <- subset(nhanesi_df, select = -c(yr.death.expanded))
