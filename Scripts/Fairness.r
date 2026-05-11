# =============================================================================
# Fairness of Machine Learning — TCGA LUAD
# =============================================================================
# Description : Evaluates fairness of logistic regression models on TCGA
#               lung data across gender and country attributes using
#               demographic parity, accuracy parity, and equalized odds.
# Author      : S.Limra Salith
# Date        : 2025-05-5
# =============================================================================


# -----------------------------------------------------------------------------
# 1. Load Libraries
# -----------------------------------------------------------------------------

library(dplyr)
library(fairness)
library(caret)
library(smotefamily)


# -----------------------------------------------------------------------------
# 2. Load Data
# -----------------------------------------------------------------------------

setwd("F:/TCGA_obj3/lung/Fair")

df_fairness <- read.csv("7_final.csv")
message("Dataset dimensions: ", paste(dim(df_fairness), collapse = " x "))


# -----------------------------------------------------------------------------
# 3. Train/Test Split — Original Data (Gender Attribute)
# -----------------------------------------------------------------------------

set.seed(123)
train_index <- createDataPartition(df_fairness$tissue_type, p = 0.7, list = FALSE)

# Gender split (excludes country column)
train_main <- df_fairness[train_index, ]  %>% dplyr::select(-country_of_residence_at_enrollment)
test_main  <- df_fairness[-train_index, ] %>% dplyr::select(-country_of_residence_at_enrollment)

# Country split (excludes gender column)
train_main1 <- df_fairness[train_index, ]  %>% dplyr::select(-gender)
test_main1  <- df_fairness[-train_index, ] %>% dplyr::select(-gender)

message("--- Original Data: Gender Split ---")
message("Train — Class:  "); print(table(train_main$tissue_type))
message("Train — Gender: "); print(table(train_main$gender))
message("Test  — Class:  "); print(table(test_main$tissue_type))
message("Test  — Gender: "); print(table(test_main$gender))


# -----------------------------------------------------------------------------
# 4. SMOTE — Balance Class Distribution
# -----------------------------------------------------------------------------

# Calculate duplication size needed to balance minority class
class_counts    <- table(df_fairness$tissue_type)
minority_class  <- names(class_counts)[which.min(class_counts)]
majority_class  <- names(class_counts)[which.max(class_counts)]
dup_size_val    <- (class_counts[majority_class] / class_counts[minority_class]) - 1

# Encode country as numeric for SMOTE
df_fair <- df_fairness %>%
  mutate(country = case_when(
    country_of_residence_at_enrollment == "Australia"     ~ 1,
    country_of_residence_at_enrollment == "Canada"        ~ 2,
    country_of_residence_at_enrollment == "Germany"       ~ 3,
    country_of_residence_at_enrollment == "Romania"       ~ 4,
    country_of_residence_at_enrollment == "Russia"        ~ 5,
    country_of_residence_at_enrollment == "Ukraine"       ~ 6,
    country_of_residence_at_enrollment == "United States" ~ 7,
    country_of_residence_at_enrollment == "Vietnam"       ~ 8,
    TRUE ~ 0
  )) %>%
  dplyr::select(-country_of_residence_at_enrollment)

# Encode gender as numeric for SMOTE
df_fair <- df_fair %>%
  mutate(genders = case_when(
    gender == "male"   ~ 1,
    gender == "female" ~ 0,
    TRUE               ~ 2
  )) %>%
  dplyr::select(-gender)

# Apply SMOTE
set.seed(123)
predictors       <- df_fair[, !(names(df_fair) %in% c("tissue_type"))]
predictors_dummy <- data.frame(model.matrix(~ . - 1, data = predictors))

smote_result <- SMOTE(
  X        = predictors_dummy,
  target   = df_fair$tissue_type,
  K        = 5,
  dup_size = dup_size_val
)

df_smote <- smote_result$data

# Decode gender back to labels
df_smote <- df_smote %>%
  mutate(gender = ifelse(genders >= 0.5, "male", "female")) %>%
  dplyr::select(-genders)

# Decode country back to labels
df_smote <- df_smote %>%
  mutate(countries = case_when(
    country == 1 ~ "Australia",
    country == 2 ~ "Canada",
    country == 3 ~ "Germany",
    country == 4 ~ "Romania",
    country == 5 ~ "Russia",
    country == 6 ~ "Ukraine",
    country == 7 ~ "UnitedStates",
    country == 8 ~ "Vietnam",
    TRUE         ~ "0"
  )) %>%
  dplyr::select(-country)

message("--- SMOTE Data: Gender Distribution ---")
print(table(df_smote$gender))


# -----------------------------------------------------------------------------
# 5. Train/Test Split — SMOTE Data (Gender Attribute)
# -----------------------------------------------------------------------------

set.seed(123)
train_index <- createDataPartition(df_smote$class, p = 0.7, list = FALSE)

# Gender split
train_bal <- df_smote[train_index, ]  %>% dplyr::select(-countries)
test_bal  <- df_smote[-train_index, ] %>% dplyr::select(-countries)

# Country split
train_bal1 <- df_smote[train_index, ]  %>% dplyr::select(-gender)
test_bal1  <- df_smote[-train_index, ] %>% dplyr::select(-gender)

message("--- SMOTE Data: Gender Split ---")
print(table(train_bal$gender)); print(table(train_bal$class))
print(table(test_bal$gender));  print(table(test_bal$class))


# -----------------------------------------------------------------------------
# 6. Equal Gender + Equal Class Sampling
# -----------------------------------------------------------------------------

# Find minimum count across all gender-class combinations
min_value <- min(
  sum(df_smote$gender == "male"   & df_smote$class == "Normal"),
  sum(df_smote$gender == "male"   & df_smote$class == "Tumor"),
  sum(df_smote$gender == "female" & df_smote$class == "Normal"),
  sum(df_smote$gender == "female" & df_smote$class == "Tumor")
)

# Sample equal counts per group
data_alt <- df_smote %>%
  group_by(gender, class) %>%
  sample_n(min_value) %>%
  ungroup()

set.seed(123)
train_index <- createDataPartition(data_alt$class, p = 0.7, list = FALSE)
train_alt   <- data_alt[train_index, ]
test_alt    <- data_alt[-train_index, ]

message("--- Equal Gender + Class Split ---")
print(table(train_alt$gender)); print(table(train_alt$class))
print(table(test_alt$gender));  print(table(test_alt$class))


# =============================================================================
# SECTION A: GENDER FAIRNESS
# =============================================================================

# -----------------------------------------------------------------------------
# 7. Train Logistic Regression Models (Gender)
# -----------------------------------------------------------------------------

# Model 1: Original data
train_main$tissue_type <- as.factor(train_main$tissue_type)
test_main$tissue_type  <- as.factor(test_main$tissue_type)

model1 <- glm(tissue_type ~ .,
              data    = train_main,
              family  = binomial(link = "logit"),
              control = glm.control(maxit = 100)
)
test_main$prob1 <- predict(model1, test_main, type = "response")

# Model 2: SMOTE class-balanced
train_bal$class <- as.factor(train_bal$class)
test_bal$class  <- as.factor(test_bal$class)

model2 <- glm(class ~ .,
              data    = train_bal,
              family  = binomial(link = "logit"),
              control = glm.control(maxit = 100)
)
test_bal$prob2 <- predict(model2, test_bal, type = "response")

# Model 3: SMOTE class + gender balanced
train_alt$class <- as.factor(train_alt$class)
test_alt$class  <- as.factor(test_alt$class)

model3 <- glm(class ~ .,
              data    = train_alt,
              family  = binomial(link = "logit"),
              control = glm.control(maxit = 100)
)
test_alt$prob3 <- predict(model3, test_alt, type = "response")


# -----------------------------------------------------------------------------
# 8. Fairness Evaluation — Demographic Parity (Gender)
# -----------------------------------------------------------------------------

res1_dem <- dem_parity(data = test_main, outcome = "tissue_type",
                       outcome_base = "Tumor", group = "gender",
                       probs = "prob1", cutoff = 0.5, base = "female")

res2_dem <- dem_parity(data = test_bal, outcome = "class",
                       outcome_base = "Tumor", group = "gender",
                       probs = "prob2", cutoff = 0.5, base = "female")

res3_dem <- dem_parity(data = test_alt, outcome = "class",
                       outcome_base = "Tumor", group = "gender",
                       probs = "prob3", cutoff = 0.5, base = "female")

message("Demographic Parity — Original:          "); print(res1_dem$Metric)
message("Demographic Parity — SMOTE Class:       "); print(res2_dem$Metric)
message("Demographic Parity — SMOTE Class+Gender:"); print(res3_dem$Metric)

res1_dem$Metric_plot
res2_dem$Metric_plot
res3_dem$Metric_plot


# -----------------------------------------------------------------------------
# 9. Fairness Evaluation — Accuracy Parity (Gender)
# -----------------------------------------------------------------------------

res1_acc <- acc_parity(data = test_main, outcome = "tissue_type",
                       outcome_base = "Tumor", group = "gender",
                       probs = "prob1", cutoff = 0.5, base = "female")

res2_acc <- acc_parity(data = test_bal, outcome = "class",
                       outcome_base = "Tumor", group = "gender",
                       probs = "prob2", cutoff = 0.5, base = "female")

res3_acc <- acc_parity(data = test_alt, outcome = "class",
                       outcome_base = "Tumor", group = "gender",
                       probs = "prob3", cutoff = 0.5, base = "female")

message("Accuracy Parity — Original:          "); print(res1_acc$Metric)
message("Accuracy Parity — SMOTE Class:       "); print(res2_acc$Metric)
message("Accuracy Parity — SMOTE Class+Gender:"); print(res3_acc$Metric)

res1_acc$Metric_plot
res2_acc$Metric_plot
res3_acc$Metric_plot


# -----------------------------------------------------------------------------
# 10. Fairness Evaluation — Equalized Odds (Gender)
# -----------------------------------------------------------------------------

res1_eq <- equal_odds(data = test_main, outcome = "tissue_type",
                      outcome_base = "Tumor", group = "gender",
                      probs = "prob1", cutoff = 0.5, base = "female")

res2_eq <- equal_odds(data = test_bal, outcome = "class",
                      outcome_base = "Tumor", group = "gender",
                      probs = "prob2", cutoff = 0.5, base = "female")

res3_eq <- equal_odds(data = test_alt, outcome = "class",
                      outcome_base = "Tumor", group = "gender",
                      probs = "prob3", cutoff = 0.5, base = "female")

message("Equalized Odds — Original:          "); print(res1_eq$Metric)
message("Equalized Odds — SMOTE Class:       "); print(res2_eq$Metric)
message("Equalized Odds — SMOTE Class+Gender:"); print(res3_eq$Metric)

res1_eq$Metric_plot
res2_eq$Metric_plot
res3_eq$Metric_plot


# =============================================================================
# SECTION B: COUNTRY FAIRNESS
# =============================================================================

# -----------------------------------------------------------------------------
# 11. Train Logistic Regression Models (Country)
# -----------------------------------------------------------------------------

# Model 4: Original data (no gender column)
train_main1$tissue_type <- as.factor(train_main1$tissue_type)
test_main1$tissue_type  <- as.factor(test_main1$tissue_type)

model4 <- glm(tissue_type ~ .,
              data    = train_main1,
              family  = binomial(link = "logit"),
              control = glm.control(maxit = 100)
)
test_main1$prob1 <- predict(model4, test_main1, type = "response")

# Model 5: SMOTE class-balanced (no gender column)
train_bal1$class <- as.factor(train_bal1$class)
test_bal1$class  <- as.factor(test_bal1$class)

model5 <- glm(class ~ .,
              data    = train_bal1,
              family  = binomial(link = "logit"),
              control = glm.control(maxit = 100)
)
test_bal1$prob2 <- predict(model5, test_bal1, type = "response")


# -----------------------------------------------------------------------------
# 12. Fairness Evaluation — Demographic Parity (Country)
# -----------------------------------------------------------------------------

res1_dem_c <- dem_parity(data = test_main1, outcome = "tissue_type",
                         outcome_base = "Tumor", group = "country_of_residence_at_enrollment",
                         probs = "prob1", cutoff = 0.5, base = "United States")

res2_dem_c <- dem_parity(data = test_bal1, outcome = "class",
                         outcome_base = "Tumor", group = "countries",
                         probs = "prob2", cutoff = 0.5, base = "UnitedStates")

message("Demographic Parity (Country) — Original:    "); print(res1_dem_c$Metric)
message("Demographic Parity (Country) — SMOTE Class: "); print(res2_dem_c$Metric)

res1_dem_c$Metric_plot
res2_dem_c$Metric_plot


# -----------------------------------------------------------------------------
# 13. Fairness Evaluation — Accuracy Parity (Country)
# -----------------------------------------------------------------------------

res1_acc_c <- acc_parity(data = test_main1, outcome = "tissue_type",
                         outcome_base = "Tumor", group = "country_of_residence_at_enrollment",
                         probs = "prob1", cutoff = 0.5, base = "United States")

res2_acc_c <- acc_parity(data = test_bal1, outcome = "class",
                         outcome_base = "Tumor", group = "countries",
                         probs = "prob2", cutoff = 0.5, base = "UnitedStates")

message("Accuracy Parity (Country) — Original:    "); print(res1_acc_c$Metric)
message("Accuracy Parity (Country) — SMOTE Class: "); print(res2_acc_c$Metric)

res1_acc_c$Metric_plot
res2_acc_c$Metric_plot


# -----------------------------------------------------------------------------
# 14. Fairness Evaluation — Equalized Odds (Country)
# -----------------------------------------------------------------------------

res1_eq_c <- equal_odds(data = test_main1, outcome = "tissue_type",
                        outcome_base = "Tumor", group = "country_of_residence_at_enrollment",
                        probs = "prob1", cutoff = 0.5, base = "United States")

res2_eq_c <- equal_odds(data = test_bal1, outcome = "class",
                        outcome_base = "Tumor", group = "countries",
                        probs = "prob2", cutoff = 0.5, base = "UnitedStates")

message("Equalized Odds (Country) — Original:    "); print(res1_eq_c$Metric)
message("Equalized Odds (Country) — SMOTE Class: "); print(res2_eq_c$Metric)

res1_eq_c$Metric_plot
res2_eq_c$Metric_plot

# =============================================================================
# End of Script
# =============================================================================