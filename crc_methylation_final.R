## Code: CRC methylation ##

# Install CRAN packages
install.packages(c(
  "dplyr",        # Data manipulation
  "readr",        # Data reading
  "tidyr",        # Data tidying
  "knitr",        # Document generation
  "ggplot2",      # Visualization
  "tableone",     # Clinical summary table
  "caret",        # Data partitioning
  "glmnet",       # Lasso cox
  "survival",     # Survival analysis
  "survminer",    # Kaplan-Meier curves
  "timeROC",      # AUC calculation
  "randomForestSRC"  # Random Survival Forest
))

# Install Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c(
  "impute",        # For imputation
  "ComplexHeatmap", # Heatmap visualization
  "circlize",      # Circular visualization
  "missMethyl",    # Methylation analysis
  "IlluminaHumanMethylation450kanno.ilmn12.hg19", # Methylation annotation
  "enrichplot",    # Pathway visualization
  "DOSE"          # Pathway analysis tools
))

# Load data manipulation and visualization packages
library(dplyr)         # Data manipulation
library(readr)         # Data reading
library(tidyr)         # Data tidying
library(ggplot2)       # Visualization
library(knitr)         # Document generation

# Load statistical and analysis packages
library(tableone)      # Clinical summary table
library(impute)        # Imputation
library(caret)         # Data partitioning
library(glmnet)        # Lasso cox

# Load survival analysis packages
library(survival)      # Survival analysis
library(survminer)     # Kaplan-Meier curves
library(timeROC)       # AUC calculation
library(randomForestSRC) # Random Survival Forest

# Load pathway analysis packages
library(ComplexHeatmap)  # Heatmap visualization
library(circlize)       # Circular visualization
library(missMethyl)     # Methylation analysis
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)  # Methylation annotation
library(enrichplot)     # Pathway visualization
library(DOSE)          # Pathway analysis tools


## Pre-processing

# Read CSV files that have been downloaded from TCGA
coad_clinical <- read.csv("tcga_coad_clinical.csv", stringsAsFactors = FALSE)
coad_methylation <- read.csv("tcga_coad_met.csv", stringsAsFactors = FALSE)
read_clinical <- read.csv("tcga_read_clinical.csv", stringsAsFactors = FALSE)
read_methylation <- read.csv("tcga_read_met.csv", stringsAsFactors = FALSE)

# Function to preprocess clinical data
preprocess_clinical <- function(data) {
  data %>%
    select_if(~!all(is.na(.)))  # Remove columns with all NAs
}

# Preprocess COAD and READ clinical data
coad_clinical_processed <- preprocess_clinical(coad_clinical)
read_clinical_processed <- preprocess_clinical(read_clinical)

# Function to preprocess methylation data
preprocess_methylation <- function(data) {
  # Transpose the data
  transposed_data <- t(data)
  
  # Convert to data frame and handle row names
  processed_data <- as.data.frame(transposed_data[-1, ]) # Remove the first row (original column names)
  colnames(processed_data) <- data[[1]] # Set column names from the first column of original data
  
  # Add bcr_patient_barcode column and format it correctly
  processed_data$bcr_patient_barcode <- gsub("\\.", "-", substr(rownames(processed_data), 1, 12))
  
  # Move bcr_patient_barcode to the front
  processed_data <- processed_data[, c("bcr_patient_barcode", setdiff(names(processed_data), c("bcr_patient_barcode")))]

  # Remove columns with all NAs
  processed_data <- processed_data %>% select_if(~!all(is.na(.)))
  
  # Convert methylation data to numeric
  methylation_columns <- setdiff(names(processed_data), "bcr_patient_barcode")
  processed_data[methylation_columns] <- lapply(processed_data[methylation_columns], as.numeric)
  
  return(processed_data)
}

# Preprocess COAD and READ methylation data
coad_methylation_processed <- preprocess_methylation(coad_methylation)
read_methylation_processed <- preprocess_methylation(read_methylation)

# Merge COAD clinical and methylation data
coad_merged <- merge(coad_clinical_processed, coad_methylation_processed, 
                     by = "bcr_patient_barcode", all = FALSE)

# Merge READ clinical and methylation data
read_merged <- merge(read_clinical_processed, read_methylation_processed, 
                     by = "bcr_patient_barcode", all = FALSE)

# Ensure that both datasets have the same columns
common_columns <- intersect(colnames(coad_merged), colnames(read_merged))

# Subset both datasets to include only the common columns
coad_subset <- coad_merged[, common_columns]
read_subset <- read_merged[, common_columns]

# Combine COAD and READ datasets
master_dataframe <- rbind(coad_subset, read_subset)

# Create the working dataframe
working_dataframe <- master_dataframe %>%
  select(
    bcr_patient_barcode,  # Keep the identifier
    gender,
    age_at_diagnosis,
    age_at_index,         # Added age_at_index
    ajcc_pathologic_stage,
    ajcc_pathologic_m,    # Metastasis column
    tissue_or_organ_of_origin,
    vital_status,         # Survival column
    days_to_death,
    days_to_last_follow_up,
    starts_with("cg")     # Only methylation data columns
  )

# Create the new columns
working_dataframe_cleaned <- working_dataframe %>%
  mutate(
    censored_status = ifelse(vital_status == "Dead", 1, 0),
    survival_time = ifelse(vital_status == "Dead", days_to_death, days_to_last_follow_up),
    survival_years = survival_time / 365.25,
    five_year_survival = ifelse(survival_time > 5*365, 1, 0),
    three_year_survival = ifelse(survival_time > 3*365, 1, 0)
  ) %>%
  filter(!is.na(survival_time))

# Reorder columns
methylation_cols <- grep("^cg", names(working_dataframe_cleaned), value = TRUE)

working_dataframe_cleaned <- working_dataframe_cleaned %>%
  select(
    bcr_patient_barcode,
    gender,
    age_at_diagnosis,
    age_at_index,
    ajcc_pathologic_stage,
    ajcc_pathologic_m,
    tissue_or_organ_of_origin,
    vital_status,
    days_to_death,
    days_to_last_follow_up,
    censored_status,
    survival_time,
    survival_years,
    five_year_survival,
    three_year_survival,
    all_of(methylation_cols)
  )

# Get summary statistics
# First, create labeled factors for these variables
working_dataframe_cleaned <- working_dataframe_cleaned %>%
  mutate(
    five_year_survival_factor = factor(five_year_survival, 
                                       levels = c(0, 1), 
                                       labels = c("Less than 5 years", "More than 5 years")),
    three_year_survival_factor = factor(three_year_survival, 
                                        levels = c(0, 1), 
                                        labels = c("Less than 3 years", "More than 3 years")),
    
    # Add age category
    age_category = factor(case_when(
      age_at_diagnosis/365.25 <= 65 ~ "<=65",
      age_at_diagnosis/365.25 > 65 ~ ">65"
    ), levels = c("<=65", ">65"))
  )

# Define the variables
vars <- c("gender", "age_category", "ajcc_pathologic_stage", "ajcc_pathologic_m", 
          "tissue_or_organ_of_origin", "vital_status", "five_year_survival_factor", "three_year_survival_factor")

# Create the table
table1 <- CreateTableOne(vars = vars, data = working_dataframe_cleaned)

# Print the table
print(table1, showAllLevels = TRUE, formatOptions = list(big.mark = ","))

# Export to a CSV file
write.csv(print(table1, showAllLevels = TRUE, printToggle = FALSE), 
          file = "clinical_characteristics_summary.csv")

# Generate jpeg image
jpeg("survival_time_histogram.jpg", width = 800, height = 600)

# Create the histogram
hist(working_dataframe_cleaned$survival_years,
     main = "Distribution of Survival Time",
     xlab = "Survival Time (Years)",
     ylab = "Frequency",
     breaks = seq(0, 5, by = 1),  # Adjust breaks to cover 0 to 5 years in one-year intervals
     col = "skyblue",
     border = "black",
     right = FALSE,  # make the intervals left-inclusive and right-exclusive
     yaxt = "n")  # Remove default y-axis

# Add custom y-axis with intervals of 10
axis(2, at = seq(0, 80, by = 10), las = 1)

# Add a vertical line at 3 years
abline(v = 3, col = "red", lty = 2)

# Add a legend
legend("topright", legend = c("3-year mark"), col = "red", lty = 2)

# Close the current graphics device
dev.off()

## Machine Learning Models
ml_data <- working_dataframe_cleaned %>%
  select(bcr_patient_barcode, survival_time, censored_status, all_of(methylation_cols))

# Function to identify and remove columns with high NA percentage
remove_high_na_cols <- function(data, threshold = 10) {
  # Separate identifier and survival variables
  survival_vars <- data %>% select(bcr_patient_barcode, survival_time, censored_status)
  
  # Select only methylation columns (predictors)
  predictors <- data %>% select(-bcr_patient_barcode, -survival_time, -censored_status)
  
  # Calculate NA percentage per column
  na_percentage <- colMeans(is.na(predictors)) * 100
  
  # Identify columns to keep (less than threshold% NA)
  cols_to_keep <- names(na_percentage[na_percentage <= threshold])
  
  print(paste("Number of columns removed due to >", threshold, "% NA values:", 
              ncol(predictors) - length(cols_to_keep)))
  
  # Create cleaned dataset
  cleaned_data <- bind_cols(
    survival_vars,
    predictors %>% select(all_of(cols_to_keep))
  )
  
  return(cleaned_data)
}

# Remove columns with high NA percentage
ml_data_cleaned <- remove_high_na_cols(ml_data)

print(paste("Original dimensions:", ncol(ml_data), "columns"))
print(paste("Cleaned dimensions:", ncol(ml_data_cleaned), "columns"))

# Perform stratified split
set.seed(123)  # for reproducibility

# Create the stratified split
split_indices <- createDataPartition(ml_data_cleaned$censored_status, p = 0.7, list = FALSE)

# Create training and testing sets
train_data <- ml_data_cleaned[split_indices, ]
test_data <- ml_data_cleaned[-split_indices, ]

# Check distribution
print("Training set distribution:")
print(table(train_data$censored_status))
print(prop.table(table(train_data$censored_status)) * 100)

print("Testing set distribution:")
print(table(test_data$censored_status))
print(prop.table(table(test_data$censored_status)) * 100)

print(paste("Training set dimensions:", dim(train_data)[1], "rows,", dim(train_data)[2], "columns"))
print(paste("Testing set dimensions:", dim(test_data)[1], "rows,", dim(test_data)[2], "columns"))

# Function to perform imputation
perform_imputation <- function(data) {
  # Separate identifier and survival variables
  survival_vars <- data %>% select(bcr_patient_barcode, survival_time, censored_status)
  
  # Select only methylation columns (predictors)
  predictors <- data %>% select(-bcr_patient_barcode, -survival_time, -censored_status)
  
  # Perform imputation
  imputed_data <- impute.knn(as.matrix(predictors))$data
  
  # Convert back to data frame
  imputed_data <- as.data.frame(imputed_data)
  
  # Recombine imputed data with survival variables
  data_imputed <- bind_cols(survival_vars, imputed_data)
  
  # Check remaining NA values
  remaining_na <- sum(is.na(data_imputed))
  print(paste("Remaining NA values after imputation:", remaining_na))
  
  return(data_imputed)
}

# Apply the function to train and test sets separately
train_data_imputed <- perform_imputation(train_data)
test_data_imputed <- perform_imputation(test_data)

# Calculate imputation percentage for training set
train_na_count <- sum(is.na(train_data))
train_total_values <- nrow(train_data) * ncol(train_data)
train_imputation_percentage <- (train_na_count/train_total_values) * 100

# Calculate imputation percentage for testing set
test_na_count <- sum(is.na(test_data))
test_total_values <- nrow(test_data) * ncol(test_data)
test_imputation_percentage <- (test_na_count/test_total_values) * 100

print(paste("Training set imputation percentage:", round(train_imputation_percentage, 2), "%"))
print(paste("Testing set imputation percentage:", round(test_imputation_percentage, 2), "%"))

print(table(train_data_imputed$censored_status))
print(table(test_data_imputed$censored_status))
print(prop.table(table(train_data_imputed$censored_status)) * 100)
print(prop.table(table(test_data_imputed$censored_status)) * 100)

# Final verification
print("Final dataset dimensions:")
print(paste("Training set:", dim(train_data_imputed)[1], "rows,", dim(train_data_imputed)[2], "columns"))
print(paste("Testing set:", dim(test_data_imputed)[1], "rows,", dim(test_data_imputed)[2], "columns"))

# Update methylation_cols
methylation_cols <- grep("^cg", names(train_data_imputed), value = TRUE)

# Prepare dataset for lasso cox
prepare_for_lasso_cox <- function(data) {
  # Ensure correct data types
  data <- data %>%
    mutate(
      survival_time = as.numeric(survival_time),
      censored_status = as.integer(censored_status)
    )
  
  # Separate predictors and survival data (keep bcr_patient_barcode)
  predictors <- data %>% select(-bcr_patient_barcode, -survival_time, -censored_status)
  survival_data <- data %>% select(bcr_patient_barcode, survival_time, censored_status)
  
  # Check range of methylation values
  if(any(predictors < 0 | predictors > 1, na.rm = TRUE)) {
    warning("Some methylation beta values are outside the 0-1 range!")
  }
  
  # Combine predictors with survival data
  final_data <- bind_cols(survival_data, predictors)
  
  # Check data
  print(str(final_data))
  print(paste("Number of predictors:", ncol(predictors)))
  
  return(final_data)
}

# Apply to both train and test sets
train_data_final_bef <- prepare_for_lasso_cox(train_data_imputed)
test_data_final <- prepare_for_lasso_cox(test_data_imputed)

# Check survival_time cos LASSO Cox regression requires event times to be positive
sum(train_data_final_bef$survival_time <= 0)
sum(test_data_final$survival_time <= 0)

# Remove the 0 cases
train_data_final_bef <- train_data_final_bef[train_data_final_bef$survival_time > 0, ]
test_data_final <- test_data_final[test_data_final$survival_time > 0, ]

print(table(train_data_final_bef$censored_status))
print(table(test_data_final$censored_status))
print(prop.table(table(train_data_final_bef$censored_status)) * 100)
print(prop.table(table(test_data_final$censored_status)) * 100)


## Perform univariate Cox filtering
# Function to perform univariate Cox regression
perform_univariate_cox <- function(data, methylation_cols) {
  # Create empty list to store results
  results <- list()
  
  # Loop through each methylation site
  for(col in methylation_cols) {
    # Fit univariate Cox model
    formula <- as.formula(paste("Surv(survival_time, censored_status) ~", col))
    cox_model <- tryCatch({
      coxph(formula, data = data)
    }, error = function(e) {
      return(NULL)
    })
    
    # If model fit successful, store results
    if(!is.null(cox_model)) {
      # Get model summary
      sum <- summary(cox_model)
      
      # Store results
      results[[col]] <- data.frame(
        methylation_site = col,
        hazard_ratio = exp(coef(cox_model)),
        p_value = sum$coefficients[5],
        z_score = sum$coefficients[4]
      )
    }
  }
  
  # Combine all results
  results_df <- do.call(rbind, results)
  
  # Sort by p-value
  results_df <- results_df %>%
    arrange(p_value)
  
  return(results_df)
}

# Perform univariate Cox regression
univariate_results <- perform_univariate_cox(train_data_final_bef, methylation_cols)

# Count significant features
significant_features <- univariate_results %>%
  filter(p_value < 0.01)
print(paste("Number of significant features:", nrow(significant_features)))

# Get names of significant methylation sites
significant_sites <- significant_features$methylation_site

# Create new dataset with only significant features
train_data_final <- train_data_final_bef %>%
  select(bcr_patient_barcode, survival_time, censored_status, all_of(significant_sites))

print(paste("Original number of features:", length(methylation_cols)))
print(paste("Number of features after univariate filtering:", length(significant_sites)))

# Update testing data 
test_data_final <- test_data_final %>%
  select(bcr_patient_barcode, survival_time, censored_status, all_of(significant_sites))

  
## Perform cross-validated LASSO Cox regression
# Prepare the data for glmnet
x_train <- as.matrix(train_data_final %>% select(-bcr_patient_barcode, -survival_time, -censored_status))
y_train <- Surv(train_data_final$survival_time, train_data_final$censored_status)

# Fit LASSO Cox model with cross-validation
set.seed(123)  # for reproducibility
lasso_cox <- cv.glmnet(x_train, y_train, family = "cox", nfolds = 10)

# Create a PNG file to save both plots
png("lasso_cox_plots.png", width = 10, height = 10, units = "in", res = 300)

# Set up a 2x1 plotting area
par(mfrow = c(2,1))

# Plot A: Coefficient Paths
plot(lasso_cox$glmnet.fit, xvar = "lambda", label = FALSE,
     xlab = "Log(Lambda)", ylab = "Coefficients")
# Add label A in top-left corner
mtext("A", side = 3, line = 2, at = par("usr")[1] - 0.07 * diff(par("usr")[1:2]), adj = 0, font = 2, cex = 1.2)

# Plot B: Cross-validation Curve
plot(lasso_cox,
     xlab = "Log(Lambda)", ylab = "Partial Likelihood Deviance")
# Add label B in top-left corner
mtext("B", side = 3, line = 2, at = par("usr")[1] - 0.07 * diff(par("usr")[1:2]), adj = 0, font = 2, cex = 1.2)

# Close the PDF device
dev.off()

# Print optimal lambda values
cat("Lambda min:", lasso_cox$lambda.min, "\n")
cat("Lambda 1se:", lasso_cox$lambda.1se, "\n")

# Print number of non-zero coefficients for each lambda
coef_min <- coef(lasso_cox, s = "lambda.min")
coef_1se <- coef(lasso_cox, s = "lambda.1se")
cat("Number of non-zero coefficients (lambda.min):", sum(coef_min != 0), "\n")
cat("Number of non-zero coefficients (lambda.1se):", sum(coef_1se != 0), "\n")

# Get the minimum lambda value
lambda_min <- lasso_cox$lambda.min

# Get the coefficients at minimum lambda
coef_min <- coef(lasso_cox, s = "lambda.min")

# Get the names of selected features (non-zero coefficients)
selected_features <- rownames(coef_min)[coef_min[,1] != 0]

# Remove the intercept if it's included
selected_features <- selected_features[selected_features != "(Intercept)"]

print(paste("Number of selected features:", length(selected_features)))

# Create a new dataset with only selected features
selected_data <- train_data_final %>%
  select(bcr_patient_barcode, survival_time, censored_status, all_of(selected_features))


## Perform multivariate Cox regression
# Get clinical variables
clinical_vars <- working_dataframe_cleaned %>%
  select(bcr_patient_barcode, ajcc_pathologic_stage, age_at_index, gender) %>%
  distinct(bcr_patient_barcode, .keep_all = TRUE) %>%
  mutate(
    gender = as.factor(gender),
    stage_simplified = factor(case_when(
      ajcc_pathologic_stage %in% c("Stage IV", "Stage IVA", "Stage IVB", "Stage IVC") ~ "Stage IV",
      ajcc_pathologic_stage %in% c("Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC") ~ "Stage III", 
      ajcc_pathologic_stage %in% c("Stage II", "Stage IIA", "Stage IIB", "Stage IIC") ~ "Stage II",
      ajcc_pathologic_stage %in% c("Stage I", "Stage IA", "Stage IB", "Stage IC") ~ "Stage I",
      TRUE ~ NA_character_
    ), levels = c("Stage I", "Stage II", "Stage III", "Stage IV")),
    age = as.numeric(age_at_index)
  ) %>%
  select(bcr_patient_barcode, stage_simplified, age, gender)

# Merge training data with clinical variables
train_with_clinical <- selected_data %>%
  left_join(clinical_vars, by = "bcr_patient_barcode") %>%
  filter(!is.na(stage_simplified)) %>%
  select(-bcr_patient_barcode)  # Remove for modeling

# Merge test data with clinical variables
test_selected <- test_data_final %>%
  select(bcr_patient_barcode, survival_time, censored_status, all_of(selected_features))

test_with_clinical <- test_selected %>%
  left_join(clinical_vars, by = "bcr_patient_barcode") %>%
  filter(!is.na(stage_simplified)) %>%
  select(-bcr_patient_barcode)  # Remove for modeling

cat("Training set final dimensions:", dim(train_with_clinical), "\n")
cat("Test set final dimensions:", dim(test_with_clinical), "\n")

# Fit Cox models on training data
# Model 1: Methylation only
cox_meth <- coxph(Surv(survival_time, censored_status) ~ ., 
                  data = train_with_clinical %>% select(survival_time, censored_status, all_of(selected_features)))

# Model 2: Clinical only
cox_clinical <- coxph(Surv(survival_time, censored_status) ~ stage_simplified + age + gender, 
                      data = train_with_clinical)

# Model 3: Combined (methylation + clinical)
formula_combined <- as.formula(paste("Surv(survival_time, censored_status) ~", 
                                     paste(c(selected_features, "stage_simplified", "age", "gender"), 
                                           collapse = " + ")))
cox_combined <- coxph(formula_combined, data = train_with_clinical)

# Calculate risk scores
# Training set
train_with_clinical$risk_meth <- predict(cox_meth, newdata = train_with_clinical, type = "lp")
train_with_clinical$risk_clinical <- predict(cox_clinical, newdata = train_with_clinical, type = "lp")
train_with_clinical$risk_combined <- predict(cox_combined, newdata = train_with_clinical, type = "lp")

# Test set
test_with_clinical$risk_meth <- predict(cox_meth, newdata = test_with_clinical, type = "lp")
test_with_clinical$risk_clinical <- predict(cox_clinical, newdata = test_with_clinical, type = "lp")
test_with_clinical$risk_combined <- predict(cox_combined, newdata = test_with_clinical, type = "lp")

# Calculate AUC
time_points <- c(365, 3*365, 5*365)

# Training AUCs
auc_train_meth <- timeROC(T = train_with_clinical$survival_time,
                          delta = train_with_clinical$censored_status,
                          marker = train_with_clinical$risk_meth,
                          cause = 1, times = time_points, iid = TRUE)

auc_train_clinical <- timeROC(T = train_with_clinical$survival_time,
                              delta = train_with_clinical$censored_status,
                              marker = train_with_clinical$risk_clinical,
                              cause = 1, times = time_points, iid = TRUE)

auc_train_combined <- timeROC(T = train_with_clinical$survival_time,
                              delta = train_with_clinical$censored_status,
                              marker = train_with_clinical$risk_combined,
                              cause = 1, times = time_points, iid = TRUE)

# Test AUCs
auc_test_meth <- timeROC(T = test_with_clinical$survival_time,
                         delta = test_with_clinical$censored_status,
                         marker = test_with_clinical$risk_meth,
                         cause = 1, times = time_points, iid = TRUE)

auc_test_clinical <- timeROC(T = test_with_clinical$survival_time,
                             delta = test_with_clinical$censored_status,
                             marker = test_with_clinical$risk_clinical,
                             cause = 1, times = time_points, iid = TRUE)

auc_test_combined <- timeROC(T = test_with_clinical$survival_time,
                             delta = test_with_clinical$censored_status,
                             marker = test_with_clinical$risk_combined,
                             cause = 1, times = time_points, iid = TRUE)

# Create AUC results table
auc_results <- data.frame(
  Time_Point = c("1 Year", "3 Years", "5 Years"),
  Train_Methylation = auc_train_meth$AUC,
  Train_Clinical = auc_train_clinical$AUC,
  Train_Combined = auc_train_combined$AUC,
  Test_Methylation = auc_test_meth$AUC,
  Test_Clinical = auc_test_clinical$AUC,
  Test_Combined = auc_test_combined$AUC
)

print("AUC Results:")
print(auc_results)
write.csv(auc_results, "final_auc_results.csv", row.names = FALSE)

# Plot ROC curves
# Training ROC
pdf("ROC_training.pdf", width = 12, height = 8)
plot(auc_train_meth, time = 3*365, col = "blue", title = "Training Set ROC at 3 Years")
plot(auc_train_clinical, time = 3*365, col = "green", add = TRUE)
plot(auc_train_combined, time = 3*365, col = "red", add = TRUE)
legend("bottomright", legend = c("Methylation", "Clinical", "Combined"), 
       col = c("blue", "green", "red"), lwd = 2)
dev.off()

# Test ROC
pdf("ROC_test.pdf", width = 12, height = 8)
plot(auc_test_meth, time = 3*365, col = "blue", title = "Test Set ROC at 3 Years")
plot(auc_test_clinical, time = 3*365, col = "green", add = TRUE)
plot(auc_test_combined, time = 3*365, col = "red", add = TRUE)
legend("bottomright", legend = c("Methylation", "Clinical", "Combined"), 
       col = c("blue", "green", "red"), lwd = 2)
dev.off()

# Create Kaplan-Meier curves
create_km <- function(data, risk_col, title, use_train_median = FALSE) {
  if (use_train_median) {
    median_risk <- median(train_with_clinical[[risk_col]], na.rm = TRUE)
  } else {
    median_risk <- median(data[[risk_col]], na.rm = TRUE)
  }
  
  data$risk_group <- ifelse(data[[risk_col]] > median_risk, "High", "Low")
  fit <- survfit(Surv(data$survival_time, data$censored_status) ~ risk_group, data = data)
  
  ggsurvplot(fit, data = data, pval = TRUE, conf.int = TRUE,
             risk.table = TRUE, risk.table.col = "strata",
             ggtheme = theme_bw(), palette = c("#E7B800", "#2E9FDF"),
             title = title)
}

# Create all KM plots
km_train_meth <- create_km(train_with_clinical, "risk_meth", "Training: Methylation Only")
km_train_clinical <- create_km(train_with_clinical, "risk_clinical", "Training: Clinical Only")
km_train_combined <- create_km(train_with_clinical, "risk_combined", "Training: Combined")
km_test_meth <- create_km(test_with_clinical, "risk_meth", "Test: Methylation Only", TRUE)
km_test_clinical <- create_km(test_with_clinical, "risk_clinical", "Test: Clinical Only", TRUE)
km_test_combined <- create_km(test_with_clinical, "risk_combined", "Test: Combined", TRUE)

# Save KM plots
ggsave("KM_train_meth.png", km_train_meth$plot, width = 10, height = 8)
ggsave("KM_train_clinical.png", km_train_clinical$plot, width = 10, height = 8)
ggsave("KM_train_combined.png", km_train_combined$plot, width = 10, height = 8)
ggsave("KM_test_meth.png", km_test_meth$plot, width = 10, height = 8)
ggsave("KM_test_clinical.png", km_test_clinical$plot, width = 10, height = 8)
ggsave("KM_test_combined.png", km_test_combined$plot, width = 10, height = 8)

# Save model coefficients for all models
# Methylation model coefficients
coeffs_meth <- coef(cox_meth)[!is.na(coef(cox_meth))]
coeff_meth_df <- data.frame(
  Model = "Methylation",
  Feature = names(coeffs_meth),
  Coefficient = coeffs_meth,
  Hazard_Ratio = exp(coeffs_meth),
  P_Value = summary(cox_meth)$coefficients[!is.na(coef(cox_meth)), "Pr(>|z|)"]
)

# Clinical model coefficients
coeffs_clinical <- coef(cox_clinical)[!is.na(coef(cox_clinical))]
coeff_clinical_df <- data.frame(
  Model = "Clinical",
  Feature = names(coeffs_clinical),
  Coefficient = coeffs_clinical,
  Hazard_Ratio = exp(coeffs_clinical),
  P_Value = summary(cox_clinical)$coefficients[!is.na(coef(cox_clinical)), "Pr(>|z|)"]
)

# Combined model coefficients
coeffs_combined <- coef(cox_combined)[!is.na(coef(cox_combined))]
coeff_combined_df <- data.frame(
  Model = "Combined",
  Feature = names(coeffs_combined),
  Coefficient = coeffs_combined,
  Hazard_Ratio = exp(coeffs_combined),
  P_Value = summary(cox_combined)$coefficients[!is.na(coef(cox_combined)), "Pr(>|z|)"]
)

# Combine all coefficients
all_coeffs_df <- rbind(coeff_meth_df, coeff_clinical_df, coeff_combined_df) %>%
  arrange(Model, P_Value)

write.csv(all_coeffs_df, "all_model_coefficients.csv", row.names = FALSE)

# Print model summaries
print("=== Methylation Model Summary ===")
print(summary(cox_meth))

print("=== Clinical Model Summary ===")
print(summary(cox_clinical))

print("=== Combined Model Summary ===")
print(summary(cox_combined))


## Pathway analysis

# Function to create annotation heatmap
create_methylation_heatmap <- function(data, significant_sites, survival_info) {
  # Extract methylation matrix
  methylation_matrix <- as.matrix(data[, significant_sites])
  
  # Scale the data
  scaled_matrix <- t(scale(methylation_matrix))
  
  # Create survival annotation with years
  ha = HeatmapAnnotation(
    survival_years = survival_info$survival_time / 365.25,
    status = survival_info$censored_status,
    col = list(
      survival_years = colorRamp2(c(min(survival_info$survival_time / 365.25), 
                                    max(survival_info$survival_time / 365.25)), 
                                  c("white", "darkblue")),
      status = c("0" = "grey", "1" = "black")
    )
  )
  
  # Create column order based on status and survival time
  column_order <- order(survival_info$censored_status, survival_info$survival_time)
  
  # Create heatmap
  heatmap <- Heatmap(scaled_matrix,
                     name = "Methylation Level",
                     show_row_names = FALSE,
                     show_column_names = FALSE,
                     clustering_distance_rows = "euclidean",
                     column_order = column_order,
                     cluster_columns = FALSE,
                     clustering_method_rows = "complete",
                     col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
                     top_annotation = ha,
                     show_row_dend = TRUE,
                     show_column_dend = FALSE,
                     row_title = "Methylation Sites",
                     column_title = "Samples")
  
  return(heatmap)
}

# Create heatmap for initial significant sites
initial_heatmap <- create_methylation_heatmap(
  train_data_final_bef,
  significant_sites,
  data.frame(
    survival_time = train_data_final_bef$survival_time,
    censored_status = train_data_final_bef$censored_status
  )
)
png("methylation_heatmap.png", width = 1200, height = 800, res = 150)
draw(initial_heatmap)
dev.off()

# Create heatmap for selected sites
final_heatmap <- Heatmap(
  # Use the data directly from train_data_final_bef
  t(scale(as.matrix(train_data_final_bef[, selected_features]))),
  name = "Methylation\nLevel",
  
  # Basic color settings
  col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
  
  # Column ordering by status and survival time
  column_order = order(train_data_final_bef$censored_status, 
                       train_data_final_bef$survival_time),
  cluster_columns = FALSE,
  
  # Simple clustering
  clustering_distance_rows = "pearson",
  
  # Basic appearance
  show_row_names = TRUE,
  show_column_names = FALSE,
  show_column_dend = FALSE,
  
  # Top annotation: survival time in years
  top_annotation = HeatmapAnnotation(
    survival_years = train_data_final_bef$survival_time / 365.25,
    status = train_data_final_bef$censored_status,
    col = list(
      survival_years = colorRamp2(c(min(train_data_final_bef$survival_time / 365.25), 
                                    max(train_data_final_bef$survival_time / 365.25)), 
                                  c("white", "darkblue")),
      status = c("0" = "grey", "1" = "black")
    )
  )
)

png("selected_sites_heatmap.png", width = 1200, height = 800, res = 150)
draw(final_heatmap)
dev.off()

# Function to create marker summary
analyze_methylation_markers <- function(methylation_data, cox_coefficients_df, survival_data) {
  # Get annotation information
  ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  
  # Calculate beta values by survival status
  beta_by_status <- lapply(cox_coefficients_df$Feature, function(cpg) {
    data.frame(
      # Mean beta for all samples
      Mean_Beta_All = mean(methylation_data[[cpg]], na.rm = TRUE),
      # Mean beta for deceased patients
      Mean_Beta_Deceased = mean(methylation_data[[cpg]][survival_data$censored_status == 1], 
                                na.rm = TRUE),
      # Mean beta for surviving patients
      Mean_Beta_Surviving = mean(methylation_data[[cpg]][survival_data$censored_status == 0], 
                                 na.rm = TRUE),
      # Wilcoxon test p-value between groups
      Wilcox_P = wilcox.test(
        methylation_data[[cpg]][survival_data$censored_status == 1],
        methylation_data[[cpg]][survival_data$censored_status == 0]
      )$p.value
    )
  })
  
  # Combine into a single data frame
  beta_stats <- do.call(rbind, beta_by_status)
  
  # Create summary dataframe
  marker_summary <- data.frame(
    CpG_Site = cox_coefficients_df$Feature,
    Coefficient = cox_coefficients_df$Coefficient,
    HR = cox_coefficients_df$Hazard_Ratio,
    beta_stats,
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      Gene_Symbol = ann450k[CpG_Site, "UCSC_RefGene_Name"],
      Chromosome = ann450k[CpG_Site, "chr"],
      Position = ann450k[CpG_Site, "pos"],
      Relation_to_Island = ann450k[CpG_Site, "Relation_to_Island"],
      # Format statistics
      HR_formatted = sprintf("%.2e", HR),
      Beta_Difference = Mean_Beta_Deceased - Mean_Beta_Surviving,
      Association = case_when(
        Coefficient > 0 ~ "Poor Survival",
        Coefficient < 0 ~ "Better Survival",
        TRUE ~ "Neutral"
      )
    )
  
  return(marker_summary)
}

# Create formatted summary table
create_summary_table <- function(marker_summary) {
  summary_table <- marker_summary %>%
    dplyr::select(
      CpG_Site,
      Gene_Symbol,
      HR_formatted,
      Mean_Beta_All,
      Mean_Beta_Deceased,
      Mean_Beta_Surviving,
      Beta_Difference,
      Wilcox_P,
      Association,
      Relation_to_Island
    ) %>%
    dplyr::mutate(
      # Format p-values and means
      Wilcox_P = format.pval(Wilcox_P, digits = 3),
      Mean_Beta_All = round(Mean_Beta_All, 3),
      Mean_Beta_Deceased = round(Mean_Beta_Deceased, 3),
      Mean_Beta_Surviving = round(Mean_Beta_Surviving, 3),
      Beta_Difference = round(Beta_Difference, 3)
    )
  
  names(summary_table) <- c(
    "CpG Site", 
    "Gene",
    "Hazard Ratio",
    "Mean β (All)",
    "Mean β (Deceased)",
    "Mean β (Surviving)",
    "β Difference",
    "Wilcoxon P",
    "Survival Association",
    "CpG Location"
  )
  
  return(summary_table)
}

# Get methylation coefficients from combined model
methylation_coeffs <- coeff_combined_df %>%
  filter(grepl("^cg", Feature))

# Create marker summary with group comparisons
marker_summary <- analyze_methylation_markers(
  methylation_data = train_data_final,
  cox_coefficients_df = methylation_coeffs,
  survival_data = data.frame(
    censored_status = train_data_final$censored_status
  )
)

# Create summary table
summary_table <- create_summary_table(marker_summary)
write.csv(summary_table, "methylation_markers_summary.csv", row.names = FALSE)

# Function to create boxplots for each CpG site
create_simple_boxplots <- function(methylation_data, feature_names, censored_status) {
  # Prepare data for plotting
  plot_data <- data.frame(
    Status = factor(censored_status, 
                    levels = c(0, 1),
                    labels = c("Surviving", "Deceased"))
  )
  
  # Add methylation data
  plot_data <- cbind(plot_data, methylation_data[, feature_names])
  
  # Convert to long format
  plot_data_long <- plot_data %>%
    pivot_longer(
      cols = -Status,
      names_to = "CpG_Site",
      values_to = "Beta_Value"
    )
  
  # Create plot
  p <- ggplot(plot_data_long, aes(x = Status, y = Beta_Value, fill = Status)) +
    geom_boxplot() +
    geom_jitter(width = 0.2, alpha = 0.2) +
    facet_wrap(~CpG_Site, ncol = 3) +
    scale_fill_manual(values = c("Surviving" = "lightblue", "Deceased" = "salmon")) +
    theme_bw() +
    labs(
      title = "Methylation Beta Values by Survival Status",
      y = "Beta Value",
      x = "Patient Status"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  return(p)
}

# Usage
boxplot <- create_simple_boxplots(
  methylation_data = train_data_final,
  feature_names = selected_features,
  censored_status = train_data_final$censored_status
)

# Save plot
ggsave("methylation_boxplots.pdf", boxplot, width = 12, height = 8)

# Create forest plot with scientific notation
create_forest_plot <- function(marker_summary) {
  plot_data <- marker_summary %>%
    mutate(
      # Remove duplicate gene names but keep unique ones
      Clean_Gene = sapply(Gene_Symbol, function(x) {
        if (is.na(x)) return(NA)
        genes <- unique(trimws(unlist(strsplit(x, "[;,]"))))
        paste(genes, collapse = " or ")
      }),
      
      # Create clean labels
      Gene_Label = ifelse(is.na(Clean_Gene), CpG_Site, 
                          paste0(CpG_Site, "\n(", Clean_Gene, ")")),
      
      log2HR = log2(HR)
    )
  
  ggplot(plot_data, aes(y = reorder(Gene_Label, log2HR))) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_point(aes(x = log2HR, color = Association), size = 3) +
    scale_color_manual(values = c("Poor Survival" = "red", "Better Survival" = "blue")) +
    theme_bw() +
    labs(x = "log2(Hazard Ratio)", y = "CpG Site", 
         title = "Forest Plot of Methylation Markers")
}

# Create and save forest plot
forest_plot <- create_forest_plot(marker_summary)
ggsave("methylation_forest_plot.pdf", forest_plot, width = 12, height = 8)


# Gene analysis and pathway enrichment
# Load required libraries
if (!require("clusterProfiler", quietly = TRUE)) {
  BiocManager::install("clusterProfiler")
  library(clusterProfiler)
}
if (!require("org.Hs.eg.db", quietly = TRUE)) {
  BiocManager::install("org.Hs.eg.db")
  library(org.Hs.eg.db)
}

# Define gene groups
high_risk_genes <- c("C14orf179", "HDAC10", "CDK10", "AKD1", "FIG4", 
                     "ATRIP", "RBP7", "ZMYND15", "CXCL16", "ITPKB")  # 10 genes

protective_genes <- c("OR5P2", "GRIP1")  # 2 genes

print(paste("High-risk genes:", length(high_risk_genes)))
print(paste("Protective genes:", length(protective_genes)))

# GO Analysis for each group
go_high_risk <- enrichGO(
  gene = high_risk_genes,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05
)

go_protective <- enrichGO(
  gene = protective_genes,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05
)

# KEGG Analysis for each group
# Convert to Entrez IDs
high_risk_entrez <- bitr(high_risk_genes, "SYMBOL", "ENTREZID", org.Hs.eg.db)
protective_entrez <- bitr(protective_genes, "SYMBOL", "ENTREZID", org.Hs.eg.db)

kegg_high_risk <- enrichKEGG(
  gene = high_risk_entrez$ENTREZID,
  organism = "hsa",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05
)

kegg_protective <- enrichKEGG(
  gene = protective_entrez$ENTREZID,
  organism = "hsa",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05
)

# Check results
print("=== RESULTS BY GROUP ===")
print(paste("High-risk GO terms:", ifelse(is.null(go_high_risk), 0, nrow(go_high_risk@result))))
print(paste("Protective GO terms:", ifelse(is.null(go_protective), 0, nrow(go_protective@result))))
print(paste("High-risk KEGG pathways:", ifelse(is.null(kegg_high_risk), 0, nrow(kegg_high_risk@result))))
print(paste("Protective KEGG pathways:", ifelse(is.null(kegg_protective), 0, nrow(kegg_protective@result))))

# Create combined plot
create_combined_plot <- function(high_risk_result, protective_result, title_prefix) {
  plot_data <- data.frame()
  
  # Add high-risk data
  if (!is.null(high_risk_result) && nrow(high_risk_result@result) > 0) {
    hr_data <- high_risk_result@result %>%
      arrange(p.adjust) %>%
      head(8) %>%
      mutate(
        Type = "High-risk",
        neg_log_p = -log10(p.adjust)
      )
    plot_data <- rbind(plot_data, hr_data)
  }
  
  # Add protective data
  if (!is.null(protective_result) && nrow(protective_result@result) > 0) {
    prot_data <- protective_result@result %>%
      arrange(p.adjust) %>%
      head(8) %>%
      mutate(
        Type = "Protective",
        neg_log_p = -log10(p.adjust)
      )
    plot_data <- rbind(plot_data, prot_data)
  }
  
  if (nrow(plot_data) == 0) {
    print(paste("No significant pathways found for", title_prefix))
    return(NULL)
  }
  
  # Create plot
  p <- ggplot(plot_data, aes(x = neg_log_p, y = reorder(Description, neg_log_p))) +
    geom_point(aes(color = Type, size = Count)) +
    scale_color_manual(values = c("High-risk" = "red", "Protective" = "blue")) +
    theme_bw() +
    labs(
      x = "-log10(Adjusted P-value)",
      y = "Pathway",
      title = paste(title_prefix, "Pathways by Risk Group"),
      color = "Gene Type"
    )
  
  return(p)
}

# Create and save plots
go_combined_plot <- create_combined_plot(go_high_risk, go_protective, "GO")
if (!is.null(go_combined_plot)) {
  ggsave("go_risk_comparison.png", go_combined_plot, width = 12, height = 8)
}

kegg_combined_plot <- create_combined_plot(kegg_high_risk, kegg_protective, "KEGG")
if (!is.null(kegg_combined_plot)) {
  ggsave("kegg_risk_comparison.png", kegg_combined_plot, width = 12, height = 8)
}

# Save all results
if (!is.null(go_high_risk)) write.csv(go_high_risk@result, "go_high_risk.csv", row.names = FALSE)
if (!is.null(go_protective)) write.csv(go_protective@result, "go_protective.csv", row.names = FALSE)
if (!is.null(kegg_high_risk)) write.csv(kegg_high_risk@result, "kegg_high_risk.csv", row.names = FALSE)
if (!is.null(kegg_protective)) write.csv(kegg_protective@result, "kegg_protective.csv", row.names = FALSE)

# Check which genes failed to map
print("High-risk genes that mapped:")
print(high_risk_entrez)

# Find which ones didn't map
failed_genes <- setdiff(high_risk_genes, high_risk_entrez$SYMBOL)
print(paste("Failed to map:", paste(failed_genes, collapse = ", ")))

