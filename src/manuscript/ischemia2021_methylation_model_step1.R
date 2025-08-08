###########################################################################
#
#                            muller_methylation_model.R
#
###########################################################################
# Author: Matthew Muller
# Date: 2025-07-21
# Script Name: muller_methylation_model.R
# Description: R version of the methylation classification script

# ======================== SETUP ========================

# Load required libraries
library(tidyverse)
library(caret)
library(randomForest)
library(e1071)
library(gbm)
library(glmnet)
library(readr)
library(jsonlite)
library(data.table)

# Set up directories
root_dir <- getwd()
experiment <- "muller_methylation_model_R"
outdir <- file.path(root_dir, "output", experiment)

# Create output directory if it doesn't exist
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}

# Set random seed for reproducibility
set.seed(420)

# Helper function to save results
save_results <- function(object, filename) {
  saveRDS(object, file.path(outdir, paste0(filename, ".rds")))
}

# ======================== DATA LOADING ========================
cat("=== LOADING DATA ===\n")

# We are going to load the data from prior processed files
data_path <- file.path("output/muller_methylation_model")

X_train <- read.csv(file.path(data_path, "X_train.csv"), row.names = 1)
X_test <- read.csv(file.path(data_path, "X_test.csv"), row.names = 1)
y_train <- read.csv(file.path(data_path, "y_train.csv"), row.names = 1)
y_test <- read.csv(file.path(data_path, "y_test.csv"), row.names = 1)

# Make some histograms of the X data to see the distributions
cat("Visualizing distributions of training features...\n")
features <- X_train %>% 
  select(where(is.numeric)) %>% 
  pivot_longer(cols = everything(), names_to = "feature", values_to = "value") %>%
  ggplot(aes(x = value, color = feature)) +
  geom_density(alpha = 0.1) +
  labs(title = "Distribution of Training Features", x = "Value", y = "Density") +
  theme_bw() +
  theme(legend.position = "none")
ggsave(file.path(outdir, "feature_distributions.png"), features, width = 10, height = 6)

# Convert labels to factors
y_train <- factor(y_train$meth_3cluster, labels = c("MS1", "MS2", "MS3"))
y_test <- factor(y_test$meth_3cluster, labels = c("MS1", "MS2", "MS3"))

cat("Training set:", nrow(X_train), "samples\n")
cat("Test set:", nrow(X_test), "samples\n")
cat("Training class distribution:\n")
print(table(y_train))
cat("Test class distribution:\n")
print(table(y_test))

# Save the feature names
feature_names <- colnames(X_train)
write_csv(data.frame(feature = feature_names), file.path(outdir, "feature_names.csv"))

# ======================== MODEL TRAINING ========================
cat("\n=== MODEL TRAINING ===\n")

# Set up cross-validation
train_control <- trainControl(
  method = "cv",
  number = 5,
  classProbs = TRUE,
  summaryFunction = multiClassSummary,
  savePredictions = "final"
)

# Create directory for saving models
model_save_dir <- file.path(outdir, "trained_models")
if (!dir.exists(model_save_dir)) {
  dir.create(model_save_dir, recursive = TRUE)
}

# Initialize results storage
model_results <- list()
cv_results <- list()
trained_models <- list()

# ======================== LOGISTIC REGRESSION (REGULARIZED) ========================
cat("Training Regularized Logistic Regression...\n")
# Use glmnet for regularized multinomial regression
model_lr <- train(
  x = X_train,
  y = y_train,
  method = "glmnet",
  trControl = train_control,
  metric = "Accuracy",
  family = "multinomial",
  tuneGrid = expand.grid(
    alpha = c(0, 0.5, 1),  # 0=ridge, 0.5=elastic net, 1=lasso
    lambda = c(0.001, 0.01, 0.1, 1)
  )
)

trained_models[["regularized_logistic"]] <- model_lr
cv_results[["regularized_logistic"]] <- model_lr$results
saveRDS(model_lr, file.path(model_save_dir, "regularized_logistic_model.rds"))
cat("✓ Saved regularized logistic regression model\n")

# ======================== SVM LINEAR ========================
cat("Training SVM (Linear Kernel)...\n")
model_svm_linear <- train(
  x = X_train,
  y = y_train,
  method = "svmLinear",
  trControl = train_control,
  metric = "Accuracy",
  tuneGrid = expand.grid(
    C = c(0.01, 0.1, 1, 10, 100)
  )
)

trained_models[["svm_linear"]] <- model_svm_linear
cv_results[["svm_linear"]] <- model_svm_linear$results
saveRDS(model_svm_linear, file.path(model_save_dir, "svm_linear_model.rds"))
cat("✓ Saved SVM linear model\n")

# ======================== SVM RBF ========================
cat("Training SVM (RBF Kernel)...\n")
model_svm_rbf <- train(
  x = X_train,
  y = y_train,
  method = "svmRadial",
  trControl = train_control,
  metric = "Accuracy",
  tuneGrid = expand.grid(
    C = c(0.1, 1, 10),
    sigma = c(0.001, 0.01, 0.1)
  )
)

trained_models[["svm_rbf"]] <- model_svm_rbf
cv_results[["svm_rbf"]] <- model_svm_rbf$results
saveRDS(model_svm_rbf, file.path(model_save_dir, "svm_rbf_model.rds"))
cat("✓ Saved SVM RBF model\n")

# ======================== RANDOM FOREST ========================
cat("Training Random Forest...\n")
model_rf <- train(
  x = X_train,
  y = y_train,
  method = "rf",
  trControl = train_control,
  metric = "Accuracy",
)

trained_models[["random_forest"]] <- model_rf
cv_results[["random_forest"]] <- model_rf$results
saveRDS(model_rf, file.path(model_save_dir, "random_forest_model.rds"))
cat("✓ Saved Random Forest model\n")

# ======================== KNN ========================
cat("Training K-Nearest Neighbors (KNN)...\n")
model_knn <- train(
  x = X_train,
  y = y_train,
  method = "knn",
  trControl = train_control,
  metric = "Accuracy",
  tuneGrid = expand.grid(k = c(3, 5, 7, 9, 11))
)

trained_models[["knn"]] <- model_knn
cv_results[["knn"]] <- model_knn$results
saveRDS(model_knn, file.path(model_save_dir, "knn_model.rds"))
cat("✓ Saved KNN model\n")

# ======================== GRADIENT BOOSTING ========================
cat("Training Gradient Boosting Machine...\n")
model_gbm <- train(
  x = X_train,
  y = y_train,
  method = "gbm",
  trControl = train_control,
  metric = "Accuracy",
  verbose = FALSE,
)

trained_models[["gradient_boosting"]] <- model_gbm
cv_results[["gradient_boosting"]] <- model_gbm$results
saveRDS(model_gbm, file.path(model_save_dir, "gradient_boosting_model.rds"))
cat("✓ Saved Gradient Boosting model\n")

# ======================== MODEL EVALUATION ========================
cat("\n=== MODEL EVALUATION ===\n")

# Function to evaluate model on test set
evaluate_model <- function(model, model_name) {
  # Make predictions
  pred_prob <- predict(model, X_test, type = "prob")
  pred_class <- predict(model, X_test)
  
  # Calculate metrics
  conf_matrix <- confusionMatrix(pred_class, y_test)
  
  # Store results
  results <- list(
    model_name = model_name,
    accuracy = conf_matrix$overall["Accuracy"],
    kappa = conf_matrix$overall["Kappa"],
    sensitivity = conf_matrix$byClass[, "Sensitivity"],
    specificity = conf_matrix$byClass[, "Specificity"],
    precision = conf_matrix$byClass[, "Precision"],
    recall = conf_matrix$byClass[, "Recall"],
    f1 = conf_matrix$byClass[, "F1"],
    confusion_matrix = conf_matrix$table,
    predictions = pred_class,
    probabilities = pred_prob
  )
  
  return(results)
}

# Evaluate all models
test_results <- list()
for (model_name in names(trained_models)) {
  cat("Evaluating", model_name, "...\n")
  test_results[[model_name]] <- evaluate_model(trained_models[[model_name]], model_name)
}

# ======================== RESULTS SUMMARY ========================
cat("\n=== RESULTS SUMMARY ===\n")

# Create summary table
summary_df <- data.frame(
  Model = character(),
  CV_Accuracy = numeric(),
  Test_Accuracy = numeric(),
  Test_Kappa = numeric(),
  stringsAsFactors = FALSE
)

for (model_name in names(trained_models)) {
  cv_acc <- max(cv_results[[model_name]]$Accuracy, na.rm = TRUE)
  test_acc <- test_results[[model_name]]$accuracy
  test_kappa <- test_results[[model_name]]$kappa
  
  summary_df <- rbind(summary_df, data.frame(
    Model = model_name,
    CV_Accuracy = round(cv_acc, 4),
    Test_Accuracy = round(test_acc, 4),
    Test_Kappa = round(test_kappa, 4)
  ))
}

# Sort by test accuracy
summary_df <- summary_df[order(summary_df$Test_Accuracy, decreasing = TRUE), ]
print(summary_df)

# Save summary
write_csv(summary_df, file.path(outdir, "model_summary.csv"))
cat("\n✓ Model summary saved to model_summary.csv\n")

# ======================== DETAILED RESULTS ========================
cat("\n=== DETAILED RESULTS ===\n")

# Print detailed results for each model
for (model_name in names(test_results)) {
  cat("\n", toupper(model_name), "RESULTS:\n")
  cat("Test Accuracy:", round(test_results[[model_name]]$accuracy, 4), "\n")
  cat("Test Kappa:", round(test_results[[model_name]]$kappa, 4), "\n")
  cat("Confusion Matrix:\n")
  print(test_results[[model_name]]$confusion_matrix)
  
  # Save detailed results as JSON
  results_json <- list(
    model_name = model_name,
    test_accuracy = as.numeric(test_results[[model_name]]$accuracy),
    test_kappa = as.numeric(test_results[[model_name]]$kappa),
    confusion_matrix = as.data.frame.matrix(test_results[[model_name]]$confusion_matrix),
    cross_validation_results = cv_results[[model_name]]
  )
  
  json_file <- file.path(model_save_dir, paste0(model_name, "_results.json"))
  write_json(results_json, json_file, pretty = TRUE)
  cat("✓ Detailed results saved to", paste0(model_name, "_results.json"), "\n")
}

# ======================== SAVE WORKSPACE ========================
cat("\n=== SAVING WORKSPACE ===\n")

# Save all important objects
workspace_objects <- list(
  X_train = X_train,
  X_test = X_test,
  y_train = y_train,
  y_test = y_test,
  trained_models = trained_models,
  cv_results = cv_results,
  test_results = test_results,
  summary_df = summary_df
)

saveRDS(workspace_objects, file.path(outdir, "workspace_objects.rds"))
cat("✓ Workspace objects saved to workspace_objects.rds\n")

# Save session info
session_info <- sessionInfo()
saveRDS(session_info, file.path(outdir, "session_info.rds"))
cat("✓ Session info saved\n")

# ======================== FINAL SUMMARY ========================
cat("\n=== FINAL SUMMARY ===\n")
cat("✓ All models trained and evaluated successfully\n")
cat("✓ Cross-validation completed for all models\n")
cat("✓ Test set evaluation completed\n")
cat("✓ Results saved to:", outdir, "\n")
cat("✓ Trained models saved to:", model_save_dir, "\n")
cat("\nBest performing model on test set:", summary_df$Model[1], 
    "with accuracy:", summary_df$Test_Accuracy[1], "\n")

cat("\nFiles created:\n")
cat("- model_summary.csv: Summary of all model performances\n")
cat("- workspace_objects.rds: All R objects for later analysis\n")
cat("- session_info.rds: R session information\n")
cat("- trained_models/: Directory with saved model objects\n")
cat("- *_results.json: Detailed results for each model\n")

cat("\nScript completed successfully!\n")
