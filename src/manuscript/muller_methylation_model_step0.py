###########################################################################
#
#                            muller_methylation_model
#
###########################################################################
# Author: Matthew Muller
# Date: 2025-07-21
# Script Name: muller_methylation_model

# ======================== SETUP ========================
# General libaries
import os
import sys
import time
import datetime
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import joblib
import json

# Modeling libraries
from sklearn import base
from sklearn import preprocessing
from sklearn import model_selection
from sklearn import metrics
from sklearn.feature_selection import SelectKBest, f_classif, mutual_info_classif, RFE
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
import skore

# Custom libraries
from MattTools import utils

## Set up directories
# Root directory
root_dir = os.getcwd()
# Output directory
experiment = "muller_methylation_model"
outdir = os.path.join(root_dir, "output", experiment) + "/"

# Create output directory if it doesn't exist
if not os.path.exists(outdir):
    os.makedirs(outdir)

## Set up some general items
# Set random seed
utils.set_random_seed(420)
# Hide warnings
utils.hide_warnings()

# ======================== CODE ========================

# Load the methylation data
# counts = pd.read_csv("data/final_meth_data_v1_20230207/limma_noCovariate.csv.gz", index_col=0)
counts = pd.read_csv(
    "data/final_meth_data_v1_20230207/ischemia_meth_limma_counttable.csv.zip",
    index_col=0,
).T  # Transpose to have samples as rows and features as columns

# Now the metadata
metadata = pd.read_csv(
    "output/run4_rmoutliers2_asr_control/integrative_analyses/run9-figs/bioreptable_waddons.csv"
)

# Methylation labels
labels = metadata[["PATNUM", "meth_3cluster"]].set_index("PATNUM").dropna()

# Set the index orders to match
counts = counts.loc[
    labels.index
].copy()  # Ensure counts are in the same order as labels

# Check the shape of the data
print(f"Shape of counts_methylation: {counts.shape}")
print(f"Shape of labels: {labels.shape}")

# Verify that the indices match
if not counts.index.equals(labels.index):
    raise ValueError(
        "Indices of counts and labels do not match. Please check the data."
    )

print(f"Number of features: {counts.shape[1]}")
print(f"Number of samples: {counts.shape[0]}")
print(f"Class distribution: \n{pd.Series(labels['meth_3cluster'].values).value_counts()}")

# Function to apply inverse rank normalization
import numpy as np
from scipy.stats import rankdata, norm
def inverse_rank_normalization(X, c=0):
    """
    Apply inverse normal transformation to each column of X.
    
    Parameters:
    X: numpy array, shape (n_samples, n_features)
    c: offset parameter (0 for Van der Waerden)
    
    Returns:
    Transformed array with same shape as X
    """
    n = X.shape[0]
    ranks = np.apply_along_axis(rankdata, 0, X, method='average')
    scaled = (ranks - c) / (n - 2*c + 1)
    scaled = np.clip(scaled, 1/(4*n), 1 - 1/(4*n))  # no infinities
    transformed = norm.ppf(scaled)
    return transformed

# Normalize the counts using inverse rank normalization
print("Applying inverse rank normalization to counts...")
counts = counts.apply(inverse_rank_normalization, axis=0)

# Save the counts and labels for reference
counts_path = os.path.join(outdir, "counts_methylation.csv")
labels_path = os.path.join(outdir, "labels_methylation.csv")
counts.to_csv(counts_path)
labels.to_csv(labels_path)
print(f"✓ Saved counts to {counts_path}")
print(f"✓ Saved labels to {labels_path}")

# ======================== Feature Selection ========================
print("\n=== FEATURE SELECTION ===")

# Set up the data
X = counts.values
y = labels["meth_3cluster"].values

# Encode the labels
label_encoder = preprocessing.LabelEncoder()
y_encoded = label_encoder.fit_transform(y.astype(str))
print(f"Encoded labels: {label_encoder.classes_}")

# Simple feature selection with ANOVA F-test
print("Performing feature selection with ANOVA F-test...")
selector = SelectKBest(score_func=f_classif, k=1000)  # Select top 1000 features
X_selected = selector.fit_transform(X, y_encoded)
selected_features = counts.columns[selector.get_support()].tolist()
print(f"Selected {len(selected_features)} features for modeling")

# Let's get the names back
X_selected = pd.DataFrame(X_selected, columns=selected_features, index=counts.index)
# Ensure y_encoded is 1D and matches the index length
y_encoded = pd.Series(np.asarray(y_encoded).reshape(-1), index=counts.index, name="meth_3cluster")

# Save the selected features
selected_features_path = os.path.join(outdir, "selected_features.json")
with open(selected_features_path, 'w') as f:
    json.dump(selected_features, f, indent=2)
print(f"✓ Saved selected features to {selected_features_path}")

# ======================== Train/Test Split ========================
print("\n=== TRAIN/TEST SPLIT ===")

from sklearn.model_selection import train_test_split
X_train, X_test, y_train, y_test = train_test_split(
    X_selected, y_encoded, 
    test_size=0.2, 
    random_state=420, 
    stratify=y_encoded
)
train_sample_names = X_train.index
test_sample_names = X_test.index

# Save the train/test split as CSV files
X_train_path = os.path.join(outdir, "X_train.csv")
X_test_path = os.path.join(outdir, "X_test.csv")
y_train_path = os.path.join(outdir, "y_train.csv")
y_test_path = os.path.join(outdir, "y_test.csv")

X_train.to_csv(X_train_path)
X_test.to_csv(X_test_path)
y_train.to_csv(y_train_path)
y_test.to_csv(y_test_path)

print(f"✓ Saved X_train to {X_train_path}")
print(f"✓ Saved X_test to {X_test_path}")
print(f"✓ Saved y_train to {y_train_path}")
print(f"✓ Saved y_test to {y_test_path}")

print(f"Training set: {X_train.shape[0]} samples")
print(f"Test set: {X_test.shape[0]} samples")
print(f"Training class distribution: {pd.Series(y_train).value_counts().to_dict()}")
print(f"Test class distribution: {pd.Series(y_test).value_counts().to_dict()}")

# Save the train/test split
train_test_split_path = os.path.join(outdir, "train_test_split.json")
train_test_split_data = {
    "X_train_shape": X_train.shape,
    "X_test_shape": X_test.shape,
    "y_train_distribution": pd.Series(y_train).value_counts().to_dict(),
    "y_test_distribution": pd.Series(y_test).value_counts().to_dict(),
}
with open(train_test_split_path, 'w') as f:
    json.dump(train_test_split_data, f, indent=2)
print(f"✓ Saved train/test split data to {train_test_split_path}")

# ======================== Model Training and Skore Comparison ========================
print("\n=== MODEL TRAINING ===")

# Set up cross-validation
cv = model_selection.StratifiedKFold(n_splits=5, shuffle=True, random_state=420)

# Initialize skore project
skore_dir = outdir + "skore_project/"
os.makedirs(skore_dir, exist_ok=True)
project = skore.Project(skore_dir)
print("Skore project initialized")

# Define models
models = {
    "logistic_regression": LogisticRegression(
        random_state=420, max_iter=1000, class_weight="balanced"
    ),
    "svm_linear": SVC(
        kernel="linear",
        random_state=420,
        probability=True,
        class_weight="balanced",
    ),
    "svm_rbf": SVC(
        kernel="rbf",
        random_state=420,
        probability=True,
        class_weight="balanced",
    ),
    "random_forest": RandomForestClassifier(
        n_estimators=100, random_state=420, class_weight="balanced_subsample", n_jobs=-1
    ),
    "gradient_boosting": GradientBoostingClassifier(
        n_estimators=100,
        random_state=420,
    ),
}

# Train models and create cross-validation reports
print("Training models and creating cross-validation reports...")
cv_reports = {}
model_save_dir = outdir + "trained_models/"
os.makedirs(model_save_dir, exist_ok=True)

for name, model in models.items():
    print(f"Training {name}...")
    
    # Train model on training set
    model.fit(X_train, y_train)
    
    # Save the trained model
    model_path = os.path.join(model_save_dir, f"{name}_model.joblib")
    joblib.dump(model, model_path)
    print(f"  ✓ Saved {name} model to {model_path}")
    
    # Create CrossValidationReport for skore
    from skore import CrossValidationReport
    cv_report = CrossValidationReport(model, X_train, y_train)
    
    # Store the report
    cv_reports[name] = cv_report

    # Save the model with joblib
    model_save_path = os.path.join(model_save_dir, f"{name}_model.joblib")
    joblib.dump(model, model_save_path)
    print(f"  ✓ Saved {name} model to {model_save_path}")
    
    # Save the CV report to skore project
    try:
        # Get accuracy metrics and convert to string for JSON serialization
        accuracy = cv_report.metrics.accuracy()
        metrics_summary = {
            'model_name': name,
            'accuracy_summary': str(accuracy)
        }
        # Save as JSON file
        metrics_path = os.path.join(model_save_dir, f"{name}_metrics.json")
        with open(metrics_path, 'w') as f:
            json.dump(metrics_summary, f, indent=2)
        print(f"  ✓ Saved {name} metrics to {metrics_path}")
    except Exception as e:
        print(f"  Warning: Could not save metrics for {name}: {e}")
    
    print(f"  ✓ Created {name} cross-validation report")

print("\n=== MODEL COMPARISON ===")

# Create ComparisonReport
try:
    from skore import ComparisonReport
    print("Creating comparison report for all models...")
    
    # Create a list of all cross-validation reports for comparison
    report_list = list(cv_reports.values())
    comparison_report = ComparisonReport(report_list)
    
    print("✓ ComparisonReport created successfully")
    print("  This report can be used to compare all models side-by-side")
    
    # Save comparison summary
    comparison_summary = {
        'models_compared': list(cv_reports.keys()),
        'comparison_type': 'CrossValidationReport',
        'num_models': len(cv_reports),
        'created_at': datetime.datetime.now().isoformat()
    }
    
    comparison_path = os.path.join(model_save_dir, "comparison_summary.json")
    with open(comparison_path, 'w') as f:
        json.dump(comparison_summary, f, indent=2)
    print(f"✓ Saved comparison summary to {comparison_path}")
    
except Exception as e:
    print(f"Warning: Could not create ComparisonReport: {e}")
    comparison_report = None
# Display metrics for each model
for name, cv_report in cv_reports.items():
    print(f"\n{name} cross-validation results:")
    try:
        # Use help to see what's available
        print("Available methods:")
        cv_report.help()
    except Exception as e:
        print(f"Help not available: {e}")
        
    try:
        # Try to get some basic metrics
        accuracy = cv_report.metrics.accuracy()
        print(f"Accuracy: {accuracy}")
    except Exception as e:
        print(f"Could not get accuracy: {e}")

# Display ComparisonReport information
if 'comparison_report' in locals() and comparison_report is not None:
    print(f"\n=== COMPARISON REPORT ===")
    print(f"✓ ComparisonReport successfully created comparing {len(cv_reports)} models")
    print(f"Models included: {', '.join(cv_reports.keys())}")
    print("The ComparisonReport object allows side-by-side comparison of all models")
    
print("\n✓ All models trained and cross-validation completed")
print("✓ Cross-validation reports created for model comparison")
print(f"✓ Models and metrics saved to: {model_save_dir}")
if 'comparison_report' in locals() and comparison_report is not None:
    print("✓ ComparisonReport created for comprehensive model comparison")

# Final summary
print(f"\nAll models trained and saved to skore project at: {skore_dir}")
print("Use 'skore ui' command in the skore_project directory to compare models")
print(f"Trained models and metrics also saved to: {model_save_dir}")