from utils import *
import os

# Configurations
input_folder_path = "input"
mdata_file_name = "1_merged_mdata.csv"
exprs_file_name = "3_deg_exprs.csv"

# Set file paths
mdata_file_path = os.path.join(input_folder_path, mdata_file_name)
exprs_file_path = os.path.join(input_folder_path, exprs_file_name)

# Set labels
labels = get_mss_statuses(mdata_file_path)

# Set features
features = get_expressions(exprs_file_path)

# Train With Naive Bayes
model, features_test, labels_test = train_with_naive_bayes(features, labels)

# Test Naive Bayes
naive_bayes_success = test_naive_bayes(model, features_test, labels_test)

print(naive_bayes_success)
