from sklearn.ensemble import RandomForestClassifier
from utils_v4 import *
import os

# Configurations
INPUT_FOLDER_PATH = "input"
METADATA_FILE_NAME = "1_merged_mdata.csv"
EXPRESSIONS_FILE_NAME = "3_deg_exprs.csv"
CLASSIFIER = RandomForestClassifier(
    min_samples_split=2,
    n_estimators=200,
)

# Set file paths
mdata_file_path = os.path.join(INPUT_FOLDER_PATH, METADATA_FILE_NAME)
exprs_file_path = os.path.join(INPUT_FOLDER_PATH, EXPRESSIONS_FILE_NAME)

# Set labels and features
print("Reading CSV files...")
labels = get_mss_statuses(mdata_file_path)
features = get_expressions(exprs_file_path)

# Seperate the selected dataset as test group
(
    features_train,
    features_test,
    labels_train,
    labels_test,
) = set_dataset_as_test_group(
    "GSE35566",
    mdata_file_path,
    features,
    labels,
)

# Train the model
print("Training...")
model = CLASSIFIER.fit(features_train, labels_train)

# Test the model
score = model.score(features_test, labels_test) * 100

print("Completed.")
print(f"Score: {score}")
