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

features_train, features_test, labels_train, labels_test = train_test_split(
    features, labels, test_size=0.25, random_state=1
)

gnb = GaussianNB()

model = gnb.fit(features_train, labels_train)

predicted_values = model.predict(features_test)
print(predicted_values)
