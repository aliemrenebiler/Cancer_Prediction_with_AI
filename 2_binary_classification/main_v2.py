from utils_v2 import *
from sklearn.naive_bayes import GaussianNB
from sklearn.ensemble import RandomForestClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.svm import SVC
from sklearn.linear_model import LogisticRegression
import os

# Configurations
input_folder_path = "input"
mdata_file_name = "1_merged_mdata.csv"
exprs_file_name = "3_deg_exprs.csv"

# Set file paths
mdata_file_path = os.path.join(input_folder_path, mdata_file_name)
exprs_file_path = os.path.join(input_folder_path, exprs_file_name)

# Set labels and features
print("Reading CSV files...")
labels = get_mss_statuses(mdata_file_path)
features = get_expressions(exprs_file_path)

# Train And Test With Naive Bayes
print("Training with Naive Bayes...")
nb_model, nb_score = train_and_test_with_classifier(
    GaussianNB(),
    features,
    labels,
)

# Train And Test With Random Forest
print("Training with Random Forest...")
rf_model, rf_score = train_and_test_with_classifier(
    RandomForestClassifier(),
    features,
    labels,
)

# Train And Test With k-NN
print("Training with k-Nearest Neighbor...")
knn_model, knn_score = train_and_test_with_classifier(
    KNeighborsClassifier(),
    features,
    labels,
)

# Train And Test With Decision Tree
print("Training with Decision Tree...")
dt_model, dt_score = train_and_test_with_classifier(
    DecisionTreeClassifier(),
    features,
    labels,
)

# Train And Test With Support Vector Machine
print("Training with Support Vector Machine...")
svm_model, svm_score = train_and_test_with_classifier(
    SVC(),
    features,
    labels,
)

# Train And Test With Logistic Regression
print("Training with Logistic Regression...")
lr_model, lr_score = train_and_test_with_classifier(
    LogisticRegression(max_iter=1000),
    features,
    labels,
)

print("Completed.")

print(f"---------------------|--------")
print(f" Classifiers         | Score  ")
print(f"---------------------|--------")
print(f" Naive Bayes         | {round(nb_score, 1)}%")
print(f" Random Forest       | {round(rf_score, 1)}%")
print(f" k-Nearest Neighbor  | {round(knn_score, 1)}%")
print(f" Decision Tree       | {round(dt_score, 1)}%")
print(f" Support Vector      | {round(svm_score, 1)}%")
print(f" Logistic Regression | {round(lr_score, 1)}%")
print(f"---------------------|--------")
