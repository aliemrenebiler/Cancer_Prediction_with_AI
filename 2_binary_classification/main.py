from utils import *
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

# Set labels
labels = get_mss_statuses(mdata_file_path)

# Set features
features = get_expressions(exprs_file_path)

# Train And Test With Naive Bayes
model, features_test, labels_test = train_with_classifier(
    GaussianNB(),
    features,
    labels,
)
naive_bayes_success = test_model(
    model,
    features_test,
    labels_test,
)

# Train And Test With Random Forest
model, features_test, labels_test = train_with_classifier(
    RandomForestClassifier(),
    features,
    labels,
)
random_forest_success = test_model(
    model,
    features_test,
    labels_test,
)

# Train And Test With k-NN
model, features_test, labels_test = train_with_classifier(
    KNeighborsClassifier(),
    features,
    labels,
)
knn_success = test_model(
    model,
    features_test,
    labels_test,
)

# Train And Test With Decision Tree
model, features_test, labels_test = train_with_classifier(
    DecisionTreeClassifier(),
    features,
    labels,
)
decision_tree_success = test_model(
    model,
    features_test,
    labels_test,
)

# Train And Test With Support Vector Machine
model, features_test, labels_test = train_with_classifier(
    SVC(),
    features,
    labels,
)
support_vector_success = test_model(
    model,
    features_test,
    labels_test,
)

# Train And Test With Logistic Regression
model, features_test, labels_test = train_with_classifier(
    LogisticRegression(max_iter=1000),
    features,
    labels,
)
logistic_regression_success = test_model(
    model,
    features_test,
    labels_test,
)

print(f"----------------------|---------")
print(f" Classifiers          | Success ")
print(f"----------------------|---------")
print(f" Naive Bayes:         | {round(naive_bayes_success, 3)}")
print(f" Random Forest:       | {round(random_forest_success, 3)}")
print(f" k-Nearest Neighbor:  | {round(knn_success, 3)}")
print(f" Decision Tree:       | {round(decision_tree_success, 3)}")
print(f" Support Vector:      | {round(support_vector_success, 3)}")
print(f" Logistic Regression: | {round(logistic_regression_success, 3)}")
print(f"----------------------|---------")
