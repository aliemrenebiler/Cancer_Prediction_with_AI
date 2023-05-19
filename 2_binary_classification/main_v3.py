from pyclbr import Class
from utils_v3 import *
from sklearn.naive_bayes import GaussianNB
from sklearn.ensemble import RandomForestClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.svm import SVC
from sklearn.linear_model import LogisticRegression
import numpy as np
import os

# Configurations
input_folder_path = "input"
mdata_file_name = "1_merged_mdata.csv"
exprs_file_name = "3_deg_exprs.csv"
classifiers = [
    Classifier(
        name="Naive Bayes",
        classifier=GaussianNB(),
        param_grid={
            "var_smoothing": [1, 0.5, 0.1, 0.05, 0.01, 0.005, 0.001],
        },
        best_params={},
        model=None,
        score=None,
    ),
    Classifier(
        name="Random Forest",
        classifier=RandomForestClassifier(),
        param_grid={
            "n_estimators": [100, 200],
            "min_samples_split": [2, 3, 4],
        },
        best_params={},
        model=None,
        score=None,
    ),
    Classifier(
        name="k-Nearest Neighbor",
        classifier=KNeighborsClassifier(),
        param_grid={
            "n_neighbors": list(range(1, 31)),
            "weights": ["uniform", "distance"],
        },
        best_params={},
        model=None,
        score=None,
    ),
    Classifier(
        name="Decision Tree",
        classifier=DecisionTreeClassifier(),
        param_grid={
            "min_samples_split": [2, 3, 4],
        },
        best_params={},
        model=None,
        score=None,
    ),
    Classifier(
        name="Support Vector",
        classifier=SVC(),
        param_grid={
            "C": [0.1, 1, 10],
            "gamma": [1, 0.1, 0.01],
        },
        best_params={},
        model=None,
        score=None,
    ),
    Classifier(
        name="Logistic Regression",
        classifier=LogisticRegression(),
        param_grid={
            "C": [0.001, 0.01, 0.1, 1, 10, 100, 1000],
            "max_iter": [1000],
            "penalty": ["l2"],
        },
        best_params={},
        model=None,
        score=None,
    ),
]

# Set file paths
mdata_file_path = os.path.join(input_folder_path, mdata_file_name)
exprs_file_path = os.path.join(input_folder_path, exprs_file_name)


# Set labels and features
print("Reading CSV files...")
labels = get_mss_statuses(mdata_file_path)
features = get_expressions(exprs_file_path)

# Create train and test groups
print("Creating Shuffle Split Cross Validator...")
shuffle_split_cv = create_shuffle_split_cross_validator(
    n_split=5,
    test_size=0.25,
    random_state=1,
)

# Train and test for each classifier
for i in range(len(classifiers)):
    print(f"Training and testing with {classifiers[i].name}...")
    (
        classifiers[i].model,
        classifiers[i].score,
        classifiers[i].best_params,
    ) = train_and_test_with_classifier(
        classifier=classifiers[i].classifier,
        param_grid=classifiers[i].param_grid,
        scoring="accuracy",
        cv=shuffle_split_cv,
        features=features,
        labels=labels,
    )

# Print the results
print_results_as_table(classifiers)
print("Completed.")
