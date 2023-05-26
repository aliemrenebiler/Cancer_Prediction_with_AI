from utils_v3 import *

from sklearn.naive_bayes import GaussianNB
from sklearn.ensemble import RandomForestClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.svm import SVC
from sklearn.linear_model import LogisticRegression
import numpy as np
import os

# Input Configurations
INPUT_FOLDER_PATH = "input"
METADATA_FILE_NAME = "1_merged_mdata.csv"
EXPRESSIONS_FILE_NAME = "3_deg_exprs.csv"

# Shuffle Split Configurations
SPLIT_AMOUNT = 5
TEST_SIZE = 0.5
RANDOM_STATE = 1

# Grid Search Configurations
SCORING = "accuracy"

# Classifier Configurations
CLASSIFIERS = [
    Classifier(
        name="Naive Bayes",
        classifier=GaussianNB(),
        param_grid={
            "var_smoothing": np.logspace(0, -4, 100),
        },
    ),
    Classifier(
        name="Support Vector",
        classifier=SVC(),
        param_grid={
            "C": [0.1, 1, 10],
        },
    ),
    Classifier(
        name="Decision Tree",
        classifier=DecisionTreeClassifier(),
        param_grid={
            "min_samples_split": [2, 3, 4],
        },
    ),
    Classifier(
        name="Random Forest",
        classifier=RandomForestClassifier(),
        param_grid={
            "n_estimators": [100, 200],
            "min_samples_split": [2, 3, 4],
        },
    ),
    Classifier(
        name="k-Nearest Neighbor",
        classifier=KNeighborsClassifier(),
        param_grid={
            "n_neighbors": list(range(1, 31)),
            "weights": ["uniform", "distance"],
        },
    ),
    Classifier(
        name="Logistic Regression",
        classifier=LogisticRegression(),
        param_grid={
            "C": [0.01, 0.1, 1, 10, 100],
            "max_iter": [1000],
            "penalty": ["l2"],
        },
    ),
]

# Set file paths
mdata_file_path = os.path.join(INPUT_FOLDER_PATH, METADATA_FILE_NAME)
exprs_file_path = os.path.join(INPUT_FOLDER_PATH, EXPRESSIONS_FILE_NAME)

# Set labels and features
print("Reading CSV files...")
labels = get_mss_statuses(mdata_file_path)
features = get_expressions(exprs_file_path)

# Create train and test groups
print("Creating Shuffle Split Cross Validator...")
shuffle_split_cv = create_shuffle_split_cross_validator(
    n_split=SPLIT_AMOUNT,
    test_size=TEST_SIZE,
    random_state=RANDOM_STATE,
)

# Train and test for each classifier
for i in range(len(CLASSIFIERS)):
    print(f"Training and testing with {CLASSIFIERS[i].name}...")
    (
        CLASSIFIERS[i].model,
        CLASSIFIERS[i].score,
        CLASSIFIERS[i].best_params,
    ) = train_and_test_with_classifier(
        classifier=CLASSIFIERS[i].classifier,
        param_grid=CLASSIFIERS[i].param_grid,
        scoring=SCORING,
        cv=shuffle_split_cv,
        features=features,
        labels=labels,
    )

# Print the results
print_results_as_table(CLASSIFIERS)
print("Completed.")
