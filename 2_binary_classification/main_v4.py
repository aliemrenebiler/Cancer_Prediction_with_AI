from utils_v4 import *

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
CONF_MAT_FOLDER_PATH = "output/conf_mat"
ROC_CURVE_FOLDER_PATH = "output/roc_curve"
METADATA_FILE_NAME = "1_merged_mdata.csv"
EXPRESSIONS_FILE_NAME = "3_deg_exprs.csv"

# Shuffle Split Configurations
SPLIT_AMOUNT = 5
TEST_SIZE = 0.25
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
            "n_estimators": [115, 130],
            "min_samples_split": [9, 10],
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

# Create cross validator
print("Creating Shuffle Split Cross Validator...")
shuffle_split_cv = create_shuffle_split_cross_validator(
    n_split=SPLIT_AMOUNT,
    test_size=TEST_SIZE,
    random_state=RANDOM_STATE,
)

# Train and test for each classifier
for i in range(len(CLASSIFIERS)):
    print(f"Training and testing with {CLASSIFIERS[i].name}...")

    # Train
    (
        CLASSIFIERS[i].model,
        CLASSIFIERS[i].best_params,
    ) = train_with_classifier(
        classifier=CLASSIFIERS[i].classifier,
        param_grid=CLASSIFIERS[i].param_grid,
        scoring=SCORING,
        cv=shuffle_split_cv,
        features=features_train,
        labels=labels_train,
    )

    # Test and predict
    (
        labels_pred,
        CLASSIFIERS[i].score,
    ) = test_the_model_and_predict(
        CLASSIFIERS[i].model,
        features_test,
        labels_test,
    )

    # Create confusion matrix and classification report
    (
        CLASSIFIERS[i].conf_mat,
        CLASSIFIERS[i].class_rep,
    ) = create_conf_mat_and_class_report(
        y_true=labels_test,
        y_pred=labels_pred,
    )

    # Save confusion matrix as PNG
    conf_mat_file_name = f"{CLASSIFIERS[i].name.replace(' ', '_').lower()}_conf_mat.png"
    conf_mat_file_path = os.path.join(CONF_MAT_FOLDER_PATH, conf_mat_file_name)
    save_confusion_matrix_plot(CLASSIFIERS[i], conf_mat_file_path)

    # Save ROC curve as PNG
    roc_curve_file_name = (
        f"{CLASSIFIERS[i].name.replace(' ', '_').lower()}_roc_curve.png"
    )
    roc_curve_file_path = os.path.join(ROC_CURVE_FOLDER_PATH, roc_curve_file_name)
    save_roc_curve_plot(CLASSIFIERS[i], labels_test, labels_pred, roc_curve_file_path)

# Print the results
print_results_as_table(CLASSIFIERS)
print("Completed.")
