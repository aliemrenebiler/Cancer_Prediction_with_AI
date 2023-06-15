from typing import List
from sklearn.model_selection import ShuffleSplit, train_test_split
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import (
    confusion_matrix,
    precision_score,
    recall_score,
    f1_score,
    accuracy_score,
    roc_auc_score,
    RocCurveDisplay,
)
from mlxtend.plotting import plot_confusion_matrix
import matplotlib.pyplot as plt
import numpy as np
import csv


class Classifier:
    def __init__(
        self,
        name,
        classifier,
        param_grid,
        model=None,
        int_score=None,
        ext_score=None,
        best_params=None,
        conf_mat=None,
        class_rep=None,
        roc_curve_auc=None,
    ):
        self.name = name
        self.classifier = classifier
        self.param_grid = param_grid
        self.best_params = best_params
        self.model = model
        self.int_score = int_score
        self.ext_score = ext_score
        self.conf_mat = conf_mat
        self.class_rep = class_rep
        self.roc_curve_auc = roc_curve_auc


# READ CSV FILES
def get_mss_statuses(mdata_file_path):
    # Open file
    file = open(mdata_file_path, "r")
    reader = csv.reader(file, delimiter=",")

    # Get columns and status column index
    headings = next(reader)
    status_column = headings.index("Status")

    # Get MSS statuses as labels
    labels = np.array([1 if row[status_column] == "MSI-H" else 0 for row in reader])

    # Close file
    file.close()

    return labels


def get_expressions(exprs_file_path):
    # Open file
    file = open(exprs_file_path, "r")
    reader = csv.reader(file, delimiter=",")

    # Read headings to pass first row
    next(reader)

    # Get MSS statuses as labels
    features = np.array(
        [
            np.array(
                [float(cell) for cell in row[1:]],
            )
            for row in reader
        ]
    )
    # Close file
    file.close()

    return features


def set_dataset_as_test_group(dataset_name, mdata_file_path, features, labels):
    # Open file
    file = open(mdata_file_path, "r")
    reader = csv.reader(file, delimiter=",")

    # Get columns and status column index
    headings = next(reader)
    dataset_column = headings.index("Dataset")

    # Get test groups' indices
    test_indices = [
        index for index, row in enumerate(reader) if row[dataset_column] == dataset_name
    ]
    train_indices = list(set([*range(len(labels))]) - set(test_indices))

    # Set the groups
    features_train = [features[index] for index in train_indices]
    features_test = [features[index] for index in test_indices]
    labels_train = [labels[index] for index in train_indices]
    labels_test = [labels[index] for index in test_indices]

    # Close file
    file.close()

    return features_train, features_test, labels_train, labels_test


# CREATE TRAIN TEST GROUPS
def create_shuffle_split_cross_validator(n_split, test_size, random_state=1):
    return ShuffleSplit(
        n_splits=n_split,
        test_size=test_size,
        random_state=random_state,
    )


def split_train_and_test_groups(features, labels):
    (
        features_train,
        features_test,
        labels_train,
        labels_test,
    ) = train_test_split(features, labels)

    return (
        features_train,
        features_test,
        labels_train,
        labels_test,
    )


# TRAIN AND TEST FUNCTIONS
def find_best_params_with_grid_search(clf, param_grid, scoring, cv, features, labels):
    # Set grid search classifier
    gs_classifier = GridSearchCV(clf, param_grid=param_grid, scoring=scoring, cv=cv)

    # Train the model
    gs_classifier.fit(features, labels)

    # Get best parameters
    return gs_classifier.best_params_


def train_with_classifier(clf, features, labels, params=None):
    # Set the hyperparameters
    if params is not None:
        clf.set_params(**params)

    # Train the model
    return clf.fit(features, labels)


def test_the_model(model, features, labels):
    return model.score(features, labels) * 100


def predict_with_model(model, features):
    return model.predict(features)


def create_conf_mat_and_class_report(y_true, y_pred):
    # Set confusion matrix dictionary
    tn, fp, fn, tp = confusion_matrix(
        y_true=y_true,
        y_pred=y_pred,
        labels=[0, 1],
    ).ravel()
    conf_mat = {"tn": tn, "fp": fp, "fn": fn, "tp": tp}

    # Set classification report dictionary
    class_rep = {
        "precision": precision_score(
            y_true=y_true,
            y_pred=y_pred,
        ),
        "recall": recall_score(
            y_true=y_true,
            y_pred=y_pred,
        ),
        "f1_score": f1_score(
            y_true=y_true,
            y_pred=y_pred,
        ),
        "accuracy": accuracy_score(
            y_true=y_true,
            y_pred=y_pred,
        ),
    }

    return conf_mat, class_rep


# PRINT TABLE
def save_confusion_matrix_plot(clf_class: Classifier, saved_file_path):
    # Set confusion matrix
    conf_mat = np.array(
        [
            [clf_class.conf_mat["tp"], clf_class.conf_mat["fp"]],
            [clf_class.conf_mat["fn"], clf_class.conf_mat["tn"]],
        ]
    )

    # Plot the matrix
    fig, ax = plot_confusion_matrix(
        conf_mat=conf_mat,
        cmap=plt.cm.Greys,
        show_absolute=True,
        show_normed=True,
        colorbar=True,
    )
    plt.xlabel("Predictions")
    plt.ylabel("Actuals")
    plt.title(f"{clf_class.name} Confusion Matrix")

    # Save matrix as PNG file
    plt.savefig(saved_file_path)
    plt.clf()


def save_roc_curve_plot(clf_class: Classifier, features, labels, saved_file_path):
    # Find model result probabilities
    y_proba = clf_class.model.predict_proba(features)

    # Plot the ROC curve
    RocCurveDisplay.from_predictions(labels, y_proba[:, 1])

    # Save matrix as PNG file
    plt.title(f"{clf_class.name} ROC Curve")
    plt.savefig(saved_file_path)
    plt.clf()

    # Calculate AUC
    return roc_auc_score(labels, y_proba[:, 1])


def print_results_as_table(classifiers: List[Classifier]):
    # Set headings
    headings = [
        "Classifier",
        "Int. Score",
        "Ext. Score",
        "Accuracy",
        "Precision",
        "Recall",
        "F1 Score",
        "AUC",
        "Best Parameters",
    ]

    # Set max lenghts of a column
    max_lengths = [
        max([len(clf.name) for clf in classifiers] + [len(headings[0])]),
        10,
        10,
        10,
        10,
        10,
        10,
        10,
        max([len(str(clf.best_params)) for clf in classifiers] + [len(headings[6])]),
    ]

    # Print row function
    def print_row(
        cell_texts,
        cell_max_lenghts,
    ):
        for i in range(len(cell_texts)):
            cell_texts[
                i
            ] = f" {cell_texts[i]}{' '*(cell_max_lenghts[i]-len(cell_texts[i]))} "
        row = "|".join(cell_texts)
        print(row)

    # Print seperator function
    def print_seperator(
        cell_max_lenghts,
    ):
        row = "|".join(
            ["-" * (cell_max_length + 2) for cell_max_length in cell_max_lenghts]
        )
        print(row)

    print_seperator(max_lengths)
    print_row(headings, max_lengths)
    print_seperator(max_lengths)
    for clf in classifiers:
        print_row(
            [
                clf.name,
                f"{str(round(clf.int_score, 6))}%",
                f"{str(round(clf.ext_score, 6))}%",
                str(round(clf.class_rep["accuracy"], 8)),
                str(round(clf.class_rep["precision"], 8)),
                str(round(clf.class_rep["recall"], 8)),
                str(round(clf.class_rep["f1_score"], 8)),
                str(round(clf.roc_curve_auc, 8)),
                str(clf.best_params),
            ],
            max_lengths,
        )
    print_seperator(max_lengths)
