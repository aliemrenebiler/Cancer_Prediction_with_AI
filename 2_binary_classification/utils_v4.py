from sklearn.model_selection import ShuffleSplit
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import (
    confusion_matrix,
    precision_score,
    recall_score,
    f1_score,
    accuracy_score,
    roc_curve,
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
        score=None,
        best_params=None,
        conf_mat=None,
        class_rep=None,
    ):
        self.name = name
        self.classifier = classifier
        self.param_grid = param_grid
        self.best_params = best_params
        self.model = model
        self.score = score
        self.conf_mat = conf_mat
        self.class_rep = class_rep


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
def create_shuffle_split_cross_validator(
    n_split,
    test_size,
    random_state=1,
):
    return ShuffleSplit(
        n_splits=n_split,
        test_size=test_size,
        random_state=random_state,
    )


# TRAIN AND TEST FUNCTIONS
def train_with_classifier(
    classifier,
    param_grid,
    scoring,
    cv,
    features,
    labels,
):
    # Set grid search classifier
    gs_classifier = GridSearchCV(
        classifier, param_grid=param_grid, scoring=scoring, cv=cv
    )

    # Train the model
    model = gs_classifier.fit(features, labels)

    # Get best parameters
    best_params = gs_classifier.best_params_

    return model, best_params


def test_the_model_and_predict(
    model,
    features,
    labels,
):
    # Predict
    predictions = model.predict(features)

    # Test and set score as precent
    score = model.score(features, labels) * 100

    return predictions, score


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
def save_confusion_matrix_plot(classifier: Classifier, saved_file_path):
    # Set confusion matrix
    conf_mat = np.array(
        [
            [classifier.conf_mat["tp"], classifier.conf_mat["fp"]],
            [classifier.conf_mat["fn"], classifier.conf_mat["tn"]],
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
    plt.title(f"{classifier.name} Confusion Matrix")

    # Save matrix as PNG file
    plt.savefig(saved_file_path)
    plt.clf()


def save_roc_curve_plot(classifier: Classifier, y_true, y_pred, saved_file_path):
    # Plot the ROC curve
    fpr, tpr, _ = roc_curve(y_true, y_pred)
    plt.plot(fpr, tpr)
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    plt.title(f"{classifier.name} ROC Curve")

    # Save matrix as PNG file
    plt.savefig(saved_file_path)
    plt.clf()


def print_results_as_table(classifiers):
    # Set headings
    headings = [
        "Classifier",
        "Score",
        "Precision",
        "Recall",
        "F1 Score",
        "Accuracy",
        "Best Parameters",
    ]

    # Set max lenghts of a column
    max_lengths = [
        max([len(clf.name) for clf in classifiers] + [len(headings[0])]),
        6,
        len(headings[2]),
        len(headings[3]),
        len(headings[4]),
        len(headings[5]),
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
                f"{str(round(clf.score, 2))}%",
                str(round(clf.class_rep["precision"], 3)),
                str(round(clf.class_rep["recall"], 3)),
                str(round(clf.class_rep["f1_score"], 3)),
                str(round(clf.class_rep["accuracy"], 3)),
                str(clf.best_params),
            ],
            max_lengths,
        )
    print_seperator(max_lengths)
