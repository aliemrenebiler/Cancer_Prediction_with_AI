from ast import List
from pyclbr import Class
from sklearn.model_selection import ShuffleSplit, train_test_split
from sklearn.model_selection import GridSearchCV
import csv
import numpy as np


class Classifier:
    def __init__(
        self,
        name,
        classifier,
        param_grid,
        model=None,
        score=None,
        best_params=None,
    ):
        self.name = name
        self.classifier = classifier
        self.param_grid = param_grid
        self.best_params = best_params
        self.model = model
        self.score = score


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
def train_and_test_with_classifier(
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

    # Test and set score as precent
    score = model.score(features, labels) * 100

    # Get best parameters
    best_params = gs_classifier.best_params_

    return model, score, best_params


# PRINT TABLE
def print_results_as_table(classifiers):
    name_title = "Classifier"
    score_title = "Score"
    params_title = "Best Parameters"
    max_name_length = max([len(clf.name) for clf in classifiers] + [len(name_title)])
    max_score_length = 6
    max_param_length = max(
        [len(str(clf.best_params)) for clf in classifiers] + [len(params_title)]
    )

    def print_row(
        cell_texts: List,
        cell_max_lenghts: List,
    ):
        for i in range(len(cell_texts)):
            cell_texts[
                i
            ] = f" {cell_texts[i]}{' '*(cell_max_lenghts[i]-len(cell_texts[i]))} "
        row = "|".join(cell_texts)
        print(row)

    def print_seperator(
        cell_max_lenghts: List,
    ):
        row = "|".join(
            ["-" * (cell_max_length + 2) for cell_max_length in cell_max_lenghts]
        )
        print(row)

    print_seperator([max_name_length, max_score_length, max_param_length])
    print_row(
        [name_title, score_title, params_title],
        [max_name_length, max_score_length, max_param_length],
    )
    print_seperator([max_name_length, max_score_length, max_param_length])
    for clf in classifiers:
        print_row(
            [
                clf.name,
                f"{str(round(clf.score, 2))}%",
                str(clf.best_params),
            ],
            [max_name_length, max_score_length, max_param_length],
        )
    print_seperator([max_name_length, max_score_length, max_param_length])
