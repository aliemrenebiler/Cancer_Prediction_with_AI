from sklearn.model_selection import train_test_split
import csv
import numpy as np


# READ CSV FILES
def get_mss_statuses(mdata_file_path):
    # Open file
    file = open(mdata_file_path, "r")
    reader = csv.reader(file, delimiter=",")

    # Get columns and status column index
    headings = next(reader)
    status_column = headings.index("Status")

    # Get MSS statuses as labels
    labels = np.array(
        ["MSI-H" if row[status_column] == "MSI-H" else "MSS" for row in reader]
    )

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


# TRAIN AND TEST FUNCTIONS
def train_with_classifier(classifier, features, labels):
    features_train, features_test, labels_train, labels_test = train_test_split(
        features, labels, test_size=0.25, random_state=1
    )

    model = classifier.fit(features_train, labels_train)

    return model, features_test, labels_test


def test_model(model, features_test, labels_test):
    return model.score(features_test, labels_test) * 100
