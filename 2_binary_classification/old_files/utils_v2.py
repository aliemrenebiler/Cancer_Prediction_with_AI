from sklearn.model_selection import ShuffleSplit
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
def train_and_test_with_classifier(classifier, features, labels):
    ss = ShuffleSplit(n_splits=5, test_size=0.25, random_state=1)

    best_model = None
    best_score = 0

    for train_index, test_index in ss.split(range(len(labels))):
        features_train = [features[i] for i in train_index]
        features_test = [features[i] for i in test_index]
        labels_train = [labels[i] for i in train_index]
        labels_test = [labels[i] for i in test_index]

        model = classifier.fit(features_train, labels_train)
        score = model.score(features_test, labels_test) * 100

        if best_model is None or (best_model is not None and best_score < score):
            best_model = model
            best_score = score

    return best_model, best_score
