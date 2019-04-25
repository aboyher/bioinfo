from scipy import stats
import numpy as np

def entropy(class_y):

    counts = np.unique(class_y, return_counts = True)[1]
    entropy = stats.entropy(counts, base = 2)
    return entropy


def partition_classes(X, y, split_attribute, split_val):

    left_idx = np.where(X[:, split_attribute] == split_val)
    right_idx = np.where(X[:, split_attribute] != split_val)

    X_left = X[left_idx]
    X_right = X[right_idx]
    y_left = y[left_idx]
    y_right = y[right_idx]
    return (X_left, X_right, y_left, y_right)


def information_gain(previous_y, current_y):

    try:
        H = entropy(previous_y)
        Hl = entropy(current_y[0])
        Hr = entropy(current_y[1])
    except:
        print(len(previous_y), len(current_y[0]), len(current_y[1]))
        exit()
    Pl = len(current_y[0]) / (len(current_y[0])+ len(current_y[1]))
    Pr = len(current_y[1]) / (len(current_y[0])+ len(current_y[1]))
    ig = H - (Hl * Pl + Hr * Pr)

    info_gain = ig

    return info_gain
