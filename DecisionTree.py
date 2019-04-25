import numpy as np
from collections import Counter
# import ast
from DT_Utils import *
import multiprocessing as mp

class DecisionTree(object):
    def __init__(self):
        self.tree = []
        self.breaks = 0
        self.pool = mp.Pool(mp.cpu_count())

    def learn(self, X, y):

        def getBestFeature(X,y):
            best_ig = 0
            split_feature = 0
            split_val = 0


            for col in range(len(X[0])):
                # print(col)
                x = X[:, col]
                values = set(X[:, col].tolist())
                for val in values:
                    # xl, xr, yl, yr = partition_classes(X, y, col, val)
                    idx = X[:, col] == val
                    lidx = [i for i,x in enumerate(idx) if x]
                    ridx = [i for i,x in enumerate(idx) if not x]
                    yl = y[lidx]
                    yr = y[ridx]
                    # if any(test) or not any(test):
                    if len(yl) < 1 or len(yr) < 1:
                        continue
                    ig = information_gain(y, [yl, yr])
                    if ig > best_ig:
                        best_ig = ig
                        split_val = val
                        split_feature = col


            return split_feature, split_val

        def build_tree(X,y):
            self.breaks += 1
            # print(self.breaks)
            if int(len(X)) == 1 or (y[0] == y).all():
                return np.array([[None, y[0], None, None]])
            elif np.all(X == X[0, :]):
                counts = Counter(y)
                most_common = counts.most_common(1)[0][0]
                return np.array([[None, most_common, None, None]])
            split_feature, split_val = getBestFeature(X,y)
            Xl, Xr, yl, yr = partition_classes(X, y, split_feature, split_val)
            if Xl.size == 0 or Xr.size == 0:
                return([[None, np.argmax(np.unique(y, return_counts = True)[1]), None, None]])
            left_tree = build_tree(Xl, yl)
            right_tree = build_tree(Xr, yr)
            root = [[split_feature, split_val, 1, len(left_tree) + 1]]
            return (np.concatenate((root, left_tree, right_tree), axis = 0))

        X = np.array(X)
        y = np.array(y)
        self.tree = build_tree(X, y)

    def classify(self, record):
        node = 0
        x = self.tree[node, 0]
        while self.tree[node, 0] is not None:
            a = record[int(self.tree[node, 0])]
            b = self.tree[node, 1]
            if record[int(self.tree[node, 0])] <= self.tree[node, 1]:
                c = self.tree[node,2]
                node += int(self.tree[node,2])
            else:
                c = self.tree[node, 3]
                node += int(self.tree[node, 3])
        pred = self.tree[node, 1]
        return pred


