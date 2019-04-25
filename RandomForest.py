from DecisionTree import DecisionTree
import pandas as pd
import numpy as np
import ast
import multiprocessing as mp

class RandomForest(object):
    num_trees = 0
    decision_trees = []
    bootstraps_datasets = []
    bootstraps_labels = []

    def __init__(self, num_trees):
        self.num_trees = num_trees
        self.decision_trees = [DecisionTree() for i in range(num_trees)]
        self.pool = mp.Pool(mp.cpu_count())

    def _bootstrapping(selfself, XX, n):
        samples = []
        labels = []
        # XX = np.array(XX)
        idx = np.random.choice(XX.shape[0] - 1, n, replace = True)
        samples = XX[idx, : -1]
        labels = XX[idx, -1]
        return (samples, labels)

    def bootstrapping(self, XX):
        for i in range(self.num_trees):
            data_sample, data_label = self._bootstrapping(XX, len(XX))
            self.bootstraps_datasets.append(data_sample)
            self.bootstraps_labels.append(data_label)

    def fitting(self):

        for i in range(self.num_trees):
            print("Fitting {}".format(i + 1))
            dataset, labels = self.bootstraps_datasets[i], self.bootstraps_labels[i]
            self.decision_trees[i].learn(self.bootstraps_datasets[i], self.bootstraps_labels[i])
            print("Bootstrap set {} fit".format(i+1))
            print("Nodes: {}".format(self.decision_trees[i].tree.shape))

    def voting(self, X):
        y = []
        for record in X:
            votes = []
            in_bag_count = 0
            for i in range(len(self.bootstraps_datasets)):
                dataset = self.bootstraps_datasets[i]
                if record not in dataset:
                    OOB_tree = self.decision_trees[i]
                    effective_vote = OOB_tree.classify(record)
                    votes.append(effective_vote)

            counts = np.bincount(votes)

            if len(counts) == 0:
                for DT in self.decision_trees:
                    votes.append(DT.classify(record))
                counts = np.bincount(votes)
                y = np.append(y, np.argmax(counts))


        return y


def main():
    data = np.array(pd.read_csv("/Users/aboyher/danforth/projects/gbs/round3/all_df.csv", index_col=0))
    rnd1 = np.array(data)[np.where(data[:,-1] == "rnd1")][:,:-1]
    X = rnd1[:,:-1]
    y = rnd1[:,-1]
    y[[i for i,x in enumerate(y=='R') if x]] = 0
    y[[i for i,x in enumerate(y=='S') if x]] = 1
    XX = rnd1
    forest_size = 10
    randomForest = RandomForest(forest_size)
    print("Bootstrapping")
    randomForest.bootstrapping(XX)
    print("Bootstrapping done\nStarting fit\n")
    randomForest.fitting()
    print("Fitting done\nPredicting")
    y_predicted = randomForest.voting(X)
    print("Prediction done\nResults\n")
    results = [prediction == truth for prediction, truth in zip(y_predicted, y)]
    accuracy = float(results.count(True)) / float(len(results))
    print("accuracy: %.4f" % accuracy)
    print("OOB estimate: %4f" % (1 - accuracy))
    #
    exit()
    print("Testing Round 3 dataset")
    rnd3 = np.array(data)[np.where(data[:, -1] == "rnd3")][:,:-1]
    X = rnd3[:,:-1]
    y = rnd3[:, -1]
    y[[i for i, x in enumerate(y == 'R') if x]] = 0
    y[[i for i, x in enumerate(y == 'S') if x]] = 1
    y_predicted = randomForest.voting(X)
    print("actual: {}".format(y))
    print("predicted: {}",y_predicted)
    results = [prediction == truth for prediction, truth in zip(y_predicted, y)]
    accuracy = float(results.count(True)) / float(len(results))
    print("accuracy: %4f" % accuracy)
    print("OOB estimate: %4f" % (1-accuracy))

if __name__ == "__main__":
    main()