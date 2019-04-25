import numpy as np
from joblib import Parallel, delayed
from sklearn import preprocessing, tree, svm
import multiprocessing

class BagLearner(object):
    def __init__(self, learner, bags):
        self.bags = bags
        self.learners = []
        self.max_threads = multiprocessing.cpu_count()
        for l in range(0, bags):
            self.learners.append(learner())

    def addEvidence(self, dataX, dataY):
        for i in range(0, self.bags):
            n_prime = np.random.choice(len(dataY), len(dataY), replace = True)
            dX = dataX[n_prime,:]
            dY = dataY[n_prime]
            self.learners[i].fit(dX, dY)

    def query(self, dataX):
        pred = []
        for i in range(0, self.bags):
            pred.append(self.learners[i].predict(dataX))
        return pred


if __name__=="__main__":
    print("nada")