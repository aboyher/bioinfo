import numpy as np

class DTLearner(object):
    def __init__(self, leaf_size=1):
        self.leaf_size = leaf_size
    def addEvidence(self, data):
        def build_tree(data):
            if int(data.shape[0]) <= self.leaf_size:
                return np.array([[np.NaN, np.mean(data[:, -1]), np.NaN, np.NAN]])
            elif (data[0, -1] == data[:, -1]).all():
                return np.array([[np.NaN, data[0, -1], np.NAN, np.NAN]])

            feat = get_corr(data)

            splitval = np.median(data[:, feat])
            if data.shape[0] == data[data[:, feat] <= splitval].shape[0]:
                return np.array([[np.NAN, np.mean(data[:, -1]), np.NAN, np.NAN]])
            left_tree = build_tree(data[data[:, feat] <= splitval])
            right_tree = build_tree(data[data[:, feat] > splitval])

            root = [[feat, splitval, 1, int(left_tree.shape[0]) + 1]]

            return (np.concatenate((root, left_tree, right_tree), axis=0))

        def get_corr(data):
            print(data)
            corrs = [np.corrcoef(data[:, i], data[:, -1])[0, 1] for i in range(0, int(data.shape[1] - 1))]
            feat = np.argmax(corrs)
            return np.argmax(corrs)

        self.model = build_tree(data)

    def query(self, data):
        pred = []
        for i in range(0, data.shape[0]):
            node = 0
            x = self.model[node, 0]
            while np.isfinite(self.model[node, 0]):
                a = data[i, int(self.model[node,0])]
                b = self.model[node,1]
                if data[i, int(self.model[node,0])] <= self.model[node,1]:
                    c = self.model[node, 2]
                    node += int(self.model[node,2])
                else:
                    c = self.model[node,3]
                    node+= int(self.model[node,3])
            val = self.model[node, 1]
            pred.append(val)
        return pred