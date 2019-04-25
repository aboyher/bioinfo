import pandas as pd
import sys
from sklearn import preprocessing, tree, svm
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score
import numpy as np
import warnings
warnings.simplefilter(action='ignore',category=FutureWarning)
# warnings.simplefilter(action='ignore',category=Conv)
import BagLearner as bl
from collections import Counter
import matplotlib.pyplot as plt


le = preprocessing.LabelEncoder()

df = pd.read_csv(sys.argv[1], index_col=0)
df2 = df[['Phenotype','Round']]
df = df.drop(['Phenotype','Round'], axis = 1)
df = df.apply(le.fit_transform)
enc = preprocessing.OneHotEncoder()
enc.fit(df)
df_OHE = enc.transform(df).toarray()

df = pd.DataFrame(data = df_OHE, index=df.index)
df = df.join(df2)


df_train = df[df['Round'] == 'rnd1']
df_trainY = df_train[['Phenotype']]
df_trainX = df_train.drop(['Phenotype','Round'], axis = 1)
df_rnd3 = df[df['Round'] == 'rnd3']
# print(df_trainX.shape, df_trainY.shape)
# exit()


prog_dic = {}
# prog_lab = {}
for prog in df_trainY.index.values:
    prog_dic[prog] = []
    # prog_lab[prog] = df_trainY[prog]



# learner = bl.BagLearner(learner = svm.LinearSVC, bags = 10)
# high_acc_prog = ["P000027","P000071","P000079","P000166","P000171","P000173","P000186","P000202","P000209","P000215","P000221","P000222","P000223","P000226","P000769","P000234","P000235","P000247","P000255","P000256","P000261","P000265","P000266","P000274","P000278","P000309","P000506","P000534","P000539","P000578","P000594","P000687","P000747","P000764","P000778","P000783","P000804","P000839","P000871","P000878","P000912","P000921","P001198","P001200","P001255","P001259","P001261","P001270"]
#
#
# trainX = df_trainX[df_trainX.index.isin(high_acc_prog)]
# trainY = df_trainY[df_trainY.index.isin(high_acc_prog)]
# testX = df_trainX[~df_trainX.index.isin(high_acc_prog)]
# testY = df_trainY[~df_trainY.index.isin(high_acc_prog)]
#
# print(trainY)
# # exit()
#
# learner.addEvidence(np.array(trainX), np.array(trainY).ravel())
# predY = learner.query(np.array(trainY))
# predY = np.transpose(predY)
#
# results = []
# results2 = []
# for i in range(len(predY)):
#     r = Counter(predY[i].tolist()).most_common(1)[0][0]
#     results2.append(r == np.array(testY)[i])
#     res = []
#     for j in range(len(predY[i])):
#         res.append(predY[i][j]==np.array(testY)[i])
#     results.append(sum(res) / float(len(res)) * 100)
#
# for prog, res in zip(testY.index.values, results):
#     prog_dic[prog].append(res)
#
# print("Iteration {}: {}%".format(z,sum(results2) / float(len(results2)) * 100))
#
# for x in prog_dic:
#     print("{}\t{}\t{}".format(x,np.average(prog_dic[x]), len(prog_dic[x])))



for z in range(10):
    trainX, testX, trainY, testY = train_test_split(df_trainX, df_trainY, test_size= 0.8)
    learner = bl.BagLearner(learner = svm.LinearSVC, bags = 10)
    learner.addEvidence(np.array(trainX), np.array(trainY).ravel())
    predY = learner.query(np.array(testX))
    predY = np.transpose(predY)
    results = []
    results2 = []

    for i in range(len(predY)):
        res = []
        for j in range(len(predY[i])):
            res.append(predY[i][j] == np.array(testY)[i])
        results.append(sum(res) / float(len(res)) * 100)

    # print(len(testY.index.values), len(results))
    for prog, res in zip(testY.index.values, results):
        prog_dic[prog].append(res)



    for i in range(len(predY)):

        r = Counter(predY[i].tolist()).most_common(1)[0][0]
        results2.append(r == np.array(testY)[i])

    print("Iteration {}: {}%".format(z,sum(results2) / float(len(results2)) * 100))

for x in prog_dic:
    print("{}\t{}\t{}".format(x,np.average(prog_dic[x]), len(prog_dic[x])))