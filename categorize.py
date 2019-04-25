import sys
from matplotlib import pyplot as plt
import pandas as pd
from sklearn import preprocessing, tree
from sklearn.feature_selection import VarianceThreshold
from sklearn.model_selection import train_test_split
import numpy as np
import warnings
from matplotlib.colors import ListedColormap
warnings.simplefilter(action='ignore',category=FutureWarning)


def OHE(data):
    df = data.apply(le.fit_transform)
    enc = preprocessing.OneHotEncoder()
    enc.fit(df)
    data_OHE = enc.transform(df).toarray()
    return data_OHE

def full_geno(df, df2, clf):

    df_OHE = OHE(df)
    print(df_OHE)
    # sel = VarianceThreshold(threshold=(.8 * (1 - .8)))
    # df_OHE2 = sel.fit_transform(df_OHE)
    results = []

    train_idx = np.where(np.array(df2)[:,1] == "rnd1")
    test_idx = np.where(np.array(df2)[:,1] == "rnd3")
    trainX = np.array(df_OHE)[train_idx][:,:]
    trainY = np.array(df2)[train_idx][:,0]
    testX = np.array(df_OHE)[test_idx][:,:]
    testY = np.array(df2)[test_idx][:,0]

    X_train, X_test, Y_train, Y_test = train_test_split(trainX, trainY, test_size = .33)
    for i in [X_train, X_test,Y_train, Y_test]:
        print(i.shape)

    clf = clf.fit(X_train, Y_train)
    pred = clf.predict(X_test)

    acc = (pred == np.array(Y_test))
    results.append(sum(acc)/float(len(pred)) * 100)

    clf = tree.DecisionTreeClassifier()
    clf = clf.fit(trainX, trainY)
    pred = clf.predict(testX)
    acc = (pred == np.array(testY))
    results.append(sum(acc)/float(len(pred)) * 100)

    X_train, X_test, Y_train, Y_test = train_test_split(testX, testY, test_size=.33)

    clf = clf.fit(X_train, Y_train)
    pred = clf.predict(X_test)
    acc = (pred == np.array(Y_test))
    results.append(sum(acc) / float(len(pred)) * 100)

    return results

def hap0(df, df2, clf):
    df3 = []
    for row in np.array(df):
        new_row = []
        for n in row:
            new_row.append(n.split('/')[0])
        df3.append(new_row)
    df3 = pd.DataFrame(data = df3, index=None, columns=None)

    df_OHE = OHE(df3)

    train_idx = np.where(np.array(df2)[:, 1] == "train")
    test_idx = np.where(np.array(df2)[:, 1] == "test")
    trainX = np.array(df_OHE)[train_idx][:, :]
    trainY = np.array(df2)[train_idx][:, 0]
    testX = np.array(df_OHE)[test_idx][:, :]
    testY = np.array(df2)[test_idx][:, 0]

    X_train, X_test, Y_train, Y_test = train_test_split(trainX, trainY, test_size=.33)

    results = []

    clf = clf.fit(X_train, Y_train)
    pred = clf.predict(X_test)
    acc = (pred == np.array(Y_test))
    results.append(sum(acc) / float(len(pred)) * 100)

    clf = tree.DecisionTreeClassifier()
    clf = clf.fit(trainX, trainY)
    pred = clf.predict(testX)
    acc = (pred == np.array(testY))
    results.append(sum(acc) / float(len(pred)) * 100)

    X_train, X_test, Y_train, Y_test = train_test_split(testX, testY, test_size=.33)

    clf = clf.fit(X_train, Y_train)
    pred = clf.predict(X_test)
    acc = (pred == np.array(Y_test))
    results.append(sum(acc) / float(len(pred)) * 100)

    return results

le = preprocessing.LabelEncoder()

df = pd.read_csv(sys.argv[1], index_col=0)

df2 = df[['Phenotype', 'Round']]
df = df.drop(['Phenotype', 'Round'], axis=1)

results_hap = []
results_full = []

clf = tree.DecisionTreeClassifier()

print(full_geno(df,df2,clf))

# for i in range(2,30):
#     # print(i)
#     clf = tree.DecisionTreeClassifier(min_samples_leaf=i, criterion='entropy')
#     # print("hap0")
#     results_hap.append(hap0(df,df2, clf))
#     # print("full")
#     results_full.append(full_geno(df, df2, clf))

# results_hap = np.transpose(results_hap)
# results_full = np.transpose(results_full)
# print(results_full)
# print(results_hap)
#
# plt.plot(results_hap[0])
# plt.plot(results_hap[1])
# plt.plot(results_full[0])
# plt.plot(results_full[1])
# plt.show()


