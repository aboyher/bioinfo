import sys
import pandas as pd
from sklearn import preprocessing, tree
from sklearn.model_selection import train_test_split
import numpy as np

















test_rati /o = np.float(sys.argv[3])

le = preprocessing.LabelEncoder()
df_train = pd.read_csv(sys.argv[1], index_col=0)
df_train = df_train.transpose()
df_train.to_csv(sys.argv[1].split(".csv")[0] + '_df_train.csv')
exit()
train_labels = pd.read_csv(sys.argv[2], sep="\t", header=None)

df_test = pd.read_csv(sys.argv[1])
df_test = df_test.transpose()
df_test.to_csv(sys.argv[1].split(".csv")[0] + '_df_test.csv')
exit()
test_labels = pd.read_csv(sys.argv[4], sep = '\t', header = None)

print(df_test.shape[0], df_train.shape[0])

df = []

exit()
Y = np.array(train_labels)[:,1]
df_2 = df.apply(le.fit_transform)

enc = preprocessing.OneHotEncoder()
enc.fit(df_2)
onehotlabels = enc.transform(df_2).toarray()
# indices = np.arange(0,onehotlabels.shape[0],1)
#
# X_train, X_test, Y_train, Y_test, i_train, i_test = train_test_split(onehotlabels,Y,indices,test_size=test_ratio)
# print(i_test,i_train)
# clf = tree.DecisionTreeClassifier()
# clf = clf.fit(X_train,Y_train)
# pred = clf.predict(X_test)

# acc = (pred == np.array(labels)[i_test][:,1])
# print(sum(acc)/np.float(len(pred))* 100.)

# testX =

clf = tree.DecisionTreeClassifier()
clf = clf.fit(onehotlabels, Y)
pred = clf.predict()
