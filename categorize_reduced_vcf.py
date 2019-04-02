import sys
import pandas as pd
from sklearn import preprocessing, tree
from sklearn.model_selection import train_test_split
import numpy as np

test_ratio = np.float(sys.argv[3])

le = preprocessing.LabelEncoder()
df = pd.read_csv(sys.argv[1])
df = df.transpose()
labels = pd.read_csv(sys.argv[2], sep="\t", header=None)
Y = np.array(labels)[:,1]
df_2 = df.apply(le.fit_transform)

enc = preprocessing.OneHotEncoder()
enc.fit(df_2)
onehotlabels = enc.transform(df_2).toarray()
indices = np.arange(0,onehotlabels.shape[0],1)

X_train, X_test, Y_train, Y_test, i_train, i_test = train_test_split(onehotlabels,Y,indices,test_size=test_ratio)
print(i_test,i_train)
clf = tree.DecisionTreeClassifier()
clf = clf.fit(X_train,Y_train)
pred = clf.predict(X_test)

acc = (pred == np.array(labels)[i_test][:,1])
print(sum(acc)/np.float(len(pred))* 100.)