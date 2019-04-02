import sys
import pandas as pd
from sklearn import preprocessing, tree
from sklearn.model_selection import train_test_split
import numpy as np


# df = pd.read_csv('/Users/aboyher/danforth/projects/gbs/round3/rnd3_reduced_allele.csv', index_col = 0)
# df = df.transpose()
# df.to_csv('/Users/aboyher/danforth/projects/gbs/round3/rnd3_reduced_allele_df_test.csv')

# df1 = pd.read_csv(sys.argv[1], index_col=0)
# df2 = pd.read_csv(sys.argv[2], index_col=0)
#
# df1 = df1.transpose()
# df2 = df2.transpose()
#
# result = pd.concat([df1,df2], axis = 1, join='inner')
# result = result.transpose()
# result.to_csv('/Users/aboyher/danforth/projects/gbs/round3/result_df.csv')

le = preprocessing.LabelEncoder()


df = pd.read_csv(sys.argv[1], index_col= 0)
df2 = df[['Phenotype','Data']]
df = df.drop(['Phenotype','Data'], axis = 1)

df = df.apply(le.fit_transform)
enc = preprocessing.OneHotEncoder()
enc.fit(df)
df_OHE = enc.transform(df).toarray()

train_idx = np.where(np.array(df2)[:,1] == "train")
test_idx = np.where(np.array(df2)[:,1] == "test")

trainX = np.array(df)[train_idx][:,:]
trainY = np.array(df2)[train_idx][:,0]
testX = np.array(df)[test_idx][:,:]
testY = np.array(df2)[test_idx][:,0]

clf = tree.DecisionTreeClassifier()
clf = clf.fit(trainX, trainY)
pred = clf.predict(testX)
acc = (pred == np.array(testY))
print(sum(acc)/float(len(pred)) * 100)