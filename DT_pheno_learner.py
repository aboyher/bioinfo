import numpy as np
import pandas as pd
import DTLearner as DT





def split_gt_fields(gt):
    return [x.split(":") for x in gt]

def split_data(data, res_idx, sus_idx):
    ran_res = np.random.choice(res_idx, int(len(res_idx)*.5))
    ran_sus = np.random.choice(sus_idx, int(len(sus_idx) * .5))
    select = ran_sus.tolist() + ran_res.tolist()
    train = data[select, :]
    test = np.delete(data,select, axis=0)
    print(test.shape)
    print(train.shape)
    return train, test


vcf = '/Users/aboyher/danforth/projects/gbs/round3/snpcall_gbsx_demux_rnd3_tme7_v0.5.3a_SS26_190326.vcf'
phenodata = '/Users/aboyher/danforth/projects/gbs/round3/progenyinfo.txt'

pdata = pd.read_table(phenodata, sep = "\t")

pdata_res = pdata[pdata['Phenotype'] == 'R']
pdata_sus = pdata[pdata['Phenotype'] == 'S']




with open(vcf) as fn:
    line = fn.readline()
    while line.startswith("##"):
        line = fn.readline()

    inds = line.rstrip("\n").split()[9:]
    inds = [i.split('.R1')[0] for i in inds]
    inds = [i.split('.2')[0] for i in inds]
    inds = [i.split('.3')[0] for i in inds]

    res_idx = np.where([i in pdata_res['Sample'].tolist() for i in inds])[0]
    sus_idx = np.where([i in pdata_sus['Sample'].tolist() for i in inds])[0]

    line = fn.readline().rstrip('\n')
    snp = []
    while line:
        snp.append([x[0] for x in split_gt_fields(line.split()[9:])])

        line = fn.readline()



pheno_line = np.array([-1] * len(snp[0]))
pheno_line[sus_idx] = 1
snp.append(pheno_line.tolist())
snp = np.array(snp).transpose()

train, test = split_data(snp, res_idx, sus_idx)
trainX = train[:,:-1]
trainY = train[:,-1]
testX = test[:, :-1]
testY = test[:,-1]

pd_trainX = pd.factorize(pd.DataFrame(trainX))
print(pd_trainX)



# learner = DT.DTLearner(leaf_size = 3)
# learner.addEvidence(train)
# pred = learner.query(testX)
# rmse = np.sqrt(((testY - pred) ** 2).sum()/testY.shape[0])
# print(rmse)
#
