import sys
import numpy as np
import pandas as pd
f = sys.argv[1]

def get_hap(geno):
    hap0 = []
    hap1 = []
    for i in geno:
        if i != ".":
            hap0.append(i[0])
            hap1.append(i[2])
        else:
            hap0.append('.')
            hap1.append('.')
    return hap0,hap1

# df = pd.read_csv(sys.argv[2], header = None, index_col= 0 )
# loc = np.array(df)



with open(f, 'r') as fn:
    for line in fn:
        if line.startswith("##"):
            continue

        if line.startswith("#"):
            header = line.split("\t")[9:]
            header = [x.split(".R1")[0] for x in header]
            print(",".join(header))
            continue

        snp = line.split()
        weight = float(snp[5])
        if weight < 100:
            continue
        # if int(snp[1]) not in loc:
        #     line = fn.readline().rstrip('\n')
        #     continue
        idx = ":".join(snp[:2])
        alleles = [snp[3]]
        alleles += snp[4].split(',')

        data = [x.split(":") for x in snp[9:]]
        geno = [x[0] for x in data]
        # eals = []
        # for g in geno:
        #     if g != '.':
        #         hap0,hap1 = g.split("/")
        #         al = alleles[int(hap0)] + "/" + alleles[int(hap1)]
        #     else:
        #         al = '.'
        #     eals.append(al)

        # print("{}".format(",".join(eals)))
        # print("{}".format(idx))
        print("{},{}".format(idx,",".join(geno)))