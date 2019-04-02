import sys
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

with open(f, 'r') as fn:
    line = fn.readline()
    while line.startswith("##"):
        line = fn.readline().rstrip('\n')
    header = line.split("\t")[9:]
    print(",".join(header))
    line = fn.readline().rstrip('\n')
    while line:
        snp = line.split()
        idx = ":".join(snp[:2])
        data = [x.split(":") for x in snp[9:]]
        geno = [x[0] for x in data]
        #print("{}".format(",".join(geno)))
        #print("{}".format(idx))
        print("{},{}".format(idx,",".join(geno)))
        line = fn.readline().rstrip('\n')