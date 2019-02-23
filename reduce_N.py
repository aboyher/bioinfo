import re
import sys


f = sys.argv[1]
N = "N"*100
n = 'n'*100
with open(f, 'r') as fn:
    line = fn.readline().rstrip("\n")
    print(line)
    line = fn.readline().rstrip("\n")
    seq = ''
    while line:
        if line.startswith(">"):
            seq = re.sub("N+",N,seq)
            # seq = re.sub('n+',n,seq)
            print(seq)
            print(line)
            seq = ""
        else:
            seq+=line
        line = fn.readline().rstrip("\n")
    seq = re.sub("N+", N,seq)
    # seq = re.sub('n+',n,seq)
    print(seq)
