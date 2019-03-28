import sys
import re

f = sys.argv[1]

with open(f, 'r') as fn:
    line = fn.readline().rstrip('\n')
    line = line.split('\t')
    scaf = line[0]
    #c = [":".join(line[1:])]
    c = line[1:]
    line = fn.readline().rstrip('\n')
    while line:
        line = line.split('\t')
        if scaf != line[0]:
            print("{} {}".format(scaf,','.join(c)))
            scaf = line[0]
            #c = [":".join(line[1:])]
            c = line[1:]
        else:
            #c+= [":".join(line[1:])]
            c += line[1:]
        line = fn.readline().rstrip('\n')
    print('{} {}'.format(scaf,','.join(c)))
