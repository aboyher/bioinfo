import sys

with open(sys.argv[1], 'r') as fn:
    line = fn.readline().rstrip('\n')
    while line.startswith("# ") or line.startswith("###"):
        print(line)
        line = fn.readline().rstrip('\n')
    good = False
    while line:
        if line.startswith("#"):
            linesplit = line.split(' ')
            chrom0,chrom1 = linesplit[6].split('&')
            if chrom0[-2:] == chrom1[-2:]:
                print(line)
                line = fn.readline().rstrip('\n')
                while not line.startswith("#"):
                    print(line)
                    line=fn.readline().rstrip('\n')
            else:
                line = fn.readline().rstrip('\n')

