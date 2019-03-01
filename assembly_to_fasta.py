import sys

f = sys.argv[1]

with open(f, 'r') as fn:
    line = fn.readline().rstrip("\n")
    prev_contig = ""
    c = 0
    while line:
        if line[0] == ">":
            if ":::" in line:
                contig = line[1:].split(":::")[0]
            else:
                contig = line[1:].split(" ")[0]
            if contig!=prev_contig:
                start = 0
                end = int(line.split(" ")[2])
                prev_contig = contig
            else:
                start = end
                end += int(line.split(" ")[2])
            print("\t".join([contig,str(start),str(end)]))
        line = fn.readline().rstrip("\n")