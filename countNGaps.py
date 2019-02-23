import sys
import re
from Bio import SeqIO

f = sys.argv[1]
total = 0
with open(f, 'r') as fn:
    for record in SeqIO.parse(fn, "fasta"):
        for match in re.finditer('N+', str(record.seq).upper()):
            total+=1
print(total)

