# go through each gene, and create a "retention score" 
# in this program, the retention score is defined as the net change in reads / the total no. of nt's passed 
# so for example, a 1000 nt gene that lost 100000 reads will have a score of -100, which is better than a 1000 nt gene that dropped 500000 - has a score of -500
import sys 
import gzip 
from matplotlib import pyplot as plt 
import numpy as np


genes = {}
with gzip.open(sys.argv[1], 'rt') as fp: 
    for line in fp: 
        words = line.split()
        gene = words[0] 
        length = float(words[1])
        reads = float(words[2])
        if gene not in genes: 
            genes[gene] = []
        genes[gene].append(length)
        genes[gene].append(reads)
    length = []
    reads = []
    for gene in genes: 
        if ((genes[gene][::-1][0] - genes[gene][1])/ genes[gene][::-1][1]) == 0: continue
        length.append(genes[gene][::-1][1])
        reads.append((genes[gene][::-1][0] - genes[gene][1])/ genes[gene][::-1][1])

for len, rds in zip(length, reads):
    print(len, rds)



x = np.array(length)
y = np.array(reads)
plt.scatter(x, y, s=0.5, color = 'black')
plt.xlabel('Gene length (nt)')
plt.ylabel('Retention Score')
plt.title(f'Net change in reads vs Gene length (nt)', fontsize = 18)
plt.savefig("len_scatter.png", dpi=300, bbox_inches='tight')
plt.show()

		
        		
                
            
