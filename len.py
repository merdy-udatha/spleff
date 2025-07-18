# go through each gene, and create a "retention score" 
# in this program, the retention score is defined as the net change in reads / the total no. of nt's passed 
# so for example, a 1000 nt gene that lost 100000 reads will have a score of -100, which is better than a 1000 nt gene that dropped 500000 - has a score of -500
import sys 
import gzip 

genes = []
scores = []
with gzip.open(sys.argv[1], 'rt') as fp: 
    for i, line in enumerate(fp): 
        words = line.split()
        gene = words[0]
        genes.append(gene)
        scores.append(words[2])
        if gene != genes[i - 1]: continue 
        else: 
            
        print(i, words[0], words[1], words[2])
                
                
            
