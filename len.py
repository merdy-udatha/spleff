# go through each gene, and create a "retention score" 
# in this program, the retention score is defined as the net change in reads / the total no. of nt's passed 
# so for example, a 1000 nt gene that lost 100000 reads will have a score of -100, which is better than a 1000 nt gene that dropped 500000 - has a score of -500
import sys 
import gzip 


genes = {}
with gzip.open(sys.argv[1], 'rt') as fp: 
    for line in fp: 
        words = line.split()
        gene = words[0]
        length = float(words[1])
        reads = float(words[2])
        if gene not in genes:
            genes[gene] = {'length': length, 'reads': reads}
        else:
            genes[gene]['length'] = length
            genes[gene]['reads']  = reads - genes[gene]['reads']
	
length = []
reads =  []
for gene in genes:
    length.append(genes[gene]['length'])
    reads.append(genes[gene]['reads'])
	
for len, read in zip (length, reads):
    print(len, read)

	
		
		
        		
                
            
