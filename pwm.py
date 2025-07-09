# pwm.py by merdy udatha 
# We want to take each splice site and create a PWM out of it 
import sys 
import gzip 

position_lists = [[], [], [], [], []]
with gzip.open(sys.argv[1], 'rt') as fp: 
    for line in fp:
        words = line.split()
        strand = str(words[0])
        if '-' in strand: continue # meaning we are doing plus strand
        site_one = list(words[1])[10:15]
        for i, nt in enumerate(site_one): 
            position_lists[i].append(nt)
    

count_matrix = [] 

for nt_list in position_lists: 
    nt_counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0} 
    for nt in nt_list: 
        if nt in nt_counts: 
            nt_counts[nt] += 1 
        else: continue 
    count_matrix.append(nt_counts) 

print()   
print("\tFrequency Matrix\t")


print("Position\tA\tC\tG\tT")

for i, nt_counts in enumerate(count_matrix): 
    a = nt_counts['A'] 
    c = nt_counts['C']
    g = nt_counts['G']
    t = nt_counts['T']
    print(f'{i + 1}\t\t{a}\t{c}\t{g}\t{t}')

prob_matrix = []

for counts in count_matrix: 
    total = counts['A'] + counts['C'] + counts['G'] + counts['T']
    probs = {
                'A': counts['A'] / total,
                'C': counts['C'] / total,
                'G': counts['G'] / total,
                'T': counts['T'] / total 
    } 
    prob_matrix.append(probs) 
    
print()
print("\tProbability Matrix\t")
print("Pos\tA\tC\tG\tT")
for i, row in enumerate(prob_matrix): 
    print(f"{i+1}\t{row['A']:.3f}\t{row['C']:.3f}\t{row['G']:.3f}\t{row['T']:.3f}")
    