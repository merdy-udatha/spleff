# pwm.py by merdy udatha 
# We want to take each splice site and create a PWM out of it 
import sys 
import gzip 
import math

def pwm_matrix(splice_site_len, nt_start, strand_no): 
	positions = [[] for i in range(splice_site_len)]
	with gzip.open(sys.argv[1], 'rt') as fp: 
		for line in fp: 
			words = line.split()
			strand = str(words[0])
			if '-' in strand: continue 
			site = list(words[1])
			
print(pwm_matrix(5))
	
"""
donor_posits = [[], [], [], [], []] # set up empty lists for each position to store a dictionary 
acc_posits = [[] , [], [], [], [], []]
with gzip.open(sys.argv[1], 'rt') as fp: 
    for line in fp:
        words = line.split()
        strand = str(words[0])
        if '-' in strand: continue # meaning we are doing plus strand
        donor_1 = list(words[1])[10:15]
        acc_1 = list(words[1])[26:32]
        for i, nt in enumerate(donor_1): 
            donor_posits[i].append(nt)
        for i, nt in enumerate(acc_1):
            acc_posits[i].append(nt)
print(acc_posits[5])
		

donor_ct_matrix = [] 

for nt_list in donor_posits: 
    nt_counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0} # stores number of nt counts per position
    for nt in nt_list: 
        if nt in nt_counts: 
            nt_counts[nt] += 1 
        else: continue 
    donor_ct_matrix.append(nt_counts) 


# printing all donor site info 
print()   
print("\tDonor Site Frequency Matrix\t")
print("Position\tA\tC\tG\tT")

for i, nt_counts in enumerate(donor_ct_matrix): 
    a = nt_counts['A'] 
    c = nt_counts['C']
    g = nt_counts['G']
    t = nt_counts['T']
    print(f'{i + 1}\t\t{a}\t{c}\t{g}\t{t}')

prob_matrix = []

for counts in donor_ct_matrix: 
    total = counts['A'] + counts['C'] + counts['G'] + counts['T']
    probs = {
                'A': counts['A'] / total,       # probability of an 'A' nt in each position
                'C': counts['C'] / total,
                'G': counts['G'] / total,
                'T': counts['T'] / total 
    } 
    prob_matrix.append(probs) 
    
print()
print("\tDonor Site Probability Matrix\t")
print("Pos\tA\tC\tG\tT")
for i, row in enumerate(prob_matrix): 
    print(f"{i+1}\t{row['A']:.3f}\t{row['C']:.3f}\t{row['G']:.3f}\t{row['T']:.3f}")
    
pwm_matrix = [] 
for probs in prob_matrix: 
    pwm_row = {}
    for base in ['A', 'C', 'G', 'T']:
        p = probs[base] 
        if p == 0: 
            pwm_row[base] = -99 
        else: pwm_row[base] = math.log2(p / 0.25)
    pwm_matrix.append(pwm_row)
    

    

scores = []
gene_exp = []
with gzip.open(sys.argv[1], 'rt') as fp: 
    for line in fp: 
        words = line.split()
        strand = str(words[0])
        if '-' in strand: continue # meaning we are doing plus strand
        site_one = list(words[1])[10:15]
        score = 0 
        for nt, row in zip(site_one, pwm_matrix):   
            if nt in row: 
                score += row[nt]
#            print(nt, row[nt])
 #       print(score, words[2])
		
print()
print("\tDonor Site Log-Odds Matrix\t")
print('Pos\tA\tC\tG\tT')
for i, row in enumerate(pwm_matrix): 
    print(f"{i+1}\t{row['A']:.2f}\t{row['C']:.2f}\t{row['G']:.2f}\t{row['T']:.2f}")
"""
