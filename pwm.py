# pwm.py by merdy udatha 
# We want to take each splice site and create a PWM out of it 
import sys 
import gzip 
import math

def pwm_matrix(splice_site_len, nt_start, strand_no, name): 
    positions = [[] for i in range(splice_site_len)]
    with gzip.open(sys.argv[1], 'rt') as fp: 
        for line in fp: 
            words = line.split()
            strand = str(words[0])
            if '-' in strand: continue 
            site = list(words[strand_no])[nt_start:nt_start + splice_site_len]
            for i, nt in enumerate(site): 
                positions[i].append(nt)
    count_matrix = [] 
    for nt_list in positions: 
        nt_counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
        for nt in nt_list: 
            if nt in nt_counts: 
                nt_counts[nt] += 1
            else: continue 
        count_matrix.append(nt_counts) 
    prob_matrix = [] 
    for counts in count_matrix: 
        total = counts['A'] + counts['C'] + counts['G'] + counts['T'] 
        probs = {'A': counts['A']/total,
                 'C': counts['C']/total,
                 'G': counts['G']/total,
                 'T': counts['T']/total}
        prob_matrix.append(probs) 
    pwm = [] 
    for probs in prob_matrix: 
        row = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
        for nt in probs: 
            if probs[nt] == 0: 
                row[nt] = -99
            else: 
                row[nt] = math.log2(probs[nt] / 0.25) 
        pwm.append(row)
    #print(f"\t{name} Log-Odds Matrix\t")
    #print('Pos\tA\tC\tG\tT')
    #for i, row in enumerate(pwm): 
    #    print(f"{i+1}\t{row['A']:.2f}\t{row['C']:.2f}\t{row['G']:.2f}\t{row['T']:.2f}")
    return pwm

scores = [] 
pwm_donor = pwm_matrix(5, 10, 1, "donor_site 1")
pwm_acc = pwm_matrix(6, 26, 1, "acc_site 1")
a = 1
b = 3
with gzip.open(sys.argv[1], 'rt') as fp: 
    for line in fp: 
        words = line.split() 
        strand = str(words[0]) 
        if '-' in strand: continue 
        donor_site_a = list(words[a])[10:15] # have to change based on positions
        acc_site_a = list(words[a])[26:32]
        score_a = 0 
        for nt, row in zip(donor_site_a, pwm_donor): 
            if nt in row: 
                score_a += row[nt] 
        for nt, row in zip(acc_site_a, pwm_acc): 
            if nt in row: 
                score_a += row[nt]
        score_b = 0 
        donor_site_b = list(words[b])[10:15] 
        acc_site_b   = list(words[b])[26:32] 
        for nt, row in zip(donor_site_b, pwm_donor): 
            if nt in row: 
                score_b += row[nt]
        for nt, row in zip(acc_site_b, pwm_acc):
            if nt in row:
                score_b += row[nt]
        print(score_a - score_b, math.log2(float(words[2]) / float(words[4]))) 
        




	



