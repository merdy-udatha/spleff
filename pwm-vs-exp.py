import argparse
import gzip
import math

import numpy as np
import scipy
from matplotlib import pyplot as plt
import korflab

def train_pwms(data):
	don = [{'A':0, 'C':0, 'G':0, 'T':0} for _ in range(10)]
	acc = [{'A':0, 'C':0, 'G':0, 'T':0} for _ in range(10)]
	for d1, a1, x1, d2, a2, x2 in data:
		for i, nt in enumerate(d1): don[i][nt] += 1
		for i, nt in enumerate(a1): acc[i][nt] += 1
		for i, nt in enumerate(d2): don[i][nt] += 1
		for i, nt in enumerate(a2): acc[i][nt] += 1
	for i in range(10):
		for nt in 'ACGT':
			don[i][nt] /= (len(data * 2))
			acc[i][nt] /= (len(data * 2))
	return don, acc

def score_pwm(pwm, seq):
	score = 1
	for i, nt in enumerate(seq):
		score *= pwm[i][nt]
	return score



parser = argparse.ArgumentParser()
parser.add_argument('datafile', help='expbias.txt.gz probably')
parser.add_argument('--xvalid', type=int, default=3,
	help='cross-validation [%(default)]')
arg = parser.parse_args()

data = []
with gzip.open(arg.datafile, 'rt') as fp:
	for line in fp:
		strand, seq1, exp1, seq2, exp2 = line.split()
		exp1 = float(exp1)
		exp2 = float(exp2)
		if strand == '-':
			seq1 = korflab.anti(seq1)
			seq2 = korflab.anti(seq2)
		d1 = seq1[10:20]
		a1 = seq1[22:32]
		d2 = seq2[10:20]
		a2 = seq2[22:32]

		# skip non-canonical (GT..AG only)
		if not d1.startswith('GT'): continue
		if not d2.startswith('GT'): continue
		if not a1.endswith('AG'): continue
		if not a2.endswith('AG'): continue

		data.append( (d1, a1, exp1, d2, a2, exp2) )

for i in range(arg.xvalid):
	train = []
	test = []
	for j, d in enumerate(data):
		if j % arg.xvalid == i: test.append(d)
		else:                   train.append(d)
	don, acc = train_pwms(train)

	scores = []
	express = []
	for d1, a1, x1, d2, a2, x2 in test:
		s1 = score_pwm(don, d1) * score_pwm(acc, a1)
		s2 = score_pwm(don, d2) * score_pwm(acc, a2)
		scores.append(math.log2(s1 / s2))
		express.append(math.log2(x1 / x2))

	x = np.array(scores)
	y = np.array(express)
	cc, pv, = scipy.stats.pearsonr(x, y)
	print(i, cc, pv, sep='\t')

	# from Merdy
	plt.scatter(x, y, s=0.5, color = 'black')
	plt.xlabel('log-odds (site a / site b)')
	plt.ylabel('SSA - SSB')
	plt.title(f'SS difference vs Gene expression', fontsize = 18)
	plt.savefig("pwm_scatter.png", dpi=300, bbox_inches='tight')
	plt.show()
