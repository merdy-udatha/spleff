# demo.py by merdy udatha 
import sys 
import gzip 

with gzip.open(sys.argv[1]) as fp: 
    for line in fp: 
        seq = line.split()
        print(seq[0], seq[1], seq[2])