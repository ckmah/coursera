import ch1functions as f1
import sys

import math

lines = sys.stdin.read().splitlines()  # read in the input from STDIN

''' String Decomposition '''
# k = int(lines[0])
# sequence = lines[1]

# result = ""
# for kmer in f1.stringComposition(k, sequence):
#     result += kmer + '\n'

# print result.strip()

''' Genome Path '''
# sequences = lines

# print f1.genomePath(sequences)


''' Overlap Graph '''
sequences = lines

for item in f1.overlapGraph(sequences):
    print item[0] + ' -> ' + item[1]
