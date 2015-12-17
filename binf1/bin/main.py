import ch1functions as f1
import ch2functions as f2
import ch3functions as f3
import sys

import math

lines = sys.stdin.read().splitlines()  # read in the input from STDIN
pattern = lines[0]

# print f3.greedyMotifSearch(k, t, DNA)
# result = ""
# for motif in f3.greedyMotifSearch(k,t,DNA):
#     print motif

print f2.patternToNumber(pattern)
