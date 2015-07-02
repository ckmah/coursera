import ch1functions as f1
import ch2functions as f2
import ch3functions as f3
import sys

import math

lines = sys.stdin.read().splitlines()  # read in the input from STDIN

# (k, d) = lines[0].split()
# DNA = lines[1:]

# result = ""
# for motif in f3.motifEnumeration(DNA, k, d):
#     result += motif + ' '

# print result.strip()

numlist = []
for line in lines:
    numlist.extend(line.split())

print numlist
entropy = 0
for num in numlist:
    if num is 'x':
        continue
    num = int(num)
    p = float(num) / 10
    if p > 0:
        entropy += p * math.log(p, 2)

print entropy
