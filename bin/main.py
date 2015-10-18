import ch1functions as f1
import ch2functions as f2
import ch3functions as f3
import sys

import math

lines = sys.stdin.read().splitlines()  # read in the input from STDIN

genome = lines[0]
(k, l, t) = lines[1].split()

# print f3.greedyMotifSearch(k, t, DNA)
print f3.createProfile(DNA)
