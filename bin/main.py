import ch1functions as f1
import ch2functions as f2
import ch3functions as f3
import sys

import math

lines = sys.stdin.read().splitlines()  # read in the input from STDIN

genome = lines[0]
k, d = lines[1].split()

print f2.approxFrequentWordsWithRevComp(genome, k, d)
