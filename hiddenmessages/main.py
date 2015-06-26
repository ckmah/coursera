import sys
import hiddenmessages as hm

lines = sys.stdin.read().splitlines()  # read in the input from STDIN

text = lines[0]
k = lines[1]
# (k, l, t) = lines[1].split(' ')

print hm.FrequentWords(text, k)
# print hm.ClumpFind(genome, k, l, t)