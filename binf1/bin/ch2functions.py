import ch1functions as f1
import os.path
import sys


def skew(sequence):
    """Calculates the difference between the total number of occurrences of
    G and the total number of occurrences of C in the first i nucleotides of
    sequence."""
    skew_diff = [0]

    C = 'C'
    G = 'G'

    # calculate diff at each index
    for i in range(len(sequence)):
        nuc = sequence[i]
        if nuc == C:
            skew_diff.append(skew_diff[i] - 1)
        elif nuc == G:
            skew_diff.append(skew_diff[i] + 1)
        else:
            skew_diff.append(skew_diff[i])

    return skew_diff


def minSkew(sequence):
    """Find a position in the given genome minimizing the skew."""
    # generate skew diff array
    skew_diff = skew(sequence)

    min_list = []
    min_val = 0

    # loop through sequence once to find minimum values; O(n) runtime
    for index in range(len(skew_diff) - 1):
        candidate = skew_diff[index]
        if candidate < min_val:
            min_list = []
            min_list.append(index)
            min_val = candidate
        elif candidate == min_val:
            min_list.append(index)

    return min_list


def hammingDistance(s1, s2):
    """Return the Hamming distance between equal-length sequences"""
    if len(s1) != len(s2):
        raise ValueError("Undefined for sequences of unequal length.")
    return sum(let1 != let2 for let1, let2 in zip(s1, s2))


def approxMatchPositions(pattern, genome, d):
    """Return all positions in which pattern occurs in genome with at most
    d mismatches."""
    indices = []

    for i in xrange(0, len(genome) - len(pattern) + 1):
        substr = genome[i:i + len(pattern)]
        if hammingDistance(substr, pattern) <= int(d):
            indices.append(i)

    return indices


def approxPatternCount(pattern, genome, d):
    """Returns the number of times a pattern occurs genome with at most d
    mismatches."""

    d = int(d)
    count = 0
    for i in range(len(genome) - len(pattern) + 1):
        subseq = genome[i:i + len(pattern)]
        if hammingDistance(pattern, subseq) <= d:
            count += 1
    return count


def neighbors(pattern, d):
    """ Generates the set of patterns that are at most d distance away
    from pattern.
    """
    d = int(d)
    nucleotides = ['A', 'C', 'G', 'T']

    # base cases
    if d == 0:
        return pattern
    if len(pattern) == 1:
        return {'A', 'C', 'G', 'T'}

    neighborhood = set()

    # dynamically generate smaller cases by taking off first letter
    suffixNeighbors = neighbors(pattern[1:], d)
    for text in suffixNeighbors:
        if hammingDistance(pattern[1:], text) < d:
            for base in nucleotides:
                neighborhood.add(base + text)
        else:
            neighborhood.add(pattern[0] + text)
    return neighborhood


def patternToNumber(pattern):
    """ Converts pattern to number using base 4 with each base assigned a number
    value.
    """
    if len(pattern) == 0:
        return 0
    symbol = pattern[-1]
    suffix = pattern[0:len(pattern) - 1]
    return 4 * patternToNumber(suffix) + baseToNumber(symbol)


def baseToNumber(base):
    """ Assigns each base to a number value of base 4.
    """
    return {'A': 0, 'C': 1, 'G': 2, 'T': 3}[base]


def numberToBase(base):
    """ Converts number back to base using base 4.
    """
    return ['A', 'C', 'G', 'T'][base]


def numberToPattern(number, k):
    """ Converts number to pattern by dividing by base 4. k is resulting
    pattern length.
    """

    number = int(number)
    k = int(k)

    # base case
    if k == 1:
        return numberToBase(number)

    # divide by 4, remainder is right most letter
    prefixIndex, remainder = divmod(number, 4)
    base = numberToBase(remainder)

    # recursively convert pattern prefix
    prefixPattern = numberToPattern(prefixIndex, k - 1)
    return prefixPattern + base


def approxFrequentWords(genome, k, d):
    """
    Return the most frequently occuring k-mers with at most d mismatches.

    Input: FrequentWords(ACGTTGCATGTCGCATGATGCATGAGAGCT, 4, 1)
    Output: GATG ATGC ATGT
    """

    k = int(k)
    d = int(d)
    freqPatterns = []  # set of most frequent words
    freqs = []  # frequency of pattern
    neighborhoods = []
    index = [] # holds patterns converted as numbers

    # generate all possible k-mers
    for i in range(len(genome) - k + 1):
        neighborhoods.extend(neighbors(genome[i:i + k], d))

    freqs = [0] * len(neighborhoods)

    # convert patterns to numbers
    for i in range(len(neighborhoods)):
        pattern = neighborhoods[i]
        index.append(patternToNumber(pattern))
        freqs[i] = 1

    # sort converted patterns
    index.sort()

    # count pattern frequencies
    for i in range(len(neighborhoods) - 1):
        if index[i] == index[i + 1]:
            freqs[i + 1] = freqs[i] + 1
    maxCount = max(freqs)

    # convert max k-mers back to pattern
    for i in range(len(neighborhoods)):
        if freqs[i] == maxCount:
            pattern = numberToPattern(index[i], k)
            freqPatterns.append(pattern)
    return freqPatterns


def approxFrequentWordsWithRevComp(genome, k, d):
    """
    Return the most frequently occuring k-mers (with mismatches and reverse
    complements) in a string.

    Input: FrequentWords(ACGTTGCATGTCGCATGATGCATGAGAGCT, 4, 1)
    Output: ATGT ACAT
    """

    k = int(k)
    d = int(d)
    freqPatterns = []  # set of most frequent words
    freqs = []  # frequency of pattern
    neighborhoods = []
    index = []

    # generate all possible k-mers
    for i in range(len(genome) - k + 1):
        neighborhoods.extend(neighbors(genome[i:i + k], d))
        neighborhoods.extend(neighbors(f1.RevComplement(genome[i:i + k]), d))

    freqs = [0] * len(neighborhoods)

    # convert patterns to numbers
    for i in range(len(neighborhoods)):
        pattern = neighborhoods[i]
        index.append(patternToNumber(pattern))
        freqs[i] = 1

    # sort converted patterns
    index.sort()

    # count pattern frequencies
    for i in range(len(neighborhoods) - 1):
        if index[i] == index[i + 1]:
            freqs[i + 1] = freqs[i] + 1
    maxCount = max(freqs)

    # convert max k-mers back to pattern
    for i in range(len(neighborhoods)):
        if freqs[i] == maxCount:
            pattern = numberToPattern(index[i], k)
            freqPatterns.append(pattern)
    return freqPatterns
