import math
import ch1functions as f1


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
    d = int(d)
    nucleotides = ['A', 'C', 'G', 'T']

    # base cases
    if d == 0:
        return pattern
    if len(pattern) == 1:
        return {'A', 'C', 'G', 'T'}

    neighborhood = set()

    # dynamically generate smaller cases by taing off first letter
    suffixNeighbors = neighbors(pattern[1:], d)
    for text in suffixNeighbors:
        if hammingDistance(pattern[1:], text) < d:
            for base in nucleotides:
                neighborhood.add(base + text)
        else:
            neighborhood.add(pattern[0] + text)
    return neighborhood


def patternToNumber(pattern):
    if len(pattern) == 0:
        return 0
    symbol = pattern[-1]
    suffix = pattern[0:len(pattern) - 1]
    return 4 * patternToNumber(suffix) + baseToNumber(symbol)


def baseToNumber(base):
    return {'A': 0, 'C': 1, 'G': 2, 'T': 3}[base]


def numberToPattern(number, k):
    values = ['A', 'C', 'G', 'T']
    number = int(number)
    k = int(k)
    pattern = ""
    while k > 1:
        quotient, remainder = divmod(number, 4)
        pattern = values[remainder] + pattern
        number = quotient
        k -= 1
    if k == 1:
        pattern = values[number] + pattern

    return pattern


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
    index = []

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

    # 2
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
    Return the most frequently occuring k-mers in genome including reverse
    complements with at most d mismatches.

    Input: FrequentWords(ACGTTGCATGTCGCATGATGCATGAGAGCT, 4, 1)
    Output: GATG ATGC ATGT
    """

    k = int(k)
    d = int(d)
    freqPatterns = []  # set of most frequent words
    freqs = []  # frequency of pattern
    neighborhoods = []
    index = []

    # generate all possible k-mers
    for i in range(len(genome) - k + 1):
        kmer = genome[i:i + k]
        neighborhoods.extend(neighbors(kmer, d))

    freqs = [0] * len(neighborhoods)

    # convert patterns to numbers
    for i in range(len(neighborhoods)):
        pattern = neighborhoods[i]
        num = patternToNumber(pattern)
        index.append(num)
        freqs[i] = 1

#############################################################
    rev_neighborhoods = []
    for pattern in neighborhoods:
        rev_neighborhoods.append(f1.RevComplement(pattern))

    # if reverse complement with at most d mismatches is in genome, + 1 freq
    # to original
    for i in range(len(rev_neighborhoods)):
        pattern = rev_neighborhoods[i]
        num = patternToNumber(pattern)
        if num in index:
            index.append(patternToNumber(neighborhoods[i]))
            freqs.append(1)

    # sort converted patterns
    index.sort()

    for i in range(len(index) - 1):
        if index[i] == index[i + 1]:
            freqs[i + 1] = freqs[i] + 1
    maxCount = max(freqs)

    # convert max k-mers back to pattern
    for i in range(len(index)):
        if freqs[i] == maxCount:
            pattern = numberToPattern(index[i], k)
            freqPatterns.append(pattern)
    return freqPatterns
