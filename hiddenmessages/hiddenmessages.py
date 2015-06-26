import copy


def PatternCount(text, pattern):
    """
    Counts number of occurences of pattern in text. Takes into account
    overlapping instances of pattern.

    Input: PatternCount(CABABAC, ABA)
    Output: 2
    """
    count = 0
    for i in xrange(0, len(text) - len(pattern) + 1):
        if text[i:i + len(pattern)] == pattern:
            count = count + 1
    return count


def FrequentWords(text, k):
    """
    Return the most frequently occuring k-mers.

    Input: FrequentWords(ACGTTGCATGTCGCATGATGCATGAGAGCT, 4)
    Output: CATG GCAT
    """
    k = int(k)
    freqPatterns = []  # set of most frequent words
    freqs = []  # frequency of pattern

    # count number of occurences of all k-mers
    for i in xrange(0, len(text) - k):
        pattern = text[i:i + k]
        freqs.append(PatternCount(text, pattern))
    maxCount = max(freqs)

    # save most frequent words
    for i in xrange(0, len(text) - k):
        if freqs[i] == maxCount:
            freqPatterns.append(text[i:i + k])

    # format output
    result = ""
    for e in set(freqPatterns):
        result += e + ' '

    return result.strip()


def RevComplement(sequence):
    """
    Given a sequence, return its reverse complement.

    Input: RevComplement(AAAACCCGGT)
    Output: ACCGGGTTTT
    """

    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

    result = ""
    for letter in sequence:
        result += complement[letter]

    return result[::-1]


def PatternPositions(pattern, genome):
    """
    Return all positions in which pattern occurs in genome.

    Input: PatternPositions(ATAT, GATATATGCATATACTT)
    Output: 1 3 9
    """
    result = ""
    for i in xrange(0, len(genome) - len(pattern) + 1):
        if genome[i:i + len(pattern)] == pattern:
            result += str(i) + " "

    return result.strip()


def ClumpFind(genome, k, l, t):
    """
    Return all distinct k-mers forming (l,t)-clumps in genome.

    Input:
    CGGACTCGACAGATGTGAAGAACGACAATGTGAAGACTCGACACGACAGAGTGAAGAGAAGAGGAAACATTGTAA
    5 50 4
    Output:CGACA GAAGA
    """
    k = int(k)
    l = int(l)
    t = int(t)

    # count occurences of each k-mer
    patterns = {}
    for i in xrange(0, len(genome) - k + 1):
        p = genome[i:i + k]
        if p in patterns:
            patterns[p].append(i)
        else:
            patterns[p] = [i]

    print len(patterns)
    print len(set(patterns))
    patterns_copy = copy.deepcopy(patterns)
    # filter
    for pattern in patterns:
        occurences = patterns[pattern]

        # delete any k-mers that don't occur at least t times
        if len(occurences) < t:
            del patterns_copy[pattern]

        # delete remaining k-mers that don't occur within l interval
        if len(occurences) > t:
            delete = True  # delete flag
            for i in xrange(0, len(occurences) - t + 1):
                if int(occurences[i + t - 1]) - int(occurences[i]) < l:
                    delete = False
                    break
            if delete:
                del patterns_copy[pattern]

        # return modified set
        patterns = patterns_copy

    return len(patterns)
    # result = ""
    # for key in patterns:
        # result += key + " "

    # return result.strip()
