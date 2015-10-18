import ch2functions as f2
import numpy as np


def motifEnumeration(DNA, k, d):
    """
    Find all k-mer motifs with at most d mismatches given a collection of
    strings.

    Input: A collection of strings Dna, and integers k and d.
    Output: All (k, d)-motifs in Dna.
    """
    k = int(k)
    d = int(d)
    patterns = []

    DNAkmers = [[] for seq in DNA]  # 2d array of all sequence kmers
    candidateMotifs = []  # contains set of generated kmer neighborhoods

    # create list of all kmers as well as kmer neighborhoods
    for dnaIndex in range(len(DNA)):
        sequence = DNA[dnaIndex]
        for seqIndex in range(len(sequence) - k + 1):
            kmer = sequence[seqIndex:seqIndex + k]
            DNAkmers[dnaIndex].append(kmer)
            candidateMotifs.extend(f2.neighbors(kmer, d))

    candidateMotifs = set(candidateMotifs)

    # exhaustive search candidates for true motifs
    for candidate in candidateMotifs:
        found = [False] * len(DNAkmers)
        for listIndex in range(len(DNAkmers)):
            for kmer in DNAkmers[listIndex]:
                if f2.hammingDistance(candidate, kmer) <= d:
                    found[listIndex] = True
                    break
        if False not in found:
            patterns.append(candidate)

    return set(patterns)


def medianString(DNA, k):
    """
    Finds a median string.

    Input: A collection of strings Dna and an integer k.
    Output: A k-mer Pattern minimizing d(Pattern, Dna) among all k-mers
    Patterns.
    """
    k = int(k)
    distance = float("inf")
    median = ""

    # generate and check all possible k-mers
    for i in range(4**k):
        pattern = f2.numberToPattern(i, k)
        d = distanceBetweenPatternAndStrings(pattern, DNA)
        if distance > d:
            distance = d
            median = pattern

    return median


def distanceBetweenPatternAndStrings(pattern, DNA):
    """
    Finds total distance between a pattern a set of seqences.

    Input: A string pattern and collection of strings DNA.
    Output: The total distance between DNA and pattern.
    """
    k = len(pattern)
    distance = 0

    # finds minimum distance by checking all kmers in each string
    for sequence in DNA:
        hammingDistance = float("inf")
        for i in range(len(sequence) - k + 1):
            kmer = sequence[i:i + k]
            d = f2.hammingDistance(pattern, kmer)
            if hammingDistance > d:
                hammingDistance = d
        distance += hammingDistance

    return distance


def profileMostProbableKmer(sequence, k, profile):
    """
    Finds a Profile-most probably k-mer in a string.

    Input: A string Text, an integer k, and a 4 x k matrix Profile.
    Output: A Profile-most probably k-mer in Text.
    """

    k = int(k)
    bases = {'A': 0, 'C': 1, 'G': 2, 'T': 3}

    freqs = []

    profile_matrix = []

    for line in profile:
        profile_matrix.append(line.split())

    for offset in range(len(sequence) - k + 1):
        frequency = 1.0
        for i in range(offset, offset + 5):
            base = sequence[i]
            frequency *= float(profile_matrix[bases[base]][i - offset])

        freqs.append(frequency)


def greedyMotifSearch(k, t, DNA):
    """
    Finds the most probable motif based on a profile generated from a collection
    of sequences.

    Input: Integers k and t, followed by a collection of strings Dna.
    Output: A collection of strings BestMotifs resulting from running
    GreedyMotifSearch(Dna, k, t). If at any step you find more than one
    Profile-most probable k-mer in a given string, select the one occurring
    first in the string.
    """

    k = int(k)
    t = int(t)

    bestMotifs = []

    motifs = []
    for index in range(len(DNA[0]) - k + 1):
        motifs.append(DNA[0][index:index + k])

    for m in motifs:
        m1 = m
        for i in range(2, t):
            profile = createProfile(DNA)


def createProfile(DNA):
    bases = {'A': 0, 'C': 1, 'G': 2, 'T': 3}

    count_matrix = [[0] * len(DNA[0])] * len(bases)
    print count_matrix
    for sequence in DNA:
        print sequence
        for base_index in range(len(sequence)):
            base = sequence[base_index]
            count_matrix[bases[base]][base_index] += 1
            print count_matrix[bases[base]]
            # print count_matrix
        # print count_matrix
    # print count_matrix
    # print np.cumsum(count_matrix, axis=0)
