import ch2functions as f2
# import numpy as np


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

    profile_matrix = profile
    # for line in profile:
    #     profile_matrix.append(line.split())

    for offset in range(len(sequence) - k + 1):
        frequency = 1.0
        for i in range(offset, offset + k):
            base = sequence[i]
            frequency *= float(profile_matrix[bases[base]][i - offset])

        freqs.append(frequency)

    max_freq = max(freqs)
    index = [i for i, j in enumerate(freqs) if j == max_freq][0]
    return sequence[index:index + k]


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

    bestMotifs = [sequence[0:k] for sequence in DNA]

    # generate k-mer motifs in first sequence from DNA
    motifs = []
    for index in range(len(DNA[0]) - k + 1):
        motifs.append(DNA[0][index:index + k])

    # iterate through generated k-mer motifs
    for m in motifs:
        probableMotifs = []
        probableMotifs.append(m)

        # creates profile based on probableMotifs, then adds next most probable
        # motif to list; repeats for each sequence in DNA
        for i in range(1, t):
            profile = createProfile(probableMotifs)
            probableMotifs.append(profileMostProbableKmer(DNA[i], k, profile))

        if scoreMotifs(probableMotifs) < scoreMotifs(bestMotifs):
            bestMotifs = probableMotifs

    return bestMotifs


def scoreMotifs(motifs):
    totalScore = 0
    for m in motifs:
        totalScore += distanceBetweenPatternAndStrings(m, motifs)

    return totalScore


def createProfile(DNA):
    """
    Generate a probability profile from a set of seqeuences.
    """
    bases = {'A': 0, 'C': 1, 'G': 2, 'T': 3}

    # 2d list for profile
    count_matrix = [[0 for i in range(len(DNA[0]))]
                    for j in range(len(bases))]

    for sequence in DNA:
        for base_index in range(len(sequence)):
            base = sequence[base_index]
            count_matrix[bases[base]][base_index] += 1 / float(len(DNA))

    return count_matrix
