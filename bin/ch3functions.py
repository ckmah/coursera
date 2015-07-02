import ch2functions as f2


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
