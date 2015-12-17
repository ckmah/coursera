def stringComposition(k, sequence):
    '''
    Given a sequence and parameter k, return the k-mer composition of
    the sequence.
    '''
    kmers = []
    for index in range(len(sequence) - k + 1):
        kmers.append(sequence[index:index + k])

    return kmers


def genomePath(sequences):
    '''
    Reconstruct a string given that all sequences given are the same length and
    that each subsequent sequence is the previous[1:] plus a new base.
    '''
    genome = ""
    for seq in sequences:
        if len(genome) is 0:
            genome = seq
        else:
            genome += seq[-1]

    return genome


def overlapGraph(sequences):
    '''
    Return an adjacency list given a set of overlapping k-mers.
    '''

    adjacencyList = []

    for index in range(len(sequences) - 1):

        element1 = sequences[index]
        # compare element with all elements after
        for element2 in sequences[index + 1:]:

            if element1[1:] == element2[:-1]:
                adjacencyList.append((element1, element2))
            if element2[1:] == element1[:-1]:
                adjacencyList.append((element2, element1))

    return adjacencyList