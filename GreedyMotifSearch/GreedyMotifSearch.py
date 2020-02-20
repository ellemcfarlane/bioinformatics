# Elle McFarlane

def greedy_motif_search(dna, k, t):
    """
    Returns collection of best motifs amongst the given strings Dna
    :param dna: List of strings
    :param k: int, size of motif to find
    :param t: int, number of strings forming Dna
    :return: List of strings for BestMotifs (most similar k-mers among Dna)
    """

    # form original motif matrix arbitrarily from first k-mers in each string

    best_motifs = []
    for sequence in dna:
        best_motifs.append(sequence[0:k])

    # for each kmer motif in the first sequence from Dna
    first_seq = dna[0]
    for i in range(len(first_seq)-k+1):
        motif = first_seq[i:i + k]
        motifs = [motif]
        for j in range(1, t):
            profile = gen_profile(motifs)
            kmer = most_prob_kmer(dna[j], k, profile)
            motifs.append(kmer)
        if score(motifs) < score(best_motifs):
            best_motifs = motifs
    return best_motifs


def gen_profile(sequences):
    """
    Generates DNA profile for list of Dna sequences
    :param sequences: list of strings
    :return: list of lists frequencies for A C G T nucleotides in each position

    Example: ["AAA", "ACG", "GCT"] -> [.66 .33 .33
                                        0  .66  0
                                       .33  0  .33
                                        0   0  .33]
    """
    profile = []
    num_seqs = len(sequences)
    if num_seqs == 0:
        return [[]]
    # for each position, calculate A C G T frequencies
    for i in range(len(sequences[0])):
        freqs = [0, 0, 0, 0]
        for seq in sequences:
            if seq[i] == 'A':
                freqs[0] = freqs[0] + 1
            elif seq[i] == 'C':
                freqs[1] = freqs[1] + 1
            elif seq[i] == 'G':
                freqs[2] = freqs[2] + 1
            elif seq[i] == 'T':
                freqs[3] = freqs[3] + 1
        # convert count to probability (divide each count by number of sequences)
        freqs[:] = [x / (num_seqs) for x in freqs]
        # append frequency count for ith position to profile
        profile.append(freqs)
    return profile


def most_prob_kmer(sequence, k, profile):
    """

    :param sequence: string of dna segment
    :param k: length of kmer
    :param profile: frequencies for A C G T in each position
    :return: most likely kmer based on profile (if multiple, returns first seen)

    Example:
    sequence = ACCAGGACTAAAGCCACA
    profile = [.66 .33 .33
                0  .66  0
               .33  0  .33
                0   0  .33]
    -> ACA
    """
    # for each possible motif in sequence, calculate probability
    # if higher than current highest probability, update
    high_p_kmer = ""
    highest_prob = -1
    for i in range(len(sequence)-k+1):
        kmer = sequence[i:i+k]
        probability = 1
        for j in range(len(kmer)):
            nucleo = kmer[j]
            if nucleo == 'A':
                probability = probability*profile[j][0]
            elif nucleo == 'C':
                probability = probability*profile[j][1]
            elif nucleo == 'G':
                probability = probability*profile[j][2]
            elif nucleo == 'T':
                probability = probability*profile[j][3]
        if probability > highest_prob:
            highest_prob = probability
            high_p_kmer = kmer
    return high_p_kmer

def score(sequences):
    """
    :param sequences: list of strings (motifs)
    :return: sum of how different each position is
    Example:
        motifs =
        "CCGACTGGCA"
        "CCGTCTTGTA"
        "TCGACTAGTA"
        "ACGACGTGCA"
         2001012020 -> 8
    """
    # for each position, count unpopular nucleotides and add to sum
    unpop_score = 0
    for i in range(len(sequences[0])):
        freqs = [0, 0, 0, 0]
        for seq in sequences:
            if seq[i] == 'A':
                freqs[0] = freqs[0] + 1
            elif seq[i] == 'C':
                freqs[1] = freqs[1] + 1
            elif seq[i] == 'G':
                freqs[2] = freqs[2] + 1
            elif seq[i] == 'T':
                freqs[3] = freqs[3] + 1

        # calculate max frequency and subtract
        # from total to find unpopular nucleotide count
        max = 0
        freq_sum = 0
        for count in freqs:
            if count > max:
                max = count
            freq_sum = freq_sum + count
        unpop_score = unpop_score + freq_sum - max
    return unpop_score


def driver(path):
    with open(path, 'r') as f:
        k, t = [int(x) for x in next(f).split()]
        sequences = []
        for line in f:
            sequences.append(line[:-1])
        motifs = (greedy_motif_search(sequences, k, t))
        output = ""
        for motif in motifs:
            output = output + motif + "\n"
        print(output)


driver('rosalind_ba2d.txt')
