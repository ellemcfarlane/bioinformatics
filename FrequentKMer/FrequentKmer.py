# Elle McFarlane
# B343
# Lab 1, pt 2

def frequent_kmer(text, k):
    """
    Returns all most frequent k-mers in text
    :param text: String representing dna
    :param k: integer representing length of dna fragment to find most frequency
    :return: list of most frequent k-mers in text
    """
    frequency_words = []
    kmers = {}

    for i in range(len(text)-k+1):
        kmer = text[i:i+k]
        if kmer in kmers:
            kmers[kmer] += 1
        else:
            kmers[kmer] = 1

    mx = -1
    for key in kmers:
        if kmers[key] > mx:
            mx = kmers[key]

    for key in kmers:
        if kmers[key] == mx:
            frequency_words.append(key)
    return frequency_words

def driver(path):
    """
    Loads test data into count_nucleotides and displays answer
    :param path: String for path of test data file
    :return: void (displays result from count_nucleotides)
    """
    k = -1
    dna = ""
    with open(path, 'r') as f:
        count = 0
        for line in f:
            if count == 0:
                dna = str(line)
            else:
                k = int(line)
            count += 1
        print(frequent_kmer(dna, k))


driver('rosalind_ba1b.txt')



