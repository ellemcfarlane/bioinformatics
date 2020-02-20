# Elle McFarlane
# B343
# Lab 1

def count_nucleotides(dna):
    """
    Returns list of number of As, Cs, Gs, and Ts in given dna string

    :param dna: string
        String of As, Cs, Gs, and Ts representing DNA sequence (ex. "AACTTTGAACG")
    :return: list
        List of ints representing number of As, Cs, Gs, and Ts in given string
    """
    counts = [0] * 4
    for char in dna:
        i = index(char)
        if i != -1:
            counts[i] += 1
    return counts

def index(char):
    """
    Returns index from 0-3 for given character to be used as index for count_nucleotides function.
    Returns -1 if character is not A, C, G, or T.

    :param char: character
        'A', 'C', 'G', or 'T'
    :return: integer from 0-3
    """
    if char == 'A':
        return 0
    elif char == 'C':
        return 1
    elif char == 'G':
        return 2
    elif char == 'T':
        return 3
    else:
        return -1


def driver(path):
    """
    Loads test data into count_nucleotides and displays answer
    :param path: String for path of test data file
    :return: void (displays result from count_nucleotides)
    """
    f = open(path, 'r')
    line = f.readlines()
    dna = str(line)
    print(count_nucleotides(dna))

#driver('rosalind_dna.txt')
