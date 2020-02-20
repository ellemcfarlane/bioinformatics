# Elle McFarlane
# B363
# 7/13/19
# Lab 3

def min_skew(dna):
    """
    Returns index at which skew is a minimum.
    :param dna: String for dna as A's, C', G's, and T's
    :return: list of integer(s) for position(s) at which skew is minimum
    """
    sz = len(dna)
    if sz == 0:
        return 0
    skew = [0]
    min_so_far = 0

# find minimum skew
    for i in range(sz):
        char = dna[i]
        part_val = 0

        if char == 'G':
            part_val = 1
        elif char == 'C':
            part_val = -1
        val = skew[i] + part_val
        skew.append(val)
        if val < min_so_far:
            min_so_far = val

# add all positions that match minimum skew to a list
    mins = []
    for i in range(len(skew)):
        if skew[i] == min_so_far:
            mins.append(i)
    return mins

def driver(path):
    """
    Loads test data into min_skew and displays answer
    :param path: String for path of test data file
    :return: void (displays result from min_skew)
    """
    dna = ""
    with open(path, 'r') as f:
        dna = str(f.readline())
        print(min_skew(dna))

# driver('rosalind_ba1f.txt')