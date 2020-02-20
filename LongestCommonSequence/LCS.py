# Elle McFarlane

def LCS_backtrack(v, w):
    """
    Finds longest string that is common to both strings via deletions and insertions
    :param v: string
    :param w: string
    :return: backtrack, matrix for paths of longest common sequence between all indices of v and w
    """
    vsz = len(v) + 1
    wsz = len(w) + 1
    # initialize sum-of-paths matrix and backtrack matrix
    s = [x[:] for x in [[0] * wsz] * vsz]
    backtrack = [x[:] for x in [[""] * wsz] * vsz]

    # initialize completely vertical and horizontal paths to 0 in sum matrix
    # and to their respective symbols in the backtrack matrix
    for i in range(vsz):
        s[i][0] = 0
        backtrack[i][0] = 'i'
    for j in range(wsz):
        s[0][j] = 0
        backtrack[0][j] = 'd'

    for i in range(1, vsz):
        for j in range(1, wsz):
            match = 0
            # check if currents characters match
            if v[i-1] == w[j-1]:
                match = 1
            # take longest path
            s[i][j] = max(s[i-1][j], s[i][j-1], s[i-1][j-1] + match)
            # determine which kind of move was the longest
            # insertion (moved down)
            if s[i][j] == s[i-1][j]:
                backtrack[i][j] = 'i'
            # deletion (moved right)
            elif s[i][j] == s[i][j - 1]:
                backtrack[i][j] = 'd'
            # match (moved diagonally)
            elif s[i][j] == s[i - 1][j - 1] + 1 and match == 1:
                backtrack[i][j] = 'm'
    return backtrack

def recur_output_LCS(backtrack, v, i, j):
    """
    :param backtrack: double matrix with characters indicating movement to get to each space in matrix
    :param v: word
    :param i: i-prefix of string v
    :param j: j-prefix of string w (other string comparing with v)
    :return: prints longest common sequence between the i-prefix of v and j-prefix of w
    """
    if i == 0 or j == 0:
        return ""
    # if previous move was down (insertion)
    if backtrack[i][j] == 'i':
        lcs = recur_output_LCS(backtrack, v, i-1, j)
    # if previous move was right (deletion)
    elif backtrack[i][j] == 'd':
        lcs = recur_output_LCS(backtrack, v, i, j-1)
    # if previous move was diagonal (match)
    else:
        lcs = recur_output_LCS(backtrack, v, i-1, j-1)
        lcs += v[i-1]
    return lcs

def iter_output_LCS(backtrack, v, i, j):
    """
    Iterative version of recur_output_LCS
    """
    lcs = ""
    while i != 0 and j != 0:
        if backtrack[i][j] == 'i':
            i -= 1
        elif backtrack[i][j] == 'd':
            j -= 1
        else:
            i -= 1
            j -= 1
            lcs = v[i] + lcs
    return lcs

def LCS(v, w):
    """
    :param v: string
    :param w: string
    :return: string for longest common subsequence between strings v and w
    """
    backtrack = LCS_backtrack(v, w)
    return iter_output_LCS(backtrack, v, len(v), len(w))

def driver(path):
    with open(path, 'r') as f:
        str1 = next(f).strip("\n")
        str2 = next(f).strip("\n")
        with open('output.txt', 'w') as f2:
            print(LCS(str1, str2), file=f2)

# writes output of test file to output.txt file
# driver('rosalind_ba5c.txt')
driver('test.txt')




