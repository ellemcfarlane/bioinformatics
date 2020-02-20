# Elle McFarlane
import pickle

def load_ob(name):
    with open(name, 'rb') as f:
        return pickle.load(f)

def affine_gap(v, w):
    """
    Constructs maximum alignment score between v and w, followed by an alignment of v and w achieving this maximum score.
    using affine gap penalties (with gap opening penalty of 11 and gap extension penalty of 1).

    *Note: Scoring matrix used for matches/mismatches is shown in BLOSUM62.txt file*

    :param v: string representing amino acids (max length 100)
    :param w: string representing amino acids (max length 100)
    :return: an int score (length of total global alignment path based on score) and matrix
    backtrack (sequence of steps for all paths)

    Example:
    affine_gap(PRTEINS, PRTWPSEIN) ->  8
                                       PRT---EINS
                                       PRTWPSEIN-
    """
    if len(v) > 100 or len(w) > 100:
        print("One or both inputs is over 100. Must be under 101.")
        return
    # E penalty for gap extension
    E = 1
    # G penalty for gap opening
    G = 7

    # load list of list holding scores for matches/mismatches
    scoring_matrix = load_ob("scoringmatrix.p")

    # Alphabet used for indexing the scoring_matrix
    Alphabet = "A  C  D  E  F  G  H  I  K  L  M  N  P  Q  R  S  T  V  W  Y".split()

    # lengths of v and w plus 1 for the empty string ""
    vsz = len(v) + 1
    wsz = len(w) + 1

    # lower matrix counts optimal path length ending with a deletion (downward)
    # upper matrix counts optimal path length ending with an insertion (rightward)
    # middle matrix counts optimal path length ending with a match/mismatch (diagonal)
    lower = [x[:] for x in [[0] * wsz] * vsz]
    upper = [x[:] for x in [[0] * wsz] * vsz]
    middle = [x[:] for x in [[0] * wsz] * vsz]

    # backtrack keeps track of the direction of the alignment step (down, right, or diagonal)
    backtrack_d = [x[:] for x in [[""] * wsz] * vsz]
    backtrack_v = [x[:] for x in [[""] * wsz] * vsz]
    backtrack_h = [x[:] for x in [[""] * wsz] * vsz]

    # initialize completely vertical paths to lengths of multiples of E (gap extension penalty)
    # initialize backtrack for diagonal and vertical paths to be all 'v' for deletions (downward movement)
    for i in range(1, vsz):
        lower[i][0] = lower[i-1][0] - E
        backtrack_d[i][0] = 'v'
        backtrack_v[i][0] = 'v'
    # initialize completely horizontal paths to lengths of multiples of E (gap extension penalty)
    # initialize backtrack for diagonal and horizontal paths to be all 'h' for insertions (rightward movement)
    for j in range(1, wsz):
        upper[0][j] = upper[0][j-1] - E
        backtrack_d[0][j] = 'h'
        backtrack_h[0][j] = 'h'

    # get optimal alignment path lengths for all motions
    for i in range(1, vsz):
        for j in range(1, wsz):
            # get index for current letter of strings v and w
            vi = Alphabet.index(v[i-1])
            wj = Alphabet.index(w[j-1])
            # compute optimal (max) alignment path length for each matrix
            score = 2 if v[i-1] == w[j-1] else -1
            # save value for lower (vertical) matrix
            lower[i][j] = max(lower[i-1][j] - E, middle[i-1][j] - G)
            # save the direction this lower value came from
            if lower[i][j] == lower[i-1][j] - E:
                backtrack_v[i][j] = 'v'  # came from vertical matrix (lower matrix)
            else:
                backtrack_v[i][j] = 'd'  # came from diagonal matrix (middle matrix)
            # save value for upper (horizontal) matrix
            upper[i][j] = max(upper[i][j-1] - E, middle[i][j-1] - G)
            # save direction this upper value came from
            if upper[i][j] == upper[i][j-1] - E:
                backtrack_h[i][j] = 'h'  # came from horizontal matrix (upper matrix)
            else:
                backtrack_h[i][j] = 'd'  # came from diagonal matrix (middle matrix)
            # save value for middle (diagonal) matrix
            middle[i][j] = max(lower[i][j],
                               middle[i-1][j-1] + score,
                               upper[i][j])
            # save direction this middle value came from
            # came from vertical matrix
            new_path_sum = middle[i][j]
            if new_path_sum == lower[i][j]:
                backtrack_d[i][j] = 'v'
            # came from horizontal matrix
            elif new_path_sum == upper[i][j]:
                backtrack_d[i][j] = 'h'
            # came from diagonal matrix
            else:
                backtrack_d[i][j] = 'd'
    # get maximum global alignment score
    score = middle[vsz-1][wsz-1]
    # get alignment from this maximum score
    al1, al2 = find_alignment(backtrack_d, backtrack_v, backtrack_h, v, w, len(v), len(w))
    printMat(lower)
    printMat(middle)
    printMat(upper)
    return score, al1, al2

def find_alignment(backtrack_d, backtrack_v, backtrack_h, v, w, i, j):
    """
    Returns the alignments for strings v and w based off of backtrack matrix

    :param backtrack_d,_v,_h: matrices of paths from each index of v to each index of w
    :param v: string (max length 100)
    :param w: string (max length 100)
    :param i: length of v used as starting position
    :param j: length of w used as starting position
    :return: str1, str2 the string representing optimal alignment of v and w
    """
    str1 = ""
    str2 = ""
    curr_matrix = 'd'
    next_matrix = ''
    while i > 0 or j > 0:
        # if in diagonal matrix
        if curr_matrix == 'd':
            # get previous direction in diagonal matrix
            next_matrix = backtrack_d[i][j]
            # if previous direction in path is still in diagonal matrix, do not switch matrices
            if next_matrix == 'd':
                str1 = v[i-1] + str1
                str2 = w[j-1] + str2
                i -= 1
                j -= 1
        # if in vertical matrix
        elif curr_matrix == 'v':
            # get previous direction in vertical matrix
            next_matrix = backtrack_v[i][j]
            str2 = '-' + str2
            str1 = v[i - 1] + str1
            i -= 1
        # if in horizontal matrix
        elif curr_matrix == 'h':
            # get previous direction in horizontal matrix
            next_matrix = backtrack_h[i][j]
            str1 = '-' + str1
            str2 = w[j - 1] + str2
            j -= 1
        # go to previous direction in path (next_matrix)
        curr_matrix = next_matrix
    return str1, str2

def printMat(mat):
    for line in mat:
        print(line)

def driver(path):
    with open(path, 'r') as f:
        v = next(f).strip()
        w = next(f).strip()
        for x in affine_gap(v, w):
            print(x)

# Run test files
driver("smalltest2")


