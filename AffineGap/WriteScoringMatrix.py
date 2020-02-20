import pickle

with open('BLOSUM62.txt') as f:
    lines = f.readlines()
    scoring_matrix = []
    for i in range(1, len(lines)):
        line = lines[i]
        scoring_matrix.append(line[1:].strip().split())

    with open('scoringmatrix.p', 'wb') as objf:
        pickle.dump(scoring_matrix, objf)

