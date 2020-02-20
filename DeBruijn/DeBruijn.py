# Elle McFarlane

import filecmp

def debruijn(k, Text):
    """
    returns de Bruijn graph of given Text in the form of an adjacency list
    :param k: integer for kmer length (length of string contained in each node)
    :param Text: string
    :return: de Bruijn graph as adjacency list (and nodes corresponding to each list)
    """

    # get set of prefix/suffix for all kmers in text
    patterns = []
    for i in range(len(Text)-k+1):
        kmer = Text[i:i+k]
        prefix = kmer[:-1]
        suffix = kmer[1:]
        if prefix not in patterns:
            patterns.append(prefix)
        if suffix not in patterns:
            patterns.append(suffix)

    # sort patterns lexicographically

    patterns.sort()
    # initialize adj_list by adding an empty list for each kmer
    adj_list = []
    for kmer in patterns:
        adj_list.append([])

    # fill adjacency list
    # (for each kmer in text, map its prefix to its suffix)
    for i in range(len(Text)-k+1):
        kmer = Text[i:i+k]
        prefix = kmer[:-1]
        suffix = kmer[1:]
        # use position in pattern as index for adjacency list
        list_id = patterns.index(prefix)
        adj_list[list_id].append(suffix)

    return adj_list, patterns

# formats output for test bank
def driver(path):
    with open(path, 'r') as f:
        k = int(next(f))
        text = next(f)[:-1]
        debgraph, nodes = debruijn(k, text)
        with open('out.txt', 'w') as f2:
            for node, lst in zip(nodes, debgraph):
                if len(lst) != 0:
                    lst.sort()
                    print(node, "->", ','.join(lst), file=f2)

# gets output section from test file
def get_out(path):
    with open(path, 'r') as f:
        with open('tstout.txt', 'w') as f2:
            for line in f.readlines()[4:]:
                print(line[:-1], file=f2)

# writes input and output sections from test file to two new files
def parse_test(path):
    with open(path, 'r') as f:
        # find output line index
        outputid = -1
        file_lines = f.readlines()
        for idx, line in enumerate(file_lines):
            if line == "Output:\n":
                outputid = idx

        # write input file
        with open('tstin.txt', 'w') as f2:
            for line in file_lines[1:outputid]:
                print(line.strip("\n"), file=f2)

        # write output file
        with open('tstout.txt', 'w') as f3:
            for line in file_lines[outputid+1:]:
                print(line.strip("\n"), file=f3)

#driver("rosalind_ba3d.txt")
parse_test('DebruijnFullTest.txt')
driver("tstin.txt")
print(filecmp.cmp('out.txt', 'tstout.txt'))