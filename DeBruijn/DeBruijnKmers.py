# Elle McFarlane

import filecmp

def debruijn(Patterns):
    """
    returns de Bruijn graph of given Patterns in the form of an adjacency list
    :param Patterns: list of string
    :return: de Bruijn graph as adjacency list (and nodes corresponding to each list)
    """

    # get set of prefix/suffix for all kmers in text
    nodes = []
    for kmer in Patterns:
        prefix = kmer[:-1]
        suffix = kmer[1:]
        if prefix not in nodes:
            nodes.append(prefix)
        if suffix not in nodes:
            nodes.append(suffix)

    # sort patterns lexicographically
    nodes.sort()

    # initialize adj_list by adding an empty list for each kmer
    adj_list = []
    for kmer in nodes:
        adj_list.append([])

    # fill adjacency list
    # (for each kmer in Patterns, map its prefix to its suffix)
    for kmer in Patterns:
        prefix = kmer[:-1]
        suffix = kmer[1:]
        # use position in pattern as index for adjacency list
        list_id = nodes.index(prefix)
        adj_list[list_id].append(suffix)

    return adj_list, nodes

# formats output for test bank
def driver(path):
    with open(path, 'r') as f:
        patterns = []
        for line in f:
            patterns.append(line.strip("\n"))
        debgraph, nodes = debruijn(patterns)
        with open('out.txt', 'w') as f2:
            for node, lst in zip(nodes, debgraph):
                if len(lst) != 0:
                    lst.sort()
                    print(node, "->", ','.join(lst), file=f2)

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

# driver("rosalind_ba3e.txt")
parse_test('DeBruijnKmersFullTest.txt')
driver('tstin.txt')
print(filecmp.cmp('out.txt', 'tstout.txt'))