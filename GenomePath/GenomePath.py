# Elle McFarlane
import filecmp
import difflib


def genome_path_to_string(kmers):
    """
    Reconstructs dna string from genome path.
    :param kmers: A sequence of k-mers Pattern1, ... , Patternn such that the last k - 1 symbols\
     of Patterni are equal to the first k - 1 symbols of Patterni+1 for i from 1 to n-1.
    :return: A string Text of length k+n-1 where the i-th k-mer in Text is equal to Patterni for all i.
    """
    if len(kmers) <= 0:
        return ""

    dna = kmers[0]
    for i in range(1, len(kmers)):
        kmer = kmers[i]
        dna = dna + kmer[-1]

    return dna

def driver(path):
    with open(path, 'r') as f:
        kmers = []
        for line in f:
            kmers.append(line.strip("\n"))
        with open('out.txt', 'w') as f2:
            print(genome_path_to_string(kmers), file=f2)

# writes input and output sections from test file to two new files
def parse_test(path):
    with open(path, 'r') as f:
        # find output line index
        outputid = -1
        file_lines = f.readlines()
        for idx, line in enumerate(file_lines):
            if line == "Output\n":
                outputid = idx

        # write input file
        with open('input.txt', 'w') as f2:
            for line in file_lines[1:outputid]:
                print(line.strip("\n"), file=f2)

        # write output file
        with open('expout.txt', 'w') as f3:
            print(file_lines[outputid+1].strip().strip("\n"), file=f3)

parse_test('bigtest.txt')
driver('input.txt')
#driver('rosalind_ba3b.txt')
print(filecmp.cmp('out.txt', 'expout.txt'))

