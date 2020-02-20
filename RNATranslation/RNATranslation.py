# Elle McFarlane

def rna_translation(Pattern):
    """
    Simulates RNA translation to Peptide
    :param Pattern: string representing RNA
    :return: string, a translation of Pattern into an amino acid string Peptide

    ex:
    RNATranslation("AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA") -> MAMAPRTEINSTRING
    """

# create codons dictionary
    codons = {
        'AUA': 'I', 'AUC': 'I', 'AUU': 'I', 'AUG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACU': 'T',
        'AAC': 'N', 'AAU': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGU': 'S', 'AGA': 'R', 'AGG': 'R',
        'CUA': 'L', 'CUC': 'L', 'CUG': 'L', 'CUU': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCU': 'P',
        'CAC': 'H', 'CAU': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGU': 'R',
        'GUA': 'V', 'GUC': 'V', 'GUG': 'V', 'GUU': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCU': 'A',
        'GAC': 'D', 'GAU': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGU': 'G',
        'UCA': 'S', 'UCC': 'S', 'UCG': 'S', 'UCU': 'S',
        'UUC': 'F', 'UUU': 'F', 'UUA': 'L', 'UUG': 'L',
        'UAC': 'Y', 'UAU': 'Y', 'UAA': '*', 'UAG': '*',
        'UGC': 'C', 'UGU': 'C', 'UGA': '*', 'UGG': 'W'}

    Peptide = ""
    for i in range(0, len(Pattern)-3+1, 3):
        codon = Pattern[i:i+3]
        amino_acid = codons[codon]
        if amino_acid:
            if amino_acid == '*':
                return Peptide
            Peptide = Peptide + amino_acid

    return Peptide

def driver(path):
    with open(path, 'r') as f:
        rna = next(f)
        with open('out.txt', 'w') as f2:
            print(rna_translation(rna), file=f2)

driver("rosalind_ba4a.txt")



