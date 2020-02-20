# Elle McFarlane
def cyclopeptide_sequencing(spectrum):
    """
    Given an ideal experimental spectrum, finds the cyclic peptide(s) whose
    theoretical spectrum matches the experimental spectrum.

    :param spectrum: a collection of (possibly repeated) integers corresponding to an ideal experimental spectrum
    :return:  every amino acid string peptide such that cyclospectrum(peptide) = spectrum (if such a string exists).

    Ex:
    cyclopeptide_sequencing([0, 113, 128, 186, 241, 299, 314, 427]) ->
    [[186,128,113], [186,113,128], [128,186,113], [128,113,186], [113,186,128] [113,128,186]]
    """
    # list containing only empty peptide
    peptides = [[]]
    # peptides that match spectrum
    winners = []
    # while peptides is nonempty
    while peptides:
        # expand peptides by one amino acid
        peptides = expand(peptides)
        # remove any peptides not consistent with spectrum
        peptides = [pep for pep in peptides if is_consistent(pep, spectrum)]
        for peptide in peptides:
            if mass(peptide) == parent_mass(spectrum):
                if cyclospectrum(peptide) == spectrum:
                    winners.append(peptide)
    winners.sort(reverse=True)
    return winners

def expand(peptides):
    """
    :param peptides: list of list of peptides where a peptide is a list of integers representing amino acid masses in Da
    :return: new list containing all possible extensions of peptides by a single amino acid mass
    Ex:
    expand([[]]) -> [[57], [71], [87], [97], [99], [101], [103], [113], [114], [115], [128], [129],\
                    [131], [137], [147], [156], [163], [186]]
    """
    if len(peptides) <= 0:
        return peptides
    new_peptides = []
    amino_acids = [57, 71, 87, 97, 99, 101, 103, 113, 114, 115, 128, 129, 131, 137, 147, 156, 163, 186]
    for peptide in peptides:
        # append each amino acid separately to each peptide
        for aa in amino_acids:
            # create copy of original peptide
            new_peptide = list(peptide)
            # expand peptide by one amino acid
            new_peptide.append(aa)
            # add to new peptide set
            new_peptides.append(new_peptide)
    return new_peptides


def mass(peptide):
    """
    :param peptide: list of integers representing mass in Da of each amino acid in peptide
    :return: sum of amino acid masses in Da

    Ex:
    Peptide VKLFPWFNQY = [99, 128, 113, 147, 97, 186, 147, 114, 128, 163]
    mass(VKLFPWFNQY) = 1322
    """
    amino_acid_sum = 0
    for aa in peptide:
        amino_acid_sum += aa
    return amino_acid_sum

def parent_mass(spectrum):
    """
    :param spectrum: collection of integers corresponding to ideal experimental spectrum
    :return: last entry, which should be the biggest value representing entire peptide

    Ex:
    parent_mass([0, 113, 128, 186, 241, 299, 314, 427]) -> 427
    """
    sz = len(spectrum)
    if sz <= 0:
        return 0
    return spectrum[sz-1]


def cyclospectrum(peptide):
    """
    :param peptide: list of integers representing mass in Da of each amino acid in cyclic peptide in natural order
    :return: sorted list representing cyclic spectrum

    Ex:
    NQEL = [114, 128, 129, 113]
    cyclospectrum(NQEL) -> [0, 113, 114, 128, 129, 227, 242, 242, 257, 355, 356, 370, 371, 484]
    """
    # create list of masses of all prefixes from peptide
    prefix_masses = [0]
    for i in range(len(peptide)):
        prefix_masses.append(prefix_masses[i] + peptide[i])
    # get peptide's mass for finding cyclic subpeptides later
    peptide_mass = mass(peptide)
    # use prex_masses to build full cyclopectrum
    cyclospec = [0]
    for i in range(len(peptide)):
        for j in range(i+1, len(peptide) + 1):
            sub_pep = prefix_masses[j]-prefix_masses[i]
            cyclospec.append(sub_pep)
            # add cyclic subpeptides if possible
            if i > 0 and j < len(peptide):
                cyclospec.append(peptide_mass-sub_pep)
    # sort spectrum
    cyclospec.sort()
    return cyclospec

def linear_spectrum(peptide):
    """
    :param peptide: list of integers representing mass in Da of each amino acid in linear peptide in natural order
    :return: sorted list representing linear spectrum

    Ex:
    NQEL = [114, 128, 129, 113]
    linear_spectrum(NQEL) -> [0, 113, 114, 128, 129, 242, 242, 257, 370, 371, 484]
    """
    # create list of masses of all prefixes from peptide
    prefix_masses = [0]
    for i in range(len(peptide)):
        prefix_masses.append(prefix_masses[i] + peptide[i])
    # use prex_masses to build full linear spectrum
    linear_spec = [0]
    for i in range(len(peptide)):
        for j in range(i+1, len(peptide) + 1):
            linear_spec.append(prefix_masses[j]-prefix_masses[i])
    # sort spectrum
    linear_spec.sort()
    return linear_spec

def is_consistent(peptide, spectrum):
    """
    :param peptide: list of integers representing mass in Da of each amino acid in peptide in natural order
    :param spectrum: a collection of (possibly repeated) integers corresponding to an ideal experimental spectrum
    :return: boolean, whether peptide's linear spectrum is a subset of spectrum

    Ex:
    tyrocidine_b1_spec = []
    VKF = [99, 128, 147]
    VKY = [99, 128, 163]
    is_consistent(VKF, tyrocidine_b1_spec) -> False
    is_consistent(VKY, tyrocidine_b1_spec) -> True
    """
    sub_spectrum = linear_spectrum(peptide)
    return all(x in spectrum for x in sub_spectrum)

def driver(path):
    with open(path, 'r') as f:
        line = next(f)
        spectrum = [int(x) for x in line.split(" ")]
        winners = cyclopeptide_sequencing(spectrum)
        for peptide in winners:
            formatted_peptide = "-".join(str(x) for x in peptide)
            print(formatted_peptide, end=" ")

# driver('rosalind_ba4e.txt')

