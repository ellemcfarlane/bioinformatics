def greedy_reversal_sort(chromosome):
    """
    Sorts chromosome until it becomes identity permutation via greedy reversal algorithm
    (sorts by position)
    Example:
    chromosome = [-3, 4, 1, 5, -2]
    greedy_reversal_sort(chromosome) =  [
                                        [-1, -4, 3, 5, -2],
                                        [1, -4, 3, 5, -2],
                                        [1, 2, -5, -3, 4],
                                        [1, 2, 3, 5, 4],
                                        [1, 2, 3, -4, -5],
                                        [1, 2, 3, 4, -5],
                                        [1, 2, 3, 4, 5]
                                        ]
    :param chromosome: list of unique integers (without any zeros)
    :return: list of permutations corresponding to steps/reversals taken to achieve the identity permutation
    """
    steps = []
    for idx, synteny_block in enumerate(chromosome):
        position = idx + 1
        if synteny_block != position:
            k_sorting_reversal(chromosome, position)
            steps.append(chromosome.copy())
            # if synteny_block is negative, flip sign
            if chromosome[idx] != abs(chromosome[idx]):
                chromosome[idx] = abs(chromosome[idx])
                steps.append(chromosome.copy())
    return steps

def k_sorting_reversal(chromosome, k):
    """
    Reverses sublist of chromosome to put the k value in the k-1 position and
    switches the sign of k if it is not positive.
    :param chromosome: list of numbers
    :param k: number to put in the k-1-th position
    :return: nothing (in-place function)
    """
    if k > len(chromosome) or k < 1:
        return chromosome
    idx = k - 1
    # find index of k (or -k)
    new_idx = chromosome.index(k) if k in chromosome else chromosome.index(k*-1)
    # do k-reversal to put k in k-1-th position
    partial_reverse(chromosome, idx, new_idx)

def partial_reverse(nums, from_, to):
    """
    Mutator method that reverses list of numbers between given indices (inclusive)
    and switches the sign of each number all in place
    :param nums: list to be reversed
    :param from_: starting index
    :param to: ending index (inclusive)
    :return: nothing
    """
    for idx in range(int((to-from_)/2)+1):
        nums[from_+idx], nums[to-idx] = -1*nums[to-idx], -1*nums[from_+idx]


def driver(path):
    """
    writes output of greedy_reversal_sort to output.txt with certain format
    :param path: string for input file path
    :return: output.txt file with greedy_reversal_sort of input
    """
    chromosome = []
    with open(path, 'r') as f:
        line = next(f).split()
        # remove parenthesis from first number
        line[0] = line[0][1:]
        # remove parenthesis from last number
        last_num_idx = len(line) - 1
        last_num_length = len(line[last_num_idx])
        line[last_num_idx] = line[last_num_idx][:last_num_length-1]
        # make into actual numbers
        for num in line:
            chromosome.append(int(num))
        # run greedy_reversal_sort on chromosome
        steps = greedy_reversal_sort(chromosome)
        for step in steps:
            print(step)
        # format output and write to output.txt
        with open('output.txt', 'w') as f:
            for step in steps:
                print('(', end="", file=f)
                for num_idx, num in enumerate(step):
                    print('+' if num == abs(num) else '-', end="", file=f)
                    if num_idx != len(step)-1:
                        print(abs(num), end=" ", file=f)
                    else:
                        print(abs(num), ")", sep="", file=f)

#driver('smalltest.txt')
driver('smalltest2.txt')
