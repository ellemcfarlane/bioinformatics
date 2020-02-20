# Elle McFarlane
from collections import defaultdict
import numpy as np

def neighbor_joining_wrapper(distances):
    """run main function neighbor_joining"""
    num_leaves = len(distances)
    new_node_id = num_leaves
    node_ids = [num for num in range(num_leaves)]
    return neighbor_joining(distances, num_leaves, new_node_id, node_ids)

def neighbor_joining(distances, n, new_node_id, node_ids):
    """
    :param distances: list of lists (nxn) of numbers representing dist between nodes at each index
    :param n: int specifying number of leaves
    :param new_node_id: int specifying id for new node (if there are originally 0,1,2 ids, then new id is 3)
    :param node_ids: list of labels for the nodes where labels can be strings or ints
    :return: adjacency list for the tree resulting from applying the neighbor-joining algorithm.
    Edge-weights are accurate to two decimal places.

    Example:
            n = 4
            distances =
            0   23  27  20
            23  0   30  28
            27  30  0   30
            20  28  30  0

            returns:

            0->4:8.00
            1->5:13.50
            2->5:16.50
            3->4:12.00
            4->5:2.00
            4->0:8.00
            4->3:12.00
            5->1:13.50
            5->2:16.50
            5->4:2.00
    """
    # when only 2 nodes are left, distance between them is known, so add to tree
    if n == 2:
        first_node = node_ids[0]
        sec_node = node_ids[1]
        dist = distances[0][1]
        # add edge from first to second node and second to first in Tree
        tree = defaultdict(dict)
        tree[first_node] = {sec_node: dist}
        tree[sec_node] = {first_node: dist}
        return tree
    # create NJ-matrix from distances
    neighbor_joining_matrix = to_neighbor_joining_matrix(distances)
    # get indices for closest pair of leaves in NJ-matrix
    min_first_node_idx, min_sec_node_idx = closest_pair(neighbor_joining_matrix)
    # get total distances for each leaf
    total_dists = total_distances(distances)
    # calculate delta for the closest pair of leaves
    delta = (total_dists[min_first_node_idx] - total_dists[min_sec_node_idx]) / (n - 2)
    # calculate limb length to both leaves in the pair for the final tree
    pair_dist = distances[min_first_node_idx][min_sec_node_idx]
    first_limb_length = (pair_dist + delta) / 2
    sec_limb_length = (pair_dist - delta) / 2
    # add new node by adding new row and column to distances
    add_row_col(distances)
    # calculate distances from new node to others (new node index is one more than previous max index n-1)
    new_node_idx = n
    for node_idx in range(n):
        distances[node_idx][new_node_idx] = (distances[node_idx][min_first_node_idx]
                                             + distances[node_idx][min_sec_node_idx]
                                             - distances[min_first_node_idx][min_sec_node_idx]) / 2
        # matrix is symmetric, so update reflection to same number
        distances[new_node_idx][node_idx] = distances[node_idx][new_node_idx]
    # remove rows and columns corresponding to the closest pairs
    # remove the larger index first so that the next index is still valid for the modified matrix
    distances = delete_row_col(delete_row_col(distances, max(min_first_node_idx, min_sec_node_idx)),
                               min(min_first_node_idx, min_sec_node_idx))
    # save ids for closest pair
    first_id = node_ids[min_first_node_idx]
    sec_id = node_ids[min_sec_node_idx]
    # delete their ids from the id list
    del node_ids[min_sec_node_idx]
    del node_ids[min_first_node_idx]
    # add new id to id list
    node_ids.append(new_node_id)
    # recursively do same for smaller distance matrix
    tree = neighbor_joining(distances, n - 1, new_node_id+1, node_ids)
    # now add known limb lengths to the tree
    tree[new_node_id][first_id] = first_limb_length
    tree[new_node_id][sec_id] = sec_limb_length
    tree[first_id][new_node_id] = first_limb_length
    tree[sec_id][new_node_id] = sec_limb_length

    return tree

def total_distances(distances):
    """
    Computes sum of distance from each node to all other nodes
    :param distances: list of lists (nxn) of numbers representing dist between nodes at each index
    :return: list of numbers
    """
    dists_to_others = []
    num_nodes = len(distances)
    # calculate total distance from each node to all other nodes
    for node_idx in range(num_nodes):
        dists_to_others.append(sum(distances[node_idx]))
    return dists_to_others


def to_neighbor_joining_matrix(distances):
    """
    Converts distances matrix to matrix representing neighbor_joining matrix (new matrix such that if
    distances is additive, the smallest element of the new matrix corresponds to a pair of neighboring leaves i and j
    in a tree representing the distances matrix).
    :param distances: list of lists of numbers representing dist between nodes at each index
    :return: list of lists of numbers
    """
    num_leaves = len(distances)
    total_dists = total_distances(distances)
    neighbor_joining_matrix = [[0] * num_leaves for leaf in range(num_leaves)]
    # fill in NJ-matrix per its definition
    for row_idx in range(num_leaves):
        for col_idx in range(num_leaves):
            if row_idx != col_idx:
                neighbor_joining_matrix[row_idx][col_idx] = (num_leaves - 2) * distances[row_idx][col_idx] \
                                                            - total_dists[row_idx] - total_dists[col_idx]
    return neighbor_joining_matrix


def closest_pair(neighbor_joining_matrix):
    """
    :param neighbor_joining_matrix: list of lists of numbers
    :return: indices of nodes that are the closest to each other in the given NJ-matrix.
    If there are more than one pairs, returns first pair seen (smallest indices).
    """
    min_first_node_idx = 0
    min_sec_node_idx = 0
    min_dist = 0
    num_leaves = len(neighbor_joining_matrix)
    for row_idx in range(num_leaves):
        for col_idx in range(num_leaves):
            # if dist is smaller than previous smallest dist, update min_dist, etc
            if row_idx != col_idx and neighbor_joining_matrix[row_idx][col_idx] < min_dist:
                min_dist = neighbor_joining_matrix[row_idx][col_idx]
                min_first_node_idx = row_idx
                min_sec_node_idx = col_idx
    return min_first_node_idx, min_sec_node_idx


def delete_row_col(distances, idx):
    """
    deletes row and column at given idx
    :param distances: list of lists (nxn) of numbers representing dist between nodes at each index
    :param idx: int for index of row and column to be deleted
    :return: list of lists (n-1xn-1) of numbers
    """
    num_rows = len(distances)
    if num_rows <= 0:
        return distances
    num_cols = len(distances[0])
    if num_rows != num_cols:
        print("Number of rows must equal numbers of columns.")
        return distances
    if idx < 0 or idx >= num_rows:
        print("Index outo of bounds")
        return distances
    # delete row at idx then column at idx
    new_dists = np.delete(np.delete(distances, idx, axis=0), idx, axis=1)
    return new_dists.tolist()


def add_row_col(distances):
    """
    adds a row of zeros and a column of zeros to distances
    :param distances: list of lists (nxn) of numbers representing dist between nodes at each index
    :return: list of lists (n+1xn+1) of numbers representing dist between nodes at each index

    Note: modifies distances input, does not return a new matrix
    """
    num_nodes = len(distances)
    # append extra column of zeros
    for row_idx in range(num_nodes):
        distances[row_idx].append(0.00)
    # append extra row of zeros
    distances.append([0.00] * (num_nodes + 1))
    return distances

def driver(path):
    with open(path, 'r') as file:
        # ignore first line for number of leaves; function calculates this already
        _ = int(next(file))
        dists = []
        for line in file:
            str_nums = line.strip().split()
            nums = [float(num) for num in str_nums]
            dists.append(nums)
        tree = neighbor_joining_wrapper(dists)
        total_nodes = len(tree)
        for node in range(total_nodes):
            adj_nodes = tree[node]
            for adj_node in adj_nodes:
                wt = tree[node][adj_node]
                print('{}->{}:{:.2f}'.format(node, adj_node, wt))

#driver("smalltest.txt")
#driver("rosalind_ba7e.txt")
