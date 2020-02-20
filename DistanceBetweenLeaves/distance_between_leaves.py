# Elle McFarlane
from collections import deque
from collections import defaultdict
import heapq
def distance_between_leaves(n, weighted_tree):
    """
    Computes the distance between leaves in a weighted tree.
    :param n: number of leaves
    :param weighted_tree: adjacency list with n leaves
    :return: A n x n matrix (di, j), where di, j is the length of the path between leaves i and j.

    Example:

    n = 4,
    weighted_tree =     0->4:11
                        1->4:2
                        2->5:6
                        3->5:7
                        4->0:11
                        4->1:2
                        4->5:4
                        5->4:4
                        5->3:7
                        5->2:6

    sample output:
                        0   13  21  22
                        13  0   12  13
                        21  12  0   13
                        22  13  13  0
    """
    distances = []
    # get distance from each leaf to all other vertices
    for leaf in range(n):
        # dists include to internal nodes, so only take the leaves dists 0-n
        distances.append(distances_from(n, leaf, weighted_tree)[0:n])
    return distances

def distances_from(n, vertex, weighted_tree):
    """
    Calculates (shortest) distances from vertex to all other vertices in weighted_tree
    :param n: number of leaves
    :param vertex: start vertex to get distances from
    :param weighted_tree: adjacency list
    :return: list of distances from vertex
    """
    # leaves plus internal nodes
    total_vertices = len(weighted_tree)
    if vertex >= n or vertex < 0:
        return None
    # initially set all dists to infinite except for starting node, which is 0 dist
    dist_from_source = [float('inf')] * total_vertices
    dist_from_source[vertex] = 0
    # add source to fringe
    fringe = [(dist_from_source[vertex], vertex)]
    # keep track of vertices whose min distance from source is already found
    min_found = set()
    # dijkstra's algorithm to find shortest distance from source to all other vertices
    while fringe:
        dist_from_parent, closest_vertex = heapq.heappop(fringe)
        # skip if min dist for node already found (
        if closest_vertex in min_found:
            continue
        children = weighted_tree[closest_vertex]
        min_found.add(closest_vertex)
        # update distances from closest_vertex to its children and add to fringe
        for child in children:
            child_vert = child[0]
            child_dist_from_parent = child[1]
            dist_calc = dist_from_parent + child_dist_from_parent
            if dist_calc < dist_from_source[child_vert]:
                dist_from_source[child_vert] = dist_calc
                heapq.heappush(fringe, (dist_from_source[child_vert], child_vert))
    return dist_from_source

def driver(path):
    """loads input file data into distance_between_leaves"""
    with open(path, 'r') as file:
        n = int(next(file))
        adj_list = defaultdict(list)
        for line in file:
            data = line.strip().split("->")
            if len(data) >= 2:
                parent = data[0]
                child_data = data[1]
                parent = int(parent)
                child_vertex, child_weight = (int(elm) for elm in child_data.split(":"))
                adj_list[parent].append((child_vertex, child_weight))

        dists = distance_between_leaves(n, adj_list)
        # space out numbers a bit and print as proper matrix
        print('\n'.join([''.join(['{:<4}'.format(item) for item in row])
                         for row in dists]))

#driver("smalltest.txt")
driver('rosalind_ba7a.txt')