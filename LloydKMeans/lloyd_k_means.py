# Elle McFarlane
import numpy as np

def lloyd_k_means(k, data):
    """
    k-means clustering heuristic that finds k clusters from data.

    :param k: int, number of centers to choose
    :param data: m-tuples (data points)
    :return: a unique list of k centers (points) resulting from applying the Lloyd algorithm
     to data and centers, where the first k points from Data are selected as the first k centers.
    """
    num_points = len(data)
    if num_points < k:
        raise Exception("Length of data cannot exceed cluster number k.")
    # choose first k points as initial centers
    centers = data[:k]
    # keep modulating centers until algorithm converges
    while True:
        # calculate new clusters based on centers
        clusters = centers_to_clusters(centers, data)
        # calculate new centers based on clusters
        new_centers = clusters_to_centers(clusters)
        # algorithm has converged if centers do not change from new clusters
        if centers == new_centers:
            return centers
        centers = new_centers

def centers_to_clusters(centers, data):
    """
    Returns k clusters for each of the k centers by assigning each point in data
    to the center it is closest to.
    :param centers: a unique list of centers (points) for each cluster
    :param k: number of clusters (should equal number of centers)
    :param data: m-tuples (data points)
    :return clusters: list of list of tuples
    """
    clusters = [[] for _ in centers]
    for point in data:
        closest_center = min(centers, key=lambda center: euclidian_dist(point, center))
        closest_center_idx = centers.index(closest_center)
        clusters[closest_center_idx].append(point)
    return clusters

def euclidian_dist(point1, point2):
    """
    Calculates euclidian distance between two points in m-dimensions.
    :param point1: m-tuple
    :param point2: m-tuple
    :return: float
    """
    m = len(point1)
    if m != len(point2):
        raise Exception("Points must have same number of dimensions.")
    square_difs_sum = 0
    for dimension in range(m):
        dif = (point1[dimension] - point2[dimension]) ** 2
        square_difs_sum += dif
    return np.sqrt(square_difs_sum)

def clusters_to_centers(clusters):
    """
    Finds the center of each cluster via its center of gravity
    :param clusters: list of tuples (points)
    :return: a unique list of centers (points) for each cluster
    """
    if not clusters:
        raise Exception("clusters must be non-empty.")
    return [center_of_gravity(cluster) for cluster in clusters]

def center_of_gravity(cluster):
    """
    Finds center of gravity point whose i-th coordinate is the average of the i-th
    coordinates in all points in cluster.
    :param cluster: m-tuples, points
    :return: m-tuple, center of gravity
    """
    if not cluster:
        raise Exception("cluster must have at least 1 point.")
    m = len(cluster[0])
    gravity_center = [0] * m
    for dimension in range(m):
        gravity_center[dimension] = np.mean([point[dimension] for point in cluster])
    return tuple(gravity_center)


def driver(path):
    with open(path, 'r') as file:
        k, m = [int(num) for num in next(file).strip().split()]
        points = []
        for line in file:
            point = tuple(float(num) for num in line.strip().split())
            points.append(point)
        centers = lloyd_k_means(k, points)
        for center in centers:
            for dim in range(m):
                print('{:.3f}'.format(center[dim]), end=" ")
            print("")
driver("rosalind_ba8c.txt")
