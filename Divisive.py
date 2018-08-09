import numpy as np
from matplotlib import pyplot as plt
from scipy.cluster.hierarchy import dendrogram
from ete3 import Tree

# function to calculate the distance between two clusters
def calculateClusterDistance(pointsDistanceMatrix, cluster1, cluster2):
    # average linkage
    sumClusterInternalDistances = 0
    counter = 0
    for item1 in cluster1:
        for item2 in cluster2:
            #print(item1,item2)
            sumClusterInternalDistances += pointsDistanceMatrix[item1][item2]
            counter+=1

    return sumClusterInternalDistances/counter

# function to reduce the clusters to the required number of clusters
# n is the required number of clusters
'''
Divisive clustering steps:
    * Obtain average dissimilarity to other clusters - Take an average of the entire row
    * Pick the cluster with max dissimilarity. This is then considered removed.
    * Find average dissimilarity among remaining clusters. From these, take the max positive value and append to previous cluster.
    * Repeat above step till all the dissimilarities are negative. The divisive step is complete.
    * Repeat these steps till desired number of clusters are obtained. Split the largest cluster each time.
    * Radius of a cluster is determined by the max dissimilarity among its points. => largest distance cell
'''

# l = np.copy(a[i:j+1, i:j+1]) => split off sub-array, and deep copy
def clustering(pointsDistanceMatrix, n, numpoints):
    counter = (2 * numpoints - 1) - 1 # max index corresponds to largest cluster
    clusters = [] # a list that maintains the current state of the clusters
    #newickList = list(range(numpoints))
    clusters.append(list(range(numpoints)))
    currentDissimilarityValues = {} # dict that maintains average dissimilarity values
    splinterDifferenceValues = {}
    currentCluster = clusters[0] # the current cluster to be split

    Z = [] # Z matrix in order to form the dendrogram

    while len(clusters) != numpoints: # run till all clusters are singletons
        # print("Cluster state: " + str(clusters))
        print("Current cluster: " + str(currentCluster))
        print("Number of clusters: " + str(len(clusters)))
        clusterLen = len(currentCluster)
        currentDissimilarityValues = {} # clear the dict

        for point in currentCluster:
            # compute average dissimilarity for each point in cluster by looping over the cluster.
            total = 0
            for i in currentCluster:
                if i != point:
                    total += pointsDistanceMatrix[point][i]

            currentDissimilarityValues[point] = total / (len(currentCluster) - 1)

        maxDissimilarPoint = max(currentDissimilarityValues, key = lambda x: currentDissimilarityValues[x])

        newCluster = [maxDissimilarPoint] # a list denoting the new cluster
        currentCluster.remove(maxDissimilarPoint) # remove from old cluster

        # obtain other splinter groups
        flag_splinter_ongoing = True
        while flag_splinter_ongoing and len(currentCluster) >= 2: # no need to check for splinters if size 1
            # keep checking if any point wants to split
            splinterDifferenceValues = {} # clear the dict
            for point in currentCluster:
                remainingDissimilarity = 0
                splinterDissimilarity = 0

                # remaining point dissimilarities
                for i in currentCluster:
                    if i != point:
                        remainingDissimilarity += pointsDistanceMatrix[point][i]
                remainingDissimilarity = remainingDissimilarity / (len(currentCluster) - 1)

                # splinter group dissimilarities
                for i in newCluster:
                    splinterDissimilarity += pointsDistanceMatrix[point][i]
                splinterDissimilarity = splinterDissimilarity / len(newCluster)

                # Difference
                splinterDifferenceValues[point] = remainingDissimilarity - splinterDissimilarity

            # at this point the splinterDifferenceValues dict is completely populated
            # obtain splinterPoint
            splinterPoint = max(splinterDifferenceValues, key = lambda x : splinterDifferenceValues[x])
            # if difference is negative, abort
            if splinterDifferenceValues[splinterPoint] < 0:
                flag_splinter_ongoing = False
            else:
                # add point to new cluster, remove from old one
                newCluster.append(splinterPoint)
                currentCluster.remove(splinterPoint)

        # all splintering complete, divisive step can be added to dendrogram
        clusters.append(newCluster)

        # index of 1st cluster being merged, index of 2nd cluster being merged, distance between clusters, number of clusters merged
        Z.insert(0, [counter-1, counter-2, calculateClusterDistance(pointsDistanceMatrix, newCluster, currentCluster), clusterLen])

        # determine next cluster to divide based on the cluster radius. Largest is taken
        largestClusterRadius = 0
        for c in clusters:
            clusterRad = 0
            for point1 in c:
                for point2 in c:
                    clusterRad = max(clusterRad, pointsDistanceMatrix[point1][point2])
            # print(str(c) + ", radius " + str(clusterRad))
            largestClusterRadius = max(largestClusterRadius, clusterRad)
            # if this is the largest cluster, update currentCluster
            if clusterRad == largestClusterRadius:
                currentCluster = c
        print("Chosen cluster radius: " + str(largestClusterRadius))
        # currentCluster for next iter has been determined.
        counter -= 2
    # WHILE LOOP END

    return Z

if __name__ == '__main__':
    # reading dataset and storing the data points
    fo = open("./genomedata.csv", "r")
    num_points = 0
    for line in fo:
        if(line[0] == '>'):
            num_points += 1
    num_points -= 1 # because there is an extra '>' at the end of the csv

    # reading the point distance matrix
    pointsDistanceMatrix = np.zeros((num_points, num_points))
    row_count = 0
    with open("./pointsDistanceMatrix.csv","r") as file_distance_matrix:
        for line in file_distance_matrix:
            row = list(map(int, line.strip().split(',')))
            pointsDistanceMatrix[row_count, :] = row
            row_count += 1

    # make all values positive
    for i in range(num_points):
        for j in range(num_points):
            if i != j: # diagonal must be 0
                pointsDistanceMatrix[i][j] = 350 - pointsDistanceMatrix[i][j]

    Z = clustering(pointsDistanceMatrix, 1, num_points)
    #print(Z)
    #newick_string = newick_string.replace('[', '(')
    #newick_string = newick_string.replace(']', ')')
    #t = Tree(newick_string+";")

    # calculate full dendrogram
    plt.figure(figsize=(25, 10))
    plt.title('Hierarchical Clustering Dendrogram')
    plt.xlabel('sample index')
    plt.ylabel('distance')
    dendrogram(
        Z,
        leaf_rotation=90.,  # rotates the x axis labels
        leaf_font_size=8.,  # font size for the x axis labels
    )
    plt.show()
