import numpy as np

if __name__ == '__main__':
    pointsDistanceMatrix = np.zeros((311, 311))
    row_count = 0
    with open("./pointsDistanceMatrix.csv","r") as file_distance_matrix:
        for line in file_distance_matrix:
            row = list(map(int, line.strip().split(',')))
            pointsDistanceMatrix[row_count] = row
            row_count += 1

    print("min: " + str(np.amin(pointsDistanceMatrix)))
    print("max: " + str(np.amax(pointsDistanceMatrix)))

    for i in range(311):
        if pointsDistanceMatrix[i][i] != 0:
            print("Problem hai at " + str(i))
