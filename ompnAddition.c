#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <time.h>
//#include "coordReader.c"
#include <omp.h>

extern int readNumOfCoords(char *filename);
extern double **readCoords(char *filename, int numOfCoords);
extern void *writeTourToFile(int *tour, int tourLength, char *filename);


// Function to calculate the Euclidean distance between all pairs of coordinates and store them in a matrix.
void matrixDistanceCalc(double *flat_coords, int numOfCoords, double **distanceMatrix) {
#pragma omp for
    for (int i = 0; i < numOfCoords; i++) {
        for (int j = i + 1; j < numOfCoords; j++) {
            double x1 = flat_coords[2 * i];     // x coordinate of the i-th point
            double y1 = flat_coords[2 * i + 1]; // y coordinate of the i-th point
            double x2 = flat_coords[2 * j];     // x coordinate of the j-th point
            double y2 = flat_coords[2 * j + 1]; // y coordinate of the j-th point

            double distance = sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
            distanceMatrix[i][j] = distance;
            distanceMatrix[j][i] = distance;
        }
        distanceMatrix[i][i] = 0;
    }
}
extern void matrixDistanceCalc2(double **coords, int numOfCoords, double **distanceMatrix){
#pragma omp parallel for
    for (int i = 0; i < numOfCoords; i++) {
        for (int j = i + 1; j < numOfCoords; j++) {
            double x1 = coords[i][0];
            double y1 = coords[i][1];
            double x2 = coords[j][0];
            double y2 = coords[j][1];

            // Calculate the Euclidean distance
            double distance = sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));

            distanceMatrix[i][j] = distance;
            distanceMatrix[j][i] = distance;
        }
        // Set diagonal elements to 0
        distanceMatrix[i][i] = 0;
    }
}
void nearestAdditionStartingFrom(double **distanceMatrix, int numOfCoords, int *tour, int startVertex) {

    int visited[numOfCoords]; // Using an integer array
//#pragma omp parallel for
    for (int i = 0; i < numOfCoords; i++) {
        visited[i] = 0; // Initialize to 0 (not visited)
    }

    tour[0] = startVertex;
    visited[startVertex] = 1; // Mark as visited by setting to 1
    int tourLength = 1;

    double nearestDist = DBL_MAX;
    int nearestVertex = -1;
    for (int i = 0; i < numOfCoords; i++) {
        if (i != startVertex && distanceMatrix[startVertex][i] < nearestDist) {
            nearestDist = distanceMatrix[startVertex][i];
            nearestVertex = i;
        }
    }

    if (nearestVertex >= 0) {
        tour[1] = nearestVertex;
        visited[nearestVertex] = 1; // Mark as visited
        tourLength++;
    }

    if (numOfCoords > 2) {
        tour[2] = startVertex;
    }

    while (tourLength < numOfCoords) {
        double globalMinAdditionCost = DBL_MAX;
        int minPosition = -1;
        int nearestNode = -1;

        double *localMinAdditionCosts = malloc(omp_get_max_threads() * sizeof(double));
        int *localMinPositions = malloc(omp_get_max_threads() * sizeof(int));
        int *localNearestNodes = malloc(omp_get_max_threads() * sizeof(int));

        for (int y = 0; y < omp_get_max_threads(); y++) {
            localMinAdditionCosts[y] = DBL_MAX;
            localMinPositions[y] = -1;
            localNearestNodes[y] = -1;
        }

        int i=0, j=0;
#pragma omp parallel for collapse(2) private(i, j)
        for (int i = 0; i < tourLength; i++) {
            for (int j = 0; j < numOfCoords; j++) {
                if (!visited[j]) {
                    int threadID = omp_get_thread_num();
                    double additionalCost = distanceMatrix[tour[i]][j];
                    if (additionalCost < localMinAdditionCosts[threadID]) {
                        localMinAdditionCosts[threadID] = additionalCost;
                        localMinPositions[threadID] = i;
                        localNearestNodes[threadID] = j;
                    }
                }
            }
        }

        // Combine the results of all threads
        for (int y = 0; y < omp_get_max_threads(); y++) {
            if (localMinAdditionCosts[y] < globalMinAdditionCost) {
                globalMinAdditionCost = localMinAdditionCosts[y];
                minPosition = localMinPositions[y];
                nearestNode = localNearestNodes[y];
            }
        }



        int indexBefore = (minPosition == 0) ? tourLength - 1 : minPosition - 1;
        int indexAfter = (minPosition == tourLength - 1) ? 0 : minPosition + 1;

// Ensure safe array access
        if (indexAfter >= numOfCoords) {
            indexAfter = numOfCoords - 1;
        }

        double distanceAfter = distanceMatrix[tour[minPosition]][nearestNode] +
                               distanceMatrix[nearestNode][tour[indexAfter]] -
                               distanceMatrix[tour[minPosition]][tour[indexAfter]];

        double distanceBefore = distanceMatrix[tour[indexBefore]][nearestNode] +
                                distanceMatrix[nearestNode][tour[minPosition]] -
                                distanceMatrix[tour[indexBefore]][tour[minPosition]];

        int insertPosition;
        if (distanceAfter < distanceBefore) {
            insertPosition = (minPosition == tourLength - 1) ? tourLength : indexAfter;
        } else {
            insertPosition = indexBefore + 1;
        }

        for (int i = tourLength; i > insertPosition; i--) {
            tour[i] = tour[i - 1];
        }
        tour[insertPosition] = nearestNode;
        visited[nearestNode] = 1; // Mark as visited
        tourLength++;
    }

    if (tourLength < numOfCoords) {
        tour[tourLength] = startVertex;
        tourLength++;
    }
    tour[numOfCoords] = startVertex; // Complete the tour by returning to the starting vertex.

}


int nearestAddition(char *argv[],double **distanceMatrix, int numOfCoords, int startVertex) {
    // Record the start time for execution time measurement.
    clock_t start_time = clock();

    argv[1] = "512_coords.coord";
    argv[2] = "natout.dat";

    // Read the number of coordinates from the file.
    //int numOfCoords = readNumOfCoords(argv[1]);

    if (numOfCoords < 0) {
        printf("Unable to open the file: %s\n", argv[1]);
        return 1;
    }

    // Read the coordinates from the file.
    double **coords = readCoords(argv[1], numOfCoords);
    // After reading coordinates from the file


    // Allocate memory for the distance matrix.
    //double **distanceMatrix = (double **)malloc(numOfCoords * sizeof(double *));

    for (int i = 0; i < numOfCoords; i++) {
        distanceMatrix[i] = (double *)malloc(numOfCoords * sizeof(double));
    }

    // Calculate the distance matrix based on the Euclidean distances between coordinates.
    matrixDistanceCalc2(coords, numOfCoords, distanceMatrix);

    // Allocate memory for the shortest tour array.
    int *shortestTour = (int *)malloc((numOfCoords + 1) * sizeof(int));
    double minTotalDistance = DBL_MAX;
#pragma omp parallel
    {
        double localMinTotalDistance = DBL_MAX;
        int *localShortestTour = (int *)malloc((numOfCoords + 1) * sizeof(int));

#pragma omp for nowait
        for (int startVertex = 0; startVertex < numOfCoords; startVertex++) {
            int *tour = (int *)malloc((numOfCoords + 1) * sizeof(int));
            nearestAdditionStartingFrom(distanceMatrix, numOfCoords, tour, startVertex);

            double totalDistance = 0.0;
            for (int i = 0; i < numOfCoords; i++) {
                totalDistance += distanceMatrix[tour[i]][tour[(i + 1) % numOfCoords]];
            }

            if (totalDistance < localMinTotalDistance) {
                localMinTotalDistance = totalDistance;
                for (int i = 0; i <= numOfCoords; i++) {
                    localShortestTour[i] = tour[i];
                }
            }
            free(tour);
        }

#pragma omp critical
        {
            if (localMinTotalDistance < minTotalDistance) {
                minTotalDistance = localMinTotalDistance;
                for (int i = 0; i <= numOfCoords; i++) {
                    shortestTour[i] = localShortestTour[i];
                }
            }
        }

        free(localShortestTour);
    }



    // Print the order of the shortest tour.
    printf("Shortest Tour Order: ");
    for (int i = 0; i <= numOfCoords; i++) {
        printf("%d ", shortestTour[i]);
    }
    printf("\n");

    // Write the shortest tour to a data file.
    writeTourToFile(shortestTour, numOfCoords , argv[2]);
    printf("Writing shortest tour to data file: %s \n", argv[2]);

    // Record the end time for execution time measurement.
    clock_t end_time = clock();
    double execution_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;
    printf("Execution time: %.6f seconds\n", execution_time);

    // Free allocated memory.
    for (int i = 0; i < numOfCoords; i++) {
        free(distanceMatrix[i]);
    }
    free(distanceMatrix);
    free(coords);
    free(shortestTour);

    return 0;
}