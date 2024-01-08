#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <omp.h>
//#include "coordReader.c"

// Declare the functions from coordReader.c
extern int readNumOfCoords(char *filename);
extern double **readCoords(char *filename, int numOfCoords);
extern void *writeTourToFile(int *tour, int tourLength, char *filename);

// Function to calculate the Euclidean distance between all pairs of coordinates and store them in a matrix.
extern void matrixDistanceCalc2(double **coords, int numOfCoords, double **distanceMatrix);
extern void matrixDistanceCalc(double *flat_coords, int numOfCoords, double **distanceMatrix);
/*  {
    // Parallelise distance calculation
//#pragma omp parallel for
   for (int i = 0; i < numOfCoords; i++) {
        for (int j = i + 1; j < numOfCoords; j++) {
            double x1 = coords[i][0];
            double y1 = coords[i][1];
            double x2 = coords[j][0];
            double y2 = coords[j][1];
            double distance = sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));

            distanceMatrix[i][j] = distance;
            distanceMatrix[j][i] = distance;
        }
        // Set diagonal elements to 0
        distanceMatrix[i][i] = 0;
    }
}
    */

// Function to find the tour using the farthest insertion heuristic for a given starting vertex.
void farthestInsertionStartingFrom(double **distanceMatrix, int numOfCoords, int *tour, int startVertex) {
    // Initialize an array to keep track of visited vertices and the tour array.
    int visited[numOfCoords];
    for (int i = 0; i < numOfCoords; i++) {
        visited[i] = 0;
        tour[i] = -1;
    }

    // Start with the specified vertex as the initial tour
    int currentVertex = startVertex;
    tour[0] = currentVertex;
    visited[currentVertex] = 1;

    // Iterate over the remaining vertices to build the tour.
    for (int currentPosition = 1; currentPosition < numOfCoords; currentPosition++) {
        int farthestVertex = -1;
        double maxDistance = -1.0;

        // Iterate over all unvisited vertices and potential positions in the current tour.
        for (int vn = 0; vn < currentPosition; vn++) {
            for (int vk = 0; vk < numOfCoords; vk++) {
                if (!visited[vk]) {
                    double dist_vn_vk = distanceMatrix[tour[vn]][vk];
                    if (dist_vn_vk > maxDistance) {
                        maxDistance = dist_vn_vk;
                        farthestVertex = vk;
                    }
                }
            }
        }

        int bestPosition = -1;
        double minDistance = DBL_MAX; // Initialize with a large value

        // Find the best position to insert the farthest vertex into the current tour.
        for (int i = 0; i < currentPosition; i++) {
            int viPlus1 = tour[(i + 1) % currentPosition];
            double insertionCost = distanceMatrix[tour[i]][farthestVertex] + distanceMatrix[farthestVertex][viPlus1] - distanceMatrix[tour[i]][viPlus1];
            {
                if (insertionCost < minDistance) {
                    minDistance = insertionCost;
                    bestPosition = i;
                }
            }
        }

        // Update the tour with the farthest vertex at the best position.
        for (int i = currentPosition; i > bestPosition; i--) {
            tour[i] = tour[i - 1];
        }
        tour[bestPosition + 1] = farthestVertex;
        visited[farthestVertex] = 1;
    }
    tour[numOfCoords] = startVertex;
}

int farthestInsertion(char *argv[],double **distanceMatrix, int numOfCoords, int startVertex ) {

    argv[1] = "512_coords.coord";
    argv[2] = "ompfOut.dat";
    double tolerance = 1e-9;

    // Check if the correct number of command-line arguments is provided
    double start_time = omp_get_wtime(); // Capture the start time

    // Read the number of coordinates from the file.
   // int numOfCoords = readNumOfCoords(argv[1]);

    if (numOfCoords < 0) {
        printf("Unable to open the file: %s\n", argv[1]);
        return 1;
    }

    // Read the coordinates from the file.
    double **coords = readCoords(argv[1], numOfCoords);

    if (coords == NULL) {
        return 1;  // Error occurred while reading coordinates
    }

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
            farthestInsertionStartingFrom(distanceMatrix, numOfCoords, tour, startVertex);

            double totalDistance = 0.0;
            for (int i = 0; i < numOfCoords; i++) {
                totalDistance += distanceMatrix[tour[i]][tour[i + 1]];
            }

            if (totalDistance < localMinTotalDistance - tolerance) {
                localMinTotalDistance = totalDistance;
                for (int i = 0; i <= numOfCoords; i++) {
                    localShortestTour[i] = tour[i];
                }
            }
            free(tour);
        }

#pragma omp critical
        {
            if (localMinTotalDistance < minTotalDistance - tolerance) {
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

    // Free allocated memory.
    for (int i = 0; i < numOfCoords; i++) {
        free(distanceMatrix[i]);
    }
    free(distanceMatrix);
    free(coords);
    free(shortestTour);

    // Record the end time for execution time measurement.
    double end_time = omp_get_wtime(); // Capture the end time
    double execution_time = end_time - start_time;
    printf("Execution time: %f seconds\n", execution_time);

    return 0;
}
