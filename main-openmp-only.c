#include <stdio.h>
#include <stdlib.h>
#include <string.h>
//#include "coordReader.c"
#include "ompcInsertion.c"
#include "ompfInsertion.c"
#include "ompnAddition.c"
#include <float.h>

extern void nearestAdditionStartingFrom(double **distanceMatrix, int numOfCoords, int *tour, int startVertex);
extern void farthestInsertionStartingFrom(double **distanceMatrix, int numOfCoords, int *tour, int startVertex);
extern void cheapestInsertionStartingFrom(double **distanceMatrix, int numOfCoords, int *tour, int startVertex);
extern void matrixDistanceCalc(double *flat_coords, int numOfCoords, double **distanceMatrix);

extern int readNumOfCoords(char *filename);
extern double **readCoords(char *filename, int numOfCoords);
extern void *writeTourToFile(int *tour, int tourLength, char *filename);

// Utility function to update the shortest tour
void updateShortestTour(int *tour, double **distanceMatrix, int numOfCoords, double *shortestDistance, int *shortestTour) {
    tour[numOfCoords] = tour[0]; // Complete the loop
    double tolerance = 1e-9;

    double totalDistance = 0.0;
    for (int i = 0; i < numOfCoords; i++) {
        totalDistance += distanceMatrix[tour[i]][tour[(i + 1) % numOfCoords]];
    }
    if (totalDistance < *shortestDistance- tolerance) {
        *shortestDistance = totalDistance;
        memcpy(shortestTour, tour, (numOfCoords + 1) * sizeof(int));
    }
}

int main(int argc, char *argv[]) {

    int numOfCoords;
    double *flat_coords = NULL;

//    argv[1]="16_coords.coord";
//    argv[2]="SerialCheapest.dat";
//    argv[3]="SerialFarthest.dat";
//    argv[4]="SerialNearest.dat";

    // Read coordinates
    numOfCoords = readNumOfCoords(argv[1]);
    printf("Number of coordinates read: %d\n", numOfCoords);

    double **coords_2d = readCoords(argv[1], numOfCoords);
    if (!coords_2d) {
        fprintf(stderr, "Failed to read coordinates\n");
        exit(EXIT_FAILURE);
    }

    flat_coords = (double *) malloc(numOfCoords * 2 * sizeof(double));
    if (!flat_coords) {
        fprintf(stderr, "Failed to allocate memory for flat_coords\n");
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < numOfCoords; i++) {
        if (coords_2d[i] != NULL) {
            flat_coords[2 * i] = coords_2d[i][0];
            flat_coords[2 * i + 1] = coords_2d[i][1];
        }
    }
    free(coords_2d);

    // Allocate memory for distance matrix
    double **distanceMatrix = (double **) malloc(numOfCoords * sizeof(double *));
    for (int i = 0; i < numOfCoords; i++) {
        distanceMatrix[i] = (double *) malloc(numOfCoords * sizeof(double));
    }

    matrixDistanceCalc(flat_coords, numOfCoords, distanceMatrix);

    double shortestDistanceCheapest = DBL_MAX;
    double shortestDistanceFarthest = DBL_MAX;
    double shortestDistanceNearest = DBL_MAX;
    int *shortestTourCheapest = (int *) malloc((numOfCoords + 1) * sizeof(int));
    int *shortestTourFarthest = (int *) malloc((numOfCoords + 1) * sizeof(int));
    int *shortestTourNearest = (int *) malloc((numOfCoords + 1) * sizeof(int));

    // Calculate tours
    for (int startVertex = 0; startVertex < numOfCoords; startVertex++) {
        int *tour = (int *) malloc((numOfCoords + 1) * sizeof(int));

        nearestAdditionStartingFrom(distanceMatrix, numOfCoords, tour, startVertex);
        updateShortestTour(tour, distanceMatrix, numOfCoords, &shortestDistanceNearest, shortestTourNearest);

        cheapestInsertionStartingFrom(distanceMatrix, numOfCoords, tour, startVertex);
        updateShortestTour(tour, distanceMatrix, numOfCoords, &shortestDistanceCheapest, shortestTourCheapest);

        farthestInsertionStartingFrom(distanceMatrix, numOfCoords, tour, startVertex);
        updateShortestTour(tour, distanceMatrix, numOfCoords, &shortestDistanceFarthest, shortestTourFarthest);

        free(tour);
    }

    // Write the shortest tour to file
    writeTourToFile(shortestTourCheapest, numOfCoords, argv[2]);
    writeTourToFile(shortestTourFarthest, numOfCoords, argv[3]);
    writeTourToFile(shortestTourNearest, numOfCoords, argv[4]);




    // Cleanup
    for (int i = 0; i < numOfCoords; i++) {
        free(distanceMatrix[i]);
    }
    free(distanceMatrix);
    free(flat_coords);
    free(shortestTourNearest);
    free(shortestTourCheapest);
    free(shortestTourFarthest);

    return 0;
}
