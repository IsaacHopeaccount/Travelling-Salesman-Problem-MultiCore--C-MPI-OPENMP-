#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
//#include "coordReader.c"
//#include "ompcInsertion.c"
//#include "ompfInsertion.c"
//#include "ompnAddition.c"
#include <omp.h>

extern void nearestAdditionStartingFrom(double **distanceMatrix, int numOfCoords, int *tour, int startVertex);
extern void farthestInsertionStartingFrom(double **distanceMatrix, int numOfCoords, int *tour, int startVertex);
extern void cheapestInsertionStartingFrom(double **distanceMatrix, int numOfCoords, int *tour, int startVertex);
extern void matrixDistanceCalc(double *flat_coords, int numOfCoords, double **distanceMatrix);

int readNumOfCoords(char *filename);
double **readCoords(char *filename, int numOfCoords);
void *writeTourToFile(int *tour, int tourLength, char *filename);

typedef struct {
    double totalDistance;
    int *tour;
} TourResult;

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    double start_time, end_time;
    start_time = MPI_Wtime();
    double tolerance = 1e-9;
    omp_set_num_threads(12);

   // argv[1] = "16_coords.coord";
  //  argv[2] = "cheapestMPI.dat";
 //   argv[3] = "farthestMPI.dat";
//    argv[4] = "nearestMPI.dat";

    int world_size, world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    int numOfCoords;
    double *flat_coords = NULL;
    double *allDistances = NULL;
    int *allTours = NULL;


    // Process 0 reads coordinates
    if (world_rank == 0) {
        numOfCoords = readNumOfCoords(argv[1]);
        printf("Number of coordinates read: %d\n", numOfCoords);

        double **coords_2d = readCoords(argv[1], numOfCoords);
        if (!coords_2d) {
            fprintf(stderr, "Failed to read coordinates\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        flat_coords = (double *) malloc(numOfCoords * 2 * sizeof(double));
        if (flat_coords == NULL) {
            fprintf(stderr, "Failed to allocate memory for flat_coords\n");
            exit(EXIT_FAILURE); // Or handle the error as needed
        }
        if (!flat_coords) {
            fprintf(stderr, "Failed to allocate memory for flat_coords\n");
            for (int i = 0; i < numOfCoords; i++) {
                free(coords_2d[i]);
            }
            free(coords_2d);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
#pragma omp parallel for
        for (int i = 0; i < numOfCoords; i++) {
            if (coords_2d[i] != NULL) {
                flat_coords[2 * i] = coords_2d[i][0];
                flat_coords[2 * i + 1] = coords_2d[i][1];
            }
        }
        free(coords_2d);
    }

    // Broadcast the number of coordinates to all processes
    MPI_Bcast(&numOfCoords, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Non-root processes allocate memory for flat_coords
    if (world_rank != 0) {
        flat_coords = (double *) malloc(numOfCoords * 2 * sizeof(double));
        if (!flat_coords) {
            fprintf(stderr, "Failed to allocate memory for flat_coords in non-root process\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }

    // Broadcast flat_coords to all processes
    MPI_Bcast(flat_coords, numOfCoords * 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Allocate memory for distance matrix
    double **distanceMatrix = (double **) malloc(numOfCoords * sizeof(double *));
    for (int i = 0; i < numOfCoords; i++) {
        distanceMatrix[i] = (double *) malloc(numOfCoords * sizeof(double));
    }


    matrixDistanceCalc(flat_coords, numOfCoords, distanceMatrix);

    int verticesPerProcess = numOfCoords / world_size;
    int remainder = numOfCoords % world_size;

// Each process gets 'verticesPerProcess' vertices plus one more if it's within the first 'remainder' processes
    int verticesForThisProcess = verticesPerProcess + (world_rank < remainder ? 1 : 0);

// Calculate the starting vertex for this process
    int startVertexForThisProcess = verticesPerProcess * world_rank + (world_rank < remainder ? world_rank : remainder);

    int endVertexForThisProcess = startVertexForThisProcess + verticesForThisProcess;


    double minDistanceCheapest = DBL_MAX;
    double minDistanceFarthest = DBL_MAX;
    double minDistanceNearest = DBL_MAX;

    int *shortestTourCheapest = (int *) malloc((numOfCoords + 1) * sizeof(int));
    int *shortestTourFarthest = (int *) malloc((numOfCoords + 1) * sizeof(int));
    int *shortestTourNearest = (int *) malloc((numOfCoords + 1) * sizeof(int));

   // printf("Process %d handling vertices from %d to %d\n", world_rank, startVertexForThisProcess,
          // endVertexForThisProcess - 1);

    TourResult *results = (TourResult *) malloc(
            (endVertexForThisProcess - startVertexForThisProcess) * sizeof(TourResult));


    // Each process calculates its part of tours
    for (int startVertex = startVertexForThisProcess; startVertex < endVertexForThisProcess; startVertex++) {
        int *tour = (int *) malloc((numOfCoords + 1) * sizeof(int));
        nearestAdditionStartingFrom(distanceMatrix, numOfCoords, tour, startVertex);
      //  printf("Process %d calculating tour starting from vertex %d\n", world_rank, startVertex);

        tour[numOfCoords] = tour[0];
        double totalDistance = 0.0;
        for (int i = 0; i < numOfCoords; i++) {
            totalDistance += distanceMatrix[tour[i]][tour[(i + 1) % numOfCoords]];
        }
       // printf("Start vertex %d has distance %lf\n", startVertex, totalDistance);

        if (totalDistance < minDistanceNearest- tolerance) {
            minDistanceNearest = totalDistance;
            memcpy(shortestTourNearest, tour, (numOfCoords + 1) * sizeof(int));
        }

        free(tour);
    }

    // Each process calculates its part of tours
    for (int startVertex = startVertexForThisProcess; startVertex < endVertexForThisProcess; startVertex++) {
        int *tour = (int *) malloc((numOfCoords + 1) * sizeof(int));
        farthestInsertionStartingFrom(distanceMatrix, numOfCoords, tour, startVertex);
       // printf("Process %d calculating tour starting from vertex %d\n", world_rank, startVertex);

        tour[numOfCoords] = tour[0];
        double totalDistance = 0.0;
        for (int i = 0; i < numOfCoords; i++) {
            totalDistance += distanceMatrix[tour[i]][tour[(i + 1) % numOfCoords]];
        }
       // printf("Start vertex %d has distance %lf\n", startVertex, totalDistance);

        if (totalDistance < minDistanceFarthest- tolerance) {
            minDistanceFarthest = totalDistance;
            memcpy(shortestTourFarthest, tour, (numOfCoords + 1) * sizeof(int));
        }

        free(tour);
    }

    // Each process calculates its part of tours
    for (int startVertex = startVertexForThisProcess; startVertex < endVertexForThisProcess; startVertex++) {
        int *tour = (int *) malloc((numOfCoords + 1) * sizeof(int));
        cheapestInsertionStartingFrom(distanceMatrix, numOfCoords, tour, startVertex);
        //printf("Process %d calculating tour starting from vertex %d\n", world_rank, startVertex);

        tour[numOfCoords] = tour[0];
        double totalDistance = 0.0;
        for (int i = 0; i < numOfCoords; i++) {
            totalDistance += distanceMatrix[tour[i]][tour[(i + 1) % numOfCoords]];
        }
        //printf("Start vertex %d has distance %lf\n", startVertex, totalDistance);

        if (totalDistance < minDistanceCheapest- tolerance) {
            minDistanceCheapest = totalDistance;
            memcpy(shortestTourCheapest, tour, (numOfCoords + 1) * sizeof(int));
        }

        free(tour);
    }

    double *allDistancesCheapest, *allDistancesFarthest, *allDistancesNearest;
    int *allToursCheapest, *allToursFarthest, *allToursNearest;
// Gather distances and tours
    if (world_rank == 0) {
        allDistancesCheapest = (double *) malloc(world_size * sizeof(double));
        allDistancesFarthest = (double *) malloc(world_size * sizeof(double));
        allDistancesNearest = (double *) malloc(world_size * sizeof(double));

        allToursCheapest = (int *) malloc(world_size * (numOfCoords + 1) * sizeof(int));
        allToursFarthest = (int *) malloc(world_size * (numOfCoords + 1) * sizeof(int));
        allToursNearest = (int *) malloc(world_size * (numOfCoords + 1) * sizeof(int));
    }

    // Gather distances and tours from all processes
    MPI_Gather(&minDistanceCheapest, 1, MPI_DOUBLE, allDistancesCheapest, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(shortestTourCheapest, numOfCoords + 1, MPI_INT, allToursCheapest, numOfCoords + 1, MPI_INT, 0, MPI_COMM_WORLD);

    MPI_Gather(&minDistanceFarthest, 1, MPI_DOUBLE, allDistancesFarthest, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(shortestTourFarthest, numOfCoords + 1, MPI_INT, allToursFarthest, numOfCoords + 1, MPI_INT, 0, MPI_COMM_WORLD);

    MPI_Gather(&minDistanceNearest, 1, MPI_DOUBLE, allDistancesNearest, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(shortestTourNearest, numOfCoords + 1, MPI_INT, allToursNearest, numOfCoords + 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (world_rank == 0) {
        // For Cheapest Insertion
        double shortestDistanceCheapest = allDistancesCheapest[0];
        int shortestTourIndexCheapest = 0;
        for (int i = 0; i < world_size; i++) {
            if (allDistancesCheapest[i] < shortestDistanceCheapest- tolerance) {
                shortestDistanceCheapest = allDistancesCheapest[i];
                shortestTourIndexCheapest = i;
            }
        }

        writeTourToFile(allToursCheapest + (shortestTourIndexCheapest * (numOfCoords + 1)), numOfCoords, argv[2]);
        printf("Cheapest Insertion tour written to file: %s\n", argv[2]);
    }
    if (world_rank == 0) {
        // For Farthest Insertion
        double shortestDistanceFarthest = allDistancesFarthest[0];
        int shortestTourIndexFarthest = 0;
        for (int i = 0; i < world_size; i++) {
            if (allDistancesFarthest[i] < shortestDistanceFarthest - tolerance) {
                shortestDistanceFarthest = allDistancesFarthest[i];
                shortestTourIndexFarthest = i;
            }
        }

        writeTourToFile(allToursFarthest + (shortestTourIndexFarthest * (numOfCoords + 1)), numOfCoords, argv[3]);
        printf("Cheapest Insertion tour written to file: %s\n", argv[3]);
    }
    if (world_rank == 0) {
        // For Nearest Insertion
        double shortestDistanceNearest = allDistancesNearest[0];
        int shortestTourIndexNearest = 0;
        for (int i = 0; i < world_size; i++) {
            if (allDistancesNearest[i] < shortestDistanceNearest - tolerance) {
                shortestDistanceNearest = allDistancesNearest[i];
                shortestTourIndexNearest = i;
            }
        }

        writeTourToFile(allToursNearest+ (shortestTourIndexNearest * (numOfCoords + 1)), numOfCoords, argv[4]);
        printf("Cheapest Insertion tour written to file: %s\n", argv[4]);
    }
    // Cleanup
    for (int i = 0; i < numOfCoords; i++) {
        free(distanceMatrix[i]);
    }
    free(distanceMatrix);
    free(flat_coords);
    if (world_rank == 0) {
        free(allDistancesCheapest);
        free(allDistancesFarthest);
        free(allDistancesNearest);
        free(allToursCheapest);
        free(allToursFarthest);
        free(allToursNearest);
    }

    free(shortestTourNearest);
    free(shortestTourCheapest);
    free(shortestTourFarthest);

    end_time = MPI_Wtime();
    printf("\nProgram took %lf seconds\n", (end_time - start_time));

    MPI_Finalize();
    return 0;
}