#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<stdbool.h>



int readNumOfCoords(char *fileName);
double **readCoords(char *filename, int numOfCoords);

/*Gets the number of the coordinates in the file. Returns as a single integer*/
int readNumOfCoords(char *filename){
    FILE *file = fopen(filename, "r");
    int numOfCoords = 0;

    if(file == NULL){
        printf("Unable to open file: %s", filename);

        return -1;
    }

    char line[100];

    while(fgets(line, sizeof(line), file) != NULL){
        numOfCoords++;
    }

    return numOfCoords;
}

/*Gets the data from the file. Returns as an array of doubles, Ignores the first integer*/
double **readCoords(char *filename, int numOfCoords){
    FILE *file = fopen(filename,"r");
    int i;

    char line[100];

    if(file == NULL) {
        printf("Unable to open file: %s", filename);
        return NULL;
    }

    double **coords = (double **)malloc(numOfCoords * sizeof(double *));

    for(i = 0; i < numOfCoords; i++){
        coords[i] = (double *) malloc(2 * sizeof(int));
        if (coords[i] == NULL){
            perror("Memory Allocation Failed");
        }
    }

    int lineNum = 0;
    while(fgets(line, sizeof(line), file) != NULL){
        double x, y;
        if (sscanf(line, "%lf,%lf", &x, &y) == 2){
            coords[lineNum][0] = x;
            coords[lineNum][1] = y;
            lineNum++;
        }
    }

    return coords;
}

void *writeTourToFile(int *tour, int tourLength, char *filename){

    FILE *file = fopen(filename, "w");
    int i;

    if(file == NULL){
        printf("Unable to open file: %s", filename);
        return NULL;
    }


    fprintf(file, "%d \n", tourLength+1);

    printf("Writing output data\n");
    for(i=0; i < tourLength+1; i++) {
        fprintf(file, "%d ", tour[i]);
    }
   // fprintf(file, "0");
    return NULL;



}