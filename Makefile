# Define compilers
GCC = gcc
MPICC = mpicc
ICC = icc
MPIICC = mpiicc
CFLAGS = -fopenmp -O2 -std=c99

# Targets
all: gomp-only gcomplete iomp-only icomplete

OBJS = coordReader.c ompcInsertion.c ompfInsertion.c ompnAddition.c

gomp-only: main-openmp-only.c coordReader.c
	gcc -fopenmp main-openmp-only.c coordReader.c  -lm -o gomp-only

iomp-only: main-openmp-only.c coordReader.c
	icc -qopenmp main-openmp-only.c coordReader.c -lm -o iomp-only

gcomplete: main-mpi.c coordReader.c
	mpicc -fopenmp main-mpi.c coordReader.c ompcInsertion.c ompfInsertion.c ompnAddition.c -lm -o gcomplete

icomplete: main-mpi.c coordReader.c
	mpiicc -qopenmp main-mpi.c coordReader.c -lm -o icomplete