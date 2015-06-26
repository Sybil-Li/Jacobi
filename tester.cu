#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "jacobi_f.h"
#include "jacobi_cpu.h"

int main(int argc, char *argv[])
{

	if(argc != 2)
	{
		printf("not enough arguments. supply matrix size.");
		exit(1);
	}

	int size = atoi(argv[1]);
	
	double* M = (double*)malloc(size*size*sizeof(double));
	int row, col;
	for (row = 0; row < size; row++)
		for (col = 0; col <= row; col++)
			*(M+row*size+col) = (double)rand()/(double)RAND_MAX*100;
	
	for (row = 0; row < size; row++)
		for (col = 0; col < (size - row - 1); col++)
			*(M+col*size+row) = *(M+row*size+col);

	printf("For matrix size %d*%d\n", size, size);

	clock_t start = clock(), diff;
	jacobi_c(M,size);
	diff = clock() - start;

	double t_in_sec = (double)diff/(double)CLOCKS_PER_SEC;
	printf("Time taken for CPU jacobi: %f seconds.\n", t_in_sec);
	
	start = clock();
	jacobi_cu(M,size);
	diff = clock() - start;

	t_in_sec = (double)diff/(double)CLOCKS_PER_SEC;
	printf("Time taken for GPU jacobi: %f seconds.\n", t_in_sec);

	return 0;
}
