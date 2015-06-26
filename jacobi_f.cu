/*
real	0m1.637s
user	0m1.142s
sys		0m0.463s 

real	0m1.528s
user	0m0.931s
sys		0m0.561s
*/

#include "jacobi_f.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>


__global__ void check (double *dev_A, int n, double tolerance, int *d_cont/*, int *d_indicator*/);
__global__ void jacobi(double *dev_A, double *dev_V, int *dev_pair, int size);

void jacobi_cu( double* M, int n)
{

	int row, col, i;	
	double tolerance = 0.000000000001;
	int *d_cont, cont = 1;
	cudaMalloc((void**) &d_cont, sizeof(int));
	cudaMemcpy(d_cont, &cont, sizeof(int), cudaMemcpyHostToDevice);

	double* A = M;

	double* V = (double*)malloc(1024*1024*sizeof(double));
	int* pair = (int*)malloc(n*sizeof(int));

	double *d_A;
	double *d_V;
	int *d_pair;
	//int *d_indicator;
	cudaMalloc( (void**) &d_A, 1024*1024*sizeof(double));
	cudaMalloc( (void**) &d_V, 1024*1024*sizeof(double));
	cudaMalloc( (void**) &d_pair, n*sizeof(int));
	//cudaMalloc( (void**) &d_indicator, 1024*1024*sizeof(int));

	

	/*copy matrix to device*/
	cudaMemcpy(d_A, A, 1024*1024*sizeof(double), cudaMemcpyHostToDevice);

	/*initializing vector matrix V */
	for (row = 0; row < n; row++) 
	{
		for (col = 0; col < n; col++) 
		{
			if (row == col)
			{
				*(V + row * n + col) = 1.0;
				//*(indicator + row * n) = 0;
			}
			else
			{
				*(V + row * n + col) = 0.0;
				//*(indicator + row * n) = 0;
			}
		}
	}

	/*copy matrix to device*/
	//cudaMemcpy(d_V, V, 1024*1024*sizeof(double), cudaMemcpyHostToDevice);
	//cudaMemcpy(d_indicator, indicator, 1024*1024*sizeof(int), cudaMemcpyHostToDevice);


	/*initializing pair matrix*/
	for (i = 0; i < n; i++)
		*(pair + i) = i;

	//for (i = 0; i < n; i++)
		//printf("%d ", *(pair + i));

	/*copy matrix to device*/
	cudaMemcpy(d_pair, pair, n*sizeof(int), cudaMemcpyHostToDevice);

	/*launch kernel here*/
	dim3 grid (1, 1, 1);
	dim3 block (n/2, 1, 1);

	int iteration = 0;
	while ((cont != 0) && (iteration <= 100000))
	{
		jacobi<<<grid, block>>>(d_A, d_V, d_pair, n);
		cont = 0;
		cudaMemcpy(d_cont, &cont, sizeof(int), cudaMemcpyHostToDevice);
		check<<<16, dim3(n/16, 1, 1)>>>(d_A, n, tolerance, d_cont/*,d_indicator*/);
		cudaMemcpy(&cont, d_cont, sizeof(int), cudaMemcpyDeviceToHost);
		/*cudaMemcpy(indicator, d_indicator, 1024*1024*sizeof(int), cudaMemcpyDeviceToHost);
		for (row = 0; row<n; row++) {
			for (col = 0; col<n; col++)
				printf("%d ", *(indicator+row*n+col));
			printf("\n");
		}*/

		iteration++;
		
	}
	
	cudaMemcpy(pair, d_pair, n*sizeof(int), cudaMemcpyDeviceToHost);
	//for (int i = 0; i < n; i++)
		//printf("%d\n", *(pair + n));

	/*write matrix back to host*/
	cudaMemcpy(A, d_A, 1024*1024*sizeof(double), cudaMemcpyDeviceToHost);
	//cudaMemcpy(V, d_V, 1024*1024*sizeof(double), cudaMemcpyDeviceToHost);
	//cudaMemcpy(indicator, d_indicator, 1024*1024*sizeof(int), cudaMemcpyDeviceToHost);
	/*for (row = 0; row<n; row++) {
		for (col = 0; col<n; col++)
			printf("%d ", *(indicator+row*n+col));
		printf("\n");
	}*/

	/*check result*/
	double* ans = (double*) malloc(n*sizeof(double));
	//double norm = 0;

	for (row = 0; row<n; row++){
		for (col = 0; col<n; col++){
			if (row==col)
			{
				//*(ans+row) = *(A+row*n+col);
				//norm += (*(ans+row))*(*(ans+col));
				printf("%lf\n", *(A+row*n+col));
			}
			
			//printf("%lf", *(A+row*n+col));
		}
		//printf("\n");
	}
	//norm = sqrt(norm);
	//printf("Norm is %lf\n", norm);

	free(A);
	free(V);
	free(pair);
	cudaFree(d_A);
	cudaFree(d_V);
	cudaFree(d_pair);

}

__global__ void jacobi(double *dev_A, double *dev_V, int *dev_pair, int size)
{
	short threadno, p, q, n, i, temp1, temp2;
	double c, s;
	threadno = threadIdx.x;
	n = size;	

	p = *(dev_pair + threadno);
	q = *(dev_pair + threadno + n/2);


	/*calculate c, s value*/
	if (*(dev_A + p * n + q) != 0)
	{
		double torque, t;
		torque = ( *(dev_A + q * n + q) - *(dev_A + p * n + p))/(2*(*(dev_A + p * n + q)));
		if (torque >= 0)
		    t = 1/(torque + sqrt(1+torque*torque));
		else
		    t = -1/(-torque + sqrt(1+torque*torque));
		
		c = 1/sqrt(1+t*t);
		s = t*c;
	}
	else
	{
		c = 1;
		s = 0;
	}

	/* A = transpose(J)*A*J */
	for (i = 0; i < n; i++)
	{
		double Api = (*(dev_A + p * n + i))*c + (*(dev_A + q * n + i))*(-s);
		double Aqi = (*(dev_A + p * n + i))*s + (*(dev_A + q * n + i))*c;
		__syncthreads();
		*(dev_A + p * n + i) = Api;
		*(dev_A + q * n + i) = Aqi;
	}


	for (i = 0; i < n; i++)
	{ 
		double Aip = (*(dev_A + i * n + p))*c + (*(dev_A + i * n + q))*(-s);
		double Aiq = (*(dev_A + i * n + p))*s + (*(dev_A + i * n + q))*c;
		__syncthreads();
		*(dev_A + i * n + p) = Aip;
		*(dev_A + i * n + q) = Aiq;
	}
	 

	/* V = V*J */
	/*for (i = 0; i < n; i++)
	{ 
		double Vpi = (*(dev_V + p * n + i))*c + (*(dev_V + q * n + i))*(-s);
		double Vqi = (*(dev_V + p * n + i))*s + (*(dev_V + q * n + i))*c;
		__syncthreads();
		*(dev_V + p * n + i) = Vpi;
		*(dev_V + q * n + i) = Vqi;
	}*/

	/* chess tournament rotate*/
	if (threadno == 0)	
	{
		temp1 = 0;
		temp2 = *(dev_pair + n/2 + 1);
	}
	else if (threadno == 1)
	{
		temp1 = *(dev_pair + n/2);
		temp2 = *(dev_pair + threadno + n/2 + 1);
	}
	else if (threadno == n/2 - 1)
	{
		temp1 = *(dev_pair + threadno - 1);
		temp2 = *(dev_pair + n/2 - 1);
	}
	else
	{
		temp1 = *(dev_pair + threadno - 1);
		temp2 = *(dev_pair + threadno + n/2 + 1);
	}

	__syncthreads();

	*(dev_pair + threadno) = temp1;
	*(dev_pair + threadno + n/2) = temp2;
	
}

__global__ void check (double *dev_A, int n, double tolerance, int *d_cont/*, int *d_indicator*/)
{
	int threadno = blockIdx.x * n/16 + threadIdx.x;
	for (int i = 0; i < n; i++)
	{	
		if (threadno != i)
		{
			if (*(dev_A + threadno * n + i) > tolerance)
			{
				//*(d_indicator + threadno * n + i) = 1;
				*d_cont = 1;
			}
			//else
				//*(d_indicator + threadno * n + i) = 0;
		} 
	}
}






		

