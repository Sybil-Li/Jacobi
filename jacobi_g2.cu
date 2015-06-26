#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <iostream>


__global__ void jacobi(double *dev_A, double *dev_V, int *dev_pair, int size, int *d_cont, int tolerance);
void check (double *A, int n, double tolerance, thrust::host_vector<int> H);


int main (void)
{

	double tolerance = 0.000000000001;
	int n = 32, cont = 1;

	thrust::host_vector<int> H(n*n);
	for (int i = 0; i < n*n; i++)
		H[i] = 0;
	thrust::device_vector<int> D = H;

	int *d_cont;
	cudaMalloc((void**) &d_cont, sizeof(int));
	cudaMemcpy(d_cont, &cont, sizeof(int), cudaMemcpyHostToDevice);

	double* A = (double*)malloc(1024*1024*sizeof(double));
	double* V = (double*)malloc(1024*1024*sizeof(double));
	int* pair = (int*)malloc(n*sizeof(int));

	double *d_A, *d_V;
	int *d_pair;
	cudaMalloc( (void**) &d_A, 1024*1024*sizeof(double));
	cudaMalloc( (void**) &d_V, 1024*1024*sizeof(double));
	cudaMalloc( (void**) &d_pair, n*sizeof(int));

	/* enter a valid matrix A*/
	int row, col, i;
	for (row=0; row<n; row++)
	{
		for (col=0; col<n-1; col++)
			scanf("%lf,", (A+row*n+col));
		scanf("%lf\n", (A+row*n+n-1));
	}

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
			}
			else
			{
				*(V + row * n + col) = 0.0;
			}
		}
	}

	/*copy matrix to device*/
	cudaMemcpy(d_V, V, 1024*1024*sizeof(double), cudaMemcpyHostToDevice);


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

	while (cont != 0)
	{
		jacobi<<<grid, block>>>(d_A, d_V, d_pair, n, d_cont, tolerance);
		cudaMemcpy(A, d_A, 1024*1024*sizeof(double), cudaMemcpyDeviceToHost);
		check(A, n, tolerance, H);
		for(int i = 0; i < H.size(); i++)
        	std::cout << "H[" << i << "] = " << H[i] << std::endl;
		cont =  thrust::reduce(H.begin(), H.end(), (int) 0, thrust::plus<int>());
		printf("%d\n", cont);
		
	}
	
	cudaMemcpy(pair, d_pair, n*sizeof(int), cudaMemcpyDeviceToHost);
	//for (int i = 0; i < n; i++)
		//printf("%d\n", *(pair + n));

	/*write matrix back to host*/
	cudaMemcpy(A, d_A, 1024*1024*sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(V, d_V, 1024*1024*sizeof(double), cudaMemcpyDeviceToHost);

	/*check result*/
	double* ans = (double*) malloc(n*sizeof(double));
	double norm = 0;

	for (row = 0; row<n; row++){
		for (col = 0; col<n; col++){
			if (row==col)
			{
				*(ans+row) = *(A+row*n+col);
				norm += (*(ans+row))*(*(ans+col));
				//printf("%lf ", *(A+row*n+col));
			}
			
			//printf("%lf", *(A+row*n+col));
		}
		//printf("\n");
	}
	norm = sqrt(norm);
	printf("Norm is %lf\n", norm);

	free(A);
	free(V);
	free(pair);
	cudaFree(d_A);
	cudaFree(d_V);
	cudaFree(d_pair);
}

__global__ void jacobi(double *dev_A, double *dev_V, int *dev_pair, int size, int *d_cont, int tolerance)
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
	for (i = 0; i < n; i++)
	{ 
		double Vpi = (*(dev_V + p * n + i))*c + (*(dev_V + q * n + i))*(-s);
		double Vqi = (*(dev_V + p * n + i))*s + (*(dev_V + q * n + i))*c;
		__syncthreads();
		*(dev_V + p * n + i) = Vpi;
		*(dev_V + q * n + i) = Vqi;
	}

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

void check (double *A, int n, double tolerance, thrust::host_vector<int> H)
{
	int i, j;
	for ( i = 0; i< n; i++)
	{
		for ( j = 0; j< n; j++)
		{
			if (i != j)
				if (*(A +i* n + j) > tolerance)	
					H[i* n + j] = 1;
				else
					H[i* n + j] = 0;
			else
				H[i * n + j] = 0;
		}
	}

}






		

