#include <stdio.h>
#include <math.h>
#include <stdlib.h>

__global__ void jacobi(double *dev_A, double *dev_V, int *dev_pair, int size, double eps);

int main (void)
{

	double tolerance = 0.0000000001;
	int n = 20;

	double* A = (double*)malloc(1024*1024*sizeof(double));
	double* V = (double*)malloc(1024*1024*sizeof(double));
	int* pair = (int*)malloc(n*sizeof(int));

	//double *d_A, *d_V;
	//int *d_pair;
	//cudaMalloc( (void**) &d_A, 1024*1024*sizeof(double));
	//cudaMalloc( (void**) &d_V, 1024*1024*sizeof(double));
	//cudaMalloc( (void**) &d_pair, n*sizeof(int));

	/* enter a valid matrix A*/
	int row, col, i;
	for (row=0; row<n; row++)
	{
		for (col=0; col<n-1; col++)
			scanf("%lf,", (A+row*n+col));
		scanf("%lf\n", (A+row*n+n-1));
	}

	//show matrix
	for (row = 0; row<n; row++){
		for (col = 0; col<n; col++){
			//if (row==col)
			//{
				printf("%lf\n", *(A+row*n+col));
			//}
		}
	}

	/*copy matrix to device*/
	//cudaMemcpy(d_A, A, 1024*1024*sizeof(double), cudaMemcpyHostToDevice);

	/*initializing vector matrix V */
	for (row = 0; row < n; row++) 
	{
		for (col = 0; col < n; col++) 
		{
			if (row == col)
				*(V + row * n + col) = 1.0;
			else
				*(V + row * n + col) = 0.0;
		}
	}

	/*copy matrix to device*/
	//cudaMemcpy(d_V, V, 1024*1024*sizeof(double), cudaMemcpyHostToDevice);

	/*initializing pair matrix*/
	for (i = 0; i < n; i++)
		*(pair + i) = i;

	/*copy matrix to device*/
	//cudaMemcpy(d_pair, pair, n*sizeof(int), cudaMemcpyHostToDevice);

	/*launch kernel here*/
	dim3 grid (1, 1, 1);
	dim3 block (10, 1, 1);
	jacobi<<<grid, block>>>(d_A, d_V, d_pair, n, tolerance);

	/*write matrix back to host*/
	cudaMemcpy(A, d_A, 1024*1024*sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(V, d_V, 1024*1024*sizeof(double), cudaMemcpyDeviceToHost);

	/*check result*/
	double* ans = (double*) malloc(n*sizeof(double));
	double norm = 0;

	for (row = 0; row<n; row++){
		for (col = 0; col<n; col++){
			//if (row==col)
			//{
				*(ans+row) = *(A+row*n+col);
				norm += (*(ans+row))*(*(ans+col));
				printf("%lf\n", *(A+row*n+col));
			//}
		}
	}
	norm = sqrt(norm);
	printf("Norm is %lf\n", norm);

	free(A);
	free(V);
	free(pair);
	//cudaFree(d_A);
	//cudaFree(d_V);
	//cudaFree(d_pair);
}

__global__ void jacobi(double *dev_A, double *dev_V, int *dev_pair, int size, double eps)
{
	short threadno, p, q, n, i, temp1, temp2, cont = 1;
	double c, s;
	threadno = threadIdx.y * blockDim.x + threadIdx.x;
	n = size;	

	__shared__ int pair[1024];
	*(pair + threadno) = *(dev_pair + threadno);
	*(pair + threadno + n/2) = *(dev_pair + threadno + n/2);
	

	while (cont == 1)
	{
		p = *(pair + threadno);
		q = *(pair + threadno + n/2);


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

		cont = 0;
		for (i = 0; i < n; i++){
			if (*(dev_A + p*n + i) > eps)
			{
				cont = 1;
				break;
			}
		}

		if (cont == 0) {
			for (i = 0; i < n; i++){
				if (*(dev_A + q*n + i) > eps)
				{
					cont = 1;
					break;
				}
			}
		}
	
		if (cont == 1)
		{
			if (threadno == 0)	
			{
				temp1 = 0;
				temp2 = *(pair + n/2 + 1);
			}
			else if (threadno == 1)
			{
				temp1 = *(pair + n/2);
				temp2 = *(pair + threadno + n/2 + 1);
			}
			else if (threadno == n/2 - 1)
			{
				temp1 = *(pair + threadno - 1);
				temp2 = *(pair + n - 1);
			}
			else
			{
				temp1 = *(pair + threadno - 1);
				temp2 = *(pair + threadno + n/2 + 1);
			}

			__syncthreads();

			*(pair + threadno) = temp1;
			*(pair + threadno + n/2) = temp2;
		}
	}
}		

