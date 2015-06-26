#include <stdio.h>
#include <math.h>
#include <stdlib.h>

__global__ void jacobi(double *dev_A, double *dev_V, int *dev_pair, int size);
__global__ void check(double *dev_A, int n, double tolerance, int *cont);

int main (void)
{

	double tolerance = 0.000000000001;
	int n = 32;

	int *cont;
	cudaMallocManaged((void**) &cont, sizeof(int), cudaMemAttachGlobal);

	double** A;
	cudaMallocManaged((void**) &A, n*sizeof(double), cudaMemAttachGlobal);
	for (int i = 0; i < n; i++)
		cudaMallocManaged((void**) &A[i], n*sizeof(double), cudaMemAttachGlobal);
	double** V;
	cudaMallocManaged((void**) &V, n*sizeof(double), cudaMemAttachGlobal);
	for (int i = 0; i < n; i++)
		cudaMallocManaged((void**) &V[i], n*sizeof(double), cudaMemAttachGlobal);
	int** pair;
	cudaMallocManaged((void**) &pair, n*sizeof(int), cudaMemAttachGlobal);

	/*double *d_A, *d_V;
	int *d_pair, *d_indicator;
	cudaMalloc( (void**) &d_A, 1024*1024*sizeof(double));
	cudaMalloc( (void**) &d_V, 1024*1024*sizeof(double));
	cudaMalloc( (void**) &d_pair, n*sizeof(int));
	cudaMalloc( (void**) &d_indicator, 1024*1024*sizeof(int));*/

	/* enter a valid matrix A*/
	int row, col, i;
	for (row=0; row<n; row++)
	{
		for (col=0; col<n-1; col++)
			scanf("%lf,", A[row][col]);
		scanf("%lf\n", A[row][n-1]);
	}


	/*initializing vector matrix V */
	for (row = 0; row < n; row++) 
	{
		for (col = 0; col < n; col++) 
		{
			if (row == col)
			{
				V [row][col] = 1.0;
			}
			else
			{
				V [row][col] = 0.0;
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
	//cudaMemcpy(d_pair, pair, n*sizeof(int), cudaMemcpyHostToDevice);

	/*launch kernel here*/
	dim3 grid (1, 1, 1);
	dim3 block (n/2, 1, 1);

	while (cont != 0)
	{
		jacobi<<<grid, block>>>(A, V, pair, n);
		cont = 0;
		check<<<1, dim3(n, 1, 1)>>>(A, n, tolerance, cont);
		/*cudaMemcpy(indicator, d_indicator, 1024*1024*sizeof(int), cudaMemcpyDeviceToHost);
		for (row = 0; row<n; row++) {
			for (col = 0; col<n; col++)
				printf("%d ", *(indicator+row*n+col));
			printf("\n");
		}*/

		printf("%d\n", *cont);
		
	}

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

	cudaFree(A);
	cudaFree(V);
	cudaFree(pair);
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

__global__ void check (double *dev_A, int n, double tolerance, int *d_cont)
{
	if (threadIdx.x != threadIdx.y)
		if (*(dev_A + threadIdx.x * n + threadIdx.y) > tolerance)
		{
			//*(d_indicator + threadIdx.x * n + threadIdx.y) = 1;
			*d_cont = 1;
		}
		//else
			//*(d_indicator + threadIdx.x * n + threadIdx.y) = 1;
}






		

