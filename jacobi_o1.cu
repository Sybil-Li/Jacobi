#include <stdio.h>
#include <math.h>
#include <stdlib.h>

__global__ void jacobi(double *dev_A, double *dev_V, int *dev_pair, int size);
__global__ void check(double *dev_A, int n, double tolerance, int *d_cont/*, int* d_indicator*/);

__device__ double get(double *mat, int n, int row, int col);
__device__ void set(double newValue, double *mat, int n, int row, int col);

int main ( int argc, char *argv[] )
{
	if ( argc != 2 ) /* argc should be 2 for correct execution */
    {
        printf( "usage: %s filename", argv[0] );
		exit(1);
    }

	double tolerance = 0.000000000001;
	int n = atoi(argv[1]), cont = 1;

	int *d_cont;
	cudaMalloc((void**) &d_cont, sizeof(int));
	cudaMemcpy(d_cont, &cont, sizeof(int), cudaMemcpyHostToDevice);

	double* A = (double*)malloc(128*128*sizeof(double));
	double* V = (double*)malloc(128*128*sizeof(double));
	int* pair = (int*)malloc(n*sizeof(int));
	//int* indicator= (int*)malloc(1024*1024*sizeof(int));

	double *d_A, *d_V;
	int *d_pair;
	//int *d_indicator;
	cudaMalloc( (void**) &d_A, 128*128*sizeof(double));
	cudaMalloc( (void**) &d_V, 128*128*sizeof(double));
	cudaMalloc( (void**) &d_pair, n*sizeof(int));
	//cudaMalloc( (void**) &d_indicator, 1024*1024*sizeof(int));


	/* enter a valid matrix A*/
	int row, col, i=0;
	double garb;
	for (row = 0; row < n; row++)
	{
		printf("row=%d\n", row);
		for (col = 0; col < n; col++)
		{
			if ( col >= row)
			{
				scanf("%lf,", A+i);
				printf("%2.0lf,", *(A+i));
				i++;
			}
			else {
				scanf("%lf,", &garb);
				printf("%2.0lf,", garb);
			}
		}
	}
	printf("scan complete\n");

	/*copy matrix to device*/
	cudaMemcpy(d_A, A, 128*128*sizeof(double), cudaMemcpyHostToDevice);

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
	cudaMemcpy(d_V, V, 1024*1024*sizeof(double), cudaMemcpyHostToDevice);
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

	while (cont != 0)
	{
		jacobi<<<grid, block>>>(d_A, d_V, d_pair, n);
		cont = 0;
		cudaMemcpy(d_cont, &cont, sizeof(int), cudaMemcpyHostToDevice);
		check<<<4, dim3(n/4, 1, 1)>>>(d_A, n, tolerance, d_cont/*, d_indicator*/);
		cudaMemcpy(&cont, d_cont, sizeof(int), cudaMemcpyDeviceToHost);
		
	}
	
	cudaMemcpy(pair, d_pair, n*sizeof(int), cudaMemcpyDeviceToHost);
	//for (int i = 0; i < n; i++)
		//printf("%d\n", *(pair + n));

	/*write matrix back to host*/
	cudaMemcpy(A, d_A, 1024*1024*sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(V, d_V, 1024*1024*sizeof(double), cudaMemcpyDeviceToHost);
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
		torque = (get(dev_A, n, q, q) - get(dev_A, n, p, p))/(2*(get(dev_A, n, p, q)));
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
		double Api = (get(dev_A, n, p, i))*c + (get(dev_A, n, q, i))*(-s);
		double Aqi = (get(dev_A, n, p, i))*s + (get(dev_A, n, q, i))*c;
		__syncthreads();
		set(Api, dev_A, n, p, i);
		set(Aqi, dev_A, n, q, i);
	}


	for (i = 0; i < n; i++)
	{ 
		double Aip = (get(dev_A, n, i, p))*c + (get(dev_A, n, i, q))*(-s);
		double Aiq = (get(dev_A, n, i, p))*s + (get(dev_A, n, i, q))*c;
		__syncthreads();
		set(Aip, dev_A, n, i, p);
		set(Aiq, dev_A, n, i, q);
	}
	 

	/* V = V*J */
	for (i = 0; i < n; i++)
	{ 
		double Vpi = (get(dev_V, n, p, i))*c + (get(dev_V, n, q, i))*(-s);
		double Vqi = (get(dev_V, n, p, i))*s + (get(dev_V, n, q, i))*c;
		__syncthreads();
		set(Vpi, dev_V, n, p, i);
		set(Vqi, dev_V, n, p, i);
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

__global__ void check (double *dev_A, int n, double tolerance, int *d_cont/*, int *d_indicator*/)
{
	int threadno = blockIdx.x * n/4 + threadIdx.x;
	for (int i = 0; i < n; i++)
	{	
		if (threadno != i)
		{
			if (get(dev_A, n, threadno, i) > tolerance)
			{
				//*(d_indicator + threadno * n + i) = 1;
				*d_cont = 1;
			}
			//else
				//*(d_indicator + threadno * n + i) = 0;
		} 
	}
}

__device__ double get(double *mat, int n, int row, int col)
{
	if (row > col)
	{
		int temp = row;
		row = col;
		col = temp;
	}
	
	int buffer = row*(row + 1)/2;
	return *(mat + row * n + col - buffer);	
}

__device__ void set(double newValue, double *mat, int n, int row, int col)
{
	if (row > col)
	{
		int temp = row;
		row = col;
		col = temp;
	}
	
	int buffer = row*(row + 1)/2;
	*(mat + row * n + col - buffer) = newValue;
}




		

