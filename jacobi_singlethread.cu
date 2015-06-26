#include <stdio.h>
#include <math.h>
#include <stdlib.h>

__global__ void jacobi(double *dev_A, double *dev_V, int *dev_pair, int size, double eps);

int main (void)
{

	double tolerance = 0.000000000001;
	int n = 32, cont = 1;

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
	int row, col, i, j;
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
				*(V + row * n + col) = 1.0;
			else
				*(V + row * n + col) = 0.0;
		}
	}

	/*copy matrix to device*/
	cudaMemcpy(d_V, V, 1024*1024*sizeof(double), cudaMemcpyHostToDevice);

	/*initializing pair matrix*/
	for (i = 0; i < n; i++)
		*(pair + i) = i;

	/*copy matrix to device*/
	cudaMemcpy(d_pair, pair, n*sizeof(int), cudaMemcpyHostToDevice);

	/*launch kernel here*/
	dim3 grid (1, 1, 1);
	dim3 block (1, 1, 1);
	jacobi<<<grid, block>>>(d_A, d_V, d_pair, n, tolerance);

	/*write matrix back to host*/
	cudaMemcpy(A, d_A, 1024*1024*sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(V, d_V, 1024*1024*sizeof(double), cudaMemcpyDeviceToHost);

	/*check result*/
	double* ans = (double *) malloc(n*sizeof(double));
	double norm = 0;

	for (i = 0; i<n; i++){
		for (j = 0; j<n; j++){
			if (i==j)
			{
				*(ans+i) = *(A+i*n+j);
				norm += (*(ans+i))*(*(ans+i));
				printf("%lf\n", *(A+i*n+j));
			}
		}
	}
	norm = sqrt(norm);
	printf("Norm is %lf\n", norm);
}


__global__ void jacobi(double *dev_A, double *dev_V, int *dev_pair, int size, double eps)
{
	short p, q, n, i, j, cont = 1;
	double c, s;
	n = size;

	while (cont != 0)
	{
		for ( i = 0; i < n/2; i++)
		{
			p = *(pair + i);
			q = *(pair + i + n/2);
			if (p > q)
			{
				int temp = p;
				p = q;
				q = temp;
			}

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
	        *(dev_A + p * n + i) = Api;
	        *(dev_A + q * n + i) = Aqi;
	    }

	    
	    for (i = 0; i < n; i++)
	    { 
	        double Aip = (*(dev_A + i * n + p))*c + (*(dev_A + i * n + q))*(-s);
	        double Aiq = (*(dev_A + i * n + p))*s + (*(dev_A + i * n + q))*c;
	        *(dev_A + i * n + p) = Aip;
	        *(dev_A + i * n + q) = Aiq;
	    }
	     
	    /* V = V*J */
	    for (i = 0; i < n; i++)
	    { 
	        double Vpi = (*(dev_V + p * n + i))*c + (*(dev_V + q * n + i))*(-s);
	        double Vqi = (*(dev_V + p * n + i))*s + (*(dev_V + q * n + i))*c;
	        *(dev_V + p * n + i) = Vpi;
	        *(dev_V + q * n + i) = Vqi;
		}
	}

		/* tournament ordering */
		int *p;
		int row1end;

		row1end = *(dev_pair + n/2 - 1);

		for (p = dev_pair + n/2 - 1; p > dev_pair; p--)
			*p = *(p - 1);

		*(dev_pair + 1) = *(dev_pair + n/2);

		for (p = dev_pair + n/2; p < dev_pair + n; p++)
			*p = *(p + 1);

		*(dev_pair + n - 1) = row1end;

		/* check entry */
		cont = 0;
		for (i = 0; i < n; i++)
		{
			for (j = 0; j < n; j++)
			{
				if (i != j)
				{
					if ((*(dev_A + i*n + j)) > eps)
					{
						cont = 1;
						break;
					}
				}
			}
		}
	}
}






