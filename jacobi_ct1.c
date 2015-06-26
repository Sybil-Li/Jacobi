/**
 * File: chess tournament ordering jacobi
 * This is a working version of the jacobi method to find eigenvalue of a symmatrix matrix
 * Each iteration applies jacobi method to a pair of rows according to the chess tournament ordering
 * Meanwhile, V matrix is also calculated.
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

void jacobi(double* A, double* V, double c, double s, int p, int q, int n);
void calcJ(double* A, double* c, double* s, int p, int q, int n);
int* CT (int* pair, int n);

int main ( int argc, char *argv[] )
{

	//command line argument takes one parameter n: size of matrix
	//use input redirection to read matrix from file
	
	if ( argc != 2 ) /* argc should be 2 for correct execution */
    {
        printf( "usage: %s filename", argv[0] );
		exit(1);
    }


	double eps = 0.0000000001; //tolerance
	int n = atoi(argv[1]), row, col, p, q;
	double max_off = 0;
	double c, s; //transformation factors calculated by jacobi
	double *cp, *sp;
	cp = &c;
	sp = &s;

	/* Allocate space for original matrix */
	double* A = (double*)malloc(1024*1024*sizeof(double));
	/* Allocate space for V matrix */
	double* V = (double*)malloc(1024*1024*sizeof(double));

	/* enter a valid matrix A*/
	int i, j;
	for (i=0; i<n; i++)
	{
		for (j=0; j<n-1; j++)
			scanf("%lf,", (A+i*n+j));
		scanf("%lf\n", (A+i*n+n-1));
	}

	//error checking code
	/*for (i = 0; i<n; i++)
		for (j=0; j<n; j++)
			printf("%lf\n", *(A+i*n+j));*/

	/*initializing pair matrix*/
	int* pair = malloc(n*sizeof(int));
	// Pair matrix is just a matrix of index 
	// corresponding to rows of the matrix
	// eg. for n = 6
	// pair = [1,2,3,4,5,6]
	i = 0;
	for (i = 0; i < n; i++)
		*(pair + i) = i;

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

	/*find largest off diagonal*/
	for (row = 0; row < n; row++) {
		for (col = 0; col < n; col++) {
			if (row != col)
			{
				if (fabsf(*(A + row * n + col)) > max_off)
				{
					max_off = fabsf(*(A + row * n + col));
					if (row < col)
					{
						p = row;
						q = col;
					}
					else 
					{
						p = col;
						q = row;
					}
				}
			}
		}
	}
	/*printf("max_off = %f\n", max_off);*/

	//stoping criterion: all off diagonal element < tolerance
	while (fabs(max_off) > eps) 
	{

		int i;
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
			calcJ(A, cp, sp, p, q, n);
			jacobi(A, V, c, s, p, q, n);
		}

		/* get new pairing from chess tournament ordering */
		pair = CT(pair, n);
		

		double* pt = A;
		/*for (; pt < A + n*n; pt++)
			printf("%f\n", *pt);
		printf("\n");*/

		max_off = 0;
		for (row = 0; row < n; row++) {
			for (col = 0; col < n; col++) {
				if (row != col)
				{
					if (fabsf(*(A + row * n + col)) > max_off)
					{
						max_off = fabsf(*(A + row * n + col));
						if (row < col)
						{
							p = row;
							q = col;
						}
						else 
						{
							p = col;
							q = row;
						}
					}
				}
			}
		}
		
		//printf("max_off = %f \n", max_off);
	}

	/*check result*/
	double* ans = malloc(n*sizeof(double));
	double norm = 0;

	for (i = 0; i<n; i++){
		for (j = 0; j<n; j++){
			if (i==j)
			{
				*(ans+i) = *(A+i*n+j);
				//norm += (*(ans+i))*(*(ans+i));
				printf("%lf\n", *(A+i*n+j));
			}
		}
	}
	//norm = sqrt(norm);
	//printf("Norm is %lf\n", norm);
}

/* calculate the transformaiton matrix J */
void calcJ(double* A, double* cp, double* sp, int p, int q, int n)
{
	if (*(A + p * n + q) != 0)
	{
		double torque, t;
        torque = ( *(A + q * n + q) - *(A + p * n + p))/(2*(*(A + p * n + q)));
        if (torque >= 0)
            t = 1/(torque + sqrt(1+torque*torque));
        else
            t = -1/(-torque + sqrt(1+torque*torque));
        
        *cp = 1/sqrt(1+t*t);
        *sp = t*(*cp);
    }
    else
    {
        *cp = 1;
        *sp = 0;
	}
}


/* performs jacobi transformation */
void jacobi(double* A, double* V, double c, double s, int p, int q, int n)
{
	int i;
	
	/* A = transpose(J)*A*J */
    for (i = 0; i < n; i++)
    {
        double Api = (*(A + p * n + i))*c + (*(A + q * n + i))*(-s);
        double Aqi = (*(A + p * n + i))*s + (*(A + q * n + i))*c;
        *(A + p * n + i) = Api;
        *(A + q * n + i) = Aqi;
    }

    
    for (i = 0; i < n; i++)
    { 
        double Aip = (*(A + i * n + p))*c + (*(A + i * n + q))*(-s);
        double Aiq = (*(A + i * n + p))*s + (*(A + i * n + q))*c;
        *(A + i * n + p) = Aip;
        *(A + i * n + q) = Aiq;
    }
     
    /* V = V*J */
    for (i = 0; i < n; i++)
    { 
        double Vpi = (*(V + p * n + i))*c + (*(V + q * n + i))*(-s);
        double Vqi = (*(V + p * n + i))*s + (*(V + q * n + i))*c;
        *(V + p * n + i) = Vpi;
        *(V + q * n + i) = Vqi;
	}
}

/* performs chess tournament ordering (rotate onece) on the pair matrix */
int* CT (int* pair, int n)
{
	int *p;
	int temp, row1end;

	row1end = *(pair + n/2 - 1);

	for ( p = pair + n/2 - 1; p > pair; p--)
		*p = *(p - 1);

	*(pair + 1) = *(pair + n/2);

	for (p = pair + n/2; p < pair + n; p++)
		*p = *(p + 1);

	*(pair + n - 1) = row1end;

	return pair;
}
