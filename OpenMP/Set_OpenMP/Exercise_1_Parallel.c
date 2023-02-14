#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <math.h>

#define threads_num 16

int main()
{
	int i, l, j, r;
	int val_1 = 0, val_2 = 0;
	double k = 0.8;
	int N = 500;
	int M = 500;
	double x_start = 0, y_start = 0, t_start = 0;
	double x_end = 2 * M_PI, y_end = 2 * M_PI, t_end = 2 * M_PI;
	double **uold, **unew;
	double dx = (x_end - x_start) / (N - 1);
	double dy = (y_end - y_start) / (M - 1);
	double dt = M_PI / 3600.0;
	double *x, *y, *t, *u;
	double ftime1 = 0, ftime2 = 0;
	ftime1 = omp_get_wtime();

	/*Memory allocation*/
	uold = (double **)malloc(N * sizeof(double *));
	unew = (double **)malloc(N * sizeof(double *));
	x = malloc(N * sizeof(double));
	y = malloc(M * sizeof(double));

	for (i = 0; i < N; i++)
	{
		uold[i] = (double *)malloc(M * sizeof(double));
		unew[i] = (double *)malloc(M * sizeof(double));
	}

#pragma omp parallel num_threads(threads_num)
{
	#pragma omp for private(i)
	for (i = 0; i < N; i++)
	{
		x[i] = x_start + i * dx;
	}

	#pragma omp for private(j)
	for (j = 0; j < M; j++)
	{
		y[j] = y_start + j * dy;
	}

	#pragma omp for collapse(2)
	for (l = 0; l < N; l++)
	{
		for (r = 0; r < M; r++)
		{
			uold[l][r] = sin(x[l]) * sin(y[r] / 2);
			unew[l][r] = sin(x[l]) * sin(y[r] / 2);
		}
	}

	/*boundary conditions*/
	#pragma omp parallel for collapse(2) private(i, j)
		for (i = 1; i < N - 1; i++)
		{
			for (j = 1; j < M - 1; j++)
			{

				unew[0][j] = 0;
				unew[N - 1][j] = 0;
				unew[i][0] = 0;
				unew[i][N - 1] = 0;
				uold[0][j] = 0;
				uold[N - 1][j] = 0;
				uold[i][0] = 0;
				uold[i][N - 1] = 0;
			}
		}


		for (int val = 0 + dt; val <= t_end; val += dt)
		{
			#pragma omp parallel for collapse(2) private(i, j)
			for (i = 1; i < N - 1; i++)
			{
				for (j = 1; j < M - 1; j++)
				{
					unew[i][j] = uold[i][j] + k * dt * ((uold[i + 1][j] - 2 * uold[i][j] + uold[i - 1][j]) / (pow(dx, 2)) + (uold[i][j + 1] - 2 * uold[i][j] + uold[i][j - 1]) / (pow(dy, 2)));
				}
			}

			#pragma omp parallel for collapse(2) private(val_1, val_2)
			for (val_1 = 0; val_1 < N; val_1++)
			{
				for (val_2 = 0; val_2 < M; val_2++)
				{
					uold[val_1][val_2] = unew[val_1][val_2];
				}
			}
		}

	}

	ftime2 = omp_get_wtime();
	printf("wall clock time= %.8f\n", ftime2 - ftime1);
	return 0;
}
