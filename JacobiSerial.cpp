#include<math.h>
#include <iostream>
#include "mpi.h"
using namespace std;
const float EPS = 0.0001;
const int max_iter = 1000;
int n=1000;
float a[10000][10000];
double x[1000];


void init(int &n, float a[10000][10000],int &i,int &j) {
	/*for (i = 1; i <= n; i++)
	{
		cout << "inserire i valori" << endl;
		for (j = 1; j <= n + 1; j++)
			cin >> a[i][j];
	
	
	   x[i] = a[i][n-1] / a[i][i];
	}*/
	for (int i = 1; i <= n; ++i) {
		for (int j = 1; j <=n+1; ++j) {
			if (i == j)
				a[i][j] =n;
			else if(j==n&&i==n-1)
				a[i][j] = n;
			else a[i][j] = 1.0;

		}
		x[i] = a[i][n-1] / a[i][i];
	}


}
int main(int argc, char* argv[])
{
	float big, temp, relerror, sum;
	int i, j, itr;
	double start_time, end_time;

	init(n, a, i, j);
	MPI_Init(&argc, &argv);

	start_time = MPI_Wtime();
	for (i = 1; i <= n; i++)
		x[i] = 0;
	for (itr = 1; itr <= max_iter; itr++)
	{
		big = 0;
		for (i = 1; i <= n; i++)
		{
			sum = 0;
			for (j = 1; j <= n; j++)
			{
				if (i != j)
					sum = sum + a[i][j] * x[j];
			}
			temp = (a[i][n + 1] - sum) / a[i][i];
			relerror = fabs((x[i] - temp) / temp);
			if (relerror>big)
				big = relerror;
			x[i] = temp;
		}
		if (big <= EPS)
		{
			end_time = MPI_Wtime();
			cout << "Iterazioni: " << itr <<endl;
			cout << "Tempo seriale: " << end_time - start_time << endl;
			MPI_Finalize();
			return 0;
			
		}

	}
	cout << "Iterazioni: " << max_iter<< endl;
	return 0;
}
