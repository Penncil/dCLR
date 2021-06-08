// This c code is for Biometrics paper (proportional likelihood ratio model)
// with more than one covariate x (xall)
// independent sparse data

#include <R.h>


void mylik_all(int *n, double *y, double *xall, double *par, double *length_par, double *result)
{
	int j, l, ini;
	result[0] = 0;
	double temp;
	for (j = 0; j < n[0] - 1; j++)
	{
		for (l = j+1; l<n[0]; l++)
		{
			temp = 0;
			for (ini=0; ini<length_par[0]; ini++){
				temp += (xall[ini*n[0]+j]-xall[ini*n[0]+l]) * par[ini];
			}
			
			result[0] += -log(1+exp(-(y[j] - y[l]) * temp));
		}
	}
}





