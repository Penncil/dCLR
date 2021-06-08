// This c code is for Biometrics paper (proportional likelihood ratio model)
// with more than one covariate x (xall)
// independent sparse data

#include <R.h>


void mylik_gradient(int *n, double *y, double *xall, double *par, double *length_par, double *result)
{
	int j, l, ini, index, row;
	double temp, temp1;
	for (index=0; index < length_par[0]; index++){
		result[index] = 0;
	}

	for(j = 0; j < n[0] - 1; j++)
	{
		for(l = j+1; l<n[0]; l++)
		{
			temp = 0;
			temp1 = 0;
			for (ini=0; ini<length_par[0]; ini++)
			{
				temp1 += (xall[ini*n[0]+j]-xall[ini*n[0]+l]) * par[ini];
			}
			temp += exp(-(y[j]-y[l])*temp1);
			
			for (row=0; row < length_par[0]; row++)
			{
					result[row] += (y[j]-y[l])*(xall[row*n[0]+j]-xall[row*n[0]+l])*(temp/(1+temp));	
			}
		}
	}
	return;
}

void mylik_gradient_second(int *n, double *y, double *xall, double *par, double *length_par, double *result)
{
	int j, l, ini, index, row, matrix_index, m, k;
	int matrix_size = (int)(length_par[0]*length_par[0]);
	double temp, temp1;
	int *arr = malloc(length_par[0]*sizeof(int));

	for (index=0; index < matrix_size; index++){
		result[index] = 0;
	}

	for(j = 0; j < n[0] - 1; j++)
	{
		for(l = j+1; l<n[0]; l++)
		{
			temp = 0;
			temp1 = 0;
			for (ini=0; ini<length_par[0]; ini++)
			{
				temp1 += (xall[ini*n[0]+j]-xall[ini*n[0]+l]) * par[ini];
			}
			temp += exp(-(y[j]-y[l])*temp1);
			

			for (row=0; row < length_par[0]; row++)
			{
				arr[row] = (xall[row*n[0]+j]-xall[row*n[0]+l]);	
			}
			
			matrix_index = 0;
			for (m = 0; m<length_par[0]; m++)
			{
				for (k = 0; k<length_par[0];k++)
				{
					result[matrix_index] += -(y[j]-y[l])*(y[j]-y[l])*arr[m]*arr[k]*(temp/((1+temp)*(1+temp)));
					matrix_index = matrix_index + 1;
				}
			}
		}
	}

	free(arr);
	return;
}





