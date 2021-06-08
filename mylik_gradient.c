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
}

void mylik_gradient_second(int *n, double *y, double *xall, double *par, double *length_par, double *result)
{
	int j, l, ini, index, row;
	double temp, temp1;
	int arr[15]; 
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
				arr[row] = (xall[row*n[0]+j]-xall[row*n[0]+l]);	
			}
			
			result[0] += -(y[j]-y[l])*(y[j]-y[l])*arr[0]*arr[0]*(temp/((1+temp)*(1+temp)));
			result[1] += -(y[j]-y[l])*(y[j]-y[l])*arr[0]*arr[1]*(temp/((1+temp)*(1+temp)));
			result[2] += -(y[j]-y[l])*(y[j]-y[l])*arr[0]*arr[2]*(temp/((1+temp)*(1+temp)));
			result[3] += -(y[j]-y[l])*(y[j]-y[l])*arr[0]*arr[3]*(temp/((1+temp)*(1+temp)));

			result[4] += -(y[j]-y[l])*(y[j]-y[l])*arr[1]*arr[0]*(temp/((1+temp)*(1+temp)));
			result[5] += -(y[j]-y[l])*(y[j]-y[l])*arr[1]*arr[1]*(temp/((1+temp)*(1+temp)));
			result[6] += -(y[j]-y[l])*(y[j]-y[l])*arr[1]*arr[2]*(temp/((1+temp)*(1+temp)));
			result[7] += -(y[j]-y[l])*(y[j]-y[l])*arr[1]*arr[3]*(temp/((1+temp)*(1+temp)));

			result[8] += -(y[j]-y[l])*(y[j]-y[l])*arr[2]*arr[0]*(temp/((1+temp)*(1+temp)));
			result[9] += -(y[j]-y[l])*(y[j]-y[l])*arr[2]*arr[1]*(temp/((1+temp)*(1+temp)));
			result[10] += -(y[j]-y[l])*(y[j]-y[l])*arr[2]*arr[2]*(temp/((1+temp)*(1+temp)));
			result[11] += -(y[j]-y[l])*(y[j]-y[l])*arr[2]*arr[3]*(temp/((1+temp)*(1+temp)));

			result[12] += -(y[j]-y[l])*(y[j]-y[l])*arr[3]*arr[0]*(temp/((1+temp)*(1+temp)));
			result[13] += -(y[j]-y[l])*(y[j]-y[l])*arr[3]*arr[1]*(temp/((1+temp)*(1+temp)));
			result[14] += -(y[j]-y[l])*(y[j]-y[l])*arr[3]*arr[2]*(temp/((1+temp)*(1+temp)));
			result[15] += -(y[j]-y[l])*(y[j]-y[l])*arr[3]*arr[3]*(temp/((1+temp)*(1+temp)));
		}
	}
}





