// calcualte the varinace of the estimator

#include <R.h>

void cal_dev_score_all(int *n, double *y, double *xall, double *par, double *length_par, double *result)
{
	int j, l, ini, index, row, col, index2;
	double temp, temp1;
	for (index=0; index < length_par[0]; index++){
		result[index] = 0;
	}
				
	for (j = 0; j < n[0] - 1; j++)
	{
		for (l = j+1; l<n[0]; l++)
		{
			temp = 0;
			temp1 = 0;
			for (ini=0; ini<length_par[0]; ini++)
			{
				temp1 += (xall[ini*n[0]+j]-xall[ini*n[0]+l]) * par[ini];
			}
			temp += exp(-(y[j]-y[l])*temp1);
			for (row=0; row < length_par[0]; row++){
				result[row] += -((y[j]-y[l])*(y[j]-y[l]))*(xall[row*n[0]+j]-xall[row*n[0]+l])*(xall[row*n[0]+j]-xall[row*n[0]+l])*(temp/((1+temp)*(1+temp)));
			}

			index2 = length_par[0];
			for (row=0; row < length_par[0]; row++){
				for (col=0; col < length_par[0]; col++){
					if (row < col){
						result[index2] += -((y[j]-y[l])*(y[j]-y[l]))*(xall[row*n[0]+j]-xall[row*n[0]+l])*(xall[col*n[0]+j]-xall[col*n[0]+l])*(temp/((1+temp)*(1+temp)));
						index2 = index2 + 1;
					}	
				}	
			}
		}
	}
}


void cal_score_all(int *n, double *y, double *xall, double *par, double *length_par, double *result)
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
