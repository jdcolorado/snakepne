#include "LLRL.h"
//devuelve el indice dentro de la matriz
double indice(MAT *mat,int indice)
{
	int i=0,cont=0,m,n,j=0;
	
	/* validate inputs (sizes and initialization) */
	if (mat == (MAT *)MNULL) 
		error(E_NULL,"indice (inputs)\n");

	m=mat->m;
	n=mat->n;

	for(i=0;i<m;i++)
	{
		for(j=0;j<m;j++)
		{
			cont++;
			if(cont==indice) break;
		}
		break;
	}
	return mat->me[i][j];
}
