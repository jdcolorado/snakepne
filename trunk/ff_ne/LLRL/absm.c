#include "LLRL.h"
//valor absoluto de una matriz
MAT *absm(MAT *mat)
{
	static	MAT *out;
	int i,j;
	int m,n;
	
	/* validate inputs (sizes and initialization) */
	if (mat == (MAT *)MNULL) 
		error(E_NULL,"absm (inputs)\n");

	m=mat->m;
	n=mat->n;

	out=m_resize(out,m,n);
	MEM_STAT_REG(out,TYPE_MAT);
	
	for(i=0;i<m;i++)
		for(j=0;j<n;j++)
		{
			if (mat->me[i][j]<0.0)
				out->me[i][j]=mat->me[i][j]*(-1.0);
			else
				out->me[i][j]=mat->me[i][j];
		}			
	return out;
}
