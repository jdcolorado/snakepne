#include "LLRL.h"

/* divide todos los elementos de dos matrices entre si */
MAT *div_m(MAT *mat,MAT *mat2,MAT *out)
{
//	static	MAT *out;
	int i,m;
	
	/* validate inputs (sizes and initialization) */
	if (mat == (MAT *)MNULL || mat2 == (MAT *)MNULL) 
		error(E_NULL,"div_m (intputs)\n");
	if (out == (MAT *)MNULL || out->m != mat->m || out->n != 1)
		out=m_resize(out,mat->m,1);

	out=m_resize(out,mat->m,1);
//	MEM_STAT_REG(out,TYPE_MAT);

	m_zero(out);
	m=mat->m;
	for (i=0 ; i < m ; i++)
		out->me[i][0]=(mat->me[i][0])/(mat2->me[i][0]);
    
	return out;
}
