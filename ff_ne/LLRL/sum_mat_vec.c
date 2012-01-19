#include "LLRL.h"
//suma una matriz con un vector fila
MAT *sum_mat_vec(MAT *mat,VEC *vt,MAT *out)
{
	int i,m,n;
//	static MAT *out;
	
	/* validate inputs (sizes and initialization) */
	if (mat == (MAT *)MNULL || vt == (VEC *)VNULL) 
		error(E_NULL,"sum_mat_vec (intputs)\n");
	if (out == (MAT *)MNULL || out->m != mat->m || out->n != 1)
		out=m_resize(out,mat->m,1);

	m=mat->m;
	n=mat->n;
	
	out=m_resize(out,m,1);
//	MEM_STAT_REG(out,TYPE_MAT);

	m_zero(out);
	for (i=0 ; i < m ; i++)
		out->me[i][0]=mat->me[i][0]+vt->ve[i];
		
	return out;
}
