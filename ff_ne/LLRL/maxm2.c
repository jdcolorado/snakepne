#include "LLRL.h"
//saca los elemntos maximos de cada posicion de dos matrices
MAT* maxm2(MAT *a,MAT *c,MAT *b)
{ 
	int i,m,n;
	//static	MAT *b=MNULL;
	
	/* validate inputs (sizes and initialization) */
	if (a == (MAT *)MNULL || c == (MAT *)MNULL) 
		error(E_NULL,"maxm2 (inputs)\n");

	m=a->m;
	n=a->n;
	b=m_resize(b,m,1);
	//MEM_STAT_REG(b,TYPE_MAT);

	m_zero(b);
	for(i=0;i<m;i++)
		if (c->me[i][0]>(a->me[i][0]))
			b->me[i][0]=c->me[i][0];
		else
			b->me[i][0]=a->me[i][0];
	return b;
}
