#include "LLRL.h"
//saca una matriz maxima entre los elementos de una matriz y un  double
MAT* maxm(MAT *a,double num,MAT *b)
{ 
	int i,m,n;
	//static	MAT *b=MNULL;
	
	/* validate inputs (sizes and initialization) */
	if (a == (MAT *)MNULL) 
		error(E_NULL,"maxm (inputs)\n");

	m=a->m;
	n=a->n;

	b=m_resize(b,m,1);
	//MEM_STAT_REG(b,TYPE_MAT);

	m_zero(b);
	for(i=0;i<m;i++)
		if (num>(a->me[i][0]))
			b->me[i][0]=num;
		else
			b->me[i][0]=a->me[i][0];
	return b;
}
