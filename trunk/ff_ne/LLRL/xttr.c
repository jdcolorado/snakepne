#include"LLRL.h"

/*************************************************************************	
	Extract a homogeneous matrix from a matrix of homogeneos contiguous 
	matrices. [T0 T1 ... Tn]

	Tout = xttr(TR, i, Tout) 
		
	TR=	matrix of homomgeneous transforms [T0 T1 ... Tn]
	i=	desired number of homogeneous matrix in TR.
	Tout=Extracted homogeneous transform.	
		  
	Copyright, Andrés Jaramillo Botero, 1999. 
***************************************************************************/
MAT	*xttr(MAT *TR,int n,MAT *Tout)
{
	/* validate inputs (sizes and initializations) */
	if (TR == (MAT *)MNULL) 
		error(E_NULL,"Null input matrix in xttr\n");
	if (Tout == (MAT *)MNULL || (Tout->m != 4 && Tout ->n != 4)) 
		Tout=m_resize(Tout,4,4);
	if (n<1 || n*4>(int)TR->n) 
	{
		printf("Matrix index must be > 1 and < %d\n",(int)TR->n/4);
		error(E_RANGE,"xttr\n");
	}
	
	m_move(TR,0,4*(n-1),4,4,Tout,0,0);
	
	return	Tout;
}
