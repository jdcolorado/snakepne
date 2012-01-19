#include "LLRL.h"
/*************************************************************************	
Compute the manipulator cartesian articulated body inertia matrix.

	M=cabi(DH,DYN,q,M) 

	DH=DH parameter matrix
	DYN=dynamic parameters matrix
	q= is an n element vector of joint state.
	
	for an n-axis manipulator returns the nxn symmetric inertia matrix 
	which relates Cartesian force/torque to Cartesian acceleration.
	  
	Copyright, Andr� Jaramillo Botero, 1999. 
***************************************************************************/
/*JCA added support for new abi, but for fixed base*/
/*01/06/2006 JCAadded support for new abi, but for fixed base*/
MAT *cabi(int Robot_type,MAT *DH,MAT *DYN,VEC *q,MAT *M) 
{
	static	MAT	*J=MNULL,*Ji=MNULL,*m1=MNULL;
	
	/* validate inputs (sizes and initialization) */
	if (DH == (MAT *)MNULL || DYN == (MAT *)MNULL || q == (VEC *)VNULL) 
		error(E_NULL,"cinertia (intputs)\n");
	if (DH->m != DYN->m && DH->m != q->dim)
		error(E_SIZES,"cinertia (inputs)\n");
	if (M == (MAT *)MNULL)	
		M=m_resize(M,q->dim,q->dim);
	else
		if (M->m != q->dim || M->n != q->dim)
			M=m_resize(M,q->dim,q->dim);

	J=m_resize(J,6,q->dim);
	Ji=m_resize(Ji,q->dim,6);
	m1=m_resize(m1,q->dim,q->dim);
	mem_stat_reg_vars(0,TYPE_MAT,&J,&Ji,&m1,NULL);

	jacob_base(DH,q,J);
	pinv(J,Ji);		/* using pseudo-inverse instead of inverse */
	abi(Robot_type,DH,DYN,q,M,MNULL);
	m_transp(Ji,J);
	m_mlt(J,M,m1);
	m_mlt(m1,Ji,M);
	
	/*	M = Ji' * M * Ji */
	
	return M;
}

