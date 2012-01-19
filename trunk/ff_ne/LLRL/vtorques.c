#include "LLRL.h"
/*************************************************************************	
Compute the manipulator velocity dependent torques matrix.

 	C=vtorque(DH,DYN,Q,Qd,C) 
	
	DH=DH parameter matrix
	DYN=dynamic parameters matrix
	Q=is an npxn element matrix of joint position states.
	Qd=is an npxn element matrix of joint velocity states.
	C=matrix of velocity dependent torques
	  
	For an n-axis manipulator returns the n element velocity torques matrix 
	each row being the corresponding joint torques (at each trajectory point) 
	at the specified pose and velocity. 
		
	Copyright, Andrés Jaramillo Botero, 1999. 
***************************************************************************/
MAT *vtorques(MAT *DH,MAT *DYN,MAT *Q,MAT *Qd, MAT *M) 
{
	static	MAT	*Q_Qd;
	static	VEC	*grav,*fext;
	
	/* validate inputs (sizes and initial values) */
	if (DH == (MAT *)MNULL || DYN == (MAT *)MNULL 
		|| Q == (MAT *)MNULL || Qd == (MAT *)MNULL) 
		error(E_NULL,"(intputs) vtorque\n");
	if (Q->n != Qd->n)
		error(E_SIZES,"vtorque (state inputs)\n");
	if (M == (MAT *)MNULL)	
		M=m_resize(M,Q->m,DH->m);
	else
		if (M->m != Q->m || M->n != DH->m)
			M=m_resize(M,Q->m,DH->m);
		
		/* allocate temporals */
		grav=v_resize(grav,3);	/* zero */
		fext=v_resize(fext,6);	/* zero */
		mem_stat_reg_vars(0,TYPE_VEC,&grav,&fext,NULL);
		
		Q_Qd=m_resize(Q_Qd,Q->m,Q->n*3);
		MEM_STAT_REG(Q_Qd,TYPE_MAT);
		m_move(Q,0,0,Q->m,Q->n,Q_Qd,0,0);
		m_move(Qd,0,0,Q->m,Q->n,Q_Qd,0,Q->n);	
		ne(DH, DYN, Q_Qd, grav, fext, M);
		
		return M;
}
