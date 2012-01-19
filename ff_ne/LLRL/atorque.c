#include "LLRL.h"
/*************************************************************************	
	Compute the manipulator inertia torque.  This corresponds to the acceleration
	dependent torques, without considering the velocity dependent torques.
	For an n-axis manipulator:
	returns the n element inertia torque vector (or a matrix each row being 
	the corresponding joint torques) at the specified pose and acceleration, 
	that is, ABI(Q)*QDD

	T = atorque(DH,DYN,Q,Qdd,T)
	
	DH=DH parameter matrix
	DYN=dynamic parameters matrix
	Q=is an npxn element matrix of joint position states.
	Qdd=is an npxn element matrix of joint acceleration states.

	for an n-axis manipulator returns the nxn symmetric inertia matrix 
	which relates Cartesian force/torque to Cartesian acceleration.

	Copyright, Andrés Jaramillo Botero, 1999. 

  07/02/2002: added validation of equal size for Q and Qdd entry arguments
***************************************************************************/
MAT *atorque(MAT *DH,MAT *DYN,MAT *Q,MAT *Qdd, MAT *M) 
{
	static	MAT	*Q_Qdd;
	static	VEC	*grav,*fext;
	
	/* validate inputs (sizes and initializations) */
	if (DH == (MAT *)MNULL || DYN == (MAT *)MNULL 
		|| Q == (MAT *)MNULL || Qdd == (MAT *)MNULL) 
		error(E_NULL,"atorque (inputs)\n");
	if (DH->m != DYN->m)
		error(E_SIZES,"atorque (parameter inputs)\n");
	if ((Q->m != Qdd->m) || (Q->n != Qdd->n))
		error(E_SIZES,"atorque (state inputs)\n");
	if (M == (MAT *)MNULL)	
		M=m_resize(M,Q->m,DH->m);
	else
		if (M->m != Q->m || M->n != DH->m)
			M=m_resize(M,Q->m,DH->m);
	
	/* allocate temporals */
	grav=v_resize(grav,3);	/* zero */
	fext=v_resize(fext,6);	/* zero */
	mem_stat_reg_vars(0,TYPE_VEC,&grav,&fext,NULL);
	
	Q_Qdd=m_resize(Q_Qdd,Q->m,Q->n*3);
	MEM_STAT_REG(Q_Qdd,TYPE_MAT);
	m_move(Q,0,0,Q->m,Q->n,Q_Qdd,0,0);
	m_move(Qdd,0,0,Q->m,Q->n,Q_Qdd,0,2*Q->n);
	ne(DH, DYN, Q_Qdd, grav, fext, M);	/* zero velocities */
	
	return M;
}

