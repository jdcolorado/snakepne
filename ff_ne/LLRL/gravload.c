#include "LLRL.h"
/*************************************************************************	
	Compute the gravity load on manipulator joints.
 
 	Tg=gravload(DH,DYN,Q,grav,Tg) 
	
	DH=DH parameter matrix
	DYN=dynamic parameters matrix
	Q=is an npxn element matrix of joint poisition states.
	Tg=matrix of gravity torques

	allows an arbitrary gravity vector to override the default of [0; 0; 9.81]	
	
	For an n-axis manipulator returns the n element friction torque matrix 
	each row being the corresponding joint torques (at each trajectory point) 
	at the specified pose and velocity required to compensate for gravity. 

	Copyright, Andrés Jaramillo Botero, 1999. 
***************************************************************************/
MAT *gravload(MAT *DH,MAT *DYN,MAT *Q,VEC *grav, MAT *Tg) 
{
	static	MAT *Q_Qd;
	static	VEC *fext;
	
	/* validate inputs (sizes and initial values) */
	if (DH == (MAT *)MNULL || DYN == (MAT *)MNULL || Q == (MAT *)MNULL) 
		error(E_NULL,"gravload (parameter inputs)\n");
	if (Q->n != DH->m)
		error(E_SIZES,"gravload (state inputs)\n");
	if (grav == (VEC *)VNULL || grav->dim != 3)
		error(E_SIZES,"gravload (gravity vector size)\n");
	if (Tg == (MAT *)MNULL) 
		Tg=m_resize(Tg,Q->m,DH->m);
	else
		if (Tg->m != Q->m || Tg->n != DH->m)
			Tg=m_resize(Tg,Q->m,DH->m);
	
	fext=v_resize(fext,6);	/* zero external forces */
	MEM_STAT_REG(fext,TYPE_VEC);
	
	Q_Qd=m_resize(Q_Qd,Q->m,Q->n*3);
	MEM_STAT_REG(Q_Qd,TYPE_MAT);
	m_move(Q,0,0,Q->m,Q->n,Q_Qd,0,0);
	ne(DH, DYN, Q_Qd, grav, fext, Tg);
	
	return Tg;
}

