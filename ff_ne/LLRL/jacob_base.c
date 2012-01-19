#include "LLRL.h"
/*************************************************************************	
	Compute manipulator Jacobian in world coordinates for the current 
	pose Q.
	
	J=jacob_base(DH,Q,J)
	
	The manipulator Jacobian matrix maps differential changes in joint space
	to differential Cartesian motion (world coord frame) of the end-effector.
	
	  dX = J dQ
	
	For an n-axis manipulator the Jacobian is a 6 x n matrix.

	Copyright, Andrés Jaramillo Botero, 1999. 
***************************************************************************/
/* 07-12-2002 AJB: changed sub_mat for move */
MAT *jacob_base(MAT *DH,VEC *q,MAT *J) 
{
	static	MAT *Tn,*Tn2,*Jn,*Q;
	int 	i,j;
	
	/* validate inputs */
	if (DH == (MAT *)MNULL || q == (VEC *)VNULL) 
		error(E_NULL, "jacob_base (inputs)\n");
	if (q->dim != DH->m) 
		error(E_SIZES, "jacob_base (inputs)\n");	
	if (J == (MAT *)MNULL || J->m != 6 || J->n != DH->m) 
		J=m_resize(J,6,DH->m);
	
	/* setup and allocate memory for temporals */
	Tn=m_resize(Tn,4,4);
	Tn2=m_resize(Tn2,3,3);
	Jn=m_resize(Jn,6,6);
	Q=m_resize(Q,1,DH->m);
	mem_stat_reg_vars(0,TYPE_MAT,&Tn,&Jn,&Q,&Tn2,NULL);
	
	/*	
	dX_tn = Jn dq where Jn = jacob_end(DH, q); 
	Jacobian from joint to wrist space
	*/

	m_zero(J);
	set_row(Q,0,q);
	fkine(DH, Q, Tn);	// end-effector transformation 
//	sub_mat(Tn,0,0,2,2,Tn2); 
	Tn2=m_move(Tn,0,0,3,3,Tn2,0,0);

	for (i=0;i<3;i++)
		for (j=0;j<3;j++)
			Jn->me[i][j]=Tn2->me[i][j];
		
	for (i=3;i<6;i++)
		for (j=3;j<6;j++)
			Jn->me[i][j]=Tn2->me[i-3][j-3];
			
	Tn=m_resize(Tn,6,(int)DH->m);
	jacob_end(DH,q,Tn);
	J=m_mlt(Jn,Tn,J);
			
	return J;
}


