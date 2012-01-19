#include "LLRL.h"
/*************************************************************************
Workspace used: 9
**************************************************************************	
	Compute the manipulator articulated body inertia matrix.

 		M = abi(DH,DYN,q,M) 

	DH=DH parameter matrix
	DYN=dynamic parameters matrix
	q= is an n element vector of joint state.

	for an n-axis manipulator returns the nxn symmetric articulated body 
	inertia matrix which relates joint torque to joint acceleration.
 
  	Copyright, Andr� Jaramillo Botero, 1999. 
*************************************************************************/
/*09-02-2006, JCA: added support for flying base inverse dynamics (ff_ne)*/
/*22/02/2006, JCA, AM: Solved problem with memory in ff_ne calls*/
/*22/02/2006 JCA  Add support for Spatial Base Position: X*/
/*23/03/2006 JCA  CHanged interface*/ 
/*30/03/2006 JCA added suport for fext, grav, energy changes in ff_ne*/
/*07/04/2006 JCA added sopport for new ff_ne  interface. and changed all null matrix for only one...*/
/*02/05/2006 JDC added support for floating-serpentine robots to obtain mass operator*/
/*07/05/2006, JDC Changed call for ff_ne function, passed MNULL for DYN_base */
/*12/05/2006* JDC, Modification for Serpentine conditions */
/*01/06/2006 JCA: removed Fbase*/

MAT *abi(int Robot_type, MAT *DH,MAT *DYN,VEC *q,MAT *M,MAT *X) 
{
	static	MAT	*Q=MNULL,*Qdd=MNULL,*XT=MNULL;
	static	VEC	*grav=VNULL,*fext=VNULL;
	int	i;

	/* validate inputs (sizes and initializations) */
	if (DH == MNULL || DYN == NULL || q == NULL) 
		error(E_NULL,"inertia (intputs)\n");
	if (M == MNULL)
		M=m_resize(M,q->dim,q->dim);

	/* allocate temporals */
	m_resize_vars(q->dim,q->dim,&Q,&Qdd,NULL);
	mem_stat_reg_vars(0,TYPE_MAT,&Q,&Qdd,NULL);
	
	grav=v_resize(grav,6);	/* zero gravity vector */
	fext=v_resize(fext,6);	/* zero external forces */
	mem_stat_reg_vars(0,TYPE_VEC,&grav,&fext,NULL);
	v_zero(grav);
	v_zero(fext);
	
	for (i=0;i<(int)DH->m;i++) 
		set_row(Q,i,q);

	Q=m_resize(Q,q->dim,q->dim*3);
	MEM_STAT_REG(Q,TYPE_MAT);
	m_ident(Qdd);
	m_move(Qdd,0,0,Qdd->m,Qdd->n,Q,0,2*q->dim);
	

	//For Flying base robots or Serpentine
	if(Robot_type==2 || Robot_type==4)
	{
		XT=m_resize(XT,6,DH->m);
		MEM_STAT_REG(XT,TYPE_MAT);
		for (i=0;i<(int)DH->m;i++) 
			m_move(X,0,0,6,1,XT,0,i);
	}	
		mem_stat_mark(9);
		ff_ne(Robot_type,DH, DYN, MNULL, Q, grav, fext, M,XT,MNULL,MNULL,MNULL);
		mem_stat_free(9);	
		
return M;
}
