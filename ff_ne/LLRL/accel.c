#include "LLRL.h"
/******************************************************************************
Workspace used: 8, 6, 9
**************************************************************************
	Compute manipulator forward dynamics (i.e. find acceleration from force)

	Qdd = accel(DH,DYN,Q_Qd,Tor,Qdd)

	DH=DH parameter matrix
	DYN=dynamic parameters matrix
	Q_Qd=is an npx2n element matrix of joint pos/vel states.
	Tor=npxn initial torque matrix
	Qdd=joint accelerations matrix

	Returns a matrix of joint accelerations that result from applying the 
	actuator TORQUE to the manipulator in state Q and QD.  It returns one 
	row entry per trajectory point (np).

	Uses Walker and Orin's method 1 to compute the forward dynamics.  This
	form is useful for simulation of manipulator dynamics, in conjunction with
	a numerical integration function (ode45.c).

	Copyright, Andr� Jaramillo Botero, 1999. 
******************************************************************************/
/*09-02-2006, JCA: added support for flying base*/
/*22/02/2006, JCA, AM: Solved problem with memory in abi and ff_ne calls*/
/*22/02/2006, JCA, JDC:Add spatial base position: X=[roll pitch yaw x y z]'*/
/*23/03/2006 JCA  CHanged interface kinetic,potential,E_Conserv*/ 
/*30/03/2006 JCA added support for fext and grav pass..*/
/*07/04/2006 JCA added sopport for new ff_ne  interface. and changed all null matrix for only one...*/
/*16/04/2006, JCA added support for friction*/
/*07/05/2006, JDC Changed call for ff_ne function, passed MNULL for DYN_base
07/05/2006, JDC added call to Kinetic and Potential energy computation functions */
/*1/06/2006: JCA removed Fbase.*/
MAT *accel(int Robot_type, MAT *DH, MAT *DYN,MAT *Q_Qd,MAT *Tor,MAT *Qdd,MAT *X, MAT *Vb, MAT *dVb,VEC *grav, VEC *fext)
{
	static	MAT	*M=MNULL,*Mi=MNULL,*Q=MNULL,*Q_Qd_Qdd=MNULL,*Qs=MNULL,*Tau=MNULL;
	static	VEC	*q=VNULL,*tau=VNULL,*tor=VNULL, *qdd=VNULL;
	int		i;
	float K,U;
	
	/* validate inputs (sizes and initial values) */
	if ((DH == MNULL) || (DYN == MNULL) || (Q_Qd == MNULL))
		error(E_NULL,"(inputs) accel\n");
	if (Tor == MNULL) 
		Tor=m_resize(Tor,(int)Q_Qd->m,(int)DH->m);	
	if (Qdd == MNULL || Qdd->n != DH->m) 
		Qdd=m_resize(Qdd,(int)Q_Qd->m,(int)DH->m);
	if(grav==VNULL)
		error(E_NULL,"NULL gravity");
	else
		if (grav->dim != 6)
			error(E_SIZES,"Wrong gravity vector");
	if(fext==VNULL)
		error(E_NULL,"NULL fext");
	else
		if (fext->dim != 6)
			error(E_SIZES,"Wrong fext vector");
	
	/* allocate temporals */
	Q=m_resize(Q,(int)Q_Qd->m,(int)DH->m);
	m_resize_vars((int)DH->m,(int)DH->m,&M,&Mi,NULL);
	Tau=m_resize(Tau,1,(int)DH->m);
	Q_Qd_Qdd=m_resize(Q_Qd_Qdd,(int)Q_Qd->m,(int)DH->m*3);
	Qs=m_resize(Qs,1,(int)DH->m*3);
	mem_stat_reg_vars(0,TYPE_MAT,&Q_Qd_Qdd,&Qs,&Q,&M,&Mi,&Tau,NULL);
	
	v_resize_vars((int)DH->m,&q,&tau,&tor,&qdd,NULL);
	mem_stat_reg_vars(0,TYPE_VEC,&q,&tau,&tor,&qdd,NULL);
	

	(void)m_move(Q_Qd,0,0,(int)Q_Qd->m,(int)DH->m,Q,0,0);
	(void)m_move(Q_Qd,0,0,(int)Q_Qd->m,(int)DH->m*2,Q_Qd_Qdd,0,0);
	
/***************************************************************************
Speed dependant Torques
***************************************************************************/
	/*printf("\nantes de speed dependant torques?");
	mem_info();*/
	
	mem_stat_mark(8);
	Tau=ff_ne(Robot_type,DH,DYN,MNULL,Q_Qd_Qdd,grav,fext,Tau,X,Vb,dVb,MNULL);
	mem_stat_free(8);
	/*printf("\ndespues de speed dependant torques?");
	mem_info();*/
	
	for (i=0;i<(int)Q_Qd->m;i++)
	{
		get_row(Q,(u_int)i,q);
/***************************************************************************
Robot mass
***************************************************************************/
		mem_stat_mark(6);
		M=abi(Robot_type,DH,DYN,q,M,X);
		mem_stat_free(6);
		
		//Compute Kinetic and Potential Energy of all the system
		/*mem_stat_mark(6);
		K=Kinetic(Q_Qd_Qdd,M);
		U=Potential(DH,DYN,Q_Qd_Qdd);
		mem_stat_free(6);*/

		get_row(Tau,(u_int)i,tau);
		get_row(Tor,(u_int)i,tor);
		
		mem_stat_mark(6);
		Mi=pinv( M, Mi);
		mem_stat_free(6);
		v_sub(tor,tau,tor);

		(void)mv_mlt(Mi,tor,qdd);
		(void)set_row(Qdd,i,qdd);
	}
	return Qdd;
}
