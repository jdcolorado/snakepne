#include "LLRL.h"
/*************************************************************************	
	Compute friction torque (that required to overcome friction).
 
 	Tau = friction(DH,DYN,Qd,Tau) 
	
	DH=DH parameter matrix
	DYN=dynamic parameters matrix
	Qd=is an npxn element matrix of joint velocity states.
	Tau=matrix of friction torques

	For an n-axis manipulator returns the n element friction torque matrix 
	each row being the corresponding joint torques (at each trajectory point) 
	at the specified pose and velocity required to overcome friction. 

	Copyright, Andrés Jaramillo Botero, 1999. 
***************************************************************************/
MAT *friction(MAT *DH, MAT *DYN, MAT *Qd,MAT *Tau) 
{
    static	VEC *b,*tcm,*tcp,*v1,*v2;
	static	MAT	*m1,*m2,*m3;
	int		i,j;
	
	/* validate inputs (sizes and initial values) */
	if (DH == (MAT *)MNULL || DYN == (MAT *)MNULL || Qd == (MAT *)MNULL) 
		error(E_NULL,"friction (intputs)\n");
	if (DYN->n < 12) 
		error(E_NULL,"No friction model present in function friction");
	if (Qd->n != DYN->m) 
		error(E_NULL,"bad data in function friction");
	if (Tau == (MAT *)MNULL) 
		Tau=m_resize(Tau,Qd->m,Qd->n);
	
	if (Tau == (MAT *)MNULL)	
		Tau=m_resize(Tau,Qd->m,Qd->n);
	else
		if (Tau->m != Qd->m || Tau->n != Qd->n)
			Tau=m_resize(Tau,Qd->m,Qd->n);
		
		/* allocate temporals */
		v_resize_vars(DYN->m,&b,&tcp,&tcm,&v1,&v2,NULL);
		mem_stat_reg_vars(0,TYPE_VEC,&b,&tcp,&tcm,&v1,&v2,NULL);
		m1=m_resize(m1,b->dim,b->dim);
		m_resize_vars(Qd->m,Qd->n,&m2,&m3,NULL);
		mem_stat_reg_vars(0,TYPE_MAT,&m1,&m2,&m3,NULL);
		
		/* get friction parameters */
		get_col(DYN,12,b);
		if (DYN->n >= 13)
		{
			get_col(DYN,13,tcp);
			if (DYN->n >= 14)
				get_col(DYN,14,tcm);
			else 
				sv_mlt(-1,tcp,tcm);
		} 
		else 
		{
			v_zero(tcp);
			v_copy(tcm,tcp);
		}
		/* refer friction values to link */
		
		get_col(DYN,11,v1);
		v_star(v1,v1,v2);
		v_star(b,v2,b);
		v_star(tcp,v1,tcp);
		v_star(tcm,v1,tcm);
		
		for (i=0;i<(int)b->dim;i++) m1->me[i][i]=b->ve[i];
		m_mlt(Qd,m1,Tau);
		
		for (i=0;i<(int)b->dim;i++) m1->me[i][i]=tcp->ve[i];
		for (i=0;i<(int)Qd->m;i++) 
			for (j=0;j<(int)Qd->n;j++)
				if (Qd->me[i][j]>0) m2->me[i][j]=1;
				else m2->me[i][j]=0;
				m_add(Tau,m_mlt(m2,m1,m3),Tau);
				
				for (i=0;i<(int)b->dim;i++) m1->me[i][i]=tcm->ve[i];
				for (i=0;i<(int)Qd->m;i++) 
					for (j=0;j<(int)Qd->n;j++)
						if (Qd->me[i][j]<0) m2->me[i][j]=1; 
						else m2->me[i][j]=0;
						m_add(Tau,m_mlt(m2,m1,m3),Tau);
						
						/*
						tcp = tcp.*dh_dyn(:,17);
						tcm = tcm.*dh_dyn(:,17);
						tau = Qd*diag(b) + (Qd>0)*diag(tcp) + (Qd<0)*diag(tcm);
						*/
						return Tau;
}

