#include"LLRL.h"
/*************************************************************************	
	Adapted from matlab's 5.3 version of the PINV function
	Matrix Computations, Golub & van Loan, Third Edition

	Copyright, Andrés Jaramillo Botero, 1999. 
***************************************************************************/
MAT	*pinv(MAT *A,MAT *PIA)
{
	static	MAT	*U,*Q,*m1,*m2,*m3;
	static	VEC	*u,*s;
	Real	tol;
	int		index,r,i;
	
	if (A == (MAT *)MNULL) error(E_NULL,"Null input in pinv\n");
	if (PIA == (MAT *)MNULL) PIA=m_resize(PIA,A->n,A->m);
	
	/* [Q u U]=svd(A) */
	U = m_resize(U,A->n,A->n);	
	Q = m_resize(Q,A->m,A->m);	
	mem_stat_reg_vars(0,TYPE_MAT,&U,&Q,NULL);
	u = v_resize(u,max(A->m,A->n));
	MEM_STAT_REG(u,TYPE_VEC);
	svd(A,Q,U,u);
	
	tol=max(A->m,A->n)*v_max(u,&index)*MACHEPS;
	/* define a projection vector from non zero elements of u */
	r=0;
	for (i = 0; i < (int)u->dim; i++)
		if (u->ve[i]>tol) r++; 
		
		/* if r=0 the pseudoinverse of A is null */
		if (r==0) {PIA=m_resize(PIA,A->n,A->m); m_zero(PIA);}
		else 
		{
			u=v_resize(u,r);
            MEM_STAT_REG(u,TYPE_VEC);
			s=v_resize(s,r);
			MEM_STAT_REG(s,TYPE_VEC);
			v_ones(s);
			v_slash(u,s,s);
			/* s = diag(ones(r,1)./s(1:r)); element by element division */
			
			Q=m_resize(Q,r,A->m);
			U=m_resize(U,r,A->n);
            mem_stat_reg_vars(0,TYPE_MAT,&U,&Q,NULL);
            
			m2=m_resize(m2,r,A->m);
			MEM_STAT_REG(m2,TYPE_MAT);
			m_resize_vars(r,r,&m1,&m3,NULL);
			mem_stat_reg_vars(0,TYPE_MAT,&m1,&m3,NULL);
			for (i = 0; i<r;i++) m1->me[i][i]=s->ve[i]; 
			m_transp(U,m2);
			m3=m_mlt(m2,m1,m3);
			PIA=m_resize(PIA,m2->m,Q->n);
			/* PIA = U(:,1:r)*m1*Q(:,1:r)'; */
			PIA=m_mlt(m3,Q,PIA);
		}
		return PIA;
}

