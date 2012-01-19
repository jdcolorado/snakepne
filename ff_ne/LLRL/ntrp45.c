#include "LLRL.h"
/**************************************************************************	
	NTRP45 Interpolation helper function for ODE45.
	YINTERP = NTRP45(TINTERP,T,Y,TNEW,YNEW,H,F) uses data computed in ODE45
	to approximate the solution at time TINTERP.
***************************************************************************/
/* 07-11-2002-AJB: correction in MEM_STAT_REG of MAT not VEC */

//funcion ntrp45 de matlab.  interpolacion para odb45
MAT *ntrp45(MAT *tinterp,double t,MAT *y,double h, MAT *f, MAT *resul)
{	
	static	MAT *BI,*s,*p2_3,*p_2,*tmp,*tmp2;
	int		i,j,mt,my;

	/* validate inputs (sizes and initializations) */
	if (tinterp == (MAT *)MNULL || y == (MAT *)MNULL || f == (MAT *)MNULL) 
		error(E_NULL,"ntrp45\n");

	my=y->m;
	mt=tinterp->m;
	BI=m_resize(BI,7,4);
	MEM_STAT_REG(BI,TYPE_MAT);
	m_zero(BI);
	s=m_resize(s,1,mt);
	MEM_STAT_REG(s,TYPE_MAT);
	m_zero(s);
	
	p2_3=m_resize(p2_3,4,mt);
	p_2=m_resize(p_2,my,mt);
	resul=m_resize(resul,my,mt);
	
	mem_stat_reg_vars(0,TYPE_MAT,&p2_3,&p_2,NULL);
	m_zero(p2_3);
	m_zero(p_2);
	m_zero(resul);

	BI->me[0][0]=1.0;	BI->me[0][1]=-183.0/64.0;	BI->me[0][2]=37.0/12.0;		BI->me[0][3]=-145.0/128.0;
	BI->me[1][0]=0.0;	BI->me[1][1]=0.0;		BI->me[1][2]=0.0;		BI->me[1][3]=0.0;
	BI->me[2][0]=0.0;	BI->me[2][1]=1500.0/371.0;	BI->me[2][2]=-1000.0/159.0;	BI->me[2][3]=1000.0/371.0;
	BI->me[3][0]=0.0;	BI->me[3][1]=-125.0/32.0;	BI->me[3][2]=125.0/12.0;	BI->me[3][3]=-375.0/64.0;
	BI->me[4][0]=0.0;	BI->me[4][1]=9477.0/3392.0;	BI->me[4][2]=-729.0/106.0;	BI->me[4][3]=25515.0/6784.0;
	BI->me[5][0]=0.0;	BI->me[5][1]=-11.0/7.0;		BI->me[5][2]=11.0/3.0;		BI->me[5][3]=-55.0/28.0;
	BI->me[6][0]=0.0;	BI->me[6][1]=3.0/2.0;		BI->me[6][2]=-4.0;		BI->me[6][3]=5.0/2.0;
	
	//y(:,ones(length(tinterp),1))
	for(i=0;i<mt;i++)
	{	
		s->me[0][i]=((tinterp->me[i][0]-t)/h);
		for(j=0;j<my;j++)
			resul->me[j][i]=y->me[j][0];
	}
    //cumprod(s(ones(4,1),:));
	for(i=0;i<mt;i++)
	{	
		p2_3->me[0][i]=s->me[0][i];
		p2_3->me[1][i]=p2_3->me[0][i]*s->me[0][i];
		p2_3->me[2][i]=p2_3->me[1][i]*s->me[0][i];
		p2_3->me[3][i]=p2_3->me[2][i]*s->me[0][i];
	}
	tmp=m_resize(tmp,BI->m,1);
	MEM_STAT_REG(tmp,TYPE_MAT);
	tmp2=m_resize(tmp2,BI->m,BI->n);
	MEM_STAT_REG(tmp2,TYPE_MAT);
	tmp2=sm_mlt(h,BI,tmp2);
	p_2=m_mlt(f,m_mlt(tmp2,p2_3,tmp),p_2);

	//yinterp = y(:,ones(length(tinterp),1)) + f*(h*BI)*cumprod(s(ones(4,1),:));
	for(i=0;i<my;i++)
		for(j=0;j<mt;j++)
			resul->me[i][j]=resul->me[i][j]+p_2->me[i][j];

	return resul;
}
