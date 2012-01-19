/*************************************************************************
Workspace used: 5
*************************************************************************
Compute base velocities respect to Mass center frame

void base_v(DH,DYN,DYN_base,Q,X)

DH: Denavit Hartenberg parameters matrix
DYN: Manipulator dynamic parameters
Q: joint Positions

This program Uses euler.c, fkine.c, tr2rpy.c

Robotics and Automation Group
	Andres Jaramillo Botero
	Juan Camilo Acosta
	Julian David Colorado
	Antonio Matta
Pontificia Universidad Javeriana, Cali.
***************************************************************************/
/*31/05/2006  JDC and JCA, Compute base velocities respect to Mass center frame */

#include "LLRL.h"
MAT *base_v(MAT *DH,MAT *DYN,MAT *DYN_base,MAT *Q_Qd, MAT *x, MAT *v)
{
	int points, df,pt,art,i;
	Real MTotal,vtemp;
	static VEC *s=VNULL,*temp=VNULL,*temp2=VNULL,*t1=VNULL,*t2=VNULL,*p=VNULL,*dQt=VNULL, *Qt=VNULL;
	static MAT *T=MNULL,*r=MNULL,*Ts=MNULL,*Tt=MNULL, *SkewP=MNULL, *identidad=MNULL, *P=MNULL, *P_trans=MNULL, *R=MNULL, *PtR=MNULL, *DHt=MNULL, *jac=MNULL;
	
	if(Q_Qd==MNULL)
		error(E_NULL,"wrong size for Q_Qd");
	points=Q_Qd->m;
	df=DH->m;	
	if(x==MNULL || x->m!=6 || x->n!=points)
		error(E_NULL,"wrong size for x in base_v");

	m_resize_vars(6,points,&v,MNULL);
	
/**********************MATRIX DECLARATION AND REGISTRATION*******************************************/
	DHt=m_resize(DHt,1,1);
	jac=m_resize(jac,1,1);
	m_resize_vars(3,3,&r,&SkewP,&identidad,NULL);
	m_resize_vars(4,4,&T,&Ts,&Tt,NULL);
	m_resize_vars(6,6,&P,&P_trans,&R,&PtR,NULL);
	mem_stat_reg_vars(0,TYPE_MAT,&T,&Ts,&Tt,&r,&identidad,&P,&P_trans,&R,&PtR,&SkewP,&jac,&DHt,NULL);
	
/**********************VECTOR DECLARATION AND REGISTRATION*******************************************/	
	v_resize_vars(6,&s,&temp,&temp2,&t1,&t2,NULL);
	Qt=v_resize(Qt,1);
	dQt=v_resize(dQt,1);
	p=v_resize(p,4);
	mem_stat_reg_vars(0,TYPE_VEC,&s,&temp,&temp2,&t1,&t2,&p,&Qt,&dQt,NULL);
	
//******************************************************************************************
//Obtain base velocities
//******************************************************************************************	
	m_ident(identidad);
	MTotal=0;
	for (pt=0;pt<points;pt++)
	{
		v_zero(temp2);
		v_zero(temp);
		for (art=0;art<df;art++)
		{
			if(pt==0)
			{	
				if(art==0)
					MTotal=MTotal+DYN_base->me[0][0]+DYN->me[art][0];
				else
					MTotal=MTotal+DYN->me[art][0];
			}
			DHt=m_resize(DHt,art+1,DH->n);
			Qt=v_resize(Qt,art+1);
		
			if(Q_Qd->m!=1)
			{
				for(i=0;i<=art;i++)
				{
					Qt->ve[i]=Q_Qd->me[pt][i+1];
				}
			}
			else
			{
				for(i=0;i<=art;i++)
				{
					Qt->ve[i]=Q_Qd->me[pt][i];
				}
			}
			m_move(DH,0,0,art+1,DH->n,DHt,0,0);
		
			//printf("antes de jacob0\n");
			//mem_info();
			mem_stat_mark(5);
			if(LLRL_ERROR==jacob0(jacobn_float,DHt,DYN,Qt,jac))
				error(E_NULL,"Problem with jacob0");
			mem_stat_free(5);
			//printf("despues de jacob0\n");
			//mem_info();
			dQt=v_resize(dQt,art+1);
	
			if(Q_Qd->m!=1)
			{
				for(i=0;i<=art;i++)
					dQt->ve[i]=Q_Qd->me[pt][df+i+1];
			}
			else
			{
				for(i=0;i<=art;i++)
					dQt->ve[i]=Q_Qd->me[pt][df+i];
			}
			
			mv_mlt(jac,dQt,temp2);
			sv_mlt(DYN->me[art][0],temp2,temp2);
			v_add(temp2,temp,temp);
		}
		
		sv_mlt(1/MTotal,temp,temp2);
	
		temp->ve[0]=temp2->ve[3];
		temp->ve[1]=temp2->ve[4];
		temp->ve[2]=temp2->ve[5];
		temp->ve[3]=temp2->ve[0];
		temp->ve[4]=temp2->ve[1];
		temp->ve[5]=temp2->ve[2];
		
		get_col(x,pt,temp2);
		euler(temp2,r,p);
	
		//Obtain the 6-dimensional P operator based on p
		skew_symetric(SkewP,p);
		m_move(identidad,0,0,3,3,P,0,0);
		m_move(identidad,0,0,3,3,P,3,3);
		m_move(SkewP,0,0,3,3,P,0,3);
		
		//Obtain the 6-dimensional R operator based on r
		sm_mlt(-1,r,r);
		m_move(r,0,0,3,3,R,0,0);
		m_move(r,0,0,3,3,R,3,3);
		
		m_transp(P,P_trans);
		m_mlt(P_trans,R,PtR);

		mv_mlt(PtR,temp,temp2);
		vtemp=temp2->ve[2];
		temp2->ve[2]=temp2->ve[0];
		temp2->ve[0]=vtemp;

		_set_col(v,pt,temp2,0);
		
	}
	return v;
}