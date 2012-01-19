/*************************************************************************
Workspace used: 5
*************************************************************************
Compute base position respect to Mass center frame

void base_x(DH,DYN,DYN_base,Q,X)

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
/*30/05/2006  JDC, Compute base position respect to Mass center frame */

#include "LLRL.h"
MAT *base_x(MAT *DH,MAT *DYN,MAT *DYN_base,MAT *Q_Qd, MAT *x)
{ 
	int points, df,pt,art,i;
	Real MTotal;
	static VEC *s=VNULL,*temp=VNULL,*temp2=VNULL,*t1=VNULL,*t2=VNULL, *p=VNULL;
	static MAT *DHtemp=MNULL,*qtemp=MNULL,*T=MNULL,*r=MNULL,*Ts=MNULL, *Tr=MNULL, *Ti=MNULL, *Tt=MNULL;
	
	if(Q_Qd==MNULL)
		error(E_NULL,"wrong size for Q_Qd");
	
	points=Q_Qd->m;
	df=DH->m;
	m_resize_vars(6,points,&x,MNULL);
	
/**********************MATRIX DECLARATION AND REGISTRATION*******************************************/
	
	m_resize_vars(3,3,&r,NULL);
	m_resize_vars(4,4,&T,&Ts,&Tt,&Ti,&Tr,NULL);
	mem_stat_reg_vars(0,TYPE_MAT,&T,&Ts,&Tt,&Ti,&Tr,&r,NULL);
	DHtemp=m_resize(DHtemp,1,1);
	qtemp=m_resize(qtemp,1,1);
	mem_stat_reg_vars(0,TYPE_MAT,&DHtemp,&qtemp,NULL);

/**********************VECTOR DECLARATION AND REGISTRATION*******************************************/	
	v_resize_vars(6,&s,&temp,&temp2,&t1,&t2,NULL);
	p=v_resize(p,4);
	mem_stat_reg_vars(0,TYPE_VEC,&s,&temp,&temp2,&t1,&t2,&p,NULL);
	
	m_ident(Ti);
	
	Tr->me[3][3]=1;

//Positions
//******************************************************************************************
//Obtain base positions
//******************************************************************************************	
	MTotal=0;
	for (pt=0;pt<points;pt++)
	{
		//v_zero(temp2);
		for (art=0;art<=df;art++)
		{
			
			//Base
			if (art==0)
			{
				//For each body of base.
				for(i=0;i<DYN_base->m;i++)
				{
					s->ve[0]=0;
					s->ve[1]=0;
					s->ve[2]=0;
					s->ve[3]=DYN_base->me[i][1];
					s->ve[4]=DYN_base->me[i][2];
					s->ve[5]=DYN_base->me[i][3];
					//mass*distance
					sv_mlt(DYN_base->me[i][0],s,temp2);
					if(pt==0)
						MTotal=MTotal+DYN_base->me[i][0];
				}
			}
			else
			{
				//mem_info();
				DHtemp=m_resize(DHtemp,art,DH->n);
				m_move(DH,0,0,art,DH->n,DHtemp,0,0);
				qtemp=m_resize(qtemp,1,art);
				
				if(Q_Qd->m!=1)
					m_move(Q_Qd,pt,0,1,art,qtemp,0,0);
				else
					m_move(Q_Qd,pt,0,1,art,qtemp,0,0);
				
				
				//running fkine from base to joint art
				mem_stat_mark(5);
				if (fkine(DHtemp,qtemp,T)==LLRL_ERROR) 
					error(E_SIZES,"fkine");
				mem_stat_free(5);
				
				s->ve[0]=0;
				s->ve[1]=0;
				s->ve[2]=0;
				s->ve[3]=DYN->me[art-1][1];  
				s->ve[4]=DYN->me[art-1][2];
				s->ve[5]=DYN->me[art-1][3];
				
				//6x1 --> 4x4 (Euler -->homogeneus)
				euler(s,r,p);
				//m_zero(Ts);
				m_move(r,0,0,3,3,Ts,0,0);
				p->ve[3]=1;
				_set_col(Ts,3,p,0);
				
				m_mlt(T,Ts,Tt);
				//temp=v_resize(temp,3);
				
				//4x4 --> 6x1 (homogeneus-->euler)
				mem_stat_mark(5);
				temp=tr2rpy(Tt,temp);
				mem_stat_free(5);
				
				//temp=v_resize(temp,6);
				temp->ve[3]=Tt->me[0][3];
				temp->ve[4]=Tt->me[1][3];
				temp->ve[5]=Tt->me[2][3];
				
				//mass*"distance"
				sv_mlt(DYN->me[art-1][0],temp,t2);
				v_add(t2,temp2,temp2);
				
				if(pt==0)
					MTotal=MTotal+DYN->me[art-1][0];
			}
		}
			
		sv_mlt(-1/MTotal,temp2,t1);
		//6x1 --> 4x4 (Euler -->homogeneus)
		
		euler(t1,r,p);
		//m_ident(T);
		
		m_move(r,0,0,3,3,Tr,0,0);
		p->ve[3]=1;
		_set_col(Ti,3,p,0);
		//Rotating...
		m_mlt(Tr,Ti,Tt);
		//temp=v_resize(temp,3);
		
		//4x4 --> 6x1 (homogeneus-->euler)
		mem_stat_mark(5); 
		temp=tr2rpy(Tt,temp);
		mem_stat_free(5);
		
		//temp=v_resize(temp,6);
		
		temp->ve[3]=Tt->me[0][3];
		temp->ve[4]=Tt->me[1][3];
		temp->ve[5]=Tt->me[2][3];
		_set_col(x,pt,temp,0);
	}
	return x;
}
