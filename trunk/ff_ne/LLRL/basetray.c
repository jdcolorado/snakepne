/*************************************************************************
Workspace used: 5
*************************************************************************
Compute base position, velocity and acceleration for floating base

void basetray(x,v,a,DH,DYN,DYN_base,Q_Qd_Qdd)

DH: Denavit Hartenberg parameters matrix
DYN: Manipulator dynamic parameters
DYN_base: Manipulator base DYnamic parameters
Q_Qd_Qdd: joint trayectory,time,x(base positions),v(base velocities),a(base accelerations)

This program Uses euler.c, fkine.c, tr2rpy.c

Robotics and Autonation Group
	Andres Jaramillo Botero
	Juan Camilo Acosta
	Julian David Colorado
	Antonio Matta
Pontificia Universidad Javeriana, Cali.
***************************************************************************/
//16/03/2006: JCA: Velocity and accelerations not working. Position: OK
//17/03/2006: Added support for base DYN parameters. Velocity and accelerations now working
//	     Resize made with Q_Qd_Qdd and Xb,Vb, Ab in order to remove first and last point.
//21/03/2006: removed x, v and a as parameters, from now included in q_qd_qdd
//		Changed time vector to first Q_Qd_Qdd column instead of last.
//21/03/2006: Changed DYN_base from VEC to MAT because the base may have more than 1 bodys
/*22/03/2006, JCA: Added support for joint limits, DH complete*/
/*25/04/2006: JCA: added suport for direct velocities: Memory OK*/
/*12/05/2006  JDC, changed interface; added Robot type */

#include "LLRL.h"
MAT *basetray(int Robot_type,MAT *DH,MAT *DYN,MAT *DYN_base,MAT *Q_Qd_Qdd)
{   
	int points, df,pt,art,i;
	Real MTotal,step=0.1,vtemp;
	static VEC *s=VNULL,*temp=VNULL,*temp2=VNULL,*t1=VNULL,*t2=VNULL,*p=VNULL,*time=VNULL, *dQt=VNULL, *Qt=VNULL;
	static MAT *DHtemp=MNULL,*qtemp=MNULL,*T=MNULL,*r=MNULL,*Ts=MNULL,*Tt=MNULL, *x=MNULL, *SkewP=MNULL, *identidad=MNULL, *P=MNULL, *P_trans=MNULL, *R=MNULL, *PtR=MNULL, *DHt=MNULL, *jac=MNULL;
	
	points=Q_Qd_Qdd->m;
	df=Q_Qd_Qdd->n/3;
	
	//m_resize_vars(6,points,&x,&v,&a,MNULL);
	
/**********************MATRIX DECLARATION AND REGISTRATION*******************************************/
	DHt=m_resize(DHt,1,1);
	x=m_resize(x,6,Q_Qd_Qdd->m);
	jac=m_resize(jac,1,1);
	m_resize_vars(3,3,&r,&SkewP,&identidad,NULL);
	m_resize_vars(4,4,&T,&Ts,&Tt,NULL);
	m_resize_vars(6,6,&P,&P_trans,&R,&PtR,NULL);
	mem_stat_reg_vars(0,TYPE_MAT,&T,&Ts,&Tt,&r,&x,&identidad,&P,&P_trans,&R,&PtR,&SkewP,&jac,&DHt,NULL);
	DHtemp=m_resize(DHtemp,1,1);
	qtemp=m_resize(qtemp,1,1);
	mem_stat_reg_vars(0,TYPE_MAT,&DHtemp,&qtemp,NULL);

/**********************VECTOR DECLARATION AND REGISTRATION*******************************************/	
	v_resize_vars(6,&s,&temp,&temp2,&t1,&t2,NULL);
	Qt=v_resize(Qt,1);
	dQt=v_resize(dQt,1);
	p=v_resize(p,4);
	time=v_resize(time,points);
	mem_stat_reg_vars(0,TYPE_VEC,&s,&temp,&temp2,&t1,&t2,&p,&time,&Qt,&dQt,NULL);
	if(Q_Qd_Qdd->m!=1)
		get_col(Q_Qd_Qdd,0,time);
	Q_Qd_Qdd=m_resize(Q_Qd_Qdd,points,Q_Qd_Qdd->n+18);
	
	//Positions
	
//******************************************************************************************
//Obtain base positions
//******************************************************************************************	
	MTotal=0;	     
	for (pt=0;pt<points;pt++)
	{
		v_zero(temp2);

		for (art=0;art<=df;art++)
		{
			v_zero(s); 
			//Base
			if (art==0)
			{
				//For each body of base.
				for(i=0;i<DYN_base->m;i++)
				{
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
				
				if(Q_Qd_Qdd->m!=1)
					m_move(Q_Qd_Qdd,pt,0+1,1,art,qtemp,0,0);
				else
					m_move(Q_Qd_Qdd,pt,0,1,art,qtemp,0,0);
				
				
				//running fkine from base to joint art
				mem_stat_mark(5);
				if (fkine(DHtemp,qtemp,T)==LLRL_ERROR) 
					error(E_SIZES,"fkine");
				mem_stat_free(5);
				
				s->ve[3]=DYN->me[art-1][1];  
				s->ve[4]=DYN->me[art-1][2];
				s->ve[5]=DYN->me[art-1][3];
				
				//6x1 --> 4x4 (Euler -->homogeneus)
				euler(s,r,p);
				m_zero(Ts);
				m_move(r,0,0,3,3,Ts,0,0);
				p->ve[3]=1;
				_set_col(Ts,3,p,0);
				Ts->me[0][0];
				m_mlt(T,Ts,Tt);
				temp=v_resize(temp,3);
				
				//4x4 --> 6x1 (homogeneus-->euler)
				mem_stat_mark(5);
				temp=tr2rpy(Tt,temp);
				mem_stat_free(5);
				
				temp=v_resize(temp,6);
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
		m_ident(T);
		m_ident(Ts);
		
		m_move(r,0,0,3,3,T,0,0);
		p->ve[3]=1;
		_set_col(Ts,3,p,0);
		//Rotating...
		m_mlt(T,Ts,Tt);
		temp=v_resize(temp,3);
		
		//4x4 --> 6x1 (homogeneus-->euler)
		mem_stat_mark(5); 
		temp=tr2rpy(Tt,temp);
		mem_stat_free(5);
		
		temp=v_resize(temp,6);
		
		temp->ve[3]=Tt->me[0][3];
		temp->ve[4]=Tt->me[1][3];
		temp->ve[5]=Tt->me[2][3];
		_set_col(x,pt,temp,0);
			
		if(Q_Qd_Qdd->m!=1)
		{
			Q_Qd_Qdd->me[pt][3*df+1]=temp->ve[0];
			Q_Qd_Qdd->me[pt][3*df+1+1]=temp->ve[1];
			Q_Qd_Qdd->me[pt][3*df+1+2]=temp->ve[2];
			Q_Qd_Qdd->me[pt][3*df+1+3]=temp->ve[3];
			Q_Qd_Qdd->me[pt][3*df+1+4]=temp->ve[4];
			Q_Qd_Qdd->me[pt][3*df+1+5]=temp->ve[5];
		}
		else
		{
			Q_Qd_Qdd->me[pt][3*df]=temp->ve[0];
			Q_Qd_Qdd->me[pt][3*df+1]=temp->ve[1];
			Q_Qd_Qdd->me[pt][3*df+2]=temp->ve[2];
			Q_Qd_Qdd->me[pt][3*df+3]=temp->ve[3];
			Q_Qd_Qdd->me[pt][3*df+4]=temp->ve[4];
			Q_Qd_Qdd->me[pt][3*df+5]=temp->ve[5];
		}
		
	}
	
	
//******************************************************************************************
//Obtain base velocities
//******************************************************************************************	
	m_ident(identidad);

	for (pt=0;pt<points;pt++)
	{
		v_zero(temp2);
		v_zero(temp);
		for (art=0;art<df;art++)
		{
			DHt=m_resize(DHt,art+1,DH->n);
			Qt=v_resize(Qt,art+1);
		
			if(Q_Qd_Qdd->m!=1)
			{
				for(i=0;i<=art;i++)
				{
					Qt->ve[i]=Q_Qd_Qdd->me[pt][i+1];
				}
			}
			else
			{
				for(i=0;i<=art;i++)
				{
					Qt->ve[i]=Q_Qd_Qdd->me[pt][i];
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
	
			if(Q_Qd_Qdd->m!=1)
			{
				for(i=0;i<=art;i++)
					dQt->ve[i]=Q_Qd_Qdd->me[pt][df+i+1];
			}
			else
			{
				for(i=0;i<=art;i++)
					dQt->ve[i]=Q_Qd_Qdd->me[pt][df+i];
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
		if(Q_Qd_Qdd->m!=1)
		{
			Q_Qd_Qdd->me[pt][3*df+1+6]=temp2->ve[0];
			Q_Qd_Qdd->me[pt][3*df+1+7]=temp2->ve[1];
			Q_Qd_Qdd->me[pt][3*df+1+8]=temp2->ve[2];
			Q_Qd_Qdd->me[pt][3*df+1+9]=temp2->ve[3];
			Q_Qd_Qdd->me[pt][3*df+1+10]=temp2->ve[4];
			Q_Qd_Qdd->me[pt][3*df+1+11]=temp2->ve[5];
		}
		else
		{
			Q_Qd_Qdd->me[pt][3*df+6]=temp2->ve[0];
			Q_Qd_Qdd->me[pt][3*df+7]=temp2->ve[1];
			Q_Qd_Qdd->me[pt][3*df+8]=temp2->ve[2];
			Q_Qd_Qdd->me[pt][3*df+9]=temp2->ve[3];
			Q_Qd_Qdd->me[pt][3*df+10]=temp2->ve[4];
			Q_Qd_Qdd->me[pt][3*df+11]=temp2->ve[5];
		}
	}
	
/*
		//Q_Qd_Qdd=m_resize(Q_Qd_Qdd,points,Q_Qd_Qdd->n+12);
		for (pt=1;pt<points-1;pt++)
		{
			step=time->ve[pt]-time->ve[pt-1];
			if (pt==1 || pt==points-2)
			{
				//v(t)=(x(t+h)-x(t-h))/(2*h)
				get_col(x,pt-1,temp);
				get_col(x,pt+1,temp2);
				v_sub(temp2,temp,temp2);
				sv_mlt(1/(2*step),temp2,temp);
				//_set_col(v,pt-1,temp,0);
				Q_Qd_Qdd->me[pt][3*df+1+6]=temp->ve[0];
				Q_Qd_Qdd->me[pt][3*df+1+7]=temp->ve[1];
				Q_Qd_Qdd->me[pt][3*df+1+8]=temp->ve[2];
				Q_Qd_Qdd->me[pt][3*df+1+9]=temp->ve[3];
				Q_Qd_Qdd->me[pt][3*df+1+10]=temp->ve[4];
				Q_Qd_Qdd->me[pt][3*df+1+11]=temp->ve[5];
			}
			else
			{
				//v(t)=(-x(t+2h)+8x(t+h)-8x(t-h)+x(t-2h))/(12*h)
				get_col(x,pt+1,temp);
				get_col(x,pt+2,temp2);
				get_col(x,pt-1,t1);
				get_col(x,pt-2,t2);
				sv_mlt(8,temp,temp);
				sv_mlt(-1,temp2,temp2);
				sv_mlt(-8,t1,t1);
				v_add(temp,temp2,temp2);
				v_add(temp2,t1,t1);
				v_add(t1,t2,t2);
				sv_mlt(1/(12*step),t2,temp);
				//_set_col(v,pt-1,temp,0);
				Q_Qd_Qdd->me[pt][3*df+1+6]=temp->ve[0];
				Q_Qd_Qdd->me[pt][3*df+1+7]=temp->ve[1];
				Q_Qd_Qdd->me[pt][3*df+1+8]=temp->ve[2];
				Q_Qd_Qdd->me[pt][3*df+1+9]=temp->ve[3];
				Q_Qd_Qdd->me[pt][3*df+1+10]=temp->ve[4];
				Q_Qd_Qdd->me[pt][3*df+1+11]=temp->ve[5];
			}
		}
	}*/
	//mem_info();
	//m_output(x);
	
	//Accelerations (differentiating.)
		
	for (pt=1;pt<points-1;pt++)
	{
		step=time->ve[pt]-time->ve[pt-1];
		if (pt==1 || pt==points-2)
		{
			//a(t)=(x(t+h)-2*x(t)+x(t-h))/(h²)
			get_col(x,pt-1,temp);

			get_col(x,pt+1,temp2);

			get_col(x,pt,t1);
			sv_mlt(-2,t1,t1);
			v_add(temp2,temp,temp2);
			v_add(temp2,t1,t2);
			sv_mlt(1/(pow(step,2)),t2,temp);

			//_set_col(a,pt-1,temp,0);
			Q_Qd_Qdd->me[pt][3*df+1+12]=temp->ve[0];
			Q_Qd_Qdd->me[pt][3*df+1+13]=temp->ve[1];
			Q_Qd_Qdd->me[pt][3*df+1+14]=temp->ve[2];
			Q_Qd_Qdd->me[pt][3*df+1+15]=temp->ve[3];
			Q_Qd_Qdd->me[pt][3*df+1+16]=temp->ve[4];
			Q_Qd_Qdd->me[pt][3*df+1+17]=temp->ve[5];
		}
		else
		{
			//a(t)=(-x(t+2h)+16x(t+h)-30x(t)+16x(t-h)-x(t-2h))/(12*h²)
			get_col(x,pt+1,temp);
			get_col(x,pt+2,temp2);
			get_col(x,pt-1,t1);
			get_col(x,pt-2,t2);
			get_col(x,pt,s);
			sv_mlt(16,temp,temp);
			sv_mlt(-1,temp2,temp2);
			sv_mlt(-1,t2,t2);
			sv_mlt(16,t1,t1);
			sv_mlt(-30,s,s);
			v_add(temp,temp2,temp2);
			v_add(temp2,t2,t2);
			v_add(t2,t1,t1);
			v_add(t1,s,t2);
			sv_mlt(1/(12*pow(step,2)),t2,temp);
			//_set_col(a,pt-1,temp,0);
			Q_Qd_Qdd->me[pt][3*df+1+12]=temp->ve[0];
			Q_Qd_Qdd->me[pt][3*df+1+13]=temp->ve[1];
			Q_Qd_Qdd->me[pt][3*df+1+14]=temp->ve[2];
			Q_Qd_Qdd->me[pt][3*df+1+15]=temp->ve[3];
			Q_Qd_Qdd->me[pt][3*df+1+16]=temp->ve[4];
			Q_Qd_Qdd->me[pt][3*df+1+17]=temp->ve[5];
		}
	}
	//removing first and last point
	if(Q_Qd_Qdd->m!=1)
	{
		m_move(Q_Qd_Qdd,1,0,points-1,Q_Qd_Qdd->n,Q_Qd_Qdd,0,0);
		Q_Qd_Qdd=m_resize(Q_Qd_Qdd,points-2,Q_Qd_Qdd->n);
	}
	return(Q_Qd_Qdd);
}

