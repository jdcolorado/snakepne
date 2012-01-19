#include"LLRL.h"
/******************************************************************************
Workspace used: 11,14,3,5,6,8,9
*******************************************************************************	
verlet solves Newton Euler differential equation with velocity verlet.

[T Q QD QDD] = verlet(func ,DH, DYN, T0, T1, Q0, dQ0,Torques, step, xb, vb, ab, grav, fext)

Integrates de forward dynamics (newton euler) using velocity verlet and returns
joint position, velocity and acceleration for each degree of freedom in Q0_Qd0.
It uses inverse dynamics algorithm (ne) to calculate mass operator 

A control torque needs to be specified by a user specified function
Torque last column is time vector

Copyright Robotics and Automation Group, Pontificia Universidad Javeriana - Cali.
		Andrés Jaramillo Botero, ajaramil@puj.edu.co
		Antonio Alejandro Matta Gómez, amatta@puj.edu.co
		Julián David Colorado, jdcolorado@puj.edu.co
		Juan Camilo Acosta, jcacosta@puj.edu.co
*************************************************************************/
/*06/02/2006, JCA: First version: fix or flying base, not support floating base*/
/*14/02/2006, JCA: velocity verlet not working, calculating velocities with euler.*/
/*16/02/2006, JCA: Improved memory management. (mem_stat_reg_vars). Added support
for differentiation in order to get velocity from positions. Still no support for floating base.*/
/*16/02/2006, JCA: changed output, it was Q_Qd_Qdd, is now: time_Q_Qd*/
/*17/02/2006, JCA: Improved memory management. No more func(,,,,NULL)*/
/*20/02/2006, JCA: Improved memory management,*/
/*22/02/2006, JCA, AM: Solved problem with memory in ff_ne calls*/
/*22/02/2006, JCA, JDC:Add spatial base position: X=[roll pitch yaw x y z]'*/
/*22/03/2006, JCA: Added support for joint limits*/
/*30/03/2006, JCA: Added pass of grav and fext. modified to work with torqfun*/
/*26/04/2006, JCA: changed interface, added Robot_type and deleted base input complete traj.
/*07/05/2006 JDC, Changed interface, added DYN_base*/
/*07/05/2006 JDC, Changed prototype call for feval, added DYN_base*/


MAT *verlet (MAT *(*fun)(), int Robot_type, MAT *DH, MAT *DYN, VEC *tspan, VEC *Q0_Qd0, MAT *Torques, MAT *Q_Qd_Qdd, double step, VEC *grav, VEC *fext,MAT *DYN_base)

{
	double pt;
	int i, j,df,points;
	static MAT *q=NULL,*dq=NULL,*d2q=NULL, *Qd_Q2d=NULL, *Qd_Q2dt=NULL,*Q0Qd0=NULL,*Ttemp=NULL;
	static VEC *dQ_temp=NULL,*time=NULL, *v_tmp=NULL, *a_tmp=NULL, *r_tmp=NULL, *inta_tmp=NULL, *intv_tmp=NULL, *dq_tmp=NULL, *q_tmp1=NULL, *q_tmp2=NULL, *q_tmp3=NULL; //,*Vbt=NULL,*dVbt=NULL, *Xt=NULL,*Tor=MNULL
	df = Q0_Qd0->dim/2; //df=degrees of freedom
	//df=Torques->n-1;
	m_resize_vars(2*df,1,&Q0Qd0,&Qd_Q2d,NULL);
	mem_stat_reg_vars(0,TYPE_MAT,&Q0Qd0,&Qd_Q2d,NULL);
	Qd_Q2dt=m_resize(Qd_Q2dt,1,2*df);
	MEM_STAT_REG(Qd_Q2dt,TYPE_MAT);
	//mem_stat_reg_vars(0,TYPE_MAT,&Q0Qd0,&Qd_Q2d,&Qd_Q2dt,NULL);
	

	if(Torques == MNULL)
	{
		//paso chambon, el que encuentre otra forma.....  dado que round no funciona ni roundf
		pt=(tspan->ve[1]-tspan->ve[0])/step+1;
		if((int)(pt)<pt)
			points=ceil((tspan->ve[1]-tspan->ve[0])/step+1);
		else
			points=floor((tspan->ve[1]-tspan->ve[0])/step+1);
		Ttemp=m_resize(Ttemp,points,df);
		MEM_STAT_REG(Ttemp,TYPE_MAT);
		time=v_resize(time,points);
		MEM_STAT_REG(time,TYPE_VEC);
		time->ve[0]=0;
		m_zero(Ttemp);
		
		for(j=0;j<points;j++)
			if(j>0)
				time->ve[j]=time->ve[j-1]+step;
	}
	else
	{
		points=Torques->m; //Trayectory points
		//gets time from Torque matrix
		time=v_resize(time,points); 
		MEM_STAT_REG(time,TYPE_VEC);
		time=get_col(Torques,df,time);
		//resize Torque to delete time column
		m_resize(Torques,points,df);
	}
	
	Q_Qd_Qdd= m_resize(Q_Qd_Qdd,points,3*df+1);
	_set_col(Q_Qd_Qdd,0,time,0);
	v_resize_vars(df,&dQ_temp,&v_tmp,&inta_tmp,&intv_tmp, &a_tmp,&r_tmp, &dq_tmp, &q_tmp1, &q_tmp2, &q_tmp3,NULL);
	mem_stat_reg_vars(0,TYPE_VEC,&dQ_temp,&v_tmp,&inta_tmp,&intv_tmp,&a_tmp, &r_tmp, &dq_tmp,&q_tmp1, &q_tmp2, &q_tmp3,NULL);
	m_resize_vars(points,df,&q,&dq,&d2q,NULL);
	mem_stat_reg_vars(0,TYPE_MAT,&q,&dq,&d2q,NULL);

	//copying initial positions and velocities
	for(i=0;i<df;i++)
	{
		q->me[0][i]=Q0_Qd0->ve[i];
		dq->me[0][i]=Q0_Qd0->ve[i+df];
	}
	//Tor=m_resize(Tor,1,df);
	//MEM_STAT_REG(Tor,TYPE_MAT);
	//Iteration for each point each trajectory.
	for(i=0;i<points;i++)
	{
		if(i>0)
		{
			//euler velocity v(t+h)=v(t)+h*a(t)
			/*get_row(dq,i-1,v_tmp);
			get_row(d2q,i-1,a_tmp);
			sv_mlt(step,a_tmp,inta_tmp);
			v_add(v_tmp,inta_tmp,dq_tmp);
			_set_row(dq,i,dq_tmp,0);*/
			//velocity with position (differentiation.)
			// v(t+h)=1/h(r(t+h)-r(t))
			if(i==1) 
			{
				get_row(q,i,q_tmp1);
				get_row(q,i-1,q_tmp2);
				v_sub(q_tmp1,q_tmp2,r_tmp);
				sv_mlt(1/step,r_tmp,v_tmp);
				_set_row(dq,i,v_tmp,0);
			}
			else
			// v(t+h)=-1/(2*h)*(-3r(t+h)+4r(t)-r(t-h))   
			{
				get_row(q,i,q_tmp2);
				sv_mlt(-3,q_tmp2,q_tmp1);
				get_row(q,i-1,q_tmp3);
				sv_mlt(4,q_tmp3,q_tmp2);
				v_add(q_tmp1,q_tmp2,q_tmp3);
				get_row(q,i-2,q_tmp1);
				v_sub(q_tmp3,q_tmp1,q_tmp1);
				sv_mlt(-1/(2*step),q_tmp1,dq_tmp);
				_set_row(dq,i,dq_tmp,0);
			}
			for(j=0;j<df;j++) 
			{
				Q0_Qd0->ve[j]=q->me[i][j];
				Q0_Qd0->ve[j+df]=dq->me[i][j];
			}
		}
		//copying Vb dVb for each tray point.
		/*if(Vb!=NULL)
		{
			m_resize_vars(6,1,&Vbt,&dVbt,&Xt,NULL);
			mem_stat_reg_vars(0,TYPE_MAT,&Vbt,&dVbt,&Xt,NULL);
			m_move(Vb,0,i,6,1,Vbt,0,0);
			m_move(dVb,0,i,6,1,dVbt,0,0);
			m_move(X,0,i,6,1,Xt,0,0);
		}*/
		//obtain time step
		if (i<(points-1))
			step=time->ve[i+1]-time->ve[i];
		
		for (j=0;j<2*df;j++)
			Q0Qd0->me[j][0]=Q0_Qd0->ve[j];
		
		//Call fot torqfun.
		//if(i==0){
		//	printf("antes de feval");mem_info();}
		if(Torques==MNULL)
		{
			mem_stat_mark(3); 
			(*fun)(torqfun,Robot_type,Ttemp,time,time->ve[i],Q0Qd0,Qd_Q2d,DH,DYN,grav,fext,DYN_base);
			mem_stat_free(3);
		}
		else
		{
			mem_stat_mark(3); 
			(*fun)(torqfun,Robot_type,Torques,time,time->ve[i],Q0Qd0,Qd_Q2d,DH,DYN,grav,fext,DYN_base);
			mem_stat_free(3);
		}
		
		//if(i==0){
		//	printf("post  feval");mem_info(); getchar();}
		//extracting accel.
		m_transp(Qd_Q2d,Qd_Q2dt);
		m_move(Qd_Q2dt,0,df,1,df,d2q,i,0);
		//velocity verlet 
		/*if(i>0)
			_set_row(dq,i,v_add(dQ_temp,sv_mlt(0.5*step,get_row(d2q,i,NULL),NULL),NULL),0);*/

		//Position with verlet taylor expantion.
		if(i<points-1)
		{
			//Mid_term velocity
			//v_add(get_row(dq,i,NULL),sv_mlt(step/2,get_row(d2q,i,NULL),NULL),dQ_temp);

			//Verlet position  r(t+h)=r(t)+h*v(t)+h*h*a(t)/2;
			get_row(q,i,r_tmp);
			get_row(dq,i,v_tmp);
			sv_mlt(step,v_tmp,intv_tmp);
			get_row(d2q,i,a_tmp);
			sv_mlt(step*step/2,a_tmp,inta_tmp);
			v_add(inta_tmp,intv_tmp,q_tmp1);
			v_add(q_tmp1,r_tmp,q_tmp2);
			
			//comparing positions obtained with joints limits
			if(Robot_type==1)
				for(j=0;j<df;j++)
				{
					if((q_tmp2->ve[j])<(DH->me[j][5]))
						q_tmp2->ve[j]=DH->me[j][5];
					if((q_tmp2->ve[j])>(DH->me[j][6]))
						q_tmp2->ve[j]=DH->me[j][6];
				}
			if(Robot_type==2)
				for(j=0;j<df-6;j++)
				{
					if((q_tmp2->ve[j])<(DH->me[j][5]))
						q_tmp2->ve[j]=DH->me[j][5];
					if((q_tmp2->ve[j])>(DH->me[j][6]))
						q_tmp2->ve[j]=DH->me[j][6];
				}	
			_set_row(q,i+1,q_tmp2,0);
			//Euler position  r(t+h)=r(t)+h*v(t)
			/*get_row(q,i,r_tmp);
			get_row(dq,i,v_tmp);
			sv_mlt(step,v_tmp,intv_tmp);
			v_add(intv_tmp,r_tmp,q_tmp1);
			_set_row(q,i+1,q_tmp1,0);*/
		}
	}
	//returning results
	
	m_move(q,0,0,points,df,Q_Qd_Qdd,0,1);
	m_move(dq,0,0,points,df,Q_Qd_Qdd,0,df+1);
	m_move(d2q,0,0,points,df,Q_Qd_Qdd,0,2*df+1);
	
	return Q_Qd_Qdd;
}
