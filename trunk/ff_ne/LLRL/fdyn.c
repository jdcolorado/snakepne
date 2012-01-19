#include "LLRL.h"
/************************************************************************** 
Workspace used: 11,14,3,4,5,6,8,9
**************************************************************************	
FDYN Integrate forward dynamics

	[T Q QD] = FDYN(DH, T0, T1, Q0, QD0, Q_Qd_Qdd, Torques, opt, step, Vb, dVb)

	Integrates the forward dynamics of manipulator described by DH over the time 
	interval T0 to T1 and returns vectors of joint position and velocity in Q0_Qd0.
	DH is a robot object and describes the manipulator dynamics and 
	kinematics (in Denavit-Hartenberg convention), and Q is an n element vector of joint state.

	step is only used with verlet, in case of torques=NULL.
	
	solver= 1: ODE45 Torques=MNULL
	solver= 2: VERLET  Torques=any
	
	A control torque needs to be specified by a user specified function

	Copyright, Andr� Jaramillo Botero, 1999.
*************************************************************************/

/*07-3-2002, AJB: Added support for initial position and velocity conditions*/
/*07-3-2002, AJB: Eliminated DYN argument NECESITA CAMBIARSE RES_RK POR UN APUNTADOR*/
/*07-09-2002, AJB: Added static declarations and memory registration of vecs and mats */
/*09-02-2006, JCA: added support for n DF, flying base, and multiple solvers*/
/*22/02/2006, JCA, AM: Solved problem with memory in ode45 and verlet calls. to andt t1 are now "double"
	     Add spatial base position: X=[roll pitch yaw x y z]'*/
/*23/02/2006, JCA, added support for real time vector, still cant free mem stat of ode45...*/
/*30/03/2006, JCA, added suppor for returning modified base trayectory for ode 45. Fixes mem stat for ode45*/
/*16/04/2006, JCA added support for friction*/
/*26/04/2006, JCA Added sopport for robot type, and flying base as it should be*/
/*07/05/2006 JDC, Changed interface prototype, added DYN_base */
/*07/05/2006 JDC, Changed prototype call for feval,Torfun,ode45 added DYN_base*/
/*12/05/2006 JDC, Modification for Serpentine conditions */
/*31/05/2006 JCA, Modificarion por base initial position and velocity*/

MAT *fdyn(int robot_type,MAT *DH,MAT *DYN, MAT *DYN_base, double t0, double t1, VEC *Q0_Qd0, MAT *Q_Qd_Qdd,MAT *Torques,MAT *Fbase,int solver,double step, MAT *Xb0_Vb0,VEC *grav, VEC *fext,MAT *xb,MAT *vb,MAT *ab)
{
	static	VEC *y0=VNULL, *tspan=VNULL;
	static 	MAT *x=MNULL, *v=MNULL, *Q_Qd=MNULL;
	double	tolerance,abs;
	int	m,i,j,points,df;
	
	m = DH->m;
	
	m_resize_vars(6,1,&x,&v,NULL);
	Q_Qd=m_resize(Q_Qd,1,2*m);
	mem_stat_reg_vars(0,TYPE_MAT,&x,&v,NULL);
	v_resize_vars(2,&y0,&tspan,NULL);
	mem_stat_reg_vars(0,TYPE_VEC,&y0,&tspan,NULL);
	
/********************************************************************************
//Validating robot type.
********************************************************************************/
	//Fixed - floating 
	if(robot_type==1 || robot_type==3)
	{
		if (Q0_Qd0 == VNULL)
			y0=v_resize(y0,2*m);
		else
		{
			if (Q0_Qd0->dim != (u_int) 2*m)
				error(E_NULL,"(wrong initial joint state) fdyn\n");
			else
				y0=v_copy(Q0_Qd0,y0);
		}
	}
	
	//Obtain spatial positions and velocities from basetray.c (Serpentine Robots)
	if(robot_type==4)
	{	
		/*printf("\nantes de base_x");
		mem_info();*/
		for(i=0;i<2*m;i++)
			Q_Qd->me[0][i]=Q0_Qd0->ve[i];
			
		mem_stat_mark(4);
		x=base_x(DH,DYN,DYN_base,Q_Qd,x);
		mem_stat_free(4);
		
		mem_stat_mark(4);
		v=base_v(DH,DYN,DYN_base,Q_Qd,x,v);
		mem_stat_free(4);
		/*printf("\ndespues de base_x");
		mem_info();	*/
		
;
		for(i=0;i<6;i++)
		{
			Xb0_Vb0->me[i][0]=Xb0_Vb0->me[i][0]-x->me[i][0];
			Xb0_Vb0->me[i][1]=Xb0_Vb0->me[i][1]-v->me[i][0];
		}
		
	}

	//Flying base - Serpentine robots 
	if(robot_type==2 || robot_type==4)
	{
		if(robot_type==2)
		{
			if(Fbase->m!=6 ||(Fbase->n!=Torques->m))
				error(E_SIZES, "Wrong Fbase sizes in flying base\n ");
			if((Xb0_Vb0->m != 6) ||(Xb0_Vb0->n != 2))
				error(E_SIZES, "Wrong initial base position and velocity sizes in flying base\n ");
		}
	
		y0=v_resize(y0,2*m+12);
		for(i=0;i<m;i++)
		{
			y0->ve[i]=Q0_Qd0->ve[i];
			y0->ve[i+6+m]=Q0_Qd0->ve[i+m];
		}
		Torques=m_resize(Torques,Torques->m,Torques->n+6);
		
	//If friction coeficientes are zero, NULL spatial base trajectory inicialization. (Serpentine Robots)
		if((robot_type==4 && DYN->me[0][10]==0 && DYN->me[0][11]==0)||(robot_type==2))
			m_zero(Xb0_Vb0);
		
		for(i=0;i<6;i++)
		{
			y0->ve[m+i]=Xb0_Vb0->me[i][0];
			y0->ve[6+2*m+i]=Xb0_Vb0->me[i][1];
			
			if(robot_type==4)
				Fbase=m_resize(Fbase,6,Torques->m);
			
			for(j=0;j<Fbase->n;j++)
			{
				if(i==0)
					Torques->me[j][m+6]=Torques->me[j][m];
				Torques->me[j][i+m]=Fbase->me[i][j];
			}
		}
	}


	/*validate inputs (sizes and initializations) */
	if (DH == MNULL) 
		error(E_NULL, "fdyn (inputs)\n");

	tspan = v_resize(tspan,2);
	tspan->ve[0]=t0;
	tspan->ve[1]=t1;
	
/*******************************************************************************/

/********************************************************************************
//ODE45.
********************************************************************************/

	if(solver==1)
	{
		tolerance=1e-3;
		abs=1e-4;
				
		/*printf("\nantes de ode45?");
		mem_info();*/
		mem_stat_mark(4);
		Q_Qd_Qdd=ode45(feval, robot_type, Torques, tspan, y0, tolerance, abs, 4, DH, DYN, Q_Qd_Qdd, grav, fext, DYN_base);
		mem_stat_free(4);
		/*printf("\ndespues de ode45?");
		mem_info();*/
	}

/********************************************************************************
//Verlet
********************************************************************************/

	if(solver==2)
	{
		mem_stat_mark(4);
		Q_Qd_Qdd=verlet(feval, robot_type, DH, DYN, tspan, y0, Torques, Q_Qd_Qdd, step, grav, fext,DYN_base);
		mem_stat_free(4);
	}


/********************************************************************************/

	if(robot_type==3)
	{
		if(solver==1)
			Q_Qd_Qdd=m_resize(Q_Qd_Qdd,Q_Qd_Qdd->m,Q_Qd_Qdd->n+(Q_Qd_Qdd->n-1)/2);
	
		mem_stat_mark(4);
		Q_Qd_Qdd=basetray(robot_type,DH,DYN,DYN_base,Q_Qd_Qdd);
		mem_stat_free(4);
	}
	
	if(robot_type!=1)
	{
		if(solver==1 && (robot_type==2 || robot_type==4))
		{
			df=(Q_Qd_Qdd->n-13)/2;
			m_resize_vars(6,Q_Qd_Qdd->m,&xb,&vb,&ab,NULL);
		}
		else
		{
			df=(Q_Qd_Qdd->n-19)/3;
			m_resize_vars(6,Q_Qd_Qdd->m,&xb,&vb,&ab,NULL);
		}	
	}
	else
	{
		if(solver==1)
			df=(Q_Qd_Qdd->n-1)/2;
		else
			df=(Q_Qd_Qdd->n-1)/3;
	}

//Flying base - Serpentine robots
	if(robot_type==2 || robot_type==4)
	{
/**********************************************************
Extracting base trajectory
**********************************************************/	
		for(i=0;i<Q_Qd_Qdd->m;i++)
		{
			xb->me[0][i]=Q_Qd_Qdd->me[i][df+1];
			xb->me[1][i]=Q_Qd_Qdd->me[i][df+1+1];
			xb->me[2][i]=Q_Qd_Qdd->me[i][df+1+2];
			xb->me[3][i]=Q_Qd_Qdd->me[i][df+1+3];
			xb->me[4][i]=Q_Qd_Qdd->me[i][df+1+4];
			xb->me[5][i]=Q_Qd_Qdd->me[i][df+1+5];
			
			vb->me[0][i]=Q_Qd_Qdd->me[i][2*df+1+6];
			vb->me[1][i]=Q_Qd_Qdd->me[i][2*df+1+1+6];
			vb->me[2][i]=Q_Qd_Qdd->me[i][2*df+1+2+6];
			vb->me[3][i]=Q_Qd_Qdd->me[i][2*df+1+3+6];
			vb->me[4][i]=Q_Qd_Qdd->me[i][2*df+1+4+6];
			vb->me[5][i]=Q_Qd_Qdd->me[i][2*df+1+5+6];
			if(solver==2)
			{
				ab->me[0][i]=Q_Qd_Qdd->me[i][3*df+1+12];
				ab->me[1][i]=Q_Qd_Qdd->me[i][3*df+1+1+12];
				ab->me[2][i]=Q_Qd_Qdd->me[i][3*df+1+2+12];
				ab->me[3][i]=Q_Qd_Qdd->me[i][3*df+1+3+12];
				ab->me[4][i]=Q_Qd_Qdd->me[i][3*df+1+4+12];
				ab->me[5][i]=Q_Qd_Qdd->me[i][3*df+1+5+12];
			}
		}
		m_move(Q_Qd_Qdd,0,(Q_Qd_Qdd->n-13)/2+7,Q_Qd_Qdd->m,(Q_Qd_Qdd->n-13)/2,Q_Qd_Qdd,0,(Q_Qd_Qdd->n-13)/2+1);
		if(solver==2)
			Q_Qd_Qdd=m_resize(Q_Qd_Qdd,Q_Qd_Qdd->m,Q_Qd_Qdd->n-18);
		else
			Q_Qd_Qdd=m_resize(Q_Qd_Qdd,Q_Qd_Qdd->m,Q_Qd_Qdd->n-12);
		
		if(robot_type==4)
		{	
			
			mem_stat_mark(4);
			Q_Qd_Qdd=basetray(robot_type,DH,DYN,DYN_base,Q_Qd_Qdd);
			mem_stat_free(4);
			
			if(solver==1)
				j=2;
			for(i=1;i<Q_Qd_Qdd->m-1;i++)
			{
				xb->me[0][i]=xb->me[0][i]+Q_Qd_Qdd->me[i-1][j*df+1];
				xb->me[1][i]=xb->me[1][i]+Q_Qd_Qdd->me[i-1][j*df+2];
				xb->me[2][i]=xb->me[2][i]+Q_Qd_Qdd->me[i-1][j*df+3];
				xb->me[3][i]=xb->me[3][i]+Q_Qd_Qdd->me[i-1][j*df+4];
				xb->me[4][i]=xb->me[4][i]+Q_Qd_Qdd->me[i-1][j*df+5];
				xb->me[5][i]=xb->me[5][i]+Q_Qd_Qdd->me[i-1][j*df+6];
				
				vb->me[0][i]=vb->me[0][i]+Q_Qd_Qdd->me[i-1][j*df+7];
				vb->me[1][i]=vb->me[1][i]+Q_Qd_Qdd->me[i-1][j*df+8];
				vb->me[2][i]=vb->me[2][i]+Q_Qd_Qdd->me[i-1][j*df+9];
				vb->me[3][i]=vb->me[3][i]+Q_Qd_Qdd->me[i-1][j*df+10];
				vb->me[4][i]=vb->me[4][i]+Q_Qd_Qdd->me[i-1][j*df+11];
				vb->me[5][i]=vb->me[5][i]+Q_Qd_Qdd->me[i-1][j*df+12];
								
			
				if(solver==2)
				{
					ab->me[0][i]=ab->me[0][i]+Q_Qd_Qdd->me[i-1][j*df+13];
					ab->me[1][i]=ab->me[1][i]+Q_Qd_Qdd->me[i-1][j*df+14];
					ab->me[2][i]=ab->me[2][i]+Q_Qd_Qdd->me[i-1][j*df+15];
					ab->me[3][i]=ab->me[3][i]+Q_Qd_Qdd->me[i-1][j*df+16];
					ab->me[4][i]=ab->me[4][i]+Q_Qd_Qdd->me[i-1][j*df+17];
					ab->me[5][i]=ab->me[5][i]+Q_Qd_Qdd->me[i-1][j*df+18];
				}
			}
			m_move(Q_Qd_Qdd,2,0,Q_Qd_Qdd->m-2,Q_Qd_Qdd->n,Q_Qd_Qdd,0,0);
			Q_Qd_Qdd=m_resize(Q_Qd_Qdd,Q_Qd_Qdd->m-4,Q_Qd_Qdd->n-18);
		}
	}
	if(robot_type==3)
	{
/**********************************************************
Extracting base trajectory
**********************************************************/	
		for(i=0;i<Q_Qd_Qdd->m;i++)
		{
			xb->me[0][i]=Q_Qd_Qdd->me[i][3*df+1];
			xb->me[1][i]=Q_Qd_Qdd->me[i][3*df+1+1];
			xb->me[2][i]=Q_Qd_Qdd->me[i][3*df+1+2];
			xb->me[3][i]=Q_Qd_Qdd->me[i][3*df+1+3];
			xb->me[4][i]=Q_Qd_Qdd->me[i][3*df+1+4];
			xb->me[5][i]=Q_Qd_Qdd->me[i][3*df+1+5];
			
			vb->me[0][i]=Q_Qd_Qdd->me[i][3*df+1+6];
			vb->me[1][i]=Q_Qd_Qdd->me[i][3*df+1+1+6];
			vb->me[2][i]=Q_Qd_Qdd->me[i][3*df+1+2+6];
			vb->me[3][i]=Q_Qd_Qdd->me[i][3*df+1+3+6];
			vb->me[4][i]=Q_Qd_Qdd->me[i][3*df+1+4+6];
			vb->me[5][i]=Q_Qd_Qdd->me[i][3*df+1+5+6];
			
			if(solver==2)
			{	
				ab->me[0][i]=Q_Qd_Qdd->me[i][3*df+1+12];
				ab->me[1][i]=Q_Qd_Qdd->me[i][3*df+1+1+12];
				ab->me[2][i]=Q_Qd_Qdd->me[i][3*df+1+2+12];
				ab->me[3][i]=Q_Qd_Qdd->me[i][3*df+1+3+12];
				ab->me[4][i]=Q_Qd_Qdd->me[i][3*df+1+4+12];
				ab->me[5][i]=Q_Qd_Qdd->me[i][3*df+1+5+12];
			}
		}
		Q_Qd_Qdd=m_resize(Q_Qd_Qdd,Q_Qd_Qdd->m,2*(Q_Qd_Qdd->n-19)/3+1);
	}
	if(solver==1)
	{
		for(i=0;i<2*df;i++)
			Q_Qd_Qdd->me[0][i+1]=Q0_Qd0->ve[i];
		
		if(robot_type!=1)
		{
			xb->me[0][0]=Xb0_Vb0->me[0][0];
			xb->me[1][0]=Xb0_Vb0->me[1][0];
			xb->me[2][0]=Xb0_Vb0->me[2][0];
			xb->me[3][0]=Xb0_Vb0->me[3][0];
			xb->me[4][0]=Xb0_Vb0->me[4][0];
			xb->me[5][0]=Xb0_Vb0->me[5][0];
		
			vb->me[0][0]=Xb0_Vb0->me[0][1];
			vb->me[1][0]=Xb0_Vb0->me[1][1];
			vb->me[2][0]=Xb0_Vb0->me[2][1];
			vb->me[3][0]=Xb0_Vb0->me[3][1];
			vb->me[4][0]=Xb0_Vb0->me[4][1];
			vb->me[5][0]=Xb0_Vb0->me[5][1];
		}
	}
	
	return Q_Qd_Qdd;	
}

