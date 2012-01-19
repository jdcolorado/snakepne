/*---------------------------------------------------------------------------------------------------------------- 
 *  inv-dyn .c --
 *
 *	   Implementation of Inverse Dynamics for Fixed/flying/floating base systems           
 *
 * Date of creation: 2006-02-05
 *
 * Requirements: inv-dyn.h
 *
 * Copyright  Robotics and Automation Group, Pontificia Universidad Javeriana - Cali.
 *		Julian David Colorado, jdcolorado@puj.edu.co	      
 *		Andres Jaramillo Botero, ajaramil@puj.edu.co
 *              Antonio Alejandro Matta Gomez, amatta@puj.edu.co
 *    	        Juan Camilo Acosta, jcacosta@puj.edu.co 
 *
 * See the file "license.terms" for information on usage and redistribution of this file, and for a
 * DISCLAIMER OF ALL WARRANTIES.
 *
 * SCCS: %Z%  %M%  %I%  %E%  %U%
 */

#include "LLRL.h"

/*************************************************************************
* Fixed,Free Flying/floating base Functions.
**************************************************************************/
//Returns the Homogeneus Transformation Matrix
//alfa ,a teta,d are the DH parameters 
//Obtain Rot:Homogeneus Transformation Matrix
/*05/2006 JCA: fixed to avoid m_zero of ROT*/
void homogeneus(float alfa,float a, float teta, float d, MAT *Rot)
{
	Rot->me[0][0]=cos(teta);
	Rot->me[1][0]=sin(teta);
	Rot->me[2][0]=0;
	
	Rot->me[0][1]=-sin(teta)*cos(alfa);
	Rot->me[1][1]=cos(alfa)*cos(teta);
	Rot->me[2][1]=sin(alfa);
	
	Rot->me[0][2]=sin(alfa)*sin(teta);
	Rot->me[1][2]=-sin(alfa)*cos(teta);
	Rot->me[2][2]=cos(alfa);
	
	//Rot->me[0][3]=a*cos(teta);
	//Rot->me[1][3]=a*sin(teta);
	//Rot->me[2][3]=d;
}

//************************************************************
//Returns a basic matrix rotation
//Rot is a homogeneus transformation matrix
//r is a basic matrix rotation

void Basic_rotation(MAT *Rot, MAT *r)
{
	int fil,col;
	for(fil=0;fil<=2;fil++)
		for(col=0;col<=2;col++)
			r->me[fil][col]=Rot->me[fil][col];
}



//*************************************************************
//Obtain a block matrixsdfsd
//fili and coli are the punters to put values into a matrix mtempo from matrix M
//M is a temporal matrix 
//Function form_matrix returns block matrix mtempo 

void form_matrix(int fili,int coli, MAT *M, MAT *mtempo)
{
	int fil,col;
	for(fil=0;fil<=2;fil++)
		for(col=0;col<=2;col++)
			mtempo->me[fil+fili][col+coli]=M->me[fil][col];
}

//************************************************************
//Returns skew symetric matrix:
//Obtain values from v_temp and Obtain M skew symetric matrix
/*05/2006 JCA: fixed to avoid m_zero of M*/
void skew_symetric(MAT *M, VEC *v_temp)
{
	M->me[0][0]=0;
	M->me[0][1]=-v_temp->ve[2];
	M->me[0][2]=v_temp->ve[1];
	M->me[1][0]=v_temp->ve[2];
	M->me[1][1]=0;
	M->me[1][2]=-v_temp->ve[0];
	M->me[2][0]=-v_temp->ve[1];
	M->me[2][1]=v_temp->ve[0];
	M->me[2][2]=0;
}

 /*****************************************************************************
 Compute inverse dynamics via recursive Newton-Euler formulation for
 a free- Fixed/flying/floating base systems.
 
 TAU = ff_ne(robot_type,DH,DYN,DYN_base,Q_Qd_Qdd,grav,fext, X, V, dV,Fbase)
 
 robot_type= 1)Fixed base 2)Flying base 3)Floating base 4)Serpentine Robots
 DH=Denavit Hartenberg parameter matrix for manipulator
 DYN=dynamic parameter matrix in the form [m,sx,sy,sz,Ixx,Iyy,Izz,Ixy,Ixz,Iyz,Jm,G,B,Tc+,Tc-]
 DYN_base=dynamic parameter of the base robot: [m,sx,sy,sz,Ixx,Iyy,Izz,Ixy,Ixz,Iyz
 Q_Qd_Qdd=[Q QD QDD] (manipulator state)
 X = base spatial position
 V = Base spatial velocity.
 dV = Base spatial acceleration.
 Fbase = Spatial base force due to Base acceleration

 Returns the joint torque required to achieve the specified joint position,
 velocity and acceleration state.
 Gravity is assumed to be acting in the -Z direction with acceleration of
 9.81m/s/s, but this may be overriden by providing a gravity acceleration
 vector grav=[gx gy gz].

 An external force/moment acting on the end of the manipulator may also be
 specified by a 6-element vector fext=[Fx Fy Fz Mx My Mz].

 where	

 Q, QD and QDD are row vectors of the manipulator state; pos, vel, and accel.
The torque computed also contains a contribution due to armature inertia.

 Copyright, Robotics and Automation Group, 2006. 
 *		Juli� David Colorado, jdcolorado@puj.edu.co
 *	 	Andr� Jaramillo Botero, ajaramil@puj.edu.co
 *  	        Antonio Alejandro Matta G�ez, amatta@puj.edu.co
 *              Juan Camilo Acosta, jcacosta@puj.edu.co 
******************************************************************************/
/*22/02/2006 JDC, Added Euler Function to Calcule matrix Rotation r and position vector p from Inertial frame to first body*/
/*22/02/2006 JCA  Add X entrace parameter: Spatial Base Position, Solved problem with robots with prismatic and rotational art.*/
/*24/02/2006 JDC, JCA added support for kinetic and potential energy computation.  */
/*01/03/2006 JDC, added support for Motor friction */
/*23/03/2006 JCA  CHanged interface kinetic,potential,E_Conserv*/ 
/*27/03/2006 JDC, added simple friction Model to support Friction with contact Surface - Serpentine Robots */
/*28/03/2006 JDC, Modification of Power Kinetics*/
/*06/04/2006 JDC, Modification of Power Kinetics, Potential, Obtain mass operator system, Interface modification: Total friction force FRIC*/
/*26/04/2006 JDC, Added support for free Flying robots: Obtained Spatial base force due to Base acceleration: Fbase*/
/*26/04/2006 JDC, Changed name of Fric for Fbase*/
/*26/04/2006 JDC, Changed interface function: Eliminated Energies from inverse_dynamics*/
/*26/04/2006 JDC, Changed interface function: added robot_type*/
/*27/04/2006 JCA, Memory OK.*/
/*29/04/2006 JDC, Modification in the simple friction model for serpentine robots*/
/*02/05/2006 JDC, Modification for fixed base systems computing gravity*/
/*07/05/2006 JDC, Changed interface prototype, added DYN_base */
/*07/05/2006 JDC, Added inertia of the base for the total mass of the system*/
/*07/05/2006 JDC, Added base's surface friction for the total friction of the systems (Serpentine) */
/*20/05/2006 JCA, improved performance, deleted form matrix, added robot type support*/
/*15/06/2006 JDC, Added Friction Torque to Friction model (Serpentine)*/
/*13/07/2006: JCA deleted redundant matrix and vectors.*/

MAT *ff_ne(int robot_type,MAT *DH,MAT *DYN,MAT *DYN_base,MAT *Q_Qd_Qdd,VEC *GRAV,VEC *fext,MAT *Torques,MAT *X, MAT *Vb, MAT *dVb,MAT *Fbase)
{
	//serial variables
	static MAT *Q=MNULL, *dQ=MNULL, *d2Q=MNULL, *r=MNULL, *R=MNULL, *P=MNULL, *RP=MNULL, *PtRt=MNULL, *sk=MNULL, *Skw=MNULL, *S=MNULL, *Ticm=MNULL, *Icm=MNULL, *I=MNULL, *Itemp=MNULL, *identidad=MNULL, *rTicm=MNULL, *Tinertia=MNULL, *Ji=MNULL, *massU=MNULL, *SIcm=MNULL, *sksk=MNULL, *vf=MNULL, *RPIant=MNULL, *C=MNULL, *MC=MNULL, *Fr1=MNULL, *T1=MNULL, *T2=MNULL, *T3=MNULL, *V=MNULL, *dV=MNULL;
	
	static VEC *RPFext=VNULL, *F=VNULL, *p=VNULL, *H=VNULL, *ws=VNULL, *s=VNULL, *ss=VNULL, *up=VNULL, *down=VNULL, *vtempo=VNULL, *v=VNULL, *term=VNULL, *rws=VNULL, *rv1=VNULL, *qddH=VNULL, *qdH=VNULL, *temp1=VNULL, *Iv=VNULL, *sotro=VNULL, *IdV=VNULL, *vxy=VNULL, *Fr6=VNULL, *Frt=VNULL;
	
	int k,n,np,i,j,cont;
	double alfa,a,teta,d,massneg;
	n=(int) DH->m;	        	/* number of joints:DOF  */
	np=(int) Q_Qd_Qdd->m;		/* number of trajectory points */

	Torques = m_resize(Torques,np,n);
	
	/* validate inputs (sizes and initializations) */
	if (DH == (MAT *)MNULL || DYN == (MAT *)MNULL)
		error(E_NULL,"ne (DH and dynamic inputs)");
	if (DH->m != DYN->m)
		error(E_SIZES,"ne (DH and dynamic inputs)\n");
	if (Q_Qd_Qdd == (MAT *)MNULL)
		error(E_NULL,"Null input (state) in newton-euler");
	if (GRAV == (VEC *)VNULL)
		error(E_NULL,"Null input (state) in Gravity");
	else if(GRAV->dim !=6)
		error(E_SIZES,"Wrong Gravity Vector size");
	if (fext == (VEC *)VNULL)
		error(E_NULL,"Null input (state) in External Force");
	else if(fext->dim !=6)
		error(E_SIZES,"Wrong External Force Vector size");		
	
	/* validate size of Q_Qd_Qddd: position, veloc, accel */
	if ((int) Q_Qd_Qdd->n == 3*n && (int) Q_Qd_Qdd->m == np)
	{
		m_resize_vars(np,n,&Q,&dQ,&d2Q,NULL);
		mem_stat_reg_vars(0,TYPE_MAT,&Q,&dQ,&d2Q,NULL);
		m_move(Q_Qd_Qdd,0,0,np,n,Q,0,0);
		m_move(Q_Qd_Qdd,0,n,np,n,dQ,0,0);
		m_move(Q_Qd_Qdd,0,2*n,np,n,d2Q,0,0);
	}
	else
		error(E_SIZES,"ne (incorrect manipulator state matrix)");

	
/**********************************************************************
		MATRIX DECLARATION AND REGISTRATION 
**********************************************************************/
	
/**************** 6x6 **************/
	
	/*R: Rotation matrix 
	P: traslation matrix 
	RP: R*P
	PtRt: m_transp(P)*m_transp(R)
	RPIant: R*P*Iant
	Skw: block 6-dimensional Skew symetric matrix
	S:  6-dimensional matrix for s.
	Icm: Inertia matrix refered to mass center frame
	I: Inertia matrix refered to body's local frame
	Sicm: temp matrix for applying parallel axes theorem to obtain the inertial operator at frame i
	Tinertia: Inertia matrix for flying base*/
	m_resize_vars(6,6, &R, &P, &RP, &PtRt, &RPIant, &Skw, &S, &Icm, &I, &SIcm, &Tinertia, NULL);
	mem_stat_reg_vars(0,TYPE_MAT, &R, &P, &RP, &PtRt, &RPIant, &Skw, &S, &Icm, &I, &SIcm, &Tinertia, NULL);
	
/**************** 3x3 **************/
	
	//r: Rotation block matrix
	//sk: skew symetric temp matrix 
	//sksk: sk*sk matrix
	//Ticm: Inertial operator referred to center of gravity of body i
	//rTicm: r* Ticm
	//Itemp, 3*3 I temp matrix
	//ji: Inertial Tensor referred to center of gravity of body i
	//massU: mass * ident
	m_resize_vars(3,3, &r, &sk, &sksk, &Ticm, &rTicm, &Itemp, &identidad, &Ji, &massU, NULL);
	mem_stat_reg_vars(0,TYPE_MAT, &r, &sk, &sksk, &Ticm, &rTicm, &Itemp, &identidad, &Ji, &massU, NULL);

/**************** 6x(n+1) **************/
	
	/*V: spatial velocities
	dV: spatial accel.*/
	m_resize_vars(6,n+1,&V,&dV,NULL);
	mem_stat_reg_vars(0,TYPE_MAT,&V,&dV,NULL);
	
/**********************************************************************
	VECTOR DECLARATION AND REGISTRATION
**********************************************************************/	

/***************** 6 ****************/
	
	/*H: proyection vector in the axes of motion
	term Accelerations temp vector for sk2 *v term
	qddH: d2Q * H
	qdH:  dQ * H
	temp1:  temp vector used for force equation and base position
	IdV: I*dV
	Fr6: 6-dimensional friction force
	Frt: Projected force at frame i, 6-dimensional
	columna: Temp matrix
	RPFext: R*P*Fext*/		
	v_resize_vars(6, &H, &term, &qddH, &qdH, &temp1, &IdV, &Fr6, &Frt, &RPFext, &F, NULL);
	mem_stat_reg_vars(0,TYPE_VEC, &H, &term, &qddH, &qdH, &temp1, &IdV, &Fr6, &Frt, &RPFext, &F, NULL);	
	
/***************** 3 ****************/
	
	//WS, up, down: temp vectors
	//ss: r* s
	v_resize_vars(3, &ws, &s, &ss, &up, &down, &vtempo, &rws, &rv1, &Iv, &sotro, &v, &p, NULL);
	mem_stat_reg_vars(0,TYPE_VEC, &ws, &s, &ss, &up, &down, &vtempo, &rws, &rv1, &Iv, &sotro, &v, &p,NULL);

/**********************************************************************
	MATRIX AND VECTORS FOR SERPENTINE FRICTION
**********************************************************************/	
	
	/*vf: Tangencial and Normal Componentes of Friction Force
	C: Friction coef.
	MC: massneg*C
	Fr1: MC*vf*/
	m_resize_vars(2,2, &vf, &C, &MC, &Fr1,NULL);
	mem_stat_reg_vars(0,TYPE_MAT, &vf, &C, &MC, &Fr1,NULL);
	
	//vxy: Projecting velocity from frame i to center of gravity
	v_resize_vars(2, &vxy, NULL);
	mem_stat_reg_vars(0,TYPE_VEC, &vxy, NULL);
	
	//m_resize_vars(1,1,&DHtemp,&qtemp,NULL);
	//mem_stat_reg_vars(0,TYPE_MAT,&DHtemp,&qtemp,NULL);
	
	//m_resize_vars(6,1,&Xcm,NULL);
	//mem_stat_reg_vars(0,TYPE_MAT,&Xcm,NULL);
	
	/**************** 4x4 **************/
	
	/*T1-T2-T3: Temp matrix.*/
	m_resize_vars(4,4,&T1,&T2,&T3,NULL);
	mem_stat_reg_vars(0,TYPE_MAT,&T1,&T2,&T3,NULL);
	
	//m_resize_vars(4,4*n,&Ttotal,NULL);
	//mem_stat_reg_vars(0,TYPE_MAT,&Ttotal,NULL);
	
	/*if(robot_type==4)
	{
		//qtemp=m_resize(qtemp,1,n);
		//m_move(Q,k,0,1,n,qtemp,0,0);
		mem_stat_mark(15);
		Xcm=base_x(DH,DYN,DYN_base,Q,Xcm);
		mem_stat_free(15);
	}*/

/**********************************************************************
	END OF MATRIX AND VECTOR DECLARATION AND REGISTRATION
**********************************************************************/

/********************CALCULATE BASE FORCE Fbase=!NULL*******************************************/	
	if (Fbase!=MNULL)
	   	m_resize_vars(6,np,&Fbase,NULL);
	
	//Set Variables values 
	//z0->ve[2]=1;
	
	m_ident(identidad);
	m_zero(P);
	m_zero(R);
	m_zero(I);
	m_zero(S);
	//v_zero(F);
	//Compute for each point of trajectory
	for(k=0;k<np;k++)
	{
/**************************************************************************  
	Forward recurrence propagating kinematics parameters
**************************************************************************/
		//Reset values
		
		for(i=0;i<n;i++)
		{
			if(DH->me[i][4]==0)
			{	
				H->ve[0]=0;
				H->ve[1]=0;
				H->ve[2]=1;               //projection Vector of rotational motion
				H->ve[3]=0;
				H->ve[4]=0;
				H->ve[5]=0;
			}
			else
			{
				H->ve[0]=0;
				H->ve[1]=0;
				H->ve[2]=0;
				H->ve[3]=0;
				H->ve[4]=0;
				H->ve[5]=1;               //projection Vector of Prismatic motion
			}

			if(i==0)
			{
				//For a Fixed Base  
				if (robot_type!=2 )
				{
					//m_ident(r);
					m_ident(T1);
					v_zero(p);
				}
				//For a Flying/floating base
				else
				{
				//Obtain 6-Dimensional Position from X   
					get_col(X,k,temp1);
					//euler(temp1,r,p);	//Obtain Rotation r and Position p
					euler(temp1,T1,p);	//Obtain Rotation r and Position p
				}
			}
			else
			{
				alfa = DH->me[i-1][0];
				a=DH->me[i-1][1];
				if(DH->me[i-1][4]==0)        //Rotational
				{
					teta=Q->me[k][i-1];
					d=DH->me[i-1][3];
				}
				else                         //Prismatic
				{
					teta=DH->me[i-1][2];
					d=Q->me[k][i-1];
				}
				homogeneus(alfa,a,teta,d,T1);	//obtain homogeneus transformation
				//m_move(Rot,0,0,3,3,r,0,0);	//obtain rotation matrix
				T1->me[0][3]=a*cos(teta);
				T1->me[1][3]=a*sin(teta);
				T1->me[2][3]=d;
				T1->me[3][3]=1;
				
				p->ve[0]=-a;			//obtain the inverse of position p
				p->ve[1]=-d*sin(alfa);
				p->ve[2]=-d*cos(alfa);
			}

			//Obtain the 6-dimensional R operator based on r
			m_move(T1,0,0,3,3,R,0,0);
			m_move(T1,0,0,3,3,R,3,3);
			
			//Obtain the 6-dimensional P operator based on p
			m_move(identidad,0,0,3,3,P,0,0);
			m_move(identidad,0,0,3,3,P,3,3);
			skew_symetric(sk,p);
			sm_mlt(-1,sk,sk);
			m_move(sk,0,0,3,3,P,0,3);
			
/*********************************************************************************
	Computing Spatial Velocities
*********************************************************************************/
			//Obtain 6-dimensional velocity: Vel
			//sv_mlt(dQ->me[k][i],H,dQH);
			//m_transp(P,P_trans);
			//m_transp(R,R_trans);
			//m_mlt(P_trans,R_trans,PtRt);
			for(cont=0;cont<=5;cont++)
				for(j=0;j<=5;j++)	
					PtRt->me[cont][j]=P->me[0][cont]*R->me[j][0] 
							+P->me[1][cont]*R->me[j][1]
							+P->me[2][cont]*R->me[j][2]
							+P->me[3][cont]*R->me[j][3]
							+P->me[4][cont]*R->me[j][4]
							+P->me[5][cont]*R->me[j][5];
			
			//spatial velocities for each joint
			if ((robot_type==2 && i==0 && Vb!=NULL) )
			{
				/*for(cont=0;cont<=5;cont++)
					V->me[cont][i]=dQH->ve[cont]
					+R_trans->me[cont][0]*Vb->me[2][k]
					+R_trans->me[cont][1]*Vb->me[1][k]
					+R_trans->me[cont][2]*Vb->me[0][k]
					+R_trans->me[cont][3]*Vb->me[3][k]
					+R_trans->me[cont][4]*Vb->me[4][k]
					+R_trans->me[cont][5]*Vb->me[5][k];*/
				
				for(cont=0;cont<=5;cont++)
					V->me[cont][i]=H->ve[cont]*dQ->me[k][i]
							+R->me[0][cont]*Vb->me[2][k]
							+R->me[1][cont]*Vb->me[1][k]
							+R->me[2][cont]*Vb->me[0][k]
							+R->me[3][cont]*Vb->me[3][k]
							+R->me[4][cont]*Vb->me[4][k]
							+R->me[5][cont]*Vb->me[5][k];
			}
			else
			{
				if(i==0)
					for(cont=0;cont<=5;cont++)
						V->me[cont][0]=H->ve[cont]*dQ->me[k][i];
						
				else
					for(cont=0;cont<=5;cont++)
						V->me[cont][i]=H->ve[cont]*dQ->me[k][i]
						+PtRt->me[cont][0]*V->me[0][i-1]
						+PtRt->me[cont][1]*V->me[1][i-1]
						+PtRt->me[cont][2]*V->me[2][i-1]
						+PtRt->me[cont][3]*V->me[3][i-1]
						+PtRt->me[cont][4]*V->me[4][i-1]
						+PtRt->me[cont][5]*V->me[5][i-1];
			}
//*********************************************************************************
// Computing Spatial Accelerations
//*********************************************************************************
			
			// Obtain the last and actual angular velocity from Vi-1 and Vi
			//v->ve[0]=V->me[0][i];
			//v->ve[1]=V->me[1][i];
			//v->ve[2]=V->me[2][i];
			if(i==0)
			{
				if(robot_type==2 && Vb!=MNULL)
				{	
					ws->ve[0]=Vb->me[2][k];
					ws->ve[1]=Vb->me[1][k];
					ws->ve[2]=Vb->me[0][k];
				}
				else
				{
					ws->ve[0]=0;
					ws->ve[1]=0;
					ws->ve[2]=0;
				}
			}
			else
			{
				ws->ve[0]=V->me[0][i-1];
				ws->ve[1]=V->me[1][i-1];
				ws->ve[2]=V->me[2][i-1];
			}

			//skew_symetric(sk,V,i);	//obtain Skew symetric wi
			sk->me[0][0]=0;
			sk->me[0][1]=-V->me[2][i];
			sk->me[0][2]=V->me[1][i];
			sk->me[1][0]=V->me[2][i];
			sk->me[1][1]=0;
			sk->me[1][2]=-V->me[0][i];
			sk->me[2][0]=-V->me[1][i];
			sk->me[2][1]=V->me[0][i];
			sk->me[2][2]=0;
			
			//Obtain 6-dimensional Sk
			m_move(sk,0,0,3,3,Skw,0,0);
			m_move(sk,0,0,3,3,Skw,3,3);
			
			//m_transp(r,r_trans);
			//mv_mlt(r_trans,ws,rws);
			for(cont=0;cont<=2;cont++)
			{
				/*rws->ve[cont]=r->me[0][cont]*ws->ve[0]
						+r->me[1][cont]*ws->ve[1]
						+r->me[2][cont]*ws->ve[2];*/
				rws->ve[cont]=T1->me[0][cont]*ws->ve[0]
						+T1->me[1][cont]*ws->ve[1]
						+T1->me[2][cont]*ws->ve[2];
			}
			
			skew_symetric(sk,rws); 	//obtain Skew symetric wi-1

			// Obtain the last and actual velocity from Vi-1 and Vi
			v->ve[0]=V->me[3][i];
			v->ve[1]=V->me[4][i];
			v->ve[2]=V->me[5][i];
			if(i==0)
			{
				if(robot_type==2 && Vb!=MNULL)
				{	
					ws->ve[0]=Vb->me[3][k];
					ws->ve[1]=Vb->me[4][k];
					ws->ve[2]=Vb->me[5][k];
				}
				else
				{
					ws->ve[0]=0;
					ws->ve[1]=0;
					ws->ve[2]=0;
				}
			}
			else
			{
				ws->ve[0]=V->me[3][i-1];
				ws->ve[1]=V->me[4][i-1];
				ws->ve[2]=V->me[5][i-1];
			}
			
			//mv_mlt(r_trans,ws,rv1);
			for(cont=0;cont<=2;cont++)
			{
				/*rv1->ve[cont]=r->me[0][cont]*ws->ve[0]
							+r->me[1][cont]*ws->ve[1]
							+r->me[2][cont]*ws->ve[2];*/
				rv1->ve[cont]=T1->me[0][cont]*ws->ve[0]
						+T1->me[1][cont]*ws->ve[1]
						+T1->me[2][cont]*ws->ve[2];
			}
			v_sub(v,rv1,v);
			mv_mlt(sk,v,vtempo);

			for(cont=0;cont<=2;cont++)
				term->ve[cont+3]=vtempo->ve[cont];

			//Obtain 6-dimensional acceleration
			sv_mlt(d2Q->me[k][i],H,qddH);
			sv_mlt(dQ->me[k][i],H,qdH);
			
			if (robot_type==2  && i==0 )	//For Floating/Flying Base || robot_type==4)
			{	
				if(dVb!=NULL)
					for(cont=0;cont<=5;cont++)
						dV->me[cont][i]=term->ve[cont]+qddH->ve[cont]
						+R->me[0][cont]*(dVb->me[2][k]+GRAV->ve[0])
						+R->me[1][cont]*(dVb->me[1][k]+GRAV->ve[1])
						+R->me[2][cont]*(dVb->me[0][k]+GRAV->ve[2])
						+R->me[3][cont]*(dVb->me[3][k]+GRAV->ve[3])
						+R->me[4][cont]*(dVb->me[4][k]+GRAV->ve[4])
						+R->me[5][cont]*(dVb->me[5][k]+GRAV->ve[5])
						+Skw->me[cont][0]*qdH->ve[0]
						+Skw->me[cont][1]*qdH->ve[1]
						+Skw->me[cont][2]*qdH->ve[2]
						+Skw->me[cont][3]*qdH->ve[3]
						+Skw->me[cont][4]*qdH->ve[4]
						+Skw->me[cont][5]*qdH->ve[5];
				else
					for(cont=0;cont<=5;cont++)
						dV->me[cont][i]=term->ve[cont]+qddH->ve[cont]
						+R->me[0][cont]*(GRAV->ve[0])
						+R->me[1][cont]*(GRAV->ve[1])
						+R->me[2][cont]*(GRAV->ve[2])
						+R->me[3][cont]*(GRAV->ve[3])
						+R->me[4][cont]*(GRAV->ve[4])
						+R->me[5][cont]*(GRAV->ve[5])
						+Skw->me[cont][0]*qdH->ve[0]
						+Skw->me[cont][1]*qdH->ve[1]
						+Skw->me[cont][2]*qdH->ve[2]
						+Skw->me[cont][3]*qdH->ve[3]
						+Skw->me[cont][4]*qdH->ve[4]
						+Skw->me[cont][5]*qdH->ve[5];
			}
			else
			{
				if(i==0)
					for(cont=0;cont<=5;cont++)
						dV->me[cont][i]=term->ve[cont]+qddH->ve[cont]
						+PtRt->me[cont][0]*(GRAV->ve[0])
						+PtRt->me[cont][1]*(GRAV->ve[1])
						+PtRt->me[cont][2]*(GRAV->ve[2])
						+PtRt->me[cont][3]*(GRAV->ve[3])
						+PtRt->me[cont][4]*(GRAV->ve[4])
						+PtRt->me[cont][5]*(GRAV->ve[5])
						+Skw->me[cont][0]*qdH->ve[0]
						+Skw->me[cont][1]*qdH->ve[1]
						+Skw->me[cont][2]*qdH->ve[2]
						+Skw->me[cont][3]*qdH->ve[3]
						+Skw->me[cont][4]*qdH->ve[4]
						+Skw->me[cont][5]*qdH->ve[5];
				else
					for(cont=0;cont<=5;cont++)
						dV->me[cont][i]=term->ve[cont]+qddH->ve[cont]
						+PtRt->me[cont][0]*dV->me[0][i-1]
						+PtRt->me[cont][1]*dV->me[1][i-1]
						+PtRt->me[cont][2]*dV->me[2][i-1]
						+PtRt->me[cont][3]*dV->me[3][i-1]
						+PtRt->me[cont][4]*dV->me[4][i-1]
						+PtRt->me[cont][5]*dV->me[5][i-1]
						+Skw->me[cont][0]*qdH->ve[0]
						+Skw->me[cont][1]*qdH->ve[1]
						+Skw->me[cont][2]*qdH->ve[2]
						+Skw->me[cont][3]*qdH->ve[3]
						+Skw->me[cont][4]*qdH->ve[4]
						+Skw->me[cont][5]*qdH->ve[5];
			}   
			
		}
		//end of kinematics propagation
		/*printf("vel");
		m_output(V);
		printf("accel");
		m_output(dV);
		getchar();*/
	
		
//****************************************************************************************
//Backward recurrence: Dynamics Propagation
//****************************************************************************************
		

//****************************************************************************************
// Computing Spatial Forces
//****************************************************************************************
	
		for(i=(n-1);i>=0;i--)
		{
			alfa = DH->me[i][0];
			a=DH->me[i][1];
			if(DH->me[i][4]==0)        //Rotational
			{
				teta=Q->me[k][i];
				d=DH->me[i][3];
			}
			else                       //Prismatic
			{
				teta=DH->me[i][2];
				d=Q->me[k][i];
			}
			homogeneus(alfa,a,teta,d,r);	//obtain homogeneus transformation
			
			p->ve[0]=-a;
			p->ve[1]=-d*sin(alfa);
			p->ve[2]=-d*cos(alfa);

			//Obtain the 6-dimensional R operator based on r
			m_move(r,0,0,3,3,R,0,0);
			m_move(r,0,0,3,3,R,3,3);

			//Obtain the 6-dimensional P operator based on p
			skew_symetric(sk,p);
			sm_mlt(-1,sk,sk);
			m_move(sk,0,0,3,3,P,0,3);
			m_move(identidad,0,0,3,3,P,0,0);
			m_move(identidad,0,0,3,3,P,3,3);

			//Obtain vector s from frame i+1 to center of mass of body i
			s->ve[0]=DYN->me[i][1];
			s->ve[1]=DYN->me[i][2];
			s->ve[2]=DYN->me[i][3];
			
			v_sub(s,p,s);
			mv_mlt(r,s,ss);
			skew_symetric(sk,ss);    //Obtain Skew symetric of ss

			//obtain 6-dimensional operator S
			m_move(identidad,0,0,3,3,S,0,0);
			m_move(identidad,0,0,3,3,S,3,3);
			m_move(sk,0,0,3,3,S,0,3);
			
			//Obtain Skew symetric of Vi
			sk->me[0][0]=0;
			sk->me[0][1]=-V->me[2][i];
			sk->me[0][2]=V->me[1][i];
			sk->me[1][0]=V->me[2][i];
			sk->me[1][1]=0;
			sk->me[1][2]=-V->me[0][i];
			sk->me[2][0]=-V->me[1][i];
			sk->me[2][1]=V->me[0][i];
			sk->me[2][2]=0;
			
			//Obtain de Inertial Tensor referred to center of gravity of body i
			Ticm->me[0][0]=DYN->me[i][4];
			Ticm->me[0][1]=DYN->me[i][7];
			Ticm->me[0][2]=DYN->me[i][9];
			Ticm->me[1][0]=DYN->me[i][7];
			Ticm->me[1][1]=DYN->me[i][5];
			Ticm->me[1][2]=DYN->me[i][8];
			Ticm->me[2][0]=DYN->me[i][9];
			Ticm->me[2][1]=DYN->me[i][8];
			Ticm->me[2][2]=DYN->me[i][6];
		
			//Obtain 6-dimensional Inertial operator in center of gravity:
			m_mlt(r,Ticm,rTicm);
			
			//m_mlt(rTicm,r_trans,Ji);
			for(cont=0;cont<=2;cont++)
				for(j=0;j<=2;j++)	
					Ji->me[cont][j]=rTicm->me[cont][0]*r->me[j][0]
							+rTicm->me[cont][1]*r->me[j][1]
							+rTicm->me[cont][2]*r->me[j][2];
			
			m_move(Ji,0,0,3,3,Icm,0,0);
			sm_mlt(DYN->me[i][0],identidad,massU);
			m_move(massU,0,0,3,3,Icm,3,3);
			
			//applying parallel axes theorem to obtain the inertial operator at frame i
			
			m_mlt(S,Icm,SIcm);
			//m_mlt(SIcm,S_trans,I);
			for(cont=0;cont<=5;cont++)
				for(j=0;j<=5;j++)	
					I->me[cont][j]=SIcm->me[cont][0]*S->me[j][0]
							+SIcm->me[cont][1]*S->me[j][1]
							+SIcm->me[cont][2]*S->me[j][2]
							+SIcm->me[cont][3]*S->me[j][3]
							+SIcm->me[cont][4]*S->me[j][4]
							+SIcm->me[cont][5]*S->me[j][5];
			
//*******************************************************************************
// Computing Force equation:
//*******************************************************************************  
			m_mlt(R,P,RP);
			
			if(i==n-1)
				mv_mlt(RP,fext,RPFext);
			else
				mv_mlt(RP,F,RPFext);
			
			//get_col(dV,i,columna);     //obtain column from spatial acceleration
			
			for(cont=0;cont<=2;cont++)
				v->ve[cont]=V->me[cont][i];
			
			m_move(I,0,0,3,3,Itemp,0,0);
			mv_mlt(Itemp,v,Iv);
			mv_mlt(sk,Iv,up);
			m_mlt(sk,sk,sksk);
			mv_mlt(sksk,ss,sotro);
			sv_mlt(DYN->me[i][0],sotro,down);

			for(cont=0;cont<=2;cont++)
			{
				temp1->ve[cont]=up->ve[cont];
				temp1->ve[cont+3]=down->ve[cont];
			}

			//Obtain 6-dimensional force
			for(cont=0;cont<=5;cont++)
			{	
				IdV->ve[cont]=I->me[cont][0]*dV->me[0][i]
					+I->me[cont][1]*dV->me[1][i]
					+I->me[cont][2]*dV->me[2][i]
					+I->me[cont][3]*dV->me[3][i]
					+I->me[cont][4]*dV->me[4][i]
					+I->me[cont][5]*dV->me[5][i];
			}
			
			v_add(RPFext,IdV,RPFext);
			v_add(RPFext,temp1,F);

			

//*************************************************************************************
//Obtain Friction Force with Contact Surface (for Serpentine Robots)
//*************************************************************************************
			if(robot_type==4)
			{
				//Tangencial and Normal Componentes of Friction Force
				vf->me[0][0]=cos(X->me[k][i]);
				vf->me[0][1]=sin(X->me[k][i]);
				vf->me[1][0]=-sin(X->me[k][i]);
				vf->me[1][1]=cos(X->me[k][i]);
					
				//Matrix of Coefficients of friction
				C->me[0][0]=DYN->me[i][10];
				C->me[0][1]=0;
				C->me[1][0]=0;
				C->me[1][1]=DYN->me[i][11];
                        	//m_transp(vf,vf_trans);
				massneg=-DYN->me[i][0];
				//sm_mlt(massneg,vf_trans,massvf);
				//m_mlt(massvf,C,MC);
				for(cont=0;cont<=1;cont++)
					for(j=0;j<=1;j++)	
						MC->me[cont][j]=massneg*(vf->me[0][cont]*C->me[0][j]
								+vf->me[1][cont]*C->me[1][j]);
				m_mlt(MC,vf,Fr1);
				
				//Projecting velocity from frame i to center of gravity
				for(cont=3;cont<=4;cont++)
					vxy->ve[cont-3]=S->me[0][cont]*V->me[0][i]
							+S->me[1][cont]*V->me[1][i]
							+S->me[2][cont]*V->me[2][i]
							+S->me[3][cont]*V->me[3][i]
							+S->me[4][cont]*V->me[4][i]
							+S->me[5][cont]*V->me[5][i];
				
				
				//6-dimensional friction force
				Fr6->ve[0]=0;
				Fr6->ve[1]=0;
				Fr6->ve[2]=-(0.3333)*DYN->me[i][11]*DYN->me[i][0]*0.05*0.05*dQ->me[k][i];
				//Fr6->ve[3]=Fr2->ve[0];
				Fr6->ve[3]=Fr1->me[0][0]*vxy->ve[0]
					+Fr1->me[0][1]*vxy->ve[1];
				//Fr6->ve[4]=Fr2->ve[1];
				Fr6->ve[4]=Fr1->me[1][0]*vxy->ve[0]
					+Fr1->me[1][1]*vxy->ve[1];
				Fr6->ve[5]=0;
				//Projecting force at frame i, 6-dimensional
				mv_mlt(S,Fr6,Frt);
				//add Friction force component to local spatial Force
				v_add(F,Frt,F);
			}
			
			//v_copy(F,Ftant);                    //reload Ftant=F 
	
			//Projecting force in the axe of motion
			if(DH->me[i][4]==0)                 //Rotational
				Torques->me[k][i]=F->ve[2];
			else                               //Prismatic
				Torques->me[k][i]=F->ve[5];

			
		}
		//end dynamics propagation
		
	}
	return Torques;

}
//end Trajectory Points
