/*
 * LLRL.h --MAT *DH, MAT *DYN,
 *
 *		This header file describes the interface for the LLRL component.
 *
 * Date : 2003-04-08
 *
 * Requirements: stdio.h math.h ../Meschach/matrix.h ../Meschach/matrix2.h
 *
 * Copyright  Robotics and Automation Group, Pontificia Universidad Javeriana - Cali.
 *		            Andr� Jaramillo Botero, ajaramil@puj.edu.co
 *
 * See the file "license.terms" for information on usage and redistribution of this file, and for a
 * DISCLAIMER OF ALL WARRANTIES.
 *
 * SCCS: %Z%  %M%  %I%  %E%  %U%
 */
/* ---------------------------------------------------------------------------------------------------------------------- */
#if !defined( __LLRL_h )
/* ----------------------------------------------------------------------------------------------------------------------*/
#define __LLRL_h
/* ----------------------------------------------------------------------------------------------------------------------*/
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "../Meschach/matrix.h"
#include "../Meschach/matrix2.h"
#include "../Meschach/matlab.h"
//#include "../CLinkage/MeschachLkg.hh"
/* ----------------------------------------------------------------------------------------------------------------------*/
#define LLRL_OK                                            1
#define LLRL_ERROR                                      0
/* --------------------------------------------------------------------------------------------------------------------*/

/*************************************************************************
* Inverse dynamic functions
**************************************************************************/
void homogeneus(float,float, float, float, MAT *Rot);
void Basic_rotation(MAT *Rot, MAT *r);
void form_matrix(int,int, MAT *M, MAT *mtempo);
void skew_symetric( MAT *M, VEC *v_temp);


/*****************************************************************************
 Compute inverse dynamics via recursive Newton-Euler formulation for
 a free- Fixed/flying/floating base systems.
 
 TAU = ff_ne(robot_type,DH,DYN,Q_Qd_Qdd,grav,fext, X, V, dV,Fbase)
 
 robot_type= 1)Fixed base 2)Flying base 3)Floating base 4)Serpentine Robots
 DH=Denavit Hartenberg parameter matrix for manipulator
 DYN=dynamic parameter matrix in the form [m,sx,sy,sz,Ixx,Iyy,Izz,Ixy,Ixz,Iyz,Jm,G,B,Tc+,Tc-]
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
/*02/05/2006 JDC, Modification for fixed base systems computing gravity*/
/*07/05/2006 JDC, Changed interface prototype, added DYN_base */
/*07/05/2006 JDC, Added inertia of the base for the total mass of the system*/
/*07/05/2006 JDC, Added base's surface friction for the total friction of the systems (Serpentine)*/
/*30/05/2006 JDC, Total mass of the system modification*/

MAT *ff_ne(int robot_type,MAT *DH,MAT *DYN,MAT *DYN_base,MAT *Q_Qd_Qdd,VEC *GRAV,VEC *fext,MAT *Torques,MAT *X, MAT *Vb, MAT *dVb,MAT *Fbase);

/*************************************************************************
* Quaternions Functions.
**************************************************************************/
/*************************************************************************	
Conversion from euler angles to quaternion 

	Q=eul2q(phi,theta,psi,Q)
  
	phi,theta,psi=scalars defining euler rotation angles (in radians)
	Q= quaternion/vector [q1,q2,q3,q4] representing basic rotation.
	
	Copyright, Andr� Jaramillo Botero, 1999. 
***************************************************************************/
VEC	*eul2q(double phi,double theta,double psi,VEC *Q);

/*************************************************************************
* Homogeneous Transformations.
**************************************************************************/
/*************************************************************************	
	Conversion from homogeneous transform (rotation only) to quaternion

	Q=tr2q(T,Q)
 
	T=4x4 homogeneous transformation matrix representing rotation.
	Q=vector defining quaternion parameters [q1,q2,q3,q4] 
	
	Copyright, Andr� Jaramillo Botero, 1999. 
***************************************************************************/
VEC *tr2q(MAT *T,VEC *Q);

/*************************************************************************
    * Homogenous Transformations.
**************************************************************************/
/*************************************************************************	
    Basic homogeneous rotation matrix about an arbitrary vector 

    ROT=rotvec[v,angle,ROT]
  
    v=arbitrary vector
    angle=scalar defining rotation angle (in radians) about v
    ROT=4x4 homogeneous transformation matrix representing basic rotation.
	
    Copyright, Andr� Jaramillo Botero, 1999. 
  ***************************************************************************/
MAT *rotvec(VEC *v,double angle,MAT *ROT);

/*************************************************************************
* Kinematics Functions.
**************************************************************************/
/*************************************************************************	
	Compute manipulator Jacobian in end-effector frame for the current 
	pose q.
	
	J=jacobn(DH,q,J)
	
	The manipulator Jacobian matrix maps differential changes in joint space
	to differential Cartesian motion of the end-effector.
	
	  dX = J dQ
	
	Uses the technique of Paul, Shimano, Mayer 
	Differential Kinematic Control Equations for Simple Manipulators
	IEEE SMC 11(6) 1981 pp. 456-460
	
	For an n-axis manipulator the Jacobian is a 6 x n matrix.

	Copyright, Andr� Jaramillo Botero, 1999. 
***************************************************************************/
int jacobn( MAT *DH, MAT *DYN, VEC *q, MAT *J );
/*************************************************************************	
	Compute manipulator Jacobian in end-effector frame for the current 
	pose q.
	
	J=jacobn_float(DH,q,J)
	
	The manipulator Jacobian matrix maps differential changes in joint space
	to differential Cartesian motion of the end-effector.
	
	  dX = J dQ
	
	Uses the technique of Paul, Shimano, Mayer 
	Differential Kinematic Control Equations for Simple Manipulators
	IEEE SMC 11(6) 1981 pp. 456-460
	
	For an n-axis manipulator the Jacobian is a 6 x n matrix.

	Copyright, Andr� Jaramillo Botero, 1999. 
***************************************************************************/
/*20/04/2006 It is the same jacobn but the last frame is referred at mass center.*/

int jacobn_float( MAT *DH, MAT *DYN, VEC *q, MAT *J );

/*************************************************************************	
	Compute manipulator Jacobian in world coordinates for the current 
	pose Q.
	
	J=jacob_base(DH,Q,J)
	
	The manipulator Jacobian matrix maps differential changes in joint space
	to differential Cartesian motion (world coord frame) of the end-effector.
	
	  dX = J dQ
	
	For an n-axis manipulator the Jacobian is a 6 x n matrix.

	Copyright, Andr� Jaramillo Botero, 1999. 
***************************************************************************/
/*20/04/2006 JCA: changed jacob0 interface.*/
int jacob0( MAT *(*fun)(),MAT *DH ,MAT *DYN,VEC *q, MAT *J );

/************************************************************************************************
  DOF -> Number of degrees of freedom
  NPT -> Number of points along a trajectory
                       Rows x  Cols
  IN  :
    mDHParameters    : DOF  x  5        - Standard Denavit Hartenberg Parameters ( Alpha, a, Theta, d, Sigma )
    mJointParameters : NPT  x  DOF      - Joint Parameters along a trajectory
  OUT :
    mHMT             :   4  x  NPT * 4  - Homogeneous matrix transformation ( 4 x 4 )
    ErrorStatus ( LLRL_OK, LLRL_ERROR )
***************************************************************************************************/
int ForwardKinematics( MAT *mDHParameters, MAT *mJointParameters, MAT *mHMT );

/*********************************************************************************************************
  DOF -> Number of degrees of freedom
  NPT -> Number of points along a trajectory
                       Rows x  Cols
  IN  :
    mDHParameters    : DOF  x  5        - Standard Denavit Hartenberg Parameters ( Alpha, a, Theta, d, Sigma )
    mHMT             :   4  x  NPT * 4  - Homogeneous matrix transformation ( 4 x 4 )
  OUT :
   mJointParameters : NPT  x  DOF      - Joint Parameters along a trajectory
    ErrorStatus ( LLRL_OK, LLRL_ERROR )
***********************************************************************************************************/
int InverseKinematics( MAT *mDHParameters, MAT *mHMT, MAT *mJointParameters );

/*********************************************************************************************************	
	Conversion from homogeneous transform difference to 6D differential vector 

	d=ttr2diff(T1,T2,d)
  
	T1/2=4x4 homogeneous transformation matrix representing basic positions and 
	orientations at different points in space.
	d=6D vector defining 3 differential translations and 3 differential rotations 
	
	Copyright, Andr� Jaramillo Botero, 1999. 
*********************************************************************************************************/
VEC	*ttr2diff(MAT *T1,MAT *T2,VEC *d);

/*************************************************************************	
	Compute the link transform from kinematic parameters
	
	  T=linktran(DH(i),q,T) 
	  
	DH(i): DH parameter table row entry corresponding to link parameters
		DH(i)=[alpha, an, theta, dn, sigma)
		alpha is the link twist angle
		an is the link length
		theta is the link rotation angle
		dn is the link offset
		sigma is 0 for a revolute joint, non-zero for prismatic
	q=joint state
	T=is a homogeneous transformation between link coordinate frames.	

	Depending on sigma the value of q is interchanged with theta or dn
	accordingly.

	Copyright, Andr� Jaramillo Botero, 1999. 
***************************************************************************/
MAT* linktran( VEC *dh_link, double q, MAT *T );

/*************************************************************************	
	Forward robot kinematics for serial link manipulator. Computes the 
	forward kinematics for each joint space point defined by Q.  DH 
	describes the manipulator kinematics in standard Denavit Hartenberg 
	notation.
	
	  T=fkine(DH,Q,T) 
	  
	DH: DH parameter table row entry corresponding to link parameters.
		DH(i)=[alpha, an, theta, dn, sigma)
		alpha is the link twist angle
		an is the link length
		theta is the link rotation angle
		dn is the link offset
		sigma is 0 for a revolute joint, non-zero for prismatic

	For an n-axis manipulator Q is an 1xn matrix or an m x n matrix. 
	The elements are interpretted as joint angle or link length according to
	the form of DH or the j'th sigma value (0 for revolute, other for 
	prismatic).
	 
	If Q is a single row matrix it is interpreted as the generalized joint coordinates, 
	and fkine(DH,Q,T) returns in T a 4x4 homogeneous transformation for the 
	final link of the manipulator.
	
	If Q is a matrix, the rows are interpretted as the generalized joint 
	coordinates for a sequence of points along a trajectory.  Q(i,j) is
	the j'th joint parameter for the i'th trajectory point.  In this case
	fkine(DH,Q,T) returns in T an m x 16 matrix with each row containing a 
	'flattened' homogeneous transform corresponding to the input joint state.  
	A row can xttr(Tin,i,Tout).

	Copyright, Andr� Jaramillo Botero, 1999. 
***************************************************************************/
int fkine(MAT *DH,MAT *Q,MAT *TT);

/*************************************************************************	
	Compute manipulator Jacobian in end-effector frame for the current 
	pose q.
	
	J=jacob_end(DH,q,J)
	
	The manipulator Jacobian matrix maps differential changes in joint space
	to differential Cartesian motion of the end-effector.
	
	  dX = J dQ
	
	Uses the technique of Paul, Shimano, Mayer 
	Differential Kinematic Control Equations for Simple Manipulators
	IEEE SMC 11(6) 1981 pp. 456-460
	
	For an n-axis manipulator the Jacobian is a 6 x n matrix.

	Copyright, Andr� Jaramillo Botero, 1999. 
***************************************************************************/
MAT *jacob_end(MAT *DH,VEC *q,MAT *J);

/*************************************************************************	
	Compute manipulator Jacobian in world coordinates for the current 
	pose Q.
	
	J=jacob_base(DH,Q,J)
	
	The manipulator Jacobian matrix maps differential changes in joint space
	to differential Cartesian motion (world coord frame) of the end-effector.
	
	  dX = J dQ
	
	For an n-axis manipulator the Jacobian is a 6 x n matrix.

	Copyright, Andr� Jaramillo Botero, 1999. 
***************************************************************************/
/* 07-12-2002 AJB: changed sub_mat for move */

MAT *jacob_base(MAT *DH,VEC *q,MAT *J);


/*************************************************************************
* Vector and Matrices Functions
**************************************************************************/
/*************************************************************************	
	v_fwrite -- saves vector V wtih name "name" in a file fp

	Copyright, Andr� Jaramillo Botero, 1999. 
***************************************************************************/
VEC *v_fwrite(FILE *fp,VEC *V,char *name);

/*************************************************************************	
	m_fwrite -- saves matrix M with name "name" in a file fp

	Copyright, Andr� Jaramillo Botero, 1999. 
***************************************************************************/
MAT *m_fwrite(FILE *fp,MAT *M,char *name);


/*************************************************************************
 * Dynamics Functions
**************************************************************************/
/*************************************************************************	
	Adapted from matlab's 5.3 version of the PINV function
	Matrix Computations, Golub & van Loan, Third Edition

	Copyright, Andr� Jaramillo Botero, 1999. 
***************************************************************************/
MAT	*pinv(MAT *A,MAT *PIA);

/*************************************************************************	
	Return the manipulator link inertia tensor matrix.  Extract the 
	manipulator link inertia for "link" from DYN and return it in matrix 
	form.

	J=linkiner(DYN,link,J)

	J=inertia tensor for "link"
	link=link number
	DYN=dynamic parameter matrix for manipulator

	Copyright, Andr� Jaramillo Botero, 1999. 
***************************************************************************/
MAT	*linkiner(MAT *DYN,u_int link,MAT *J);

/*************************************************************************	
	Extract a homogeneous matrix from a matrix of homogeneos contiguous 
	matrices. [T0 T1 ... Tn]

	Tout = xttr(TR, i, Tout) 
		
	TR=	matrix of homomgeneous transforms [T0 T1 ... Tn]
	i=	desired number of homogeneous matrix in TR.
	Tout=Extracted homogeneous transform.	
		  
	Copyright, Andr� Jaramillo Botero, 1999. 
***************************************************************************/
MAT	*xttr(MAT *TR,int n,MAT *Tout);

/*****************************************************************************
	Compute inverse dynamics via recursive Newton-Euler formulation
 
 	TAU = ne(DH,DYN,Q_Qd_Qdd,grav,fext)

	DH=Denavit Hartenberg parameter matrix for manipulator
	DYN=dynamic parameter matrix in the form [m,sx,sy,sz,Ixx,Iyy,Izz,Ixy,Ixz,Iyz]
	Q_Qd_Qdd=[Q QD QDD] (manipulator state)

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
 

	Copyright, Andr� Jaramillo Botero, 1999. 
******************************************************************************/
MAT *ne(MAT *DH,MAT *DYN,MAT *Q_Qd_Qdd,VEC *gravity,VEC *fext,MAT *Fout);

/*************************************************************************	
Compute the manipulator cartesian articulated body inertia matrix.

	M=cabi(DH,DYN,q,M) 

	DH=DH parameter matrix
	DYN=dynamic parameters matrix
    q= is an n element vector of joint state.
	
	for an n-axis manipulator returns the nxn symmetric inertia matrix 
	which relates Cartesian force/torque to Cartesian acceleration.
	  
	Copyright, Andr� Jaramillo Botero, 1999. 
***************************************************************************/

MAT *cabi(int Robot_type,MAT *DH,MAT *DYN,VEC *q,MAT *M) ;

/************************************************************************	
	Compute manipulator Jacobian in end-effector frame for the current 
	pose q.

	J=jacob_end(DH,q,J)

	The manipulator Jacobian matrix maps differential changes in joint space
	to differential Cartesian motion of the end-effector.
	
	  dX = J dQ
	
	Uses the technique of Paul, Shimano, Mayer 
	Differential Kinematic Control Equations for Simple Manipulators
	IEEE SMC 11(6) 1981 pp. 456-460
	
	For an n-axis manipulator the Jacobian is a 6 x n matrix.

	Copyright, Andr� Jaramillo Botero, 1999. 
***************************************************************************/
MAT *jacob_end(MAT *DH,VEC *q,MAT *J);

/*************************************************************************	
	Compute the manipulator inertia torque.  This corresponds to the acceleration
	dependent torques, without considering the velocity dependent torques.
	For an n-axis manipulator:
	returns the n element inertia torque vector (or a matrix each row being 
	the corresponding joint torques) at the specified pose and acceleration, 
	that is, ABI(Q)*QDD

	T = atorque(DH,DYN,Q,Qdd,T)
	
	DH=DH parameter matrix
	DYN=dynamic parameters matrix
	Q=is an npxn element matrix of joint position states.
	Qdd=is an npxn element matrix of joint acceleration states.

	for an n-axis manipulator returns the nxn symmetric inertia matrix 
	which relates Cartesian force/torque to Cartesian acceleration.

	Copyright, Andr� Jaramillo Botero, 1999. 

  07/02/2002: added validation of equal size for Q and Qdd entry arguments
***************************************************************************/
MAT *atorque(MAT *DH,MAT *DYN,MAT *Q,MAT *Qdd, MAT *M) ;

/*************************************************************************
Workspace used: 9
**************************************************************************	
	Compute the manipulator articulated body inertia matrix.

 		M = abi(DH,DYN,q,M) 

	DH=DH parameter matrix
	DYN=dynamic parameters matrix
	q= is an n element vector of joint state.

	for an n-axis manipulator returns the nxn symmetric articulated body 
	inertia matrix which relates joint torque to joint acceleration.
 
  	Copyright, Andr� Jaramillo Botero, 1999. 
*************************************************************************/
/*09-02-2006, JCA: added support for flying base inverse dynamics (ff_ne)*/
/*22/02/2006, JCA, AM: Solved problem with memory in ff_ne calls*/
/*22/02/2006 JCA  Add support for Spatial Base Position: X*/
/*23/03/2006 JCA  CHanged interface*/ 
/*30/03/2006 JCA added suport for fext, grav, energy changes in ff_ne*/
/*07/04/2006 JCA added sopport for new ff_ne  interface. and changed all null matrix for only one...*/
/*02/05/2006 JDC added support for floating-serpentine robots to obtain mass operator*/
/*07/05/2006, JDC Changed call for ff_ne function, passed MNULL for DYN_base */
/*12/05/2006* JDC, Modification for Serpentine conditions */
/*1/06/2006: JCA removed Fbase.*/

MAT *abi(int Robot_type, MAT *DH,MAT *DYN,VEC *q,MAT *M,MAT *X);

/*************************************************************************	
Compute the manipulator velocity dependent torques matrix.

 	C=vtorque(DH,DYN,Q,Qd,C) 
	
	DH=DH parameter matrix
	DYN=dynamic parameters matrix
	Q=is an npxn element matrix of joint position states.
	Qd=is an npxn element matrix of joint velocity states.
	C=matrix of velocity dependent torques
	  
	For an n-axis manipulator returns the n element velocity torques matrix 
	each row being the corresponding joint torques (at each trajectory point) 
	at the specified pose and velocity. 
		
	Copyright, Andr� Jaramillo Botero, 1999. 
***************************************************************************/
MAT *vtorques(MAT *DH,MAT *DYN,MAT *Q,MAT *Qd, MAT *M) ;

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

	Copyright, Andr� Jaramillo Botero, 1999. 
***************************************************************************/
MAT *friction(MAT *DH, MAT *DYN, MAT *Qd,MAT *Tau) ;

/*************************************************************************	
	Compute the gravity load on manipulator joints.
 
 	Tg=gravload(DH,DYN,Q,grav,Tg) 
	
	DH=DH parameter matrix
	DYN=dynamic parameters matrix
	Q=is an npxn element matrix of joint poisition states.
	Tg=matrix of gravity torques

	allows an arbitrary gravity vector to override the default of [0; 0; 9.81]	
	
	For an n-axis manipulator returns the n element friction torque matrix 
	each row being the corresponding joint torques (at each trajectory point) 
	at the specified pose and velocity required to compensate for gravity. 

	Copyright, Andr� Jaramillo Botero, 1999. 
***************************************************************************/
MAT *gravload(MAT *DH,MAT *DYN,MAT *Q,VEC *grav, MAT *Tg) ;

/******************************************************************************
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
/*07/05/2006, JDC Changed call for ff_ne function, passed MNULL for DYN_base */
/*07/05/2006, JDC added call to Kinetic and Potential energy computation functions */
/*1/06/2006: JCA removed Fbase.*/

MAT *accel(int Robot_type, MAT *DH, MAT *DYN,MAT *Q_Qd,MAT *Tor,MAT *Qdd,MAT *X, MAT *Vb, MAT *dVb,VEC *grav, VEC *fext);

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
/*12/05/2006* JDC, Modification for Serpentine conditions */

MAT *fdyn(int robot_type, MAT *DH, MAT *DYN, MAT *DYN_base, double t0, double t1, VEC *Q0_Qd0, MAT *Q_Qd_Qdd, MAT *Torques, MAT *Fbase, int solver, double step, MAT *Xb0_Vb0, VEC *grav, VEC *fext, MAT *outXb, MAT *outVb, MAT *outAb);

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
/*26/04/2006, JCA: changed interface, added Robot_type and deleted base input complete traj.*/
/*07/05/2006 JDC, Changed interface, added DYN_base*/
/*07/05/2006 JDC, Changed prototype call for feval, added DYN_base*/

MAT *verlet (MAT *(*fun)(), int Robot_type, MAT *DH, MAT *DYN, VEC *tspan, VEC *Q0_Qd0, MAT *Torques, MAT *Q_Qd_Qdd, double step, VEC *grav, VEC *fext,MAT *DYN_base);

/*************************************************************************
Workspace used: 11, 5, 6, 8, 9
**************************************************************************	
	The feval function is implemented in C by passing a function in the 
	function header:
	return type (* function)(argument types)
	Ex: MAT *(*fun)()
	Esta notaci� equivale a un apuntador a un apuntador de una funcion

	FEVAL Execute the specified function
	
	FEVAL(F,x1,...,xn) evaluates the function specified by a function
	handle or function name, F, at the given arguments, x1,...,xn.
	For example, if F = @foo, FEVAL(F,9.64) is the same as foo(9.64).
	
	Copyright, Andr� Jaramillo Botero, 1999. 
**************************************************************************/
/* 07-12-2002-AJB: changed assignment for block move in t_res_feval */
/* 07-12-2002-AJB: eliminated mem_stat_marks, added torque argument */
/* 09-23-2004-AJB: el llamado a la funci� abi ten� un argumento de m�
   para soportar la funci� de torque, se elimina temp_x por ahora */
/* 09-02-2006-JCA: fixed problem with 6 df, added suport for n DF and flying base
  where Vb and dVb are base velocity and acceleration. Note: Still works with fixed base*/
/*22/02/2006, JCA, AM: Solved problem with memory in accel calls*/
/*22/02/2006, JCA, JDC:Add spatial base position: X=[roll pitch yaw x y z]'*/
/*23/02/2006, JCA, deleted redundant DH and DYN matrix.*/
/*24/02/2006 JCA Solved problems with tor reg.*/
/*28/02/2006 JCA Added torque function support*/
/*16/04/2006, JCA added support for friction*/
/*26/04/2006, JCA Added sopport for robot type, and flying base as it should be*/
/*07/05/2006 JDC, Added some additional conditions*/
/*07/05/2006 JDC, Changed interface, added DYN_base*/
/*07/05/2006 JDC, Changed prototype call for AcelBase, added DYN_base*/
/*12/05/2006* JDC, Modification for Serpentine conditions */

MAT *feval(MAT *(*fun)(),int Robot_type,MAT *tau,VEC *time_vector,double time, MAT *x, MAT *res_feval, MAT *DH, MAT *DYN,VEC *grav, VEC *fext,MAT *DYN_base);

/**************************************************************************	
	NTRP45 Interpolation helper function for ODE45.
	YINTERP = NTRP45(TINTERP,T,Y,TNEW,YNEW,H,F) uses data computed in ODE45
	to approximate the solution at time TINTERP.
***************************************************************************/
/* 07-11-2002-AJB: correction in MEM_STAT_REG of MAT not VEC */

//funcion ntrp45 de matlab.  interpolacion para odb45
MAT *ntrp45(MAT *tinterp,double t,MAT *y,double h, MAT *f, MAT *resul);

/*************************************************************************
Workspace used: 11,14,3, 5, 6, 8, 9
**************************************************************************		
	ODE45  Solve non-stiff differential equations, medium order method.
    
	[T,Y] = ODE45(ODEFUN,TSPAN,Y0) 
	
	with TSPAN = [T0 TFINAL] integrates the system of differential equations 
	y' = f(t,y) from time T0 to TFINAL with initial conditions Y0. 
	Function ODEFUN(T,Y) must return a column vector
    corresponding to f(t,y). Each row in the solution array Y corresponds to
    a time returned in the column vector T. To obtain solutions at specific
    times T0,T1,...,TFINAL (all increasing or all decreasing), use 
    TSPAN = [T0 T1 ... TFINAL].     
		  
    ODE45 can solve problems M(t,y)*y' = f(t,y) with mass matrix M that is
    nonsingular. 

	Derived from Matlab's ode45 function
	Copyright, Marcial Qui�nez y Andr� Jaramillo Botero, 1999. 
***************************************************************************/

/* 07-5-2002, AJB: Eliminated function number -n_funcion- argument in ODE45 */
/* 07-09-2002, AJB: Added static declarations and memory management functions for feval calls */
/* 07-11-2003, AJB: Need to check fevals arguments ... int */
/* 09-02-2006, JCA: Added support for n df but the problem was with feval, not with ode45.*/
/*22-02-2006, JCA, JDC:Add spatial base position: X=[roll pitch yaw x y z]'*/
/*23-02-2006, JCA, mem stat reg and meme stat free still not working...*/
/*30/03/2006, JCA, added support for passing fext and grav and;
		added support for torques!=NULL  and flying base!!
		added support for mem_stat_mark and mem_stat_free in *fun calls.*/
/*16/04/2006, JCA added support for base force*/
/*16/04/2006 JCA Solved memory problems with maxm and maxm2*/
/*26/04/2006, JCA Added sopport for robot type, and flying base as it should be*/
/*07/05/2006 JDC, Changed interface, added DYN_base*/
/*07/05/2006 JDC, Changed prototype call for feval, added DYN_base*/

//funcion ode45 de matlab,  recibe un vector con tiempos,  uno con puntos de evaluacion, una tolerancia,  otra tolerancia
// un indice de refinamiento.

MAT	*ode45 (MAT *(*fun)(),int Robot_type,MAT *Torques,VEC *tspan,VEC *y0, double  rtol, double atol,int refine,MAT* input_dh, MAT *input_dyn, MAT *rK, VEC *grav, VEC *fext,MAT *DYN_base);

MAT	*PA10(MAT *pa10);
MAT	*PUMA560(MAT *p560);
MAT	*PUMA560_DH(MAT *p560_DH);

//devuelve el indice dentro de la matriz
double indice(MAT *mat,int indice);

/* divide todos los elementos de dos matrices entre si */
MAT *div_m(MAT *mat,MAT *mat2,MAT *out);

//valor absoluto de una matriz
MAT *absm(MAT *mat);

/********************************************************************************************************
    DOF -> Number of degrees of freedom
    NPT -> Number of points along a trajectory
    Rows x  Cols
    IN  :
      mDHParameters    : DOF  x  5    - Standard Denavit Hartenberg Parameters ( Alpha, a, Theta, d, Sigma )
      mHMT             :   4  x  NPT * 4  - Homogeneous matrix transformation ( 4 x 4 )
    OUT :
      mJointParameters : NPT  x  DOF      - Joint Parameters along a trajectory
      ErrorStatus ( LLRL_OK, LLRL_ERROR )
**********************************************************************************************************/
int ikine( MAT *mDHParameters, MAT *mHMT, MAT *mJointParameters );
//saca una matriz maxima entre los elementos de una matriz y un  double
MAT* maxm(MAT *a,double num,MAT *b);

//saca los elemntos maximos de cada posicion de dos matrices
MAT* maxm2(MAT *a,MAT *c,MAT *b);
//suma una matriz con un vector fila
MAT *sum_mat_vec(MAT *mat,VEC *vt,MAT *out);
/*************************************************************************	
	Conversion from homogeneous transform difference to 6D differential vector 

	d=ttr2diff(T1,T2,d)
  
	T1/2=4x4 homogeneous transformation matrix representing basic positions and 
	orientations at different points in space.
	d=6D vector defining 3 differential translations and 3 differential rotations 
	
	Copyright, Andr� Jaramillo Botero, 1999. 
***************************************************************************/
VEC	*ttr2diff(MAT *T1,MAT *T2,VEC *d);
/*************************************************************************	
    v_cross
	
    Copyright, Andr� Jaramillo Botero, 1999. 
***************************************************************************/
VEC* v_cross( VEC *v1, VEC *v2, VEC *v3 );

/*************************************************************************	
	Conversion from homogeneous transform to euler angles (y-convention)

	E=tr2rpy(T,E)
  
	T=4x4 homogeneous transformation matrix representing rotation.
	E=vector defining [roll,pitch,jaw] rotation angles (in radians)
	
	Copyright, Andr� Jaramillo Botero, 1999. 
***************************************************************************/
VEC *tr2rpy(MAT *T,VEC *E);

/*************************************************************************
Compute base position, velocity and acceleration for floating base

void basetray(x,v,a,DH,DYN,Q_Qd_Qdd

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
/*22/03/2006, JCA: Added support for joint limits*/
/*25/04/2006: JCA: added suport for direct velocities: Memory OK*/
/*12/05/2006  JDC, changed interface; added Robot type */

MAT *basetray(int Robot_type,MAT *DH,MAT *DYN,MAT *DYN_base,MAT *Q_Qd_Qdd);

/*
 *  serpenoid .c --
 *
 *	   This program compute inverse kinematics to obtain a Serpenoid Curve         
 *
 * Date of creation: 22/11/2005
 *
 *  
 * Copyright  Robotics and Automation Group, Pontificia Universidad Javeriana - Cali.
 *	          Andr� Jaramillo Botero, ajaramil@puj.edu.co
 *            	Juan Manuel Florez,      jmflorez@puj.edu.co  
 *            	Juli� David Colorado, jdcolorado@puj.edu.co
 *            	Juan Camilo Acosta, jcacosta@puj.edu.co 
 *
 * See the file "license.terms" for information on usage and redistribution of this file, and for a
 * DISCLAIMER OF ALL WARRANTIES.
 *
%DH: Denavit Hartenberg parameters matrix
%DYN: Manipulator dynamic parameters
%DYN_base: Manipulator base DYnamic parameters
%Q_Qd_Qdd: joint trayectory,time vector,x(base positions),v(base velocities),a(base accelerations) 
%n=# of bodies of the Serpentine
%a=ondulation degree of the Curve
%b=#periods per lenght unit
%c=Serpentine motion bias
%w=velocity propagation along all the bodies
%step= step trajectory time
%end_t=Final Simulation time
*/
// Fecha de creaci�: 22/11/2005
//16/03/2006 JCA: Added suport for base trayecory. basetray.c
//17/03/2006: Added support for base DYN parameters and corrected all identation problems.
//21/03/2006: removed x, v and a as parameters, from now included in q_qd_qdd
//		Changed time vector to first Q_Qd_Qdd column instead of last.
//21/03/2006: Changed DYN_base from VEC to MAT because the base may have more than 1 bodys
/*12/05/2006  JDC, changed interface basetray prototype called; added Robot type */
/*24/05/2006: JCA: removed base tray call. Serpenoid curve complete.*/
/*31/05/2006: JCA: Arrenged time vector*/
MAT *serpenoid( MAT *DH, MAT *DYN, MAT *DYN_base, MAT *Q_Qd_Qdd, Real n, Real a, Real b, Real c, Real w, Real step, Real end_t);

//**************************************************************************/
//Returns the Homogeneus Transformation Matrix for the base
//X: Spatial Position of the Base: [roll pitch yaw X  Y  Z]' 
//Obtain [r p]:Matrix Rotation, Position Vector
//original euler, vec*p and vec *columna instead of mat*p and mat*columna 
void euler(VEC *columna, MAT *r, VEC *p);

/***********************************************************************************************************
torqfun.c

		torqfun make a cubic spline of inTAU matrix in order to estimate a torque at any time.
		it uses spline_2 to calculate the second part of the spline for base trayectory.

Copyright 	Robotic and Automation Group, Pontificia Universidad Javeriana, Cali
	
		Andrés Jaramillo Botero
		Juan Camilo Acosta
		Julian David Colorado
		Antonio Matta.

%MAT *inTAU: input torques. Size points x df(Degrees of Freedom)
%VEC *time_vector: input time vector. Dim: points
%double time: estimation torque time input.
%MAT *outTAU: output matrix with torque evaluated at time. Size: 1 x df
**************************************************************************************************************/

/*28/03/2006: JCA: First version, doing straight line between 2 points.*/
/*29/03/2006: JCA: Improved with cubic splines.*/
/*02/05/2006: JCA: Removed base trayectory*/

MAT *torqfun(MAT *DH, MAT *DYN,MAT *inTAU,VEC *time_vector, double time, MAT *outTAU);

/*************************************************************************
Workspaces used: none
*************************************************************************
Compute base Acceleration  due to the Spatial Force of the base

void AcelBase(DH,DYN,Q_Qd_Qdd,Fric,Ab)

DH: Denavit Hartenberg parameters matrix
DYN: Manipulator dynamic parameters
Q_Qd_Qdd: joint trayectory,time,x(base positions),v(base velocities),a(base accelerations)
Fbase: Spatial base Force

Robotics and Automation Group
	Andres Jaramillo Botero
	Juan Camilo Acosta
	Julian David Colorado
	Antonio Matta

Pontificia Universidad Javeriana, Cali.

09/04/06   JDC: Compute base Acceleration  due to the Spatial Force of the base
26/04/06   JDC: Changed name of Fric for Fbase
26/04/06   JDC: Added robot type: 1)Fixed base 2)Flying base 3)Floating base 4)Serpentine Robots
02/05/06   JDC: None gravity in base acceleration
07/05/2006 JDC, Changed interface prototype, added DYN_base 
07/05/2006 JDC, Added inertia of the base for the total mass of the system
30/05/2006 JDC, Orienting MASS operator to mass center frame
***************************************************************************/
//Compute base Acceleration  due to the Spatial Force of the base
MAT *AcelBase(int robot_type,MAT *DH, MAT *DYN, MAT *DYN_base,VEC *grav, MAT *Q_Qd_Qdd,MAT* Fbase, MAT* Ab);

/*************************************************************************
Workspaces used: none
*************************************************************************
Compute the friction force with the surface referred at base frame

MAT* Friction_surface(MAT *DH, MAT *DYN,MAT *Q_Qd_Qdd,MAT* Vb)

DH: Denavit Hartenberg parameters matrix
DYN: Manipulator dynamic parameters
Q_Qd_Qdd: joint trayectory,time,x(base positions),v(base velocities),a(base accelerations)
Vb: Spatial base velocity

Robotics and Automation Group
	Andres Jaramillo Botero
	Juan Camilo Acosta
	Julian David Colorado
	Antonio Matta

Pontificia Universidad Javeriana, Cali.

29/04/06   JDC: Compute the friction force with the surface referred at base frame*/
/*07/05/2006 JDC, Changed interface prototype, added DYN_base */
/*07/05/2006 JDC, Added base's surface friction for the total friction of the systems (Serpentine)
/*19/05/2006 JDC, Added base velocity for initial conditions (Compute spatial vel)*/
/*19/05/2006 JDC, Elimination of base friction*/
/*30/05/2006 JDC, Interface Prototype has changed: Elimination of MAT Xb */
/*30/05/2006 JDC, New Call for base_x.c to obtain spatial base position Xb*/
/*30/05/2006 JDC, Orienting Friction to mass center frame*/
/*2/06/2006 JDC,  Interface has changed; added Xb*/
/****************************************************************************/
MAT *Friction_surface(MAT *DH, MAT *DYN,MAT *DYN_base,MAT *Q_Qd_Qdd,MAT *Xb, MAT *Vb, MAT *Fric);

/*************************************************************************
Compute Kinetic Energy of all the system

float Kinetic(VEC*dQ,MAT*M)

dQ: 1xn joint velocities, where n=#DOF
M: nxn Mass matrix of the system

Robotics and Automation Group
	Andres Jaramillo Botero
	Juan Camilo Acosta
	Julian David Colorado
	Antonio Matta

Pontificia Universidad Javeriana, Cali.

26/04/2006 JDC  Compute Kinetic Energy of all the system
***************************************************************************/
//Compute Kinetic Energy of all the system
float Kinetic(MAT *dQ,MAT *M);

/*************************************************************************
Compute Potential Energy of all the system

float Potential(MAT *DH, MAT *DYN, VEC *Q)

DH: Denavit Hartenberg parameters matrix
DYN: Manipulator dynamic parameters
Q: 1xn joint velocities, where n=#DOF

Robotics and Automation Group
	Andres Jaramillo Botero
	Juan Camilo Acosta
	Julian David Colorado
	Antonio Matta

Pontificia Universidad Javeriana, Cali.

26/04/2006 JDC  Compute Potential Energy of all the system
***************************************************************************/
//Compute Potential Energy of all the system
float Potential(MAT *DH, MAT *DYN, MAT* Q);

/*************************************************************************
Workspace used: 5
*************************************************************************
Compute base position respect to Mass center frame

void base_x(DH,DYN,Q,X)

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
MAT *base_x(MAT *DH,MAT *DYN,MAT *DYN_base,MAT *Q_Qd, MAT *x);

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

MAT *base_v(MAT *DH,MAT *DYN,MAT *DYN_base,MAT *Q_Qd, MAT *x, MAT *v);


#endif

/* ---------------------------------------------------------------------------------------------------------------------- */
