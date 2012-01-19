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

#include "LLRL.h"

//Compute Kinetic Energy of all the system
float Kinetic(MAT *Q_Qd_Qdd,MAT *M)
{
	static MAT *dQM=MNULL,*dQ=MNULL;
	float K;
	int n=M->m,np;	        	/* number of joints:DOF  */
	np=(int) Q_Qd_Qdd->m;
	m_resize_vars(1,n,&dQ,NULL);
	mem_stat_reg_vars(0,TYPE_MAT,&dQ,NULL);
	m_move(Q_Qd_Qdd,0,n,np,n,dQ,0,0);
	dQM=m_resize(dQM,1,n);
	MEM_STAT_REG(dQM,TYPE_MAT);
	m_mlt(dQ,M,dQM);
	K=0.5*((dQM->me[0][0]*dQ->me[0][0])+(dQM->me[0][1]*dQ->me[0][1])+(dQM->me[0][2]*dQ->me[0][2])+
		(dQM->me[0][3]*dQ->me[0][3])+(dQM->me[0][4]*dQ->me[0][4])+(dQM->me[0][5]*dQ->me[0][5]));
	
	return K;
}
