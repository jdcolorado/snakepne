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

#include "LLRL.h"

//Compute Potential Energy of all the system
float Potential(MAT *DH, MAT *DYN, MAT* Q_Qd_Qdd)
{
	static MAT *DHtemp=MNULL,*qtemp=MNULL,*T=MNULL,*Q=MNULL;
	float U;
	int n=(int) DH->m;	/* number of joints:DOF  */
	int cont,np; 
	
	np=(int) Q_Qd_Qdd->m;
	m_resize_vars(1,n,&Q,NULL);
	mem_stat_reg_vars(0,TYPE_MAT,&Q,NULL);
	m_move(Q_Qd_Qdd,0,0,np,n,Q,0,0);
	m_resize_vars(4,4,&T,NULL);
	mem_stat_reg_vars(0,TYPE_MAT,&T,NULL);    	
	
	U=0;
	for (cont=0;cont<n;cont++)
	{
		DHtemp=m_resize(DHtemp,cont,DH->n);
		qtemp=m_resize(qtemp,1,cont);
		mem_stat_reg_vars(0,TYPE_MAT,&DHtemp,&qtemp,NULL);
		
		m_move(DH,0,0,cont,DH->n,DHtemp,0,0);
		m_move(Q,0,0,1,cont,qtemp,0,0);

		fkine(DHtemp,qtemp,T);

		U=(DYN->me[cont][0]*9.81*T->me[2][3])+U;
	}		
	
	return U;
}
