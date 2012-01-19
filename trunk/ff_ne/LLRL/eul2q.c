#include "LLRL.h"
/*************************************************************************	
Conversion from euler angles to quaternion 

	Q=eul2q(phi,theta,psi,Q)
  
	phi,theta,psi=scalars defining euler rotation angles (in radians)
	Q= quaternion/vector [q1,q2,q3,q4] representing basic rotation.
	
	Copyright, Andrés Jaramillo Botero, 1999. 
***************************************************************************/

VEC	*eul2q(double phi,double theta,double psi,VEC *Q)
{
  static double	ct,st;
	
  /* validate inputs (sizes and initializations) */
  if (Q == (VEC *)VNULL || Q->dim != 4) Q=v_resize(Q,4);
	
  ct=cos(theta/2);
  st=sin(theta/2);
	
  Q->ve[0]=cos((phi+psi)/2)*ct;
  Q->ve[1]=cos((phi-psi)/2)*st;
  Q->ve[2]=sin((phi-psi)/2)*st;
  Q->ve[3]=sin((phi+psi)/2)*ct;
	
  return	Q;
}

