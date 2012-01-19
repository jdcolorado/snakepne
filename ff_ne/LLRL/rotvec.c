#include "LLRL.h"
/*************************************************************************	
    Basic homogeneous rotation matrix about an arbitrary vector 

    ROT=rotvec[v,angle,ROT]
  
    v=arbitrary vector
    angle=scalar defining rotation angle (in radians) about v
    ROT=4x4 homogeneous transformation matrix representing basic rotation.
	
    Copyright, Andrés Jaramillo Botero, 1999. 
  ***************************************************************************/

MAT *rotvec(VEC *v,double angle,MAT *ROT)
{
  static	double	ct,st,vt;
  static	VEC	*v1=VNULL,*v2=VNULL;
  static	MAT	*m1=MNULL,*m2=MNULL,*m3=MNULL;
	
  /* validate inputs (sizes and initializations) */
  if (v == (VEC *)VNULL) 
    error(E_NULL, "rotvec (input vector)\n");
  if (v->dim != 3) 
    error(E_SIZES, "rotvec (input vector)\n");
  if (ROT == (MAT *)MNULL || ROT->m != 4 || ROT->n != 4) 
    ROT=m_resize(ROT,4,4);
	
  /* allocate memory */
  v_resize_vars(3,&v1,&v2,NULL);
  mem_stat_reg_vars(0,TYPE_VEC,&v1,&v2,NULL);
  m_resize_vars(4,4,&m1,&m2,&m3,NULL);
  mem_stat_reg_vars(0,TYPE_MAT,&m1,&m2,&m3,NULL);
  /* end of memory allocation */
	
  m_ident(ROT);
  ct=cos(angle);
  st=sin(angle);
  vt=1-ct;
	
  ROT->me[0][0]=ct;
  ROT->me[0][1]=-v->ve[2]*st;
  ROT->me[0][2]=v->ve[1]*st;
  ROT->me[1][0]=v->ve[2]*st;
  ROT->me[1][1]=ct;
  ROT->me[1][2]=-v->ve[0]*st;
  ROT->me[2][0]=-v->ve[1]*st;
  ROT->me[2][1]=v->ve[0]*st;
  ROT->me[2][2]=ct;
	
  set_col(m1,0,v);
  m_transp(m1,m2);
  m3=m_mlt(m1,m2,m3);
  sm_mlt(vt,m3,m3);
  m_add(m3,ROT,ROT);
	
	/* ROT=[vx^2*vt+ct		vx*vy*vt-st*vz	vx*vz*vt+st*vy	0
  vx*vy*vt+st*vz	vy^2*vt+ct		vy*vz*vt-st*vx	0
  vx*vz*vt-st*vy	vy*vz*vt+st*vx	vz^2*vt+ct		0
  0				0				0				1]
    */
  V_FREE(v1);
  V_FREE(v2);
  M_FREE(m1);
  M_FREE(m2);
  M_FREE(m3);
  return ROT;
}
