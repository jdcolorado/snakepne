#include "LLRL.h"
/*************************************************************************
Workspace used: 7
*************************************************************************
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
int fkine(MAT *DH,MAT *Q,MAT *TT)
{
	static	VEC	*link=VNULL,*q=VNULL;
	static	MAT	*m1=MNULL,*m2=MNULL,*T=MNULL;
	int		i,j;
	
	/* validate inputs (sizes and initialization) */
	if (DH == (MAT *)MNULL || Q == (MAT *)MNULL)
        return( LLRL_ERROR );
	if (TT == (MAT *)MNULL || TT->m != 4) 
		TT=m_resize(TT,4,4); 
	
	m_zero(TT);
	/* setup temporals and allocate memory */
	v_resize_vars(DH->m,&link,&q,NULL);
	mem_stat_reg_vars(0,TYPE_VEC,&link,&q,NULL);
	m_resize_vars(4,4,&m1,&m2,&T,NULL);
	mem_stat_reg_vars(0,TYPE_MAT,&m1,&m2,&T,NULL);
	
	if ((int) Q->n == (int) DH->m) 
	{
		for (j=0;j<(int)Q->m;j++)
		{
			m_ident(T);
			get_row(Q,j,q);
			for (i=0;i<(int) DH->m;i++)
			{
				get_row(DH,i,link);
				mem_stat_mark(7);
				linktran(link,q->ve[i],m1);
				mem_stat_free(7);
				m_mlt(T,m1,m2);
				m_copy(m2,T);
			}
			if (j>0) TT=m_resize(TT,4,4*j);
			m_move(T,0,0,4,4,TT,0,4*j);
		}
	}
    else { return( LLRL_ERROR );}
	
    return( LLRL_OK );
}

