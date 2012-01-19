#include "LLRL.h"
/*************************************************************************
Workspace used: 7
*************************************************************************	
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

int jacobn_float( MAT *DH, MAT *DYN, VEC *q, MAT *J )
{
	static  MAT *T = MNULL, *Ts = MNULL, *Tt = MNULL, *r=MNULL, *m1 = MNULL, *m2 = MNULL;
	static  VEC *link = VNULL, *s=VNULL, *p=VNULL;
	int i;
	// validate inputs (sizes and initialization)
	if( DH == ( MAT* ) MNULL || q == ( VEC* ) VNULL ) 
		return( LLRL_ERROR );
	
	if( q->dim != DH->m )
		return( LLRL_ERROR );
	
	if( J == ( MAT* ) MNULL || J->m != 6 || J->n != DH->m )
		J = m_resize( J, 6, DH->m );
	
	// allocate temporals
	m_resize_vars( 4, 4, &T, &Ts, &Tt, &m1, &m2, NULL );
	r=m_resize(r,3,3);
	mem_stat_reg_vars( 0, TYPE_MAT, &T, &Ts, &Tt, &r, &m1, &m2, NULL );
	p=v_resize(p,4);
	MEM_STAT_REG(p,TYPE_VEC);
	s=v_resize(s,6);
	MEM_STAT_REG(s,TYPE_VEC);
	v_zero(s);
	link = v_resize( link, DH->n );
	MEM_STAT_REG( link, TYPE_VEC );
	T = m_ident( T );
	for( i = ( int ) DH->m - 1; i >= 0; i-- )
	{
		get_row( DH, i, link );
		mem_stat_mark(7);
		m1 = linktran( link, q->ve[i], m1 );
		mem_stat_mark(7);
		
		m2 = m_mlt( m1, T, m2 );
		T = m_copy( m2, T );
		if(i==DH->m-1)
		{
			s->ve[3]=DYN->me[DH->m-1][1];
			s->ve[4]=DYN->me[DH->m-1][2];
			s->ve[5]=DYN->me[DH->m-1][3];
		
			//6x1 --> 4x4 (Euler -->homogeneus)
			euler(s,r,p);
			m_zero(Ts);
			m_move(r,0,0,3,3,Ts,0,0);
			p->ve[3]=1;
			_set_col(Ts,3,p,0);
			Ts->me[0][0];
			m_mlt(T,Ts,Tt);
			m_copy(Tt,T);
		}
	
		
		if( DH->me[ i ][ 4 ] == 0 )
		{
			// rotational joint
			J->me[ 0 ][ i ] = -T->me[ 0 ][ 0 ] * T->me[ 1 ][ 3 ] + T->me[ 1 ][ 0 ] * T->me[ 0 ][ 3 ];
			J->me[ 1 ][ i ] = -T->me[ 0 ][ 1 ] * T->me[ 1 ][ 3 ] + T->me[ 1 ][ 1 ] * T->me[ 0 ][ 3 ];
			J->me[ 2 ][ i ] = -T->me[ 0 ][ 2 ] * T->me[ 1 ][ 3 ] + T->me[ 1 ][ 2 ] * T->me[ 0 ][ 3 ];
			J->me[ 3 ][ i ] = T->me[ 2 ][ 0 ];  // nz
			J->me[ 4 ][ i ] = T->me[ 2 ][ 1 ];  // oz
			J->me[ 5 ][ i ] = T->me[ 2 ][ 2 ];  // az
		}
		else
		{
			// prismatic joint
			J->me[ 0 ][ i ] = T->me[ 2 ][ 0 ];  // nz
			J->me[ 1 ][ i ] = T->me[ 2 ][ 1 ];  // oz
			J->me[ 2 ][ i ] = T->me[ 2 ][ 2 ];  // az
			J->me[ 3 ][ i ] = 0;
			J->me[ 4 ][ i ] = 0;
			J->me[ 5 ][ i ] = 0;
		}
	}
	return( LLRL_OK );
}
