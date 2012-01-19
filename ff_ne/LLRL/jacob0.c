#include "LLRL.h"
/*************************************************************************
Workspace used: 10,7
*************************************************************************	
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
/*26/04/2006 JCA: changed identation.*/
int jacob0( MAT *(*fun)(),MAT *DH ,MAT *DYN,VEC *q, MAT *J )
{
	static  MAT *Tn = MNULL, *Jn = MNULL, *Q = MNULL;
	int i, j, ErrorStatus;
	
	// validate inputs
	if( DH == ( MAT* ) MNULL || q == ( VEC* )VNULL )
		return( LLRL_ERROR );
	
	if( q->dim != DH->m )
		return( LLRL_ERROR );
	
	if( J == ( MAT* ) MNULL || J->m != 6 || J->n != DH->m )
		J = m_resize( J, 6, DH->m );
	
	// setup and allocate memory for temporals
	Tn = m_resize( Tn, 4, 4);
	Jn = m_resize( Jn, 6, 6 );
	Q = m_resize( Q, 1, DH->m );
	mem_stat_reg_vars( 0, TYPE_MAT, &Tn, &Jn, &Q, NULL );
	m_zero( J );
	set_row( Q, 0, q );
	// end-effector transformation
	mem_stat_mark(10);
	ErrorStatus = fkine( DH, Q, Tn );
	mem_stat_free(10);
	if( ErrorStatus == LLRL_ERROR )
		return( LLRL_ERROR );
	sub_mat( Tn, 0, 0, 2, 2, Tn );
	
	for( i = 0; i < 3; i++ )
		for( j = 0; j < 3; j++ )
			Jn->me[ i ][ j ] = Tn->me[ i ][ j ];
	
	for( i = 3; i < 6; i++ )
		for( j = 3; j < 6; j++ )
			Jn->me[ i ][ j ] = Tn->me[ i - 3 ][ j - 3 ];
	
	Tn = m_resize( Tn, 6, DH->m );
	mem_stat_mark(10);
	ErrorStatus = (*fun)( DH, DYN, q, Tn );
	mem_stat_free(10);

	if( ErrorStatus == LLRL_ERROR )
		return( LLRL_ERROR );
	
	J = m_mlt( Jn, Tn, J );
	return( LLRL_OK );
}

