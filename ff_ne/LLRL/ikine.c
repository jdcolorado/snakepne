#include "LLRL.h"
/********************************************************************************************************
Workspace used: 3,7
*************************************************************************
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
/*20/04/2006 JCA: changed jacob0 interface.*/
/*26/04/2006 JCA: changed identation.*/
int ikine( MAT *mDHParameters, MAT *mHMT, MAT *mJointParameters )
{
	static  VEC *dq = VNULL, *e = VNULL, *q = VNULL;
	static  MAT *m1 = MNULL, *m2 = MNULL, *m3 = MNULL, *DYN=MNULL;
	static  MAT *T = MNULL, *Qi = MNULL, *DQ = MNULL;
	int ilimit, count, tcount, i, ErrorStatus;
	double stol, nm;
	// validate inputs (sizes and initial values)
	if( mDHParameters == ( MAT* ) MNULL || mHMT == ( MAT* ) MNULL )
		return( LLRL_ERROR );
	if( mHMT->m != 4 )
		return( LLRL_ERROR );
	// Initial guess value for mJointParameters = 0
	if( mJointParameters == ( MAT* ) MNULL )
	{
		mJointParameters = m_resize( mJointParameters, mHMT->n / 4, mDHParameters->m );
		m_zero( mJointParameters );
	}
	// set iteration limit
	ilimit = 2000;
	// set convergence criteria
	stol = 1e-8;
	// Setup and allocate space for temporals
	m1 = m_resize( m1, 6, mDHParameters->m );
	m_resize_vars( 4, 4, &m2, &T, NULL );
	m_resize_vars( 1, mDHParameters->m, &Qi, &DQ, NULL );
	m3 = m_resize( m3, m1->n, m1->m );
	e = v_resize( e, 6 );
	v_resize_vars( m1->n, &dq, &q, NULL );
	//v_resize_vars( 6, &e, NULL );
	mem_stat_reg_vars( 0, TYPE_MAT, &m1, &m2, &m3, &T, &Qi, &DQ, NULL );
	mem_stat_reg_vars( 0, TYPE_VEC, &e, &dq, &q, NULL );
	// initialize point trajectory counter
	tcount = 0;
	if( mHMT->m == 4 && mHMT->n == 4 )
	{
		nm = 1;
		count = 0;
		q = get_row( mJointParameters, 0, q );
		while( nm > stol )
		{
			mem_stat_mark( 3 );
			ErrorStatus = fkine( mDHParameters, mJointParameters, m2 );
			mem_stat_free( 3 );
			if( ErrorStatus == LLRL_ERROR )
				return( LLRL_ERROR );
			
			mem_stat_mark( 3 );
			ttr2diff( m2, mHMT, e );
			mem_stat_free( 3 );
			
			mem_stat_mark( 3 );
			ErrorStatus = jacob0(jacobn, mDHParameters,DYN, q, m1 );
			mem_stat_free( 3 );
			if( ErrorStatus == LLRL_ERROR )
				return( LLRL_ERROR );
			mem_stat_mark( 3 );
			m3 = pinv( m1, m3);
			mem_stat_free( 3 );
			dq = mv_mlt( m3,e, dq );
			set_row( DQ, 0, dq );
			v_add(q,dq,q);
			m_add( mJointParameters, DQ, mJointParameters );
			//nm = v_norm2( dq );
			nm = v_norm2( dq );
			
			count++;
			if( count > ilimit )
				return( LLRL_ERROR );
		}
	}
	else
	{
		mJointParameters = m_resize( mJointParameters, mHMT->n/4, mDHParameters->m );
		m_move( mJointParameters, 0, 0, 1, mDHParameters->m, Qi, 0, 0 );
		q = get_row( Qi, 0, q );
		for( i = 0; i < ( int ) mHMT->n / 4; i++ )
		{
			nm = 1;
			T = xttr( mHMT, i + 1, T);
			count = 0;
			while( nm > stol )
			{
				mem_stat_mark( 3 );
				ErrorStatus = fkine( mDHParameters, Qi, m2 );
				mem_stat_free( 3 );
				if( ErrorStatus == LLRL_ERROR )
					return( LLRL_ERROR );
				mem_stat_mark( 3 );
				ttr2diff( m2, T, e );
				mem_stat_free( 3 );
				mem_stat_mark( 3 );
				ErrorStatus = jacob0( jacobn,mDHParameters, DYN,q, m1 );
				mem_stat_free( 3 );
				if( ErrorStatus == LLRL_ERROR )
					return( LLRL_ERROR );
				mem_stat_mark( 3 );
				m3 = pinv( m1, m3 );
				mem_stat_free( 3 );
				dq = mv_mlt( m3, e, dq );
				set_row( DQ, 0, dq );
				Qi = m_add( Qi, DQ, Qi );
				nm = v_norm2( dq );
				
				count++;
				if( count > ilimit )
					return( LLRL_ERROR );
				q = get_row( Qi, 0, q );
			}
			m_move( Qi, 0, 0, 1, mDHParameters->m, mJointParameters, i, 0 );
			tcount = tcount + count;
		}
	}
	return( LLRL_OK );
}
