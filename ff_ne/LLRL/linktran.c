#include "LLRL.h"
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

	Copyright, Andrés Jaramillo Botero, 1999. 
***************************************************************************/
MAT* linktran( VEC *dh_link, double q, MAT *T )
{
  double theta, dn, ca, ct, sa, st;
  // validate inputs (sizes and initialization)
  if( dh_link == ( VEC* ) VNULL )
  {
    return( MNULL );
  }
  if( dh_link->dim < 4 )
  {
    return( MNULL );
  }
  if( ( T == ( MAT* ) MNULL ) || ( T->m != 4 ) || ( T->n != 4 ) )
  {
    T = m_resize( T, 4, 4 );
  }
  // set variable parameter according to sigma
  if( dh_link->dim > 4 )
  {
    if( dh_link->ve[ 4 ] == 0 )
    {
      theta = q;
      dn = dh_link->ve[ 3 ];
    }
    else
    {
      theta = dh_link->ve[ 2 ];
      dn = q;
    }
  } 
  else 
  {
    theta = q;
    dn = dh_link->ve[ 3 ];   
  }
  sa = sin( dh_link->ve[ 0 ] ); 
  ca = cos( dh_link->ve[ 0 ] );
  st = sin( theta ); 
  ct = cos( theta );
  // set up basic homogeneous transform matrix
  m_zero( T );
  T->me[ 0 ][ 0 ] = ct;
  T->me[ 0 ][ 1 ] = -st * ca;
  T->me[ 0 ][ 2 ] = st * sa;
  T->me[ 0 ][ 3 ] = dh_link->ve[ 1 ] * ct;
  T->me[ 1 ][ 0 ] = st;
  T->me[ 1 ][ 1 ] = ct * ca;
  T->me[ 1 ][ 2 ] = -ct * sa;
  T->me[ 1 ][ 3 ] = dh_link->ve[ 1 ] * st;
  T->me[ 2 ][ 1 ] = sa;
  T->me[ 2 ][ 2 ] = ca;
  T->me[ 2 ][ 3 ] = dn;
  T->me[ 3 ][ 3 ] = 1;
  /*
     T = [
           ct -st*ca  st*sa  an*ct
           st  ct*ca -ct*sa  an*st
            0     sa     ca     dn
            0      0      0      1
         ];
  */
  return( T );
}

