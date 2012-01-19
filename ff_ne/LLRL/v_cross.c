#include "LLRL.h"
/*************************************************************************	
    v_cross
	
    Copyright, Andrés Jaramillo Botero, 1999. 
***************************************************************************/
VEC* v_cross( VEC *v1, VEC *v2, VEC *v3 )
{
  /*
    (    0 -v1z +v1y )   ( v2x )   ( v1y * v2z ) - ( v1z * v2y )
    ( +v1z    0 -v1x ) * ( v2y ) = ( v1z * v2x ) - ( v1x * v2z )
    ( -v1y +v1x   0  )   ( v2z )   ( v1x * v2y ) - ( v1y * v2x ) */

  v3->ve[ 0 ] = ( v1->ve[ 1 ] * v2->ve[ 2 ] ) - ( v1->ve[ 2 ] * v2->ve[ 1 ] );
  v3->ve[ 1 ] = ( v1->ve[ 2 ] * v2->ve[ 0 ] ) - ( v1->ve[ 0 ] * v2->ve[ 2 ] );
  v3->ve[ 2 ] = ( v1->ve[ 0 ] * v2->ve[ 1 ] ) - ( v1->ve[ 1 ] * v2->ve[ 0 ] );
  return( v3 );
}
