/*************************************************************************	
	Conversion from homogeneous transform to euler angles (y-convention)

	E=tr2rpy(T,E)
  
	T=4x4 homogeneous transformation matrix representing rotation.
	E=vector defining [roll,pitch,jaw] rotation angles (in radians)
	
	Copyright, Andr� Jaramillo Botero, 1999. 
***************************************************************************/

//#include <math.h>
//#include "matrix2.h"
#include "LLRL.h"
VEC *tr2rpy(MAT *T,VEC *E)
{
	/* validate inputs (sizes and initializations) */
	if (T == (MAT *)MNULL) 
		error(E_NULL, "tr2rpy (input matrix)\n");
	if (T->m != 4 || T->n != 4)	
		error(E_SIZES,"tr2rpy (input matrix)\n");
	if (E == (VEC *)VNULL || E->dim < 3) E=v_resize(E,3);
	
	if (fabs(T->me[0][0]) < D_MACHEPS && fabs(T->me[1][0]) < D_MACHEPS)
    {
		E->ve[0]=0;									/* roll */
		E->ve[1]=atan2(-T->me[2][0],T->me[0][0]);	/* pitch */
		E->ve[2]=atan2(-T->me[1][2],T->me[1][1]);	/* jaw */
    }
	else
    {
		E->ve[0]=atan2(T->me[1][0],T->me[0][0]);	/* roll */
		E->ve[1]=atan2(-T->me[2][0],cos(E->ve[0])*T->me[0][0]+
			sin(E->ve[0])*T->me[1][0]);				/* pitch*/
		E->ve[2]=atan2(sin(E->ve[0])*T->me[0][2]-cos(E->ve[0])*T->me[1][2],
			-sin(E->ve[0])*T->me[0][1]+cos(E->ve[0])*T->me[1][1]); /*jaw */
    }
	
	return E;
}

