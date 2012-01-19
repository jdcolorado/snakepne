#include "LLRL.h"
/*************************************************************************	
	Return the manipulator link inertia tensor matrix.  Extract the 
	manipulator link inertia for "link" from DYN and return it in matrix 
	form.

	J=linkiner(DYN,link,J)

	J=inertia tensor for "link"
	link=link number
	DYN=dynamic parameter matrix for manipulator

	Copyright, Andrés Jaramillo Botero, 1999. 
***************************************************************************/
MAT	*linkiner(MAT *DYN,u_int link,MAT *J)
{
	if (DYN == (MAT *) MNULL) 
		error(E_NULL,"linkiner (dynamic parameters)");
	if (J == (MAT *) MNULL || (J->m != 3 && J->n != 3)) 
		J=m_resize(J,3,3);

	/* moments of inertia */
	J->me[0][0]=DYN->me[link][4];
	J->me[1][1]=DYN->me[link][5];;
	J->me[2][2]=DYN->me[link][6];;

	/* products of inertia */
	J->me[0][1]=DYN->me[link][7];
	J->me[0][2]=DYN->me[link][9];
	J->me[1][0]=DYN->me[link][7];
	J->me[1][2]=DYN->me[link][8];
	J->me[2][0]=DYN->me[link][9];
	J->me[2][1]=DYN->me[link][8];

	/*
	J =[Ixx	Ixy	Ixz;
		Ixy	Iyy Iyz;
		Ixz Iyz Izz]
	*/
	return J;
}
