#include "LLRL.h"

//**************************************************************************/
//Returns the Homogeneus Transformation Matrix for the base
//X: Spatial Position of the Base: [roll pitch yaw X  Y  Z]'
//Obtain [r p]:Matrix Rotation, Position Vector
/*21/04/2006 JCA: check for input sizes*/
void euler(VEC *columna, MAT *r, VEC *p)
{
	float roll,pitch,yaw;
	if ((columna->dim) <3 || (r->n<3 && r->m<3) )
			error(E_SIZES,"EULER: wrong input sizes");
	roll=columna->ve[0];
	pitch=columna->ve[1];
	yaw=columna->ve[2];
	//Orientation Matrix
	r->me[0][0]=cos(roll)*cos(pitch);
	r->me[0][1]=cos(roll)*sin(pitch)*sin(yaw)-sin(roll)*cos(yaw);
	r->me[0][2]=cos(roll)*sin(pitch)*cos(yaw)+sin(roll)*sin(yaw);
	r->me[1][0]=sin(roll)*cos(pitch);
	r->me[1][1]=sin(roll)*sin(pitch)*sin(yaw)+cos(roll)*cos(yaw);
	r->me[1][2]=sin(roll)*sin(pitch)*cos(yaw)-cos(roll)*sin(yaw);
	r->me[2][0]=-sin(pitch);
	r->me[2][1]=cos(pitch)*sin(yaw);
	r->me[2][2]=cos(pitch)*cos(yaw);
	//Position Vector
	p->ve[0]=columna->ve[3];
	p->ve[1]=columna->ve[4];
	p->ve[2]=columna->ve[5];
}
