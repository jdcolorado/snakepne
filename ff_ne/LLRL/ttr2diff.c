#include"LLRL.h"
/*************************************************************************	
	Conversion from homogeneous transform difference to 6D differential vector 

	d=ttr2diff(T1,T2,d)
  
	T1/2=4x4 homogeneous transformation matrix representing basic positions and 
	orientations at different points in space.
	d=6D vector defining 3 differential translations and 3 differential rotations 
	
	Copyright, Andrés Jaramillo Botero, 1999. 
***************************************************************************/
VEC	*ttr2diff(MAT *T1,MAT *T2,VEC *d)
{
	/* setup temporals */
	static	VEC	*nc,*oc,*ac,*v1,*v2;
	
	/* validate inputs (sizes and initializations) */
	if (T1 == (MAT *)MNULL || T2 == (MAT *)MNULL) 
		error(E_NULL, "ttr2diff (input matrices)\n");
	if (T1->m != 4 || T1->n != 4 || T2->m != 4 || T2->n != 4)	
		error(E_SIZES,"ttr2diff (input matrices)\n");
	if (d == (VEC *)VNULL || d->dim != 6) 
		d=v_resize(d,6);
	
	/* allocate memory */
	v_resize_vars(4,&v1,&v2,&nc,&oc,&ac,NULL);
	mem_stat_reg_vars(0,TYPE_VEC,&v1,&v2,&nc,&oc,&ac,NULL);
	
	/* set translational component */
	d->ve[0]=T2->me[0][3]-T1->me[0][3];
	d->ve[1]=T2->me[1][3]-T1->me[1][3];
	d->ve[2]=T2->me[2][3]-T1->me[2][3];
	
	get_col(T2,0,v1);
	get_col(T1,0,v2);
	v_resize_vars(3,&v1,&v2,NULL);
	v_cross(v2,v1,nc);
	
	v_resize_vars(4,&v1,&v2,NULL);
	get_col(T2,1,v1);
	get_col(T1,1,v2);
	v_resize_vars(3,&v1,&v2,NULL);
	v_cross(v2,v1,oc);
	
	v_resize_vars(4,&v1,&v2,NULL);
	get_col(T2,2,v1);
	get_col(T1,2,v2);
	v_resize_vars(3,&v1,&v2,NULL);
	v_cross(v2,v1,ac);
	
	v_add(nc,oc,v1);
	v_add(v1,ac,v2);
	sv_mlt(0.5,v2,v2);
	
	/* set rotational component */
	d->ve[3]=v2->ve[0];
	d->ve[4]=v2->ve[1];
	d->ve[5]=v2->ve[2];
	
	return d;
}

