#include "LLRL.h"
/*************************************************************************	
	Conversion from homogeneous transform (rotation only) to quaternion

	Q=tr2q(T,Q)
  
	T=4x4 homogeneous transformation matrix representing rotation.
	Q=vector defining quaternion parameters [q1,q2,q3,q4] 
	
	Copyright, Andrés Jaramillo Botero, 1999. 
***************************************************************************/
VEC *tr2q(MAT *T,VEC *Q)
{
  static	double	kx, ky, kz, kx1, ky1, kz1, nm, s;
  static	int	add;
	
  /* validate inputs (sizes and initializations) */
  if (T == (MAT *)MNULL) 
    error(E_NULL, "tr2q (input matrix)\n");
  if (T->m != 4 || T->n != 4)	
    error(E_SIZES,"tr2q (input matrix)\n");
  if (Q == (VEC *)VNULL || Q->dim != 4) Q=v_resize(Q,4);
	
  Q->ve[0]=sqrt(T->me[0][0]+T->me[1][1]+T->me[2][2]+T->me[3][3])/2;
  kx=T->me[2][1]-T->me[1][2];		/* Oz-Ay */
  ky=T->me[0][2]-T->me[2][0];		/* Ax-Nz */
  kz=T->me[1][0]-T->me[0][1];		/* Ny-Ox */
	
  if (T->me[0][0] >= T->me[1][1] && T->me[0][0] >= T->me[2][2])
  {
    kx1=T->me[0][0]-T->me[1][1]-T->me[2][2]+1;	/* Nx - Oy - Az + 1 */
    ky1=T->me[1][0]+T->me[0][1];	/* Ny + Ox */
    kz1=T->me[2][0]+T->me[0][2];	/* Nz + Ax */
		
    if (kx >= 0) add=1;		
    else add=0;
  }
  else
  {
    if (T->me[1][1] >= T->me[2][2])
    {
      kx1=T->me[1][0]+T->me[0][1];	/* Ny + Ox */
      ky1=T->me[1][1]-T->me[0][0]-T->me[2][2]+1;	/* Oy - Nx - Az + 1 */
      kz1=T->me[2][1]+T->me[1][2];	/* Oz + Ay */
			
      if (ky >= 0) add=1;
      else add=0;
    }
    else
    {
      kx1=T->me[2][0]+T->me[0][2];	/* Nz + Ax */
      ky1=T->me[2][1]+T->me[1][2];	/* Oz + Ay */
      kz1=T->me[2][2]-T->me[0][0]-T->me[1][1]+1;	/* Az - Nx - Oy + 1 */
			
      if (kz >= 0) add=1;
      else add=0;
    }
  }
  if (add == 1)
  {
    kx=kx+kx1;
    ky=ky+ky1;
    kz=kz+kz1;
  }
  else
  {
    kx=kx-kx1;
    ky=ky-ky1;
    kz=kz-kz1;
  }
  nm=sqrt(pow(kx,2)+pow(ky,2)+pow(kz,2));
  if (nm==0) Q->ve[0]=1;
  else
  {
    s=sqrt(1-pow(Q->ve[0],2))/nm;
    Q->ve[1]=s*kx;
    Q->ve[2]=s*ky;
    Q->ve[3]=s*kz;
  }
	
  return Q;
}
