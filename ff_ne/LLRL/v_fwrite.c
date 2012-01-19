#include "LLRL.h"
/*************************************************************************	
	v_fwrite -- saves vector V wtih name "name" in a file fp

	Copyright, Andrés Jaramillo Botero, 1999. 
***************************************************************************/
VEC *v_fwrite(FILE *fp,VEC *V,char *name)
{
  static char    *format = "%14.8g";
  int     i;
	
  if ( ! V )
    error(E_NULL,"v_fwrite");
	
  fprintf(fp,"%s=[",name);
  /* write actual data */
	
  for ( i = 0; i < (int) V->dim; i++)
    fprintf(fp,format,V->ve[i]);
  fprintf(fp,"]\n");
	
  return V;
}

