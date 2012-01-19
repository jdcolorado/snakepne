#include "LLRL.h"
/*************************************************************************	
	m_fwrite -- saves matrix M with name "name" in a file fp

	Copyright, Andrés Jaramillo Botero, 1999. 
***************************************************************************/
MAT *m_fwrite(FILE *fp,MAT *M,char *name)
{
	static char    *format = "%14.8g";
	int     i,j;
	
	if ( ! M )
		error(E_NULL,"m_fwrite");
	
	fprintf(fp,"%s=[",name);
	/* write actual data */
	for ( j = 0; j < (int) M->m; j++ )
	{
		for ( i = 0; i < (int) M->n; i++)
		{
			fprintf(fp,format,M->me[j][i]);
			if (i < (int)M->n-1) putc(' ',fp);
		}
		if ( i==(int)M->n && j < (int)M->m-1 ) 
		    fprintf(fp,";\n");  //Aquí iba un ,fp 
	}
	fprintf(fp,"]\n");
	
	return M;
}


