MAT *conv(MAT *Q_Qd_Qddi,MAT *Q_Qd_Qddf, MAT *Xi, MAT *Vi, MAT *Ai, MAT *Xf, MAT *Vf, MAT *Af, double stepi, double stepf)
{
	int rel;
	int i,j;
	rel=int(stepf/stepi);
	Q_Qd_Qddf=m_resize(Q_Qd_Qddf,Q_Qd_Qddi->m/rel,Q_Qd_Qddi->n);
	
	for(j=0;j<Q_Qd_Qddi->n;j++)
		Q_Qd_Qddf->me[0][j]=Q_Qd_Qddf->me[0][j];
		
	for(i=1;i<Q_Qd_Qddi->m;i++)
	{
		if(i%rel=0)
		{
			for(j=0;j<Q_Qd_Qddi->n;j++)
				Q_Qd_Qddf->me[int(i/rel)][j]=Q_Qd_Qddf->me[i][j];
		}
	}


}