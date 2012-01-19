#include "../LLRL/LLRL.h"

MAT *serp_DH(MAT *DH,int df);
MAT *serp_DYN(MAT *DYN,int df);
MAT *serp_DYN_base(MAT *DYN_base);

MAT *serp_DH(MAT *DH,int df)
{
	int i;
	DH=m_resize(DH,df,7);
	//MEM_STAT_REG(DH,TYPE_MAT);
	for (i=0;i<df;i++)
	{
		DH->me[i][0]=0; DH->me[i][1]=-0.1; DH->me[i][2]=0; DH->me[i][3]=0; DH->me[i][4]=0; DH->me[i][5]=-1.5707963268; DH->me[i][6]=1.5707963268;
	}
	return(DH);
}


MAT *serp_DYN(MAT *DYN,int df)
{
	int i;
	DYN=m_resize(DYN,df,12);
	//MEM_STAT_REG(DYN,TYPE_MAT);
	for (i=0;i<df;i++)
	{
		DYN->me[i][0]=0.3; DYN->me[i][1]=0.03; DYN->me[i][2]=1e-8; DYN->me[i][3]=0; DYN->me[i][4]=1.618e-4; DYN->me[i][5]=1.584e-4; DYN->me[i][6]=9.38e-5; DYN->me[i][7]=0; DYN->me[i][8]=0; DYN->me[i][9]=-3e-8;  DYN->me[i][10]=0.001; DYN->me[i][11]=0.1;
	}
	return(DYN);
}

MAT *serp_DYN_base(MAT *DYN_base)
{
	DYN_base=m_resize(DYN_base,1,4);
	//MEM_STAT_REG(DYN_base,TYPE_MAT);
	DYN_base->me[0][0]=0.3; DYN_base->me[0][1]=0.03; DYN_base->me[0][2]=0; DYN_base->me[0][3]=0;
	return(DYN_base);
}

int main()
{
	int df,pt,k,j,i,cont;
	Real a,n,b,c,w,step,end_t,pi;
	double alfa,teta,d;
	static MAT *Q_Qd_Qdd=MNULL,*DH=MNULL, *DYN=MNULL, *DYN_base=MNULL, *xb=MNULL, *vb=MNULL, *ab=MNULL, *Xcm=MNULL, *T1=MNULL, *Ttotal=MNULL, *orient=MNULL, *T2=MNULL;
	static VEC *Fr6=VNULL;

	FILE *fp;
	df=9;

	pi=3.14159265359;
	n=df+1;
	//a=pi/3;
	a=1.0471;
	//b=2*pi;
	b=6.2832;
	c=0;
	w=1.0471;
	//w=pi/3;

	step=0.001;
	end_t=2.002;
//	step=0.1;
//	end_t=2;
	v_resize_vars(6,&Fr6,NULL);
	mem_stat_reg_vars(0,TYPE_VEC,&Fr6,NULL);
	
	m_resize_vars(1,1,&Ttotal,&orient,NULL);
	mem_stat_reg_vars(0,TYPE_MAT,&Ttotal,&orient,NULL);
	
	m_resize_vars(1,1,&xb,&vb,&ab,&Xcm,NULL);
	mem_stat_reg_vars(0,TYPE_MAT,&xb,&vb,&ab,&Xcm,NULL);
	/**************** 4x4 **************/
	/*T1-T2-T2: Temp matrix.*/
	m_resize_vars(4,4,&T1,&T2,NULL);
	mem_stat_reg_vars(0,TYPE_MAT,&T1,&T2,NULL);
	

	DH=serp_DH(DH,df);
	DYN=serp_DYN(DYN,df);
	DYN_base=serp_DYN_base(DYN_base);


	mem_stat_mark(2);
	Q_Qd_Qdd=serpenoid(DH,DYN,DYN_base,Q_Qd_Qdd, n, a, b, c, w, step, end_t);
	mem_stat_free(2);
	
	
	
	#ifdef print_traj
	if(my_id==0)
	{
		printf("\n-------------------------------Trajectory:----------------------------------\n\n");
		m_output(Q_Qd_Qdd);
	}
	#endif
	m_resize_vars(6,Q_Qd_Qdd->m,&xb,&vb,&ab,&Xcm,NULL);
	m_resize_vars(4*Q_Qd_Qdd->m,4*df,&Ttotal,NULL);

	for(pt=0;pt<Q_Qd_Qdd->m;pt++)
	{
		xb->me[0][pt]=Q_Qd_Qdd->me[pt][3*df+1];
		xb->me[1][pt]=Q_Qd_Qdd->me[pt][3*df+1+1];
		xb->me[2][pt]=Q_Qd_Qdd->me[pt][3*df+1+2];
		xb->me[3][pt]=Q_Qd_Qdd->me[pt][3*df+1+3];
		xb->me[4][pt]=Q_Qd_Qdd->me[pt][3*df+1+4];
		xb->me[5][pt]=Q_Qd_Qdd->me[pt][3*df+1+5];
		
		vb->me[0][pt]=Q_Qd_Qdd->me[pt][3*df+1+6];
		vb->me[1][pt]=Q_Qd_Qdd->me[pt][3*df+1+1+6];
		vb->me[2][pt]=Q_Qd_Qdd->me[pt][3*df+1+2+6];
		vb->me[3][pt]=Q_Qd_Qdd->me[pt][3*df+1+3+6];
		vb->me[4][pt]=Q_Qd_Qdd->me[pt][3*df+1+4+6];
		vb->me[5][pt]=Q_Qd_Qdd->me[pt][3*df+1+5+6];
		
		ab->me[0][pt]=Q_Qd_Qdd->me[pt][3*df+1+12];
		ab->me[1][pt]=Q_Qd_Qdd->me[pt][3*df+1+1+12];
		ab->me[2][pt]=Q_Qd_Qdd->me[pt][3*df+1+2+12];
		ab->me[3][pt]=Q_Qd_Qdd->me[pt][3*df+1+3+12];
		ab->me[4][pt]=Q_Qd_Qdd->me[pt][3*df+1+4+12];
		ab->me[5][pt]=Q_Qd_Qdd->me[pt][3*df+1+5+12];
	}
	
	m_move(Q_Qd_Qdd,0,1,Q_Qd_Qdd->m,3*df,Q_Qd_Qdd,0,0);
	Q_Qd_Qdd=m_resize(Q_Qd_Qdd,Q_Qd_Qdd->m,3*df);
	m_resize_vars(Q_Qd_Qdd->m,df,&orient,NULL);
	
	mem_stat_mark(2);
	Xcm=base_x(DH,DYN,DYN_base,Q_Qd_Qdd,Xcm);
	mem_stat_free(2);
	
	
	for(k=0;k<Q_Qd_Qdd->m;k++)
	{
		for(i=1;i<=df;i++)
		{
			alfa = DH->me[i-1][0];
			a=DH->me[i-1][1];
			if(DH->me[i-1][4]==0)        //Rotational
			{
				teta=Q_Qd_Qdd->me[k][i-1];
				d=DH->me[i-1][3];
			}
			else                         //Prismatic
			{
				teta=DH->me[i-1][2];
				d=Q_Qd_Qdd->me[k][i-1];
			}
			homogeneus(alfa,a,teta,d,T1);	//obtain homogeneus transformation
				//m_move(Rot,0,0,3,3,r,0,0);	//obtain rotation matrix
			T1->me[0][3]=a*cos(teta);
			T1->me[1][3]=a*sin(teta);
			T1->me[2][3]=d;
			T1->me[3][3]=1;
			if(i==1)
			{
				m_move(T1,0,0,4,4,Ttotal,4*k,0);
			}
			else
			{
				for(cont=0;cont<=3;cont++)
					for(j=0;j<=3;j++)	
						Ttotal->me[cont+4*k][j+4*(i-1)]=Ttotal->me[cont+4*k][0+4*(i-2)]*T1->me[0][j] 
								+Ttotal->me[cont+4*k][1+4*(i-2)]*T1->me[1][j]
								+Ttotal->me[cont+4*k][2+4*(i-2)]*T1->me[2][j]
								+Ttotal->me[cont+4*k][3+4*(i-2)]*T1->me[3][j];
			}
		}
	}
	for(k=0;k<Q_Qd_Qdd->m;k++)
	{
		get_col(Xcm,k,Fr6);
		euler(Fr6,T2,Fr6);	//Obtain base homogeneus transformation
					
		T1->me[0][0]=T2->me[0][0];
		T1->me[1][0]=T2->me[1][0];
		T1->me[2][0]=T2->me[2][0];
		T1->me[3][0]=0;
		
		T1->me[0][1]=T2->me[0][1];
		T1->me[1][1]=T2->me[1][1];
		T1->me[2][1]=T2->me[2][1];
		T1->me[3][1]=0;
		
		T1->me[0][2]=T2->me[0][2];
		T1->me[1][2]=T2->me[1][2];
		T1->me[2][2]=T2->me[2][2];
		T1->me[3][2]=0;
		
		T1->me[0][3]=Fr6->ve[0];
		T1->me[1][3]=Fr6->ve[1];
		T1->me[2][3]=Fr6->ve[2];
		T1->me[3][3]=1;
		
		for(i=0;i<df;i++)
		{
			for(cont=0;cont<=3;cont++)
				for(j=0;j<=3;j++)	
					T2->me[cont][j]=T1->me[cont][0]*Ttotal->me[0+4*k][j+4*i] 
							+T1->me[cont][1]*Ttotal->me[1+4*k][j+4*i]
							+T1->me[cont][2]*Ttotal->me[2+4*k][j+4*i]
							+T1->me[cont][3]*Ttotal->me[3+4*k][j+4*i];
				
					//m_mlt(T1,T2,T2);
			
			mem_stat_mark(15);
			Fr6=tr2rpy(T2,Fr6);
			mem_stat_free(15);
			orient->me[k][i]=Fr6->ve[0];
		}
	}
	if ( (fp=fopen("results.txt","wr+")) != NULL )
	{	
		
		m_foutput(fp,DH);
		m_foutput(fp,DYN);
		m_foutput(fp,DYN_base);
		m_foutput(fp,Q_Qd_Qdd);
		m_foutput(fp,orient);
		//m_foutput(fp,Xcm);
		//m_foutput(fp,Ttotal);
		//m_foutput(fp,xb);
		//m_foutput(fp,vb);
		//m_foutput(fp,ab);
		
		fclose(fp);	
	}
	
	return 0;
}
