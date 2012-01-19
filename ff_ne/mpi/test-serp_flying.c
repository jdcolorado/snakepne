/******************************************************************************
Program for testing inverse and forward dynamics for any serpentine robot fix or flying base.

	Copyright Robotics and Automation Group, Pontificia Universidad Javeriana - Cali.
		Andr� Jaramillo Botero, ajaramil@puj.edu.
		Antonio Alejandro Matta G�ez, amatta@puj.edu.
		Juli� David Colorado, jdcolorado@puj.edu.
		Juan Camilo Acosta, jcacosta@puj.edu.co
*************************************************************************
/*06/02/2006, JCA: First version: fix or flying base, not support floating base
  17/04/2006  JDC: Added support for Fbase model */

//#include<stdio.h>
//#include<stdlib.h>
#include "../LLRL/LLRL.h"

MAT *serp_DH(MAT *DH,int df);
MAT *serp_DYN(MAT *DYN,int df);
MAT *serp_DYN_base(MAT *DYN_base);

MAT *serp_DH(MAT *DH,int df)
{
	int i;
	DH=m_resize(DH,df,7);
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
	for (i=0;i<df;i++)
	{
		DYN->me[i][0]=0.3; DYN->me[i][1]=0.03; DYN->me[i][2]=1e-8; DYN->me[i][3]=0; DYN->me[i][4]=1.618e-4; DYN->me[i][5]=1.584e-4; DYN->me[i][6]=9.38e-5; DYN->me[i][7]=0; DYN->me[i][8]=0; DYN->me[i][9]=-3e-8;  DYN->me[i][10]=0.001; DYN->me[i][11]=0.1;
	}
	return(DYN);
}

MAT *serp_DYN_base(MAT *DYN_base)
{
	DYN_base=m_resize(DYN_base,1,4);
	DYN_base->me[0][0]=0.3; DYN_base->me[0][1]=0.03; DYN_base->me[0][2]=0; DYN_base->me[0][3]=0;
	return(DYN_base);
}

int main()
{
	mem_info_on(TRUE);
//	mem_stat_mark(1);
	//*******************************************************************************************
	static MAT *DH=MNULL,*DYN=MNULL,*DYN_base=MNULL, *Q_Qd_Qdd=MNULL, *Torques=MNULL,*vb=MNULL, *ab=MNULL, *Xcm=MNULL, *Fbase=MNULL;
	
	static VEC *fext=VNULL,*grav=VNULL;
	int robot_type;
	Real pi,n;
	clock_t start, middle, end;
	
	//*******************************************************************************************
	FILE *fp;
	grav=v_resize(grav,6);	/* zero gravity vector */
	fext=v_resize(fext,6);	/* zero external forces */
	mem_stat_reg_vars(0,TYPE_VEC,&grav,&fext,NULL);
	v_zero(grav);
	v_zero(fext);
	//grav->ve[5]=9.81;
	m_resize_vars(1,1,&Q_Qd_Qdd,NULL);//&DH,&DYN,&DYN_base,NULL);
	mem_stat_reg_vars(0,TYPE_MAT,&Q_Qd_Qdd,NULL);//&DH,&DYN,&DYN_base,NULL);
	
	
	Torques=m_resize(Torques,1,1);
	Fbase=m_resize(Fbase,1,1);
	mem_stat_reg_vars(0,TYPE_MAT,&Torques,&Fbase,NULL);
	
	
	pi=3.14159265359;
//Robot_type
	robot_type=4;
//Serpenoid parameters
	start=clock();

	if ( (fp=fopen("results.txt","r+")) != NULL )
	{	
		
		DH=m_finput(fp,MNULL);
		//MEM_STAT_REG(DH,TYPE_MAT);
		DYN=m_finput(fp,MNULL);
		//MEM_STAT_REG(DYN,TYPE_MAT);
		DYN_base=m_finput(fp,MNULL);
		//MEM_STAT_REG(DYN_base,TYPE_MAT);
		Q_Qd_Qdd=m_finput(fp,MNULL);
		//MEM_STAT_REG(Q_Qd_Qdd,TYPE_MAT);
		Xcm=m_finput(fp,MNULL);
		
		//xb=m_finput(fp,MNULL);
		//MEM_STAT_REG(xb,TYPE_MAT);
		//vb=m_finput(fp,MNULL);
		//MEM_STAT_REG(vb,TYPE_MAT);
		//ab=m_finput(fp,MNULL);
		//MEM_STAT_REG(ab,TYPE_MAT);
		
		fclose(fp);	
	}
	
	n=DH->m;

	
/**********************************************************
Inverse Dynamics
**********************************************************/	
	middle=clock();
	mem_stat_mark(2); 
	Torques = ff_ne(robot_type,DH,DYN,DYN_base,Q_Qd_Qdd,grav,fext,Torques,Xcm,vb,ab,Fbase);
	mem_stat_free(2); 
	end=clock();
	m_output(Torques);
	printf("\nTtotal= %f\nTinverse=%f\nTaccesodisco=%f\n",((double)(end-start))/CLOCKS_PER_SEC,((double)(end-middle))/CLOCKS_PER_SEC,((double)(middle-start))/CLOCKS_PER_SEC);


//	mem_stat_free(1);
	//mem_info();
	return 0;
}
