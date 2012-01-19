#include"LLRL.h"

MAT *spline_2(MAT *inMAT, VEC *b, VEC *X, VEC *h, VEC *time_vector, MAT *A_inv, int k, int points, double time, MAT *outMAT);

MAT *spline_2(MAT *inMAT, VEC *b, VEC *X, VEC *h, VEC *time_vector, MAT *A_inv, int k, int points, double time, MAT *outMAT)
{
	int pos,i;
	double aj,bj,cj,dj;
	for(pos=0;pos<6;pos++)
	{
		//Generating spline A matrix and b vector.
		for(i=1;i<=2;i++)
		{
			if(k==0)
			{	
				
				b->ve[i]=3*(inMAT->me[pos][k+1+i]-inMAT->me[pos][k+i])/h->ve[i]-3*(inMAT->me[pos][k+i]-inMAT->me[pos][k-1+i])/h->ve[i-1];
			}
			else
			{
				if(k==points-2)
					b->ve[i]=3*(inMAT->me[pos][k-1+i]-inMAT->me[pos][k-2+i])/h->ve[i]-3*(inMAT->me[pos][k-2+i]-inMAT->me[pos][k-3+i])/h->ve[i-1];
				else
					b->ve[i]=3*(inMAT->me[pos][k+i]-inMAT->me[pos][k-1+i])/h->ve[i]-3*(inMAT->me[pos][k-1+i]-inMAT->me[pos][k-2+i])/h->ve[i-1];
			}
		}
					
		//A.X=b --->> X=(A^-1)*b
		mv_mlt(A_inv,b,X);
		
		if(k==0)
		{
			aj=inMAT->me[pos][k];
			cj=X->ve[0];
			bj=(1/h->ve[0])*(inMAT->me[pos][k+1]-aj)-(h->ve[0]/3)*(2*cj+X->ve[1]);
			dj=(X->ve[1]-cj)/(3*h->ve[0]);
		}
		else
		{
			if(k==points-2)
			{
				aj=inMAT->me[pos][k];
				cj=X->ve[2];
				bj=(1/h->ve[2])*(inMAT->me[pos][k+1]-aj)-(h->ve[2]/3)*(2*cj+X->ve[3]);
				dj=(X->ve[3]-cj)/(3*h->ve[2]);
			}
			else
			{
				aj=inMAT->me[pos][k];
				cj=X->ve[1];
				bj=(1/h->ve[1])*(inMAT->me[pos][k+1]-aj)-(h->ve[1]/3)*(2*cj+X->ve[2]);
				dj=(X->ve[2]-cj)/(3*h->ve[1]);
			}
		}
		outMAT->me[pos][0]=aj+bj*(time-time_vector->ve[k])+cj*pow((time-time_vector->ve[k]),2)+dj*pow((time-time_vector->ve[k]),3);
	}
}

/***********************************************************************************************************
torqfun.c

		torqfun make a cubic spline of inTAU matrix in order to estimate a torque at any time.
		it uses spline_2 to calculate the second part of the spline for base trayectory.

Copyright 	Robotic and Automation Group, Pontificia Universidad Javeriana, Cali
	
		Andrés Jaramillo Botero
		Juan Camilo Acosta
		Julian David Colorado
		Antonio Matta.

%MAT *inTAU: input torques. Size points x df(Degrees of Freedom)
%VEC *time_vector: input time vector. Dim: points
%double time: estimation torque time input.
%MAT *outTAU: output matrix with torque evaluated at time. Size: 1 x df
**************************************************************************************************************/

/*28/03/2006: JCA: First version, doing straight line between 2 points.*/
/*29/03/2006: JCA: Improved with cubic splines.*/
/*02/05/2006: JCA: Removed base trayectory*/

MAT *torqfun(MAT *DH, MAT *DYN,MAT *inTAU,VEC *time_vector, double time, MAT *outTAU)
{
	int i,j,df,k,points,art,pos;
	double delta,delta2,aj,bj,cj,dj;
	static VEC  *b=VNULL, *X=VNULL, *h=VNULL;//*temp=VNULL,
	static MAT *A=MNULL, *A_inv=MNULL, *t=MNULL; 
	
	h=v_resize(h,3);
	v_resize_vars(4,&b,&X,NULL);
	
	mem_stat_reg_vars(0,TYPE_VEC,&h,&b,&X,NULL);
	
	t=m_resize(t,6,1);
	m_resize_vars(4,4,&A,&A_inv,NULL);
	mem_stat_reg_vars(0,TYPE_MAT,&A,&A_inv,&t,NULL);
	
	if(inTAU!=NULL)
		if(inTAU->m != time_vector->dim)
			error(E_SIZES,"torqfun TAU vs time-vector rows");
		
	points=time_vector->dim;
	if(time>(time_vector->ve[points-1]))
	{
		printf("WARNING: time(%f) out of range, replaced with last time(%f).\n",time,(time_vector->ve[points-1]));
		time=(time_vector->ve[points-1]);
	}
	
	
	//solucionar problema de correr el ciclo hasta el final siempre....
	for(i=0;i<(points-1);i++)
	{
		if(time_vector->ve[i]<=time && time_vector->ve[i+1]>=time)
			k=i;
	}

	//Generating spline step vector h.
	for(i=0;i<3;i++)
	{
		if(k==0)
		{
			h->ve[0]=time_vector->ve[k+1]-time_vector->ve[k];
			h->ve[1]=time_vector->ve[k+2]-time_vector->ve[k+1];
			h->ve[2]=time_vector->ve[k+3]-time_vector->ve[k+2];
			
		}
		else
		{
			if(k==points-2)
			{
				h->ve[0]=time_vector->ve[k-1]-time_vector->ve[k-2];
				h->ve[1]=time_vector->ve[k]-time_vector->ve[k-1];
				h->ve[2]=time_vector->ve[k+1]-time_vector->ve[k];	
			}
			else
			{
				h->ve[0]=time_vector->ve[k]-time_vector->ve[k-1];
				h->ve[1]=time_vector->ve[k+1]-time_vector->ve[k];
				h->ve[2]=time_vector->ve[k+2]-time_vector->ve[k+1];
				
			}
		}
	}
	A->me[0][0]=1;
	A->me[3][3]=1;
	for(i=1;i<=2;i++)
		{
			A->me[i][i-1]=h->ve[i-1];
			A->me[i][i]=2*(h->ve[i-1]+h->ve[i]);
			A->me[i][i+1]=h->ve[i];
		}
	
	A_inv=pinv(A,A_inv);
	if(inTAU!=NULL)
	{
		df=inTAU->n;
		//loop for each joint torque.
		for(art=0;art<df;art++)
		{
			//Generating spline A matrix and b vector.
			for(i=1;i<=2;i++)
			{
				if(k==0)
				{
					b->ve[i]=3*(inTAU->me[k+1+i][art]-inTAU->me[k+i][art])/h->ve[i]-3*(inTAU->me[k+i][art]-inTAU->me[k-1+i][art])/h->ve[i-1];
				}
				else
				{
					if(k==points-2)
						b->ve[i]=3*(inTAU->me[k-1+i][art]-inTAU->me[k-2+i][art])/h->ve[i]-3*(inTAU->me[k-2+i][art]-inTAU->me[k-3+i][art])/h->ve[i-1];
					else
						b->ve[i]=3*(inTAU->me[k+i][art]-inTAU->me[k-1+i][art])/h->ve[i]-3*(inTAU->me[k-1+i][art]-inTAU->me[k-2+i][art])/h->ve[i-1];
				}
			}
			
			//A.X=b --->> X=(A^-1)*b
			mv_mlt(A_inv,b,X);
			
			if(k==0)
			{
				aj=inTAU->me[k][art];
				cj=X->ve[0];
				bj=(1/h->ve[0])*(inTAU->me[k+1][art]-aj)-(h->ve[0]/3)*(2*cj+X->ve[1]);
				dj=(X->ve[1]-cj)/(3*h->ve[0]);
			}
			else
			{
				if(k==points-2)
				{
					aj=inTAU->me[k][art];
					cj=X->ve[2];
					bj=(1/h->ve[2])*(inTAU->me[k+1][art]-aj)-(h->ve[2]/3)*(2*cj+X->ve[3]);
					dj=(X->ve[3]-cj)/(3*h->ve[2]);
				}
				else
				{
					aj=inTAU->me[k][art];
					cj=X->ve[1];
					bj=(1/h->ve[1])*(inTAU->me[k+1][art]-aj)-(h->ve[1]/3)*(2*cj+X->ve[2]);
					dj=(X->ve[2]-cj)/(3*h->ve[1]);
				}
			}
			outTAU->me[0][art]=aj+bj*(time-time_vector->ve[k])+cj*pow((time-time_vector->ve[k]),2)+dj*pow((time-time_vector->ve[k]),3);
		}
	}
		
	//doing a straight line	instead of a cubic spline 
	/*delta=time_vector->ve[i+1]-time_vector->ve[i];
	delta2=time-time_vector->ve[i];

	for(j=0;j<df;j++)
	{
		row1->ve[j]=inTAU->me[k][j];
		row2->ve[j]=inTAU->me[k+1][j];
	}
	
	v_sub(row2,row1,temp);
	sv_mlt(delta2/delta,temp,temp);
	v_add(temp,row1,temp);
	for(j=0;j<df;j++)
		outTAU->me[0][j]=temp->ve[j];*/
	
	return outTAU;
}
