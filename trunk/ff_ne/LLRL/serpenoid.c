/*************************************************************************
Workspace used.: none
*************************************************************************
  serpenoid .c --

	   This program compute inverse kinematics to obtain a Serpenoid Curve         

 Date of creation: 22/11/2005

  
 Copyright  Robotics and Automation Group, Pontificia Universidad Javeriana - Cali.
	         Andr� Jaramillo Botero, ajaramil@puj.edu.co
            	Juan Manuel Florez,      jmflorez@puj.edu.co  
            	Juli� David Colorado, jdcolorado@puj.edu.co
            	Juan Camilo Acosta, jcacosta@puj.edu.co 

 See the file "license.terms" for information on usage and redistribution of this file, and for a
 DISCLAIMER OF ALL WARRANTIES.

DH: Denavit Hartenberg parameters matrix
DYN: Manipulator dynamic parameters
DYN_base: Manipulator base DYnamic parameters
Q_Qd_Qdd: joint trayectory,time vector,x(base positions),v(base velocities),a(base accelerations) 
n=# of bodies of the Serpentine
a=ondulation degree of the Curve
b=#periods per lenght unit
c=Serpentine motion bias
w=velocity propagation along all the bodies
step= step trajectory time
end_t=Final Simulation time

*************************************************************************/
// Fecha de creaci�: 22/11/2005
//16/03/2006 JCA: Added suport for base trayecory. basetray.c
//17/03/2006: Added support for base DYN parameters and corrected all identation problems.
//21/03/2006: removed x, v and a as parameters, from now included in q_qd_qdd
//		Changed time vector to first Q_Qd_Qdd column instead of last.
//21/03/2006: Changed DYN_base from VEC to MAT because the base may have more than 1 bodys
/*25/04/2006: Memory perfect. Removed Q_Qd_Qdd time vector values change*/
/*12/05/2006:  JDC, changed interface basetray prototype called; added Robot type */
/*24/05/2006: JCA: removed base tray call. Serpenoid curve complete.*/
/*31/05/2006: JCA: Arrenged time vector*/

#include "LLRL.h"
MAT *serpenoid( MAT *DH, MAT *DYN, MAT *DYN_base, MAT *Q_Qd_Qdd, Real n, Real a, Real b, Real c, Real w, Real step, Real end_t)
{
	int j,z,pt;
	Real i,k,beta,gama,alfa,points;
	Real longi=0.1,pi,n2,T,rel,tempx,tempy;
	static VEC *s=VNULL,*temp=VNULL,*temp2=VNULL,*t1=VNULL,*t2=VNULL;
	
	v_resize_vars(6,&s,&temp,&temp2,&t1,&t2,NULL);
	mem_stat_reg_vars(0,TYPE_VEC,&s,&temp,&temp2,&t1,&t2,NULL);
	//Matrix Inicialization
	//paso chambon, el que encuentre otra forma.....  dado que round no funciona ni roundf
	points=end_t/step+1;
	pi=3.14159265359;
	
	if((int)(points)<points)
		pt=ceil(end_t/step+1);
	else
		pt=floor(end_t/step+1);	
		
//******************************************************************************************
//Obtain joint positions,velocities and Accelerations
//******************************************************************************************	

	beta=b/n;
	alfa=a*sin(beta/2);
	gama=-c/n;
	T=2*pi/w;
	n2=T/step;
	rel=n2/n;
	Q_Qd_Qdd= m_resize(Q_Qd_Qdd,pt,(int)(3*(n-1))+1+18);
	
	k=0;
	for(z=0;z<pt;z++)
	{
		j=0;
		for(i=1;i<n;i++)
		{
			//Joint positions
			Q_Qd_Qdd->me[z][j+1]=2*alfa*sin(w*k+(i-1)*beta)+gama;
			//Joint velocities
			Q_Qd_Qdd->me[z][j+(int)(n-1)+1]=2*alfa*w*cos(w*k+(i-1)*beta);
			//Joint Accelerations
			Q_Qd_Qdd->me[z][j+(int)(2*(n-1))+1]=-2*alfa*pow(w,2)*sin(w*k+(i-1)*beta);
			j++;
		}
		
		//Obtain Time Vector
		if(z>0)
			Q_Qd_Qdd->me[z][0]=Q_Qd_Qdd->me[z-1][0]+step;
		else
			Q_Qd_Qdd->me[z][0]=0;
		
		k=k+step;
	}
	
	k=0;
	//Spatial Trajectory for the Serpentine: Flying Base
	for (i=0;i<=pt+rel-1;i++)
	{
		if(i<pt)
		{
			for (k=1;k<=i+1;k++)
			{
			
				//Spatial positions in Function of Time
				Q_Qd_Qdd->me[(int)i][(int)(3*(n-1)+1+3)]=-(1/n2)*cos(a*cos(k*b/n2)+(k*c)/n2)+Q_Qd_Qdd->me[(int)i][(int)(3*(n-1)+1+3)];
				
				Q_Qd_Qdd->me[(int)i][(int)(3*(n-1)+1+4)]=(1/n2)*sin(a*cos(k*b/n2)+(k*c)/n2)+Q_Qd_Qdd->me[(int)i][(int)(3*(n-1)+1+4)];
			}
			//orientation
			if(i>=rel)
				Q_Qd_Qdd->me[(int)(i-rel)][(int)(3*(n-1)+1)]=-Q_Qd_Qdd->me[(int)(i-rel)][1]+
				-atan2(
				Q_Qd_Qdd->me[(int)i][(int)(3*(n-1)+1+4)]-Q_Qd_Qdd->me[(int)(i-rel)][(int)(3*(n-1)+1+4)],
				-Q_Qd_Qdd->me[(int)i][(int)(3*(n-1)+1+3)]+Q_Qd_Qdd->me[(int)(i-rel)][(int)(3*(n-1)+1+3)]);
		}
		
		//printf("oh oh %f %f %d %f %f %d\n",n,n2,pt,step,end_t,(int)(3*(n-1)+1+3));
		tempx=0;
		tempy=0;
		if(i>=pt)
		{
			//Spatial positions in Function of Time
			for (k=1;k<=i+1;k++)
			{
				tempx=-(1/n2)*cos(a*cos(k*b/n2)+(k*c)/n2)+tempx;
				tempy=(1/n2)*sin(a*cos(k*b/n2)+(k*c)/n2)+tempy;
			}
			//orientation
			if(i>=rel)
				Q_Qd_Qdd->me[(int)(i-rel)][(int)(3*(n-1)+1)]=-Q_Qd_Qdd->me[(int)(i-rel)][1]+-atan2(
								tempy-Q_Qd_Qdd->me[(int)(i-rel)][(int)(3*(n-1)+1+4)],
								-tempx+Q_Qd_Qdd->me[(int)(i-rel)][(int)(3*(n-1)+1+3)]);
		}	
	}
	//Velocities, differentiating
	for (j=1;j<pt-1;j++)
	{
		//step=time->ve[j]-time->ve[j-1];
		if (j==1 || j==pt-2)
		{
			//v(t)=(x(t+h)-x(t-h))/(2*h)
			
			//get_col(x,j-1,temp);
			temp->ve[0]=Q_Qd_Qdd->me[j-1][(int)(3*(n-1)+1)];
			temp->ve[1]=0;
			temp->ve[2]=0;
			temp->ve[3]=Q_Qd_Qdd->me[j-1][(int)(3*(n-1)+1+3)];
			temp->ve[4]=Q_Qd_Qdd->me[j-1][(int)(3*(n-1)+1+4)];
			temp->ve[5]=0;
			
			//get_col(x,j+1,temp2);
			temp2->ve[0]=Q_Qd_Qdd->me[j+1][(int)(3*(n-1)+1)];
			temp2->ve[1]=0;
			temp2->ve[2]=0;
			temp2->ve[3]=Q_Qd_Qdd->me[j+1][(int)(3*(n-1)+1+3)];
			temp2->ve[4]=Q_Qd_Qdd->me[j+1][(int)(3*(n-1)+1+4)];
			temp2->ve[5]=0;
			
			v_sub(temp2,temp,temp2);
			sv_mlt(1/(2*step),temp2,temp);
			//_set_col(v,j-1,temp,0);
			Q_Qd_Qdd->me[j][(int)(3*(n-1)+1+6)]=temp->ve[0];
			Q_Qd_Qdd->me[j][(int)(3*(n-1)+1+7)]=temp->ve[1];
			Q_Qd_Qdd->me[j][(int)(3*(n-1)+1+8)]=temp->ve[2];
			Q_Qd_Qdd->me[j][(int)(3*(n-1)+1+9)]=temp->ve[3];
			Q_Qd_Qdd->me[j][(int)(3*(n-1)+1+10)]=temp->ve[4];
			Q_Qd_Qdd->me[j][(int)(3*(n-1)+1+11)]=temp->ve[5];
		}
		else
		{
			//v(t)=(-x(t+2h)+8x(t+h)-8x(t-h)+x(t-2h))/(12*h)
			
			//get_col(x,j+1,temp);
			temp->ve[0]=Q_Qd_Qdd->me[j+1][(int)(3*(n-1)+1)];
			temp->ve[1]=0;
			temp->ve[2]=0;
			temp->ve[3]=Q_Qd_Qdd->me[j+1][(int)(3*(n-1)+1+3)];
			temp->ve[4]=Q_Qd_Qdd->me[j+1][(int)(3*(n-1)+1+4)];
			temp->ve[5]=0;
			
			//get_col(x,j+2,temp2);
			temp2->ve[0]=Q_Qd_Qdd->me[j+2][(int)(3*(n-1)+1)];
			temp2->ve[1]=0;
			temp2->ve[2]=0;
			temp2->ve[3]=Q_Qd_Qdd->me[j+2][(int)(3*(n-1)+1+3)];
			temp2->ve[4]=Q_Qd_Qdd->me[j+2][(int)(3*(n-1)+1+4)];
			temp2->ve[5]=0;
			
			//get_col(x,j-1,t1);
			t1->ve[0]=Q_Qd_Qdd->me[j-1][(int)(3*(n-1)+1)];
			t1->ve[1]=0;
			t1->ve[2]=0;
			t1->ve[3]=Q_Qd_Qdd->me[j-1][(int)(3*(n-1)+1+3)];
			t1->ve[4]=Q_Qd_Qdd->me[j-1][(int)(3*(n-1)+1+4)];
			t1->ve[5]=0;
			
			//get_col(x,j-2,t2);
			t2->ve[0]=Q_Qd_Qdd->me[j-2][(int)(3*(n-1)+1)];
			t2->ve[1]=0;
			t2->ve[2]=0;
			t2->ve[3]=Q_Qd_Qdd->me[j-2][(int)(3*(n-1)+1+3)];
			t2->ve[4]=Q_Qd_Qdd->me[j-2][(int)(3*(n-1)+1+4)];
			t2->ve[5]=0;
			
			sv_mlt(8,temp,temp);
			sv_mlt(-1,temp2,temp2);
			sv_mlt(-8,t1,t1);
			v_add(temp,temp2,temp2);
			v_add(temp2,t1,t1);
			v_add(t1,t2,t2);
			sv_mlt(1/(12*step),t2,temp);
			
			//_set_col(v,j-1,temp,0);
			Q_Qd_Qdd->me[j][(int)(3*(n-1)+1+6)]=temp->ve[0];
			Q_Qd_Qdd->me[j][(int)(3*(n-1)+1+7)]=temp->ve[1];
			Q_Qd_Qdd->me[j][(int)(3*(n-1)+1+8)]=temp->ve[2];
			Q_Qd_Qdd->me[j][(int)(3*(n-1)+1+9)]=temp->ve[3];
			Q_Qd_Qdd->me[j][(int)(3*(n-1)+1+10)]=temp->ve[4];
			Q_Qd_Qdd->me[j][(int)(3*(n-1)+1+11)]=temp->ve[5];
		}
	}
	//Accelerations (differentiating.)
	for (j=1;j<pt-1;j++)
	{
		//step=time->ve[j]-time->ve[j-1];
		if (j==1 || j==pt-2)
		{
			//a(t)=(x(t+h)-2*x(t)+x(t-h))/(h²)
			
			//get_col(x,j-1,temp);
			temp->ve[0]=Q_Qd_Qdd->me[j-1][(int)(3*(n-1)+1)];
			temp->ve[1]=0;
			temp->ve[2]=0;
			temp->ve[3]=Q_Qd_Qdd->me[j-1][(int)(3*(n-1)+1+3)];
			temp->ve[4]=Q_Qd_Qdd->me[j-1][(int)(3*(n-1)+1+4)];
			temp->ve[5]=0;
			
			//get_col(x,j+1,temp2);
			temp2->ve[0]=Q_Qd_Qdd->me[j+1][(int)(3*(n-1)+1)];
			temp2->ve[1]=0;
			temp2->ve[2]=0;
			temp2->ve[3]=Q_Qd_Qdd->me[j+1][(int)(3*(n-1)+1+3)];
			temp2->ve[4]=Q_Qd_Qdd->me[j+1][(int)(3*(n-1)+1+4)];
			temp2->ve[5]=0;
			
			//get_col(x,j,t1);
			t1->ve[0]=Q_Qd_Qdd->me[j][(int)(3*(n-1)+1)];
			t1->ve[1]=0;
			t1->ve[2]=0;
			t1->ve[3]=Q_Qd_Qdd->me[j][(int)(3*(n-1)+1+3)];
			t1->ve[4]=Q_Qd_Qdd->me[j][(int)(3*(n-1)+1+4)];
			t1->ve[5]=0;
				
			sv_mlt(-2,t1,t1);
			v_add(temp2,temp,temp2);
			v_add(temp2,t1,t2);
			sv_mlt(1/(pow(step,2)),t2,temp);

			//_set_col(a,j-1,temp,0);
			Q_Qd_Qdd->me[j][(int)(3*(n-1)+1+12)]=temp->ve[0];
			Q_Qd_Qdd->me[j][(int)(3*(n-1)+1+13)]=temp->ve[1];
			Q_Qd_Qdd->me[j][(int)(3*(n-1)+1+14)]=temp->ve[2];
			Q_Qd_Qdd->me[j][(int)(3*(n-1)+1+15)]=temp->ve[3];
			Q_Qd_Qdd->me[j][(int)(3*(n-1)+1+16)]=temp->ve[4];
			Q_Qd_Qdd->me[j][(int)(3*(n-1)+1+17)]=temp->ve[5];
		}
		else
		{
			//a(t)=(-x(t+2h)+16x(t+h)-30x(t)+16x(t-h)-x(t-2h))/(12*h²)
			//get_col(x,j+1,temp);
			temp->ve[0]=Q_Qd_Qdd->me[j+1][(int)(3*(n-1)+1)];
			temp->ve[1]=0;
			temp->ve[2]=0;
			temp->ve[3]=Q_Qd_Qdd->me[j+1][(int)(3*(n-1)+1+3)];
			temp->ve[4]=Q_Qd_Qdd->me[j+1][(int)(3*(n-1)+1+4)];
			temp->ve[5]=0;
			
			//get_col(x,j+2,temp2);
			temp2->ve[0]=Q_Qd_Qdd->me[j+2][(int)(3*(n-1)+1)];
			temp2->ve[1]=0;
			temp2->ve[2]=0;
			temp2->ve[3]=Q_Qd_Qdd->me[j+2][(int)(3*(n-1)+1+3)];
			temp2->ve[4]=Q_Qd_Qdd->me[j+2][(int)(3*(n-1)+1+4)];
			temp2->ve[5]=0;
			
			//get_col(x,j-1,t1);
			t1->ve[0]=Q_Qd_Qdd->me[j-1][(int)(3*(n-1)+1)];
			t1->ve[1]=0;
			t1->ve[2]=0;
			t1->ve[3]=Q_Qd_Qdd->me[j-1][(int)(3*(n-1)+1+3)];
			t1->ve[4]=Q_Qd_Qdd->me[j-1][(int)(3*(n-1)+1+4)];
			t1->ve[5]=0;
			
			//get_col(x,j-2,t2);
			t2->ve[0]=Q_Qd_Qdd->me[j-1][(int)(3*(n-1)+1)];
			t2->ve[1]=0;
			t2->ve[2]=0;
			t2->ve[3]=Q_Qd_Qdd->me[j-1][(int)(3*(n-1)+1+3)];
			t2->ve[4]=Q_Qd_Qdd->me[j-1][(int)(3*(n-1)+1+4)];
			t2->ve[5]=0;
			
			//get_col(x,j,s);
			s->ve[0]=Q_Qd_Qdd->me[j][(int)(3*(n-1)+1)];
			s->ve[1]=0;
			s->ve[2]=0;
			s->ve[3]=Q_Qd_Qdd->me[j][(int)(3*(n-1)+1+3)];
			s->ve[4]=Q_Qd_Qdd->me[j][(int)(3*(n-1)+1+4)];
			s->ve[5]=0;

			sv_mlt(16,temp,temp);
			sv_mlt(-1,temp2,temp2);
			sv_mlt(-1,t2,t2);
			sv_mlt(16,t1,t1);
			sv_mlt(-30,s,s);
			v_add(temp,temp2,temp2);
			v_add(temp2,t2,t2);
			v_add(t2,t1,t1);
			v_add(t1,s,t2);
			sv_mlt(1/(12*pow(step,2)),t2,temp);
			
			//_set_col(a,j-1,temp,0);
			Q_Qd_Qdd->me[j][(int)(3*(n-1)+1+12)]=temp->ve[0];
			Q_Qd_Qdd->me[j][(int)(3*(n-1)+1+13)]=temp->ve[1];
			Q_Qd_Qdd->me[j][(int)(3*(n-1)+1+14)]=temp->ve[2];
			Q_Qd_Qdd->me[j][(int)(3*(n-1)+1+15)]=temp->ve[3];
			Q_Qd_Qdd->me[j][(int)(3*(n-1)+1+16)]=temp->ve[4];
			Q_Qd_Qdd->me[j][(int)(3*(n-1)+1+17)]=temp->ve[5];
		}
	}
	//removing first and last point
	if(Q_Qd_Qdd->m!=1)
	{
		m_move(Q_Qd_Qdd,1,0,points-1,Q_Qd_Qdd->n,Q_Qd_Qdd,0,0);
		Q_Qd_Qdd=m_resize(Q_Qd_Qdd,points-2,Q_Qd_Qdd->n);
	}
	
	//arrenging time vector due to first and last point removal.	
	for (z=points-3;z>=0;z--)
		Q_Qd_Qdd->me[z][0]=Q_Qd_Qdd->me[z][0]-Q_Qd_Qdd->me[0][0];
	return(Q_Qd_Qdd);
}
