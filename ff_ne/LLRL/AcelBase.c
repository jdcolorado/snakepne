#include "LLRL.h"
/*************************************************************************
Workspaces used: none
*************************************************************************
Compute base Acceleration  due to the Spatial Force of the base

void AcelBase(DH,DYN,Q_Qd_Qdd,Fric,Ab)

DH: Denavit Hartenberg parameters matrix
DYN: Manipulator dynamic parameters
Q_Qd_Qdd: joint trayectory,time,x(base positions),v(base velocities),a(base accelerations)
Fbase: Spatial base Force

Robotics and Automation Group
	Andres Jaramillo Botero
	Juan Camilo Acosta
	Julian David Colorado
	Antonio Matta

Pontificia Universidad Javeriana, Cali.

09/04/06   JDC: Compute base Acceleration  due to the Spatial Force of the base
26/04/06   JDC: Changed name of Fric for Fbase
26/04/06   JDC: Added robot type: 1)Fixed base 2)Flying base 3)Floating base 4)Serpentine Robots
02/05/06   JDC: None gravity in base acceleration
07/05/2006 JDC, Changed interface prototype, added DYN_base 
07/05/2006 JDC, Added inertia of the base for the total mass of the system
30/05/2006 JDC, Orienting MASS operator to mass center frame
***************************************************************************/

//Compute base Acceleration  due to the Spatial Force of the base
MAT *AcelBase(int robot_type,MAT *DH, MAT *DYN,MAT *DYN_base,VEC *grav, MAT *Q_Qd_Qdd,MAT* Fbase, MAT* Ab)
{
  
	static MAT *Q=MNULL, *r=MNULL, *Rot=MNULL, *R=MNULL, *P=MNULL, *SkewP=MNULL, *sk=MNULL, *sks=MNULL, *S=MNULL, *Ticm=MNULL, *Icm=MNULL, *I=MNULL, *identidad=MNULL, *P_trans=MNULL, *R_trans=MNULL, *PtRt=MNULL, *r_trans=MNULL, *rTicm=MNULL, *Ji=MNULL, *massU=MNULL, *S_trans=MNULL, *SIcm=MNULL, *RP=MNULL, *Tinertia=MNULL, *Iant=MNULL, *RPIant=MNULL, *Iup=MNULL, *invmass=MNULL, *X=MNULL, *RT=MNULL,*Tmass=MNULL,*RPTinertia=MNULL,*Ine=MNULL;
	
	static VEC *p=VNULL, *s=VNULL, *ss=VNULL, *columna=VNULL, *resul=VNULL;

	int k,n,np,i,cont;
	float alfa,a,teta,d;
	n=(int) DH->m;	        		/* number of joints:DOF  */
	np=(int) Q_Qd_Qdd->m;			/* number of trajectory points */
	Ab = m_resize(Ab,6,np);
	Q=m_resize(Q,np,n);
	MEM_STAT_REG(Q,TYPE_MAT);
	m_move(Q_Qd_Qdd,0,0,np,n,Q,0,0);

	m_resize_vars(6,6,&R,&P,&S,&Icm,&I,&P_trans,&R_trans,&PtRt,&SIcm,&RP,&S_trans,&Tinertia,&Iant,&RPIant,&Iup,&invmass,&RT,&Tmass,&RPTinertia,&Ine,NULL);
	mem_stat_reg_vars(0,TYPE_MAT,&R,&P,&S,&Icm,&I,&P_trans,&R_trans,&PtRt,&SIcm,&RP,&S_trans,&Tinertia,&Iant,&RPIant,&Iup,&invmass,&RT,&Tmass,&RPTinertia,&Ine,NULL);
	m_resize_vars(3,3,&r,&SkewP,&sk,&sks,&Ticm,&identidad,&r_trans,&rTicm,&Ji,&massU,NULL);
	mem_stat_reg_vars(0,TYPE_MAT,&r,&SkewP,&sk,&sks,&Ticm,&identidad,&r_trans,&rTicm,&Ji,&massU,NULL);
	v_resize_vars(3,&p,&s,&ss,NULL);
	mem_stat_reg_vars(0,TYPE_VEC,&p,&s,&ss,NULL);
	v_resize_vars(6,&columna,&resul,NULL);
	mem_stat_reg_vars(0,TYPE_VEC,&columna,&resul,NULL);
	Rot=m_resize(Rot,4,4);
	MEM_STAT_REG(Rot,TYPE_MAT);
	m_ident(identidad);

	//Obtain Spatial base Position respect to Mass Center frame
	if(robot_type==4)
	{
		X=m_resize(X,6,1);
		mem_stat_reg_vars(0,TYPE_MAT,&X,NULL);
		X=base_x(DH,DYN,DYN_base,Q_Qd_Qdd,X);	
	}
	
	k=0;
	//for(k=0;k<np;k++)
	//{
	//Backward recurrence propagating Inertia parameters
		for(i=(n-1);i>=0;i--)
		{
			m_zero(R);	//Reset values
			m_zero(P);
			m_zero(r);
			v_zero(p);
			m_zero(Rot);
			m_zero(I);
			m_zero(S);
			alfa = DH->me[i][0];
			a=DH->me[i][1];
			if(DH->me[i][4]==0)        //Rotational
			{
				teta=Q->me[k][i];
				d=DH->me[i][3];
			}
			else                       //Prismatic
			{
				teta=DH->me[i][2];
				d=Q->me[k][i];
			}
			homogeneus(alfa,a,teta,d,Rot);	//obtain homogeneus transformation
			Basic_rotation(Rot,r);		//Obtain basic matrix Rotation
			p->ve[0]=-a;
			p->ve[1]=-d*sin(alfa);
			p->ve[2]=-d*cos(alfa);

			//Obtain the 6-dimensional R operator based on r
			form_matrix(0,0,r,R);
			form_matrix(3,3,r,R);
					
			//Obtain the 6-dimensional P operator based on p
			skew_symetric(SkewP,p);
			form_matrix(0,0,identidad,P);
			form_matrix(3,3,identidad,P);
			sm_mlt(-1,SkewP,SkewP);
			form_matrix(0,3,SkewP,P);  
	
			//Obtain vector s from frame i+1 to center of mass of body i
			s->ve[0]=DYN->me[i][1];
			s->ve[1]=DYN->me[i][2];
			s->ve[2]=DYN->me[i][3];
			v_sub(s,p,s);
			mv_mlt(r,s,ss);
			skew_symetric(sks,ss);    //Obtain Skew symetric of ss

			//obtain 6-dimensional operator S
			form_matrix(0,0,identidad,S);
			form_matrix(3,3,identidad,S);
			form_matrix(0,3,sks,S);
	
			//Obtain de Inertial Tensor referred to center of gravity of body i
			Ticm->me[0][0]=DYN->me[i][4];
			Ticm->me[0][1]=DYN->me[i][7];
			Ticm->me[0][2]=DYN->me[i][9];
			Ticm->me[1][0]=DYN->me[i][7];
			Ticm->me[1][1]=DYN->me[i][5];
			Ticm->me[1][2]=DYN->me[i][8];
			Ticm->me[2][0]=DYN->me[i][9];
			Ticm->me[2][1]=DYN->me[i][8];
			Ticm->me[2][2]=DYN->me[i][6];
		
			//Obtain 6-dimensional Inertial operator in center of gravity:
			m_mlt(r,Ticm,rTicm);
			m_transp(r,r_trans);
			m_mlt(rTicm,r_trans,Ji);
			form_matrix(0,0,Ji,Icm);
			sm_mlt(DYN->me[i][0],identidad,massU);
			form_matrix(3,3,massU,Icm);
			
			//applying parallel axes theorem to obtain the inertial operator at frame i
			m_transp(S,S_trans);
			m_mlt(S,Icm,SIcm);
			m_mlt(SIcm,S_trans,I);
			//************************************************************************************
			//Add Inertia of all the system at frame i
			//************************************************************************************
			if (i==(n-1))
			{
				m_copy(I,Tinertia);
				m_copy(Tinertia,Iant);
					
			}
			else
			{
				m_mlt(R,P,RP);
				m_transp(R,R_trans);
				m_transp(P,P_trans);
				m_mlt(RP,Iant,RPIant);
				m_mlt(P_trans,R_trans,PtRt);
				m_mlt(RPIant,PtRt,Iup);
				m_add(I,Iup,Tinertia);
				m_copy(Tinertia,Iant);
				
				if (i==0)
				{
					//Obtain de Inertial Tensor of base
					Ticm->me[0][0]=DYN_base->me[0][4];
					Ticm->me[0][1]=DYN_base->me[0][7];
					Ticm->me[0][2]=DYN_base->me[0][9];
					Ticm->me[1][0]=DYN_base->me[0][7];
					Ticm->me[1][1]=DYN_base->me[0][5];
					Ticm->me[1][2]=DYN_base->me[0][8];
					Ticm->me[2][0]=DYN_base->me[0][9];
					Ticm->me[2][1]=DYN_base->me[0][8];
					Ticm->me[2][2]=DYN_base->me[0][6];
					form_matrix(0,0,Ticm,Icm);
					sm_mlt(DYN_base->me[0][0],identidad,massU);
					form_matrix(3,3,massU,Icm);
					//Obtain vector s 
					s->ve[0]=DYN_base->me[0][1];
					s->ve[1]=DYN_base->me[0][2];
					s->ve[2]=DYN_base->me[0][3];
					skew_symetric(sks,s);    //Obtain Skew symetric of s
					//obtain 6-dimensional operator S
					form_matrix(0,0,identidad,S);
					form_matrix(3,3,identidad,S);
					form_matrix(0,3,sks,S);
					m_transp(S,S_trans);
					m_mlt(S,Icm,SIcm);
					m_mlt(SIcm,S_trans,I);
					//add inertia of base
					m_add(Tinertia,I,Tinertia);
					//Inverting Mass operator
					pinv(Tinertia,invmass);

					//For Serpentine Robots
					if(robot_type==4)
					{
						//Referred at mass center frame
						get_col(X,k,columna);
						euler(columna,r,p);
						//Obtain the 6-dimensional R operator based on r
						form_matrix(0,0,r,R);
						form_matrix(3,3,r,R);
						//Obtain the 6-dimensional P operator based on p
						skew_symetric(SkewP,p);
						form_matrix(0,0,identidad,P);
						form_matrix(3,3,identidad,P);
						sm_mlt(-1,SkewP,SkewP);
						form_matrix(0,3,SkewP,P);
						m_mlt(R,P,RP);
						m_transp(R,R_trans);
						m_transp(P,P_trans);
						m_mlt(RP,Tinertia,RPTinertia);
						m_mlt(P_trans,R_trans,PtRt);
						m_mlt(RPTinertia,PtRt,Ine);
						//Orienting to mass center frame
						sm_mlt(-1,r,r);
						form_matrix(0,0,r,R);
						form_matrix(3,3,r,R);
						m_transp(R,R_trans);
						m_mlt(R,Tinertia,RT);
						m_mlt(RT,R_trans,Tmass);
						//Inverting Mass operator
						pinv(Tmass,invmass);
					}
					for(cont=0;cont<=5;cont++)
						columna->ve[cont]=Fbase->me[cont][k];
					mv_mlt(invmass,columna,resul);
					Ab->me[0][k]=resul->ve[2];
					Ab->me[1][k]=resul->ve[1];
					Ab->me[2][k]=resul->ve[0];
					Ab->me[3][k]=resul->ve[3];
					Ab->me[4][k]=resul->ve[4];
					Ab->me[5][k]=resul->ve[5];
					//for Serpentine robots
					if (robot_type==4)
						Ab->me[5][k]=0;
										
				}
			}
		}
	//}
	return Ab;
}
