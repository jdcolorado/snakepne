#include "LLRL.h"
/*************************************************************************
Workspaces used: none
*************************************************************************
Compute the friction force with the surface referred at base frame

MAT* Friction_surface(MAT *DH, MAT *DYN,MAT *Q_Qd_Qdd,MAT* Vb)

DH: Denavit Hartenberg parameters matrix
DYN: Manipulator dynamic parameters
Q_Qd_Qdd: joint trayectory,time,x(base positions),v(base velocities),a(base accelerations)
Vb: Spatial base velocity

Robotics and Automation Group
	Andres Jaramillo Botero
	Juan Camilo Acosta
	Julian David Colorado
	Antonio Matta

Pontificia Universidad Javeriana, Cali.

29/04/06   JDC: Compute the friction force with the surface referred at base frame
/*07/05/2006 JDC, Changed interface prototype, added DYN_base */
/*07/05/2006 JDC, Added base's surface friction for the total friction of the systems (Serpentine) 
/*19/05/2006 JDC, Added base velocity for initial conditions (Compute spatial vel)*/
/*19/05/2006 JDC, Elimination of base friction*/
/*30/05/2006 JDC, Interface Prototype has changed: Elimination of MAT Xb */
/*30/05/2006 JDC, New Call for base_x.c to obtain spatial base position Xb*/
/*30/05/2006 JDC, Orienting Friction to mass center frame*/
/*2/06/2006 JDC,  Interface has changed; added Xb*/

/***************************************************************************/
MAT *Friction_surface(MAT *DH, MAT *DYN,MAT *DYN_base,MAT *Q_Qd_Qdd,MAT *Xb, MAT *Vb, MAT *Fric)
{
	static MAT *Q=MNULL, *dQ=MNULL, *d2Q=MNULL, *V=MNULL, *dV=MNULL, *r=MNULL, *Rot=MNULL, *R=MNULL, *P=MNULL, *SkewP=MNULL, *sks=MNULL, *S=MNULL,*identidad=MNULL, *P_trans=MNULL, *R_trans=MNULL, *PtRt=MNULL, *r_trans=MNULL, *S_trans=MNULL, *RP=MNULL, *vf=MNULL, *C=MNULL, *vf_trans=MNULL, *massvf=MNULL, *MC=MNULL, *Fr1=MNULL, *Xcm=MNULL, *Vcm=MNULL, *X=MNULL,*DHtemp=MNULL,*qtemp=MNULL,*T=MNULL,*Ts=MNULL;

	static VEC *H=VNULL, *p=VNULL, *s=VNULL, *ss=VNULL, *Vant=VNULL, *columna=VNULL, *dQH=VNULL,*RPFriction=VNULL,*Ang=VNULL, *PtRtVant=VNULL,*Vicm=VNULL,*vxy=VNULL,*Fr2=VNULL,*Fr6=VNULL,*Frt=VNULL,*Friction=VNULL,*Frtant=VNULL,*RPFrtant=VNULL, *Frtotal=VNULL;
	
	int k,n,np,i,cont;
	float alfa,a,teta,d,massneg,tauf;
	n=(int) DH->m;	        		/* number of joints:DOF  */
	np=(int) Q_Qd_Qdd->m;			/* number of trajectory points */
	Fric = m_resize(Fric,6,1);
	m_resize_vars(np,n,&Q,&dQ,&d2Q,NULL);
	mem_stat_reg_vars(0,TYPE_MAT,&Q,&dQ,&d2Q,NULL);
	m_move(Q_Qd_Qdd,0,0,np,n,Q,0,0);
	m_move(Q_Qd_Qdd,0,n,np,n,dQ,0,0);
	
/**********************MATRIX DECLARATION AND REGISTRATION*******************************************/	
	m_resize_vars(6,n+1,&V,&dV,NULL);
	mem_stat_reg_vars(0,TYPE_MAT,&V,&dV,NULL);
	m_resize_vars(2,2,&vf,&C,&vf_trans,&massvf,&MC,&Fr1,NULL);
	mem_stat_reg_vars(0,TYPE_MAT,&vf,&C,&vf_trans,&massvf,&MC,&Fr1,NULL);
	m_resize_vars(4,4,&Rot,&T,&Ts,NULL);
	mem_stat_reg_vars(0,TYPE_MAT,&Rot,&T,&Ts,NULL);
	
	m_resize_vars(6,6,&R,&P,&S,&P_trans,&R_trans,&PtRt,&RP,&S_trans,NULL);
	mem_stat_reg_vars(0,TYPE_MAT,&R,&P,&S,&P_trans,&R_trans,&PtRt,&RP,&S_trans,NULL);
	
	m_resize_vars(3,3,&r,&SkewP,&sks,&identidad,&r_trans,NULL);
	mem_stat_reg_vars(0,TYPE_MAT,&r,&SkewP,&sks,&identidad,&r_trans,NULL);
/**********************VECTOR DECLARATION AND REGISTRATION*******************************************/	
	v_resize_vars(3,&p,&s,&ss,NULL);
	mem_stat_reg_vars(0,TYPE_VEC,&p,&s,&ss,NULL);
	
	v_resize_vars(2,&vxy,&Fr2,NULL);
	mem_stat_reg_vars(0,TYPE_VEC,&vxy,&Fr2,NULL);
	
	v_resize_vars(6,&H,&Vant,&columna,&dQH,&PtRtVant,&Vicm,&Fr6,&Frt,&Friction,&Frtant,&RPFrtant,&Frtotal,&RPFriction,&Ang,NULL);
	mem_stat_reg_vars(0,TYPE_VEC,&H,&Vant,&columna,&dQH,&PtRtVant,&Vicm,&Fr6,&Frt,&Friction,&Frtant,&RPFrtant,&Frtotal,&RPFriction,&Ang,NULL);
		
	k=0;

	//Obtain Spatial base Position respect to Mass Center frame
	Xcm=m_resize(Xcm,6,1);
	X=m_resize(X,6,1);
	mem_stat_reg_vars(0,TYPE_MAT,&Xcm,&X,NULL);
	Xcm=base_x(DH,DYN,DYN_base,Q_Qd_Qdd,Xcm);	

	Vcm=m_resize(Vcm,6,1);
	mem_stat_reg_vars(0,TYPE_MAT,&Vcm,NULL);
	Vcm=base_v(DH,DYN,DYN_base,Q_Qd_Qdd,Xcm,Vcm);

	m_add(Xb,Xcm,X);
	m_add(Vb,Vcm,Vb);
	//obtain column from Vb-dVb (Base trajectory)          
	Vant->ve[0]=Vb->me[2][k];
	Vant->ve[1]=Vb->me[1][k];
	Vant->ve[2]=Vb->me[0][k];
	Vant->ve[3]=Vb->me[3][k];
	Vant->ve[4]=Vb->me[4][k];
	Vant->ve[5]=Vb->me[5][k];
	
	for(cont=0;cont<=5;cont++)
		Vb->me[cont][k]=Vant->ve[cont];
	m_ident(identidad);
	
	//Computing spatial velocities of each link	
	for(i=0;i<n;i++)
	{
		m_zero(R);	//Reset values
		m_zero(P);
		m_zero(r);
		v_zero(p);
		m_zero(Rot);
		v_zero(H);
		if(DH->me[i][4]==0)
			H->ve[2]=1;               //projection Vector of rotational motion
		else
			H->ve[5]=1;               //projection Vector of Prismatic motion
		if(i==0)
		{
			//Obtain 6-Dimensional Position from X   
			get_col(X,k,columna);
			euler(columna,r,p);	//Obtain Rotation r and Position p	
		}
		else
		{
			alfa = DH->me[i-1][0];
			a=DH->me[i-1][1];
			if(DH->me[i-1][4]==0)        //Rotational
			{
				teta=Q->me[k][i-1];
				d=DH->me[i-1][3];
			}
			else                         //Prismatic
			{
				teta=DH->me[i-1][2];
				d=Q->me[k][i-1];
			}
			
			homogeneus(alfa,a,teta,d,Rot);	//obtain homogeneus transformation
			Basic_rotation(Rot,r);	  	//obtain rotation matrix
			p->ve[0]=-a;			//obtain the inverse of position p
			p->ve[1]=-d*sin(alfa);
			p->ve[2]=-d*cos(alfa);
		}
		//Obtain the 6-dimensional R operator based on r
   	
		form_matrix(0,0,r,R);
		form_matrix(3,3,r,R);
//*********************************************************************************
// Computing Spatial Velocities
//*********************************************************************************
		//Obtain the 6-dimensional P operator based on p
		skew_symetric(SkewP,p);
		form_matrix(0,0,identidad,P);
		form_matrix(3,3,identidad,P);
		sm_mlt(-1,SkewP,SkewP);
		form_matrix(0,3,SkewP,P);

		//Obtain 6-dimensional velocity: Vel
		sv_mlt(dQ->me[k][i],H,dQH);
		m_transp(P,P_trans);
		m_transp(R,R_trans);
		m_mlt(P_trans,R_trans,PtRt);
		
	
		if (i==0)	
		   mv_mlt(R_trans,Vant,PtRtVant);
		else
		    mv_mlt(PtRt,Vant,PtRtVant);		
		
	 	v_add(dQH,PtRtVant,Vant);
	
		//spatial velocities for each joint
		_set_col(V,i,Vant,0);
	}

	//Computing Friction Surface	
	for(i=(n-1);i>=0;i--)
	{
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
		//Tangencial and Normal Componentes of Friction Force
		get_col(X,k,columna);
		euler(columna,r,p);	//Obtain base homogeneus transformation
		Rot->me[0][0]=r->me[0][0];
		Rot->me[1][0]=r->me[1][0];
		Rot->me[2][0]=r->me[2][0];
		Rot->me[3][0]=0;
	
		Rot->me[0][1]=r->me[0][1];
		Rot->me[1][1]=r->me[1][1];
		Rot->me[2][1]=r->me[2][1];
		Rot->me[3][1]=0;
	
		Rot->me[0][2]=r->me[0][2];
		Rot->me[1][2]=r->me[1][2];
		Rot->me[2][2]=r->me[2][2];
		Rot->me[3][2]=0;
	
		Rot->me[0][3]=p->ve[0];
		Rot->me[1][3]=p->ve[1];
		Rot->me[2][3]=p->ve[2];
		Rot->me[3][3]=1;
		
		DHtemp=m_resize(DHtemp,i,DH->n);
		qtemp=m_resize(qtemp,1,i);
		mem_stat_reg_vars(0,TYPE_MAT,&DHtemp,&qtemp,NULL);
		m_move(DH,0,0,i,DH->n,DHtemp,0,0);
		m_move(Q,0,0,1,i,qtemp,0,0);
		fkine(DHtemp,qtemp,T);
		m_mlt(Rot,T,Ts);
		Ang=tr2rpy(Ts,Ang);
		vf->me[0][0]=cos(Ang->ve[2]);		
		vf->me[0][1]=sin(Ang->ve[2]);
		vf->me[1][0]=-sin(Ang->ve[2]);
		vf->me[1][1]=cos(Ang->ve[2]);
		
		//obtain homogeneus transformation		
		homogeneus(alfa,a,teta,d,Rot);	
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
//*************************************************************************************
//Obtain Friction Force with Contact Surface (for Serpentine Robots)
//*************************************************************************************
		//Matrix of Coefficients of friction
		C->me[0][0]=DYN->me[i][10];
		C->me[0][1]=0;
		C->me[1][0]=0;
		C->me[1][1]=DYN->me[i][11];
		m_transp(vf,vf_trans);
		massneg=-DYN->me[i][0];
		sm_mlt(massneg,vf_trans,massvf);
		m_mlt(massvf,C,MC);
		m_mlt(MC,vf,Fr1);
		
		//Projecting velocity from frame i to center of gravity
		get_col(V,i,columna);   	  //obtain column from spatial Velocity
		m_transp(S,S_trans);
		mv_mlt(S_trans,columna,Vicm);
		vxy->ve[0]=Vicm->ve[3];
		vxy->ve[1]=Vicm->ve[4];
		mv_mlt(Fr1,vxy,Fr2);	
		//adding friction Torque
		tauf=-(0.3333)*DYN->me[i][11]*DYN->me[i][0]*0.05*0.05*dQ->me[k][i];	
		//6-dimensional friction force
		Fr6->ve[0]=0;
		Fr6->ve[1]=0;
		Fr6->ve[2]=tauf;
		Fr6->ve[3]=Fr2->ve[0];
		Fr6->ve[4]=Fr2->ve[1];
		Fr6->ve[5]=0;	 	
		
		//Projecting force at frame i, 6-dimensional
		mv_mlt(S,Fr6,Frt);
		//Projecting and adding Friction force of each link to first frame 		
		if(i==(n-1))
		{
			v_copy(Frt,Friction);
			v_copy(Friction,Frtant);
		}
		else
		{	
			m_mlt(R,P,RP);
			mv_mlt(RP,Frtant,RPFrtant);
			v_add(Frt,RPFrtant,Friction);
			v_copy(Friction,Frtant);
			//add friction component due to Serpentine base
			if (i==0)
			{
				vf->me[0][0]=cos(X->me[0][k]);		
				vf->me[0][1]=sin(X->me[0][k]);
				vf->me[1][0]=-sin(X->me[0][k]);
				vf->me[1][1]=cos(X->me[0][k]);
				C->me[0][0]=DYN_base->me[0][10];
				C->me[0][1]=0;
				C->me[1][0]=0;
				C->me[1][1]=DYN_base->me[0][11];
                        	m_transp(vf,vf_trans);
				massneg=-DYN_base->me[0][0];
				sm_mlt(massneg,vf_trans,massvf);
				m_mlt(massvf,C,MC);
				m_mlt(MC,vf,Fr1);
				//Obtain vector s from frame i+1 to center of mass of body i
				s->ve[0]=DYN_base->me[i][1];
				s->ve[1]=DYN_base->me[i][2];
				s->ve[2]=DYN_base->me[i][3];
				skew_symetric(sks,s);   	 //Obtain Skew symetric of s
				//obtain 6-dimensional operator S
				form_matrix(0,0,identidad,S);
				form_matrix(3,3,identidad,S);
				form_matrix(0,3,sks,S);
				get_col(Vb,k,columna);   	 //obtain column from spatial Velocity
				m_transp(S,S_trans);
				mv_mlt(S_trans,columna,Vicm);
				vxy->ve[0]=Vicm->ve[3];
				vxy->ve[1]=Vicm->ve[4];
				mv_mlt(Fr1,vxy,Fr2);	
				//adding friction Torque of base
				tauf=-(0.3333)*DYN_base->me[i][11]*DYN_base->me[i][0]*0.05*0.05*Vicm->ve[2];
				Fr6->ve[0]=0;
				Fr6->ve[1]=0;
				Fr6->ve[2]=tauf;
				Fr6->ve[3]=Fr2->ve[0];
				Fr6->ve[4]=Fr2->ve[1];
				Fr6->ve[5]=0;	
				mv_mlt(S,Fr6,Frtotal);
				//add Friction of base
				v_add(Friction,Frtotal,Friction);
				//Referred at mass center frame
				get_col(Xcm,k,columna);
				euler(columna,r,p);	
				//Obtain the 6-dimensional R operator based on r
				form_matrix(0,0,r,R);
				form_matrix(3,3,r,R);
				//Obtain the 6-dimensional P operator based on p
				skew_symetric(SkewP,p);
				form_matrix(0,0,identidad,P);
				form_matrix(3,3,identidad,P);
				//sm_mlt(-1,SkewP,SkewP);
				form_matrix(0,3,SkewP,P);
				m_mlt(R,P,RP);
				mv_mlt(RP,Friction,RPFriction);
				//Orienting to mass center frame
				sm_mlt(-1,r,r);
				form_matrix(0,0,r,R);
				form_matrix(3,3,r,R);
				mv_mlt(R,RPFriction,columna);
				for(cont=0;cont<=5;cont++)
					Fric->me[cont][0]=columna->ve[cont];
			}

					
		}
	}
return Fric;
}
