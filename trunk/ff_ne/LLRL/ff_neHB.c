
 /*****************************************************************************
 Compute inverse dynamics via recursive Newton-Euler formulation for
 a free- Fixed/flying/floating base systems.
 
 TAU = ff_ne(robot_type,DH,DYN,DYN_base,Q_Qd_Qdd,grav,fext, X, V, dV,Fbase)
 
 robot_type= 1)Fixed base 2)Flying base 3)Floating base 4)Serpentine Robots
 DH=Denavit Hartenberg parameter matrix for manipulator
 DYN=dynamic parameter matrix in the form [m,sx,sy,sz,Ixx,Iyy,Izz,Ixy,Ixz,Iyz,Jm,G,B,Tc+,Tc-]
 DYN_base=dynamic parameter of the base robot: [m,sx,sy,sz,Ixx,Iyy,Izz,Ixy,Ixz,Iyz
 Q_Qd_Qdd=[Q QD QDD] (manipulator state)
 X = base spatial position
 V = Base spatial velocity.
 dV = Base spatial acceleration.
 Fbase = Spatial base force due to Base acceleration

 Returns the joint torque required to achieve the specified joint position,
 velocity and acceleration state.
 Gravity is assumed to be acting in the -Z direction with acceleration of
 9.81m/s/s, but this may be overriden by providing a gravity acceleration
 vector grav=[gx gy gz].

 An external force/moment acting on the end of the manipulator may also be
 specified by a 6-element vector fext=[Fx Fy Fz Mx My Mz].

 where	

 Q, QD and QDD are row vectors of the manipulator state; pos, vel, and accel.
 The torque computed also contains a contribution due to armature inertia.

 Copyright, Robotics and Automation Group, 2006. 
 *				Juli� David Colorado, jdcolorado@puj.edu.co
 *		 		Andr� Jaramillo Botero, ajaramil@puj.edu.co
 *  	        Antonio Alejandro Matta G�ez, amatta@puj.edu.co
 *              Juan Camilo Acosta, jcacosta@puj.edu.co 
 ******************************************************************************/
MAT *ff_neHB(int robot_type,MAT *DH,MAT *DYN,MAT *DYN_base,MAT *Q_Qd_Qdd,VEC *GRAV,VEC *fext,MAT *Torques,MAT *X, MAT *Vb, MAT *dVb,MAT *Fbase, int son, MAT *intcon)
{
	MAT *DHtmp=NULL, *DYNtmp=NULL, *Q_Qd_Qddtmp=NULL;
	int lrama, j;
	
//Searching branch limits
	for(i=son;i<df;i++)
	{	k=0;
		for(j=0;j<df;j++)
			if(intcon->me[j][i]==1)
				k++;
		if(k>=2)
		{	
			lrama=i;
			i=df;
		}
	}

//Compute for each point of trajectory
	for(k=0;k<np;k++)
	{
//**************************************************************************
//Forward recurrence propagating kinematics parameters
//**************************************************************************
		m_zero(R);	//Reset values
		m_zero(Rot);
		for(i=son;i<lrama;i++)
		{
			if(DH->me[i][4]==0)
			{	
				H->ve[0]=0;
				H->ve[1]=0;
				H->ve[2]=1;               //projection Vector of rotational motion
				H->ve[3]=0;
				H->ve[4]=0;
				H->ve[5]=0;
			}
			else
			{
				H->ve[0]=0;
				H->ve[1]=0;
				H->ve[2]=0;
				H->ve[3]=0;
				H->ve[4]=0;
				H->ve[5]=1;               //projection Vector of Prismatic motion
			}

			if(i==son)
			{
			//For a Fixed Base  
				if (robot_type!=2 )
				{
					m_ident(r);
					v_zero(p);
				}
			//For a Flying/floating base
				else
				{
					//Obtain 6-Dimensional Position from X
					get_col(X,k,columna);
					euler(columna,r,p);	//Obtain Rotation r and Position p
				}
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
				m_move(Rot,0,0,3,3,r,0,0);	//obtain rotation matrix
				p->ve[0]=-a;			//obtain the inverse of position p
				p->ve[1]=-d*sin(alfa);
				p->ve[2]=-d*cos(alfa);
			}

			//Obtain the 6-dimensional R operator based on r
			m_move(r,0,0,3,3,R,0,0);
			m_move(r,0,0,3,3,R,3,3);

//*********************************************************************************
// Computing Spatial Velocities
//*********************************************************************************
			
			//Obtain the 6-dimensional P operator based on p
			skew_symetric(SkewP,p);
			m_move(identidad,0,0,3,3,P,0,0);
			m_move(identidad,0,0,3,3,P,3,3);
			
			sm_mlt(-1,SkewP,SkewP);
			m_move(SkewP,0,0,3,3,P,0,3);
			
			//Obtain 6-dimensional velocity: Vel
			sv_mlt(dQ->me[k][i],H,dQH);
			m_transp(P,P_trans);
			m_transp(R,R_trans);
			m_mlt(P_trans,R_trans,PtRt);

			//spatial velocities for each joint
			if ((robot_type==2 && i==son && Vb!=NULL) )//|| (robot_type==4 && i==0 ))
			{
				if(son==0)
				{
					for(cont=0;cont<=5;cont++)
						V->me[cont+6*k][i]=dQH->ve[cont]
							+R_trans->me[cont][0]*Vb->me[2][k]
							+R_trans->me[cont][1]*Vb->me[1][k]
							+R_trans->me[cont][2]*Vb->me[0][k]
							+R_trans->me[cont][3]*Vb->me[3][k]
							+R_trans->me[cont][4]*Vb->me[4][k]
							+R_trans->me[cont][5]*Vb->me[5][k];
				}
				else
				{
					for(cont=0;cont<=5;cont++)
						V->me[cont+6*k][i]=dQH->ve[cont]
								+PtRt->me[cont][0]*Vb->me[2][k]
								+PtRt->me[cont][1]*Vb->me[1][k]
								+PtRt->me[cont][2]*Vb->me[0][k]
								+PtRt->me[cont][3]*Vb->me[3][k]
								+PtRt->me[cont][4]*Vb->me[4][k]
								+PtRt->me[cont][5]*Vb->me[5][k];
				}	
			}
			else
			{
				if(i==son)
					for(cont=0;cont<=5;cont++)
						V->me[cont+6*k][i]=dQH->ve[cont];
						
				else
					for(cont=0;cont<=5;cont++)
						V->me[cont+6*k][i]=dQH->ve[cont]
								+PtRt->me[cont][0]*V->me[0+6*k][i-1]
								+PtRt->me[cont][1]*V->me[1+6*k][i-1]
								+PtRt->me[cont][2]*V->me[2+6*k][i-1]
								+PtRt->me[cont][3]*V->me[3+6*k][i-1]
								+PtRt->me[cont][4]*V->me[4+6*k][i-1]
								+PtRt->me[cont][5]*V->me[5+6*k][i-1];
			}
//*********************************************************************************
// Computing Spatial Accelerations
//*********************************************************************************
			
			// Obtain the last and actual angular velocity from Vi-1 and Vi
			for(cont=0;cont<=2;cont++)
			{
				v->ve[cont]=V->me[cont+6*k][i];
				ws->ve[cont]=V->me[cont+6*k][i-1];
			}

			skew_symetric(sk,v);	//obtain Skew symetric wi

			//Obtain 6-dimensional Sk
			m_move(sk,0,0,3,3,Skw,0,0);
			m_move(sk,0,0,3,3,Skw,3,3);
			
			m_transp(r,r_trans);
			mv_mlt(r_trans,ws,rws);
			skew_symetric(sk2,rws); 	//obtain Skew symetric wi-1

			// Obtain the last and actual velocity from Vi-1 and Vi
			for(cont=0;cont<=2;cont++)
			{
				v->ve[cont]= V->me[cont+3+6*k][i];
				ws->ve[cont]= V->me[cont+3+6*k][i-1];
			}
			mv_mlt(r_trans,ws,rv1);
			v_sub(v,rv1,v);
			mv_mlt(sk2,v,vtempo);

			for(cont=0;cont<=2;cont++)
				term->ve[cont+3]=vtempo->ve[cont];

			//Obtain 6-dimensional acceleration
			sv_mlt(d2Q->me[k][i],H,qddH);
			sv_mlt(dQ->me[k][i],H,qdH);
			
			if (robot_type==2  && i==0 )	//For Floating/Flying Base || robot_type==4)
			{	
				if(dVb!=NULL)
				{
					if(son==0)
					{
						for(cont=0;cont<=5;cont++)
							dV->me[cont+6*k][i]=term->ve[cont]+qddH->ve[cont]
								+R_trans->me[cont][0]*(dVb->me[2][k]+GRAV->ve[0])
								+R_trans->me[cont][1]*(dVb->me[1][k]+GRAV->ve[1])
								+R_trans->me[cont][2]*(dVb->me[0][k]+GRAV->ve[2])
								+R_trans->me[cont][3]*(dVb->me[3][k]+GRAV->ve[3])
								+R_trans->me[cont][4]*(dVb->me[4][k]+GRAV->ve[4])
								+R_trans->me[cont][5]*(dVb->me[5][k]+GRAV->ve[5])
								+Skw->me[cont][0]*qdH->ve[0]
								+Skw->me[cont][1]*qdH->ve[1]
								+Skw->me[cont][2]*qdH->ve[2]
								+Skw->me[cont][3]*qdH->ve[3]
								+Skw->me[cont][4]*qdH->ve[4]
								+Skw->me[cont][5]*qdH->ve[5];
					}
					else
					{
						for(cont=0;cont<=5;cont++)
							dV->me[cont+6*k][i]=term->ve[cont]+qddH->ve[cont]
									+PtRt->me[cont][0]*(dVb->me[2][k])
									+PtRt->me[cont][1]*(dVb->me[1][k])
									+PtRt->me[cont][2]*(dVb->me[0][k])
									+PtRt->me[cont][3]*(dVb->me[3][k])
									+PtRt->me[cont][4]*(dVb->me[4][k])
									+PtRt->me[cont][5]*(dVb->me[5][k])
									+Skw->me[cont][0]*qdH->ve[0]
									+Skw->me[cont][1]*qdH->ve[1]
									+Skw->me[cont][2]*qdH->ve[2]
									+Skw->me[cont][3]*qdH->ve[3]
									+Skw->me[cont][4]*qdH->ve[4]
									+Skw->me[cont][5]*qdH->ve[5];	
					}
				}
				else
					for(cont=0;cont<=5;cont++)
						dV->me[cont+6*k][i]=term->ve[cont]+qddH->ve[cont]
								+R_trans->me[cont][0]*(GRAV->ve[0])
								+R_trans->me[cont][1]*(GRAV->ve[1])
								+R_trans->me[cont][2]*(GRAV->ve[2])
								+R_trans->me[cont][3]*(GRAV->ve[3])
								+R_trans->me[cont][4]*(GRAV->ve[4])
								+R_trans->me[cont][5]*(GRAV->ve[5])
								+Skw->me[cont][0]*qdH->ve[0]
								+Skw->me[cont][1]*qdH->ve[1]
								+Skw->me[cont][2]*qdH->ve[2]
								+Skw->me[cont][3]*qdH->ve[3]
								+Skw->me[cont][4]*qdH->ve[4]
								+Skw->me[cont][5]*qdH->ve[5];
			}
			else
			{
				if(i==son)
					for(cont=0;cont<=5;cont++)
						dV->me[cont+6*k][i]=term->ve[cont]+qddH->ve[cont]
								+PtRt->me[cont][0]*(GRAV->ve[0])
								+PtRt->me[cont][1]*(GRAV->ve[1])
								+PtRt->me[cont][2]*(GRAV->ve[2])
								+PtRt->me[cont][3]*(GRAV->ve[3])
								+PtRt->me[cont][4]*(GRAV->ve[4])
								+PtRt->me[cont][5]*(GRAV->ve[5])
								+Skw->me[cont][0]*qdH->ve[0]
								+Skw->me[cont][1]*qdH->ve[1]
								+Skw->me[cont][2]*qdH->ve[2]
								+Skw->me[cont][3]*qdH->ve[3]
								+Skw->me[cont][4]*qdH->ve[4]
								+Skw->me[cont][5]*qdH->ve[5];
				else
					for(cont=0;cont<=5;cont++)
						dV->me[cont+6*k][i]=term->ve[cont]+qddH->ve[cont]
								+PtRt->me[cont][0]*dV->me[0+6*k][i-1]
								+PtRt->me[cont][1]*dV->me[1+6*k][i-1]
								+PtRt->me[cont][2]*dV->me[2+6*k][i-1]
								+PtRt->me[cont][3]*dV->me[3+6*k][i-1]
								+PtRt->me[cont][4]*dV->me[4+6*k][i-1]
								+PtRt->me[cont][5]*dV->me[5+6*k][i-1]
								+Skw->me[cont][0]*qdH->ve[0]
								+Skw->me[cont][1]*qdH->ve[1]
								+Skw->me[cont][2]*qdH->ve[2]
								+Skw->me[cont][3]*qdH->ve[3]
								+Skw->me[cont][4]*qdH->ve[4]
								+Skw->me[cont][5]*qdH->ve[5];
			}
		}
		//end of kinematics propagation
	}	
		"puco, cuadrar fext y cuadrar lo Xb.vb, dvb"
		
	for(i==lrama;i<df;i++)
	{
		if(intcon->me[i][lrama]==1)
		{
			Fbase=ff_neHB(2,DH,DYN,DYN_base,Q_Qd_Qdd,GRAV,fext,MAT *Torques,Xb,Vb,dVb,Fbase,i,intcon);
			m_add(Fbase,Fext,Fext);
		}
	}
//****************************************************************************************
//Backward recurrence: Dynamics Propagation
//****************************************************************************************

//v_copy(fext,Ftant);    //initial  value for external force get from fext
//****************************************************************************************
// Computing Spatial Forces
//****************************************************************************************
	//Compute for each point of trajectory
	for(k=0;k<np;k++)
	{
		m_zero(P);
		m_zero(R);
		
		for(i=lrama;i>=son;i--)
		{
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
			
			m_move(Rot,0,0,3,3,r,0,0);	//Obtain basic matrix Rotation
			
			p->ve[0]=-a;
			p->ve[1]=-d*sin(alfa);
			p->ve[2]=-d*cos(alfa);

			//Obtain the 6-dimensional R operator based on r
			m_move(r,0,0,3,3,R,0,0);
			m_move(r,0,0,3,3,R,3,3);

			//Obtain the 6-dimensional P operator based on p
			skew_symetric(SkewP,p);
			
			m_move(identidad,0,0,3,3,P,0,0);
			m_move(identidad,0,0,3,3,P,3,3);
			
			sm_mlt(-1,SkewP,SkewP);
			
			//form_matrix(0,3,SkewP,P);
			m_move(SkewP,0,0,3,3,P,0,3);

			//Obtain vector s from frame i+1 to center of mass of body i
			s->ve[0]=DYN->me[i][1];
			s->ve[1]=DYN->me[i][2];
			s->ve[2]=DYN->me[i][3];
			
			v_sub(s,p,s);
			mv_mlt(r,s,ss);
			skew_symetric(sks,ss);    //Obtain Skew symetric of ss

			//obtain 6-dimensional operator S
			m_move(identidad,0,0,3,3,S,0,0);
			m_move(identidad,0,0,3,3,S,3,3);
			m_move(sks,0,0,3,3,S,0,3);
			

			//Obtain angular velocity from Vi
			for(cont=0;cont<=2;cont++)
				v->ve[cont]=V->me[cont+6*k][i];

			skew_symetric(sk,v);   //Obtain Skew symetric of v

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
			m_move(Ji,0,0,3,3,Icm,0,0);
			sm_mlt(DYN->me[i][0],identidad,massU);
			m_move(massU,0,0,3,3,Icm,3,3);
			
			//applying parallel axes theorem to obtain the inertial operator at frame i
			m_transp(S,S_trans);
			m_mlt(S,Icm,SIcm);
			m_mlt(SIcm,S_trans,I);
			
//*******************************************************************************
// Computing Force equation:
//*******************************************************************************  
			m_mlt(R,P,RP);
			
			if(i==n-1)
				mv_mlt(RP,fext,RPFext);
			else
				mv_mlt(RP,F,RPFext);
			
			//get_col(dV,i,columna);     //obtain column from spatial acceleration
			m_move(I,0,0,3,3,Itemp,0,0);

			for(cont=0;cont<=2;cont++)
				v->ve[cont]=V->me[cont+6*k][i];

			mv_mlt(Itemp,v,Iv);
			mv_mlt(sk,Iv,up);
			m_mlt(sk,sk,sksk);
			mv_mlt(sksk,ss,sotro);
			sv_mlt(DYN->me[i][0],sotro,down);

			for(cont=0;cont<=2;cont++)
			{
				termino->ve[cont]=up->ve[cont];
				termino->ve[cont+3]=down->ve[cont];
			}

			//Obtain 6-dimensional force
			for(cont=0;cont<=5;cont++)
			{	
				IdV->ve[cont]=I->me[cont][0]*dV->me[0+6*k][i]
						+I->me[cont][1]*dV->me[1+6*k][i]
						+I->me[cont][2]*dV->me[2+6*k][i]
						+I->me[cont][3]*dV->me[3+6*k][i]
						+I->me[cont][4]*dV->me[4+6*k][i]
						+I->me[cont][5]*dV->me[5+6*k][i];
			}
			v_add(RPFext,IdV,RPFext);
			v_add(RPFext,termino,F);
			if(i==son)	
				_set_col(Fall,k,F,0);
//***************************************************************************			
//Add 6-dimensional mass operator of all system referred at base frame
			/*if(robot_type==2)
			{
				if (i==(n-1))
				{
					m_copy(I,Tinertia);
				}
				else
				{
					m_transp(R,R_trans);
					m_transp(P,P_trans);
					m_mlt(RP,Tinertia,RPIant);
					m_mlt(R_trans,P_trans,PtRt);
					m_mlt(RPIant,PtRt,Iup);
					m_add(I,Iup,Tinertia);
				
					if (i==0 && DYN_base!=MNULL)
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
						
						m_move(Ticm,0,0,3,3,Icm,0,0);
						sm_mlt(DYN_base->me[0][0],identidad,massU);
						m_move(massU,0,0,3,3,Icm,3,3);
					
						//Obtain vector s 
						s->ve[0]=DYN_base->me[0][1];
						s->ve[1]=DYN_base->me[0][2];
						s->ve[2]=DYN_base->me[0][3];
						skew_symetric(sks,s);    //Obtain Skew symetric of s
					
					
						//obtain 6-dimensional operator S
						m_move(identidad,0,0,3,3,Sb,0,0);
						m_move(identidad,0,0,3,3,Sb,3,3);
						m_move(sks,0,0,3,3,Sb,0,3);
						
						m_transp(Sb,S_transb);
						m_mlt(Sb,Icm,SIcm);
						m_mlt(SIcm,S_transb,I);
						//add inertia of base
						m_add(Tinertia,I,Tinertia);
					}
				}
			}*/
//*************************************************************************************
//Obtain Spatial Base Force for free Flying Robots
//*************************************************************************************
/*			if(robot_type==2 && dVb!=NULL && Fbase!=NULL)
			{
			//Obtain spatial force of the base: F=ma	
				if (i==0)
				{
					//mv_mlt(Tinertia,Ace,resul);
					for(cont=0;cont<=5;cont++)
						//Fbase->me[cont][k]=resul->ve[cont];
						Fbase->me[cont][k]=Tinertia->me[cont][2]*dVb->me[1][k]
								+Tinertia->me[cont][1]*dVb->me[1][k]
								+Tinertia->me[cont][2]*dVb->me[0][k]
								+Tinertia->me[cont][3]*dVb->me[3][k]
								+Tinertia->me[cont][4]*dVb->me[4][k]
								+Tinertia->me[cont][5]*dVb->me[5][k];
				}
			}*/
	
			//Projecting force in the axe of motion
			if(DH->me[i][4]==0)                 //Rotational
				Torques->me[k][i]=F->ve[2];
			else                               //Prismatic
				Torques->me[k][i]=F->ve[5];

		}
		//end dynamics propagation
	}
		
	"cuadrar F"
	return F(:,son);

}