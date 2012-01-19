#include "LLRL.h"
/*****************************************************************************
	Compute inverse dynamics via recursive Newton-Euler formulation
 
 	TAU = ne(DH,DYN,Q_Qd_Qdd,grav,fext)

	DH=Denavit Hartenberg parameter matrix for manipulator
	DYN=dynamic parameter matrix in the form [m,sx,sy,sz,Ixx,Iyy,Izz,Ixy,Ixz,Iyz]
	Q_Qd_Qdd=[Q QD QDD] (manipulator state)

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
 

	Copyright, Andrés Jaramillo Botero, 1999. 
******************************************************************************/
MAT *ne(MAT *DH,MAT *DYN,MAT *Q_Qd_Qdd,VEC *gravity,VEC *fext,MAT *Fout) 
{
	static	VEC	*z0,*v1,*v2,*v3,*v4,*vhat,*m,*grav;    
	static	VEC *q,*qd,*qdd,*F,*N,*f,*nn,*pstar,*w,*wd,*v,*vd;
	static	MAT *Q,*Qd,*Qdd,*Fm,*Nm,*Jm,*J,*Rm,*R,*pstarm,*r,*mt;    
	int		n,np,j,p;
	double	alpha,A,theta,D;
	double	sa,ca,st,ct;
	//FILE *fp;
    
	/****************** INVERSE DYNAMICS *************************/
	
	n=(int) DH->m;				/* number of joints  */
	np=(int) Q_Qd_Qdd->m;		/* number of trajectory points */
	Fout=m_resize(Fout,np,n);
	v_resize_vars(3,&v1,&z0,NULL);
	mem_stat_reg_vars(0,TYPE_VEC,&v1,&z0,NULL);
	z0->ve[2]=1;
	
	/* validate inputs (sizes and initializations) */
	if (DH == (MAT *)MNULL || DYN == (MAT *)MNULL) 
		error(E_NULL,"ne (DH and dynamic inputs)");
	if (DH->m != DYN->m) 
		error(E_SIZES,"ne (DH and dynamic inputs)\n");
	if (Q_Qd_Qdd == (MAT *)MNULL) 
		error(E_NULL,"Null input (state) in newton-euler");
	if (gravity == (VEC *)VNULL) 
	{
		/* Define the gravity vector */
		v1=sv_mlt(9.81, z0, v1);
//		grav=v_resize(grav,3);
		grav=v_copy(v1,grav);
	}
	else grav=v_copy(gravity,grav);
	MEM_STAT_REG(grav,TYPE_VEC);
	if (fext == (VEC *)VNULL) fext=v_resize(fext,6);
	
	/* allocate temporals */
	v_resize_vars(3,&v2,&v3,&v4,&vhat,&F,&N,&f,&nn,NULL);
	mem_stat_reg_vars(0,TYPE_VEC,&v2,&v3,&v4,&vhat,&F,&N,&f,&nn,NULL);
	
	/* validate size of Q_Qd_Qddd: position, veloc, accel */
	if ((int) Q_Qd_Qdd->n == 3*n && (int) Q_Qd_Qdd->m == np)
	{
		m_resize_vars(np,n,&Q,&Qd,&Qdd,NULL);
		mem_stat_reg_vars(0,TYPE_MAT,&Q,&Qd,&Qdd,NULL);
		m_move(Q_Qd_Qdd,0,0,np,n,Q,0,0);
		m_move(Q_Qd_Qdd,0,n,np,n,Qd,0,0);
		m_move(Q_Qd_Qdd,0,2*n,np,n,Qdd,0,0);
	}
	else
		error(E_SIZES,"ne (incorrect manipulator state matrix)");
	
		/* Example of input format for Q,Qd,Qdd joint
		Qd=[Qd11	Qd21	Qd31	.	Qdn1	timestep
		Qd12	Qd22	Qd32	.	Qdn2
		.		.		.		.	.
		Qd1np	Qd2np	Qd3np	.	Qdnnp
		]
	*/
	
	/* link masses */
	m=v_resize(m,n);		
	MEM_STAT_REG(m,TYPE_VEC);
	get_col(DYN,0,m);
	
	/* Center of Mass */
	r=m_resize(r,n,3);		
	m_move(DYN,0,1,n,3,r,0,0);
	
	/* Initialize temporal mt and inertia tensor J, for first joint */
	m_resize_vars(3,3,&mt,&Jm,&J,NULL);
	mem_stat_reg_vars(0,TYPE_MAT,&r,&mt,&Jm,&J,NULL);
	
	linkiner(DYN,0,Jm);
	/* calculate inertia tensor for each link (2..n) and put in one 3 X n*3 matrix */
	for (j=1;j<n;j++)
	{
		linkiner(DYN,j,mt);
		Jm=m_resize(Jm,3,3*(j+1));
		m_move(mt,0,0,3,3,Jm,0,3*j);
	}
	v_resize_vars(n,&q,&qd,&qdd,NULL);
	mem_stat_reg_vars(0,TYPE_VEC,&q,&qd,&qdd,NULL);
	m_resize_vars(3,3,&Rm,&R,NULL);
	m_resize_vars(3,n,&pstarm,&Fm,&Nm,NULL);
	mem_stat_reg_vars(0,TYPE_MAT,&Rm,&R,&pstarm,&Fm,&Nm,NULL);
	
	v_resize_vars(3,&w,&wd,&v,&vd,&pstar,NULL);
	mem_stat_reg_vars(0,TYPE_VEC,&w,&wd,&v,&vd,&pstar,NULL);
	
	Rm=m_resize(Rm,3,3*n);
	for (p=0;p<np;p++)
	{
		get_row(Q,p,q);
		get_row(Qd,p,qd);
		get_row(Qdd,p,qdd);
		
		/* initial acceleration = grav */
		vd=v_copy(grav,vd);		
		v_zero(w);
		v_zero(wd);
		v_zero(v);
		
		/* Define variable joint and compute rotation matrices 
		and translation vectors */
		for (j=0;j<n;j++)
		{
			alpha = DH->me[j][0];
			A = DH->me[j][1];
			if (DH->me[j][4] == 0)
			{
				theta = q->ve[j];
				D = DH->me[j][3];
				
      } else 
			{
				theta = DH->me[j][3];  // No debería ser theta = DH->me[j]["2"] ?
				D = q->ve[j];
      }
			
			sa = sin(alpha); 
			ca = cos(alpha); 
			st = sin(theta); 
			ct = cos(theta);
			
			/*	R = [ct		-st*ca		  st*sa
							 st    ct*ca		 -ct*sa
								0       sa			  ca ]; */
			
			R->me[0][0]=ct;
			R->me[0][1]=-st*ca;
			R->me[0][2]=st*sa;
			R->me[1][0]=st;
			R->me[1][1]=ct*ca;
			R->me[1][2]=-ct*sa;
			R->me[2][0]=0;
			R->me[2][1]=sa;
			R->me[2][2]=ca;
			
			m_move(R,0,0,3,3,Rm,0,3*j); 
			pstar->ve[0]=A; 
			pstar->ve[1]=D*sa; 
			pstar->ve[2]=D*ca;
			pstarm=m_resize(pstarm,3,j+1);
			set_col(pstarm,j,pstar);
        }
		
		/* the forward recursion */
		
		for (j=0;j<n;j++) 
		{
			m_move(Rm,0,3*j,3,3,mt,0,0);
			m_transp(mt,R);
			get_col(pstarm,j,pstar);
			
			/* statement order is important here */
			if (DH->me[j][4] == 0)
			{
				/* revolute axis */
				sv_mlt(qd->ve[j],z0,v1);
				v_cross(w,v1,v2);
				sv_mlt(qdd->ve[j],z0,v1);
				v_add(v1,v2,v3);
				v_add(wd,v3,v1);
				mv_mlt(R,v1,wd);
				
				/*
				wd=mv_mlt(R,v_add(wd,
				v_add(v1,v_cross(w,sv_mlt(qd->ve[j],z0,v1),v2),v2),
				sv_mlt(qdd->ve[j],z0,v1)),wd);
				*/				
				mv_mlt(R,v_add(w,sv_mlt(qd->ve[j],z0,v1),v2),w);
				
				v_cross(wd,pstar,v1);
				v_cross(w,pstar,v2);
				v_cross(w,v2,v3);
				v_add(v1,v3,v2);
				mv_mlt(R,vd,v1);
				v_add(v1,v2,vd);
				
				/*
				wd = R*(wd + z0*qdd(j) + cross(w,z0*qd(j)));
				w = R*(w + z0*qd(j));
				vd = cross(wd,pstar) + cross(w, cross(w,pstar)) +R*vd;
				*/
			}
			else
			{
				/* prismatic axis */
				mv_mlt(R,w,v1);
				v_copy(v1,w);			// w
				mv_mlt(R,wd,v2);		
				v_copy(v2,wd);			//wd
				
				mv_mlt(R,v_add(sv_mlt(qdd->ve[j],z0,v1),vd,v2),v3);
				v_add(v3,v_cross(wd,pstar,v1),v2);		// vd temp 1 in v2
				sv_mlt(2,v_cross(w,mv_mlt(R,sv_mlt(qd->ve[j],z0,v1),v3),v1),v1); // vd temp 2 in v1 
				v_add(v1,v2,v3);	// vd temp 1a in v3
				v_add(v3,mv_mlt(R,vd,v1),v2);	//vd temp 2a in v2
				v_add(v2,v_cross(w,v_cross(w,pstar,v1),v3),vd);
				/*			
				w = R*w;
				wd = R*wd;
				v = R*(z0*qd(j) + v) + cross(w,pstar);
				vd = R*(z0*qdd(j)+vd) + cross(wd,pstar) + 2*cross(w,R*z0*qd(j)) ...
				+ cross(w, cross(w,pstar)) +R*vd;
				*/
			}		
			/*
			J = Jm(:,3*j-2:3*j);
			vhat = cross(wd,r(j,:)') + cross(w,cross(w,r(j,:)')) + vd;
			F = m(j)*vhat;
			N = J*wd + cross(w,J*w);
			Fm = [Fm F];
			Nm = [Nm N];
			*/
			linkiner(DYN,j,J);
			v_cross(wd,get_row(r,j,v1),v2);
			v_cross(w,get_row(r,j,v1),v3);
			v_cross(w,v3,v1);
			v_add(v2,v1,v3);
			v_add(v3,vd,vhat);
			sv_mlt(m->ve[j],vhat,F);
			v_add(mv_mlt(J,wd,v1),v_cross(w,mv_mlt(J,w,v2),v3),N);
			set_col(Fm,j,F);
			set_col(Nm,j,N);
		}
		
		/* the backward recursion */
		/* force/moments on end of arm */
		f->ve[0] = fext->ve[0];	
		f->ve[1] = fext->ve[1];	
		f->ve[2] = fext->ve[2];	
		nn->ve[0] = fext->ve[3];	
		nn->ve[1] = fext->ve[4];	
		nn->ve[2] = fext->ve[5];	
		
		for (j=n-1;j>-1;j--)
		{
			get_col(pstarm,j,pstar);
			/* order of these statements is important, 
			since both nn and f are functions of previous f */
			if (j == n-1)
				m_ident(R);
			else
				m_move(Rm,0,3*j+3,3,3,R,0,0);
			
			m_transp(R,mt);
			mv_mlt(mt,pstar,v1);
			v_cross(v1,f,v2);
			v_add(v2,nn,v1);
			mv_mlt(R,v1,v3);		//R*(nn + cross(R'*pstar,f))
			get_row(r,j,v1);
			v_add(pstar,v1,v2);
			get_col(Fm,j,v4);
			v_cross(v2,v4,nn);
			v_add(nn,v3,v2);
			get_col(Nm,j,v1);
			
			v_add(v2,v1,nn);	/* cross(pstar+r(j,:)',Fm(:,j)) + Nm(:,j) */
			v_add(mv_mlt(R,f,v1),get_col(Fm,j,v2),f);
			m_move(Rm,0,3*j,3,3,R,0,0);
			
			/*
			nn = R*(nn + cross(R'*pstar,f)) + cross(pstar+r(j,:)',Fm(:,j)) + Nm(:,j);
			f = R*f + Fm(:,j);
			R = Rm(:,3*j-2:3*j);
			*/
			
			/* revolute */
			if (DH->me[j][4] == 0) 
			{
				mv_mlt(m_transp(R,mt),z0,v2);
				Fout->me[p][j]=in_prod(nn,v2)+
					DYN->me[j][11]*DYN->me[j][11]*(
					DYN->me[j][10]*qdd->ve[j]+
					DYN->me[j][12]*qd->ve[j]/*+
											(qd->ve[j]>0 ? DYN->me[j][13] : 0.0) +
											(qd->ve[j]<0 ? DYN->me[j][14] : 0.0) 
											*/					);
											/* los terminos de friccion no permiten la convergencia de ode45 */
			}
			else					
				/* prismatic */
			{
				mv_mlt(m_transp(R,mt),z0,v2);
				Fout->me[p][j]=in_prod(f,v2)+
					DYN->me[j][11]*DYN->me[j][11]*(
					DYN->me[j][10]*qdd->ve[j]+
					DYN->me[j][12]*qd->ve[j]+
					(qd->ve[j]>0 ? DYN->me[j][13] : 0.0) +
					(qd->ve[j]<0 ? DYN->me[j][14] : 0.0) 
					);
				//DYN->me[j][11]*DYN->me[j][11]*(DYN->me[j][10]*qdd->ve[j]+DYN->me[j][12]*qd->ve[j]+coulomb);
			}
		}
	}
	
    /*if((fp=fopen("forces.m","a+"))==NULL)
		printf("\n can not open file");
	else
	{
		for(i=0;i<Fout->n;i++)

		m_fwrite(fp,Fout,"fuerza");
		fclose(fp);
	}
    */return Fout;
}
