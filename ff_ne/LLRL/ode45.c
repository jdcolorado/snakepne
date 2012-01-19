#include "LLRL.h"
/*************************************************************************
Workspace used: 11,14,3, 5, 6, 8, 9
**************************************************************************		
	ODE45  Solve non-stiff differential equations, medium order method.
    
	[T,Y] = ODE45(ODEFUN,TSPAN,Y0) 
	
	with TSPAN = [T0 TFINAL] integrates the system of differential equations 
	y' = f(t,y) from time T0 to TFINAL with initial conditions Y0. 
	Function ODEFUN(T,Y) must return a column vector
    corresponding to f(t,y). Each row in the solution array Y corresponds to
    a time returned in the column vector T. To obtain solutions at specific
    times T0,T1,...,TFINAL (all increasing or all decreasing), use 
    TSPAN = [T0 T1 ... TFINAL].     
		  
    ODE45 can solve problems M(t,y)*y' = f(t,y) with mass matrix M that is
    nonsingular. 

	Derived from Matlab's ode45 function
	Copyright, Marcial Qui�nez y Andr� Jaramillo Botero, 1999. 
***************************************************************************/

/* 07-5-2002, AJB: Eliminated function number -n_funcion- argument in ODE45 */
/* 07-09-2002, AJB: Added static declarations and memory management functions for feval calls */
/* 07-11-2003, AJB: Need to check fevals arguments ... int */
/* 09-02-2006, JCA: Added support for n df but the problem was with feval, not with ode45.*/
/*22-02-2006, JCA, JDC:Add spatial base position: X=[roll pitch yaw x y z]'*/
/*23-02-2006, JCA, mem stat reg and meme stat free still not working...*/
/*30/03/2006, JCA, added support for passing fext and grav and;
		added support for torques!=NULL  and flying base!!
		added support for mem_stat_mark and mem_stat_free in *fun calls.*/
/*16/04/2006, JCA added support for base force*/
/*16/04/2006 JCA Solved memory problems with maxm and maxm2*/
/*26/04/2006, JCA Added sopport for robot type, and flying base as it should be*/
/*07/05/2006 JDC, Changed interface, added DYN_base*/
/*07/05/2006 JDC, Changed prototype call for feval, added DYN_base*/

//funcion ode45 de matlab,  recibe un vector con tiempos,  uno con puntos de evaluacion, una tolerancia,  otra tolerancia
// un indice de refinamiento.

MAT	*ode45 (MAT *(*fun)(),int Robot_type,MAT *Torques,VEC *tspan,VEC *y0, double  rtol, double atol,int refine,MAT* input_dh, MAT *input_dyn, MAT *rK, VEC *grav, VEC *fext,MAT *DYN_base)
{
	double	t0,tfinal,tdir,t,t_send,hmax,pow2,hmin,absh,threshold,h,rh,tnew=0,err=0,temp;
	int	ntspan,neq,outflag,i,chunk=1,nout,next,done,nofailed,oldnout,m,n,j,tmp59, siempre=1,points,df;	
	static	MAT *y=MNULL,*S=MNULL,*tout=MNULL,*yout=MNULL,*A=MNULL,*B=MNULL;
	static	MAT *tmp=MNULL,*sol=MNULL,*E=MNULL,*f=MNULL,*f0=MNULL,*hA=MNULL;
	static	MAT *hB=MNULL,*ynew=MNULL,*iii=MNULL,*mtmp=MNULL,*mtmp2=MNULL, *tempm=MNULL;
	static	VEC *vtmp=VNULL,*vtmp2=VNULL,*time=VNULL;
	static	IVEC *ii=IVNULL;

	/* NEED TO CHECK ERRORS IN INPUT/OUTPUT ARGS */
	if(Torques!=NULL)
	{
		df=(Torques->n-1);
		points=Torques->m; //Trayectory points
		//gets time from Torque matrix
		time=v_resize(time,points); 
		MEM_STAT_REG(time,TYPE_VEC);
		time=get_col(Torques,df,time);
		//resize Torque to delete time column
		m_resize(Torques,points,df);
	}
	tempm = m_resize(tempm,1,1);
	mtmp = m_resize(mtmp,1,1);
	mtmp2 = m_resize(mtmp2,1,1);
	//inicializacion de las matrices  a evaluar
	A = m_resize(A,6,1);
	B = m_resize(B,7,6);
	E = m_resize(E,7,1);
	f = m_resize(f,y0->dim,7);
	f0= m_resize(f0,y0->dim,1);
	y=m_resize(y,y0->dim,1);
	ynew=m_resize(ynew,y0->dim,1);
	S=m_resize(S,refine-1,1);
	hA = m_resize(hA,6,1);
	hB = m_resize(hB,7,6);
	ntspan=tspan->dim;
	iii=m_resize(iii,1,1);
	tmp=m_resize(tmp,1,1);
	sol=m_resize(sol,1,1);
		
	mem_stat_reg_vars(0, TYPE_MAT, &y, &S, &A, &B, &tmp, &sol, &E, &f, &f0, &hA, &hB, &ynew, &iii, &mtmp, &mtmp2, &tempm, NULL);//&tor, &tout, &yout,NULL);
	vtmp=v_resize(vtmp,1);
	vtmp2=v_resize(vtmp2,1);
	mem_stat_reg_vars(0, TYPE_VEC, &vtmp, &vtmp2,NULL);
	ii=iv_resize(ii,1);
	MEM_STAT_REG(ii,TYPE_IVEC);
	
	t0=tspan->ve[0]*1.0;
	next=1;
	tfinal=tspan->ve[ntspan-1];
	
	if(tfinal - t0>0)
		tdir=1;
	else
		tdir=-1;
	if(t0<0)
		t=t0*(-1.0);
	else
		t=t0*1.0;
	
	neq=y->m;
	
	for(i=0;i<neq;i++)
		y->me[i][0]=y0->ve[i];

	if (rtol<100*2.2204e-016)
		rtol=100*2.2204e-016;
	
	threshold = (double) atol/rtol;
	hmax=fabs((tfinal-t)*0.1);
	
	//ntspan es el numero de datos en el vector de tiempo
	if (ntspan >2)
		outflag=1;
	else
	{	
		outflag=3;
		for(i=1;i<=refine-1;i++)
			S->me[i-1][0]=i*1.0/refine*1.0;
	}

	//dimensionamiento inicial de las matrices de salidas de los datos
	if (ntspan >2)
	{
		tout=m_resize(tout,ntspan,1);
		yout=m_resize(yout,ntspan,neq);
		
	}
	else
	{
		chunk=(int)max(ceil(128/neq),refine);
		//chunk=min(max(100,50*refine),floor((2^13)/neq));
		tout=m_resize(tout,chunk,1);
		yout=m_resize(yout,chunk,neq);
	}
	mem_stat_reg_vars(0,TYPE_MAT,&yout,&tout,NULL);	
	m_zero(yout);
	m_zero(tout);
	
	//asignacion de valores iniciales a las matrices de salida
	nout=1;
	tout->me[nout-1][0]=t;
	vtmp=get_row(y,0,vtmp);
	set_row(yout,0,vtmp);

	pow2 = 1.0 / 5.0;
	m_zero(A);
	m_zero(B);
	m_zero(E);
	m_zero(f);
	
	//inicializacion de los parametros del metodo
	A ->me[0][0]=1.0/5.0; A ->me[1][0]= 3.0/10.0; A ->me[2][0]=4.0/5.0; A ->me[3][0]= 8.0/9.0;  A ->me[4][0]=1.0;
	A ->me[5][0] =1.0; 
	B->me[0][0]=1.0/5.0;
	B->me[0][1]=3.0/40.0;       B->me[1][1]=9.0/40.0;
	B->me[0][2]=44.0/45.0;      B->me[1][2]=-56.0/15.0;      B->me[2][2]=32.0/9.0;
	B->me[0][3]=19372.0/6561.0; B->me[1][3]=-25360.0/2187.0; B->me[2][3]=64448.0/6561.0; B->me[3][3]=-212.0/729.0;
	B->me[0][4]=9017.0/3168.0;  B->me[1][4]=-355.0/33.0;     B->me[2][4]=46732.0/5247.0; B->me[3][4]=49.0/176.0; B->me[4][4]=-5103.0/18656.0;
	B->me[0][5]=35.0/384.0;     B->me[1][5]=0.0;             B->me[2][5]=500.0/1113.0;   B->me[3][5]=125.0/192.0; B->me[4][5]=-2187.0/6784.0; B->me[5][5]=11.0/84.0;
	
	E ->me[0][0]=71.0/57600.0; E ->me[1][0]= 0.0; E ->me[2][0]=-71.0/16695.0; E ->me[3][0]= 71.0/1920.0;  E ->me[4][0]=-17253.0/339200.0;
	E ->me[5][0] =22.0/525.0; E ->me[6][0]= -1.0/40.0;

	//evaluacion de la funcion con los valores iniciales

mem_stat_mark(3);
	(*fun)(torqfun,Robot_type,Torques,time,t,y,f0,input_dh,input_dyn,grav,fext,DYN_base);
mem_stat_free(3);

	vtmp=v_resize(vtmp,f0->n);
	vtmp=get_col(f0,0,vtmp);
	set_col(f,0,vtmp);

	hmin=16*2.2204e-016*fabs(t);
	absh=min(hmax,fabs((tspan->ve[next])-t));
mem_stat_mark(3);
	tempm=maxm(absm(y),threshold,tempm);
mem_stat_free(3);
	
mem_stat_mark(3);
	mtmp=div_m(f0,tempm,mtmp);
mem_stat_free(3);
	
	vtmp=get_col(mtmp,0,vtmp);
	rh=v_norm_inf(vtmp)/(0.8*pow(rtol,pow2));

	if (absh*rh>1)
		absh=1/rh;
	
	absh=max(absh,hmin);
	
	done=1;
	while (done==1)
	{
		hmin=16*2.2204e-016*fabs(t);

		absh=min(hmax,max(hmin,absh));
		h=tdir*absh;

		if(1.1*absh>=fabs(tfinal-t))
		{
			h=tfinal-t;
			absh=fabs(h);
			done=0;
		}
		nofailed=1;
		
		while (siempre==1)
		{
			m_zero(f0);
			sm_mlt(h,A,hA);
			sm_mlt(h,B,hB);
			vtmp=v_resize(vtmp,hB->m);
			vtmp2=v_resize(vtmp2,f->n);
		
			t_send=t+hA->me[0][0];
mem_stat_mark(3);	
			vtmp=get_col(hB,0,vtmp);
		
			vtmp2=mv_mlt(f,vtmp,vtmp2);
			mtmp=sum_mat_vec(y,vtmp2,mtmp);
			(*fun)(torqfun,Robot_type,Torques,time,t_send,mtmp,f0,input_dh,input_dyn,grav,fext,DYN_base);
mem_stat_free(3);
			vtmp=v_resize(vtmp,f0->m);
			vtmp=get_col(f0,0,vtmp);
			set_col(f,1,vtmp);
			m_zero(f0);
			vtmp=v_resize(vtmp,hB->m);
			t_send=t+hA->me[1][0];
mem_stat_mark(3);
			vtmp=get_col(hB,1,vtmp);
			mtmp=sum_mat_vec(y,mv_mlt(f,vtmp,vtmp2),mtmp);
			(*fun)(torqfun,Robot_type,Torques,time,t_send,mtmp,f0,input_dh,input_dyn,grav,fext,DYN_base);
mem_stat_free(3);
			vtmp=v_resize(vtmp,f0->m);
			vtmp=get_col(f0,0,vtmp);
			set_col(f,2,vtmp);
			m_zero(f0);
			vtmp=v_resize(vtmp,hB->m);
			t_send=t+hA->me[2][0];
mem_stat_mark(3);
			vtmp=get_col(hB,2,vtmp);
			mtmp=sum_mat_vec(y,mv_mlt(f,vtmp,vtmp2),mtmp);
			(*fun)(torqfun,Robot_type,Torques,time,t_send,mtmp,f0,input_dh,input_dyn,grav,fext,DYN_base);
mem_stat_free(3);
			vtmp=v_resize(vtmp,f0->m);
			vtmp=get_col(f0,0,vtmp);
			set_col(f,3,vtmp);
			m_zero(f0);
			vtmp=v_resize(vtmp,hB->m);
			t_send=t+hA->me[3][0];
mem_stat_mark(3);
			vtmp=get_col(hB,3,vtmp);
			mtmp=sum_mat_vec(y,mv_mlt(f,vtmp,vtmp2),mtmp);
			(*fun)(torqfun,Robot_type,Torques,time,t_send,mtmp,f0,input_dh,input_dyn,grav,fext,DYN_base);
mem_stat_free(3);
			vtmp=v_resize(vtmp,f0->m);
			vtmp=get_col(f0,0,vtmp);
			set_col(f,4,vtmp);
			m_zero(f0);
			vtmp=v_resize(vtmp,hB->m);
			t_send=t+hA->me[4][0];
mem_stat_mark(3);
			vtmp=get_col(hB,4,vtmp);
			mtmp=sum_mat_vec(y,mv_mlt(f,vtmp,vtmp2),mtmp);
			(*fun)(torqfun,Robot_type,Torques,time,t_send,mtmp,f0,input_dh,input_dyn,grav,fext,DYN_base);
//CHECK MNULL AS IN ARG and size of mtmp
mem_stat_free(3);
			vtmp=v_resize(vtmp,f0->m);
			vtmp=get_col(f0,0,vtmp);
			set_col(f,5,vtmp);
			m_zero(f0);
			vtmp=v_resize(vtmp,hB->m);
			vtmp=get_col(hB,5,vtmp);
			ynew=sum_mat_vec(y,mv_mlt(f,vtmp,vtmp2),ynew);
			
			m_zero(f0);
			tnew=t+indice(hA,5);
			t_send=t+hA->me[5][0];
			
mem_stat_mark(3);
			(*fun)(torqfun,Robot_type,Torques,time,t_send,ynew,f0,input_dh,input_dyn,grav,fext,DYN_base);
mem_stat_free(3);
			vtmp=v_resize(vtmp,f0->m);
			vtmp=get_col(f0,0,vtmp);
			set_col(f,6,vtmp);
			
			m=f->n;
			n=E->n;
			
			vtmp=v_resize(vtmp,f->m);
			tmp=m_resize(tmp,m,n);
			
mem_stat_mark(3);
			tempm=maxm2(absm(y),absm(ynew),tempm);
mem_stat_free(3);
mem_stat_mark(3);
			mtmp2=maxm(tempm,threshold,mtmp2);
mem_stat_free(3);
			mtmp=div_m(m_mlt(f,E,tmp),mtmp2,mtmp);
			vtmp=get_col(mtmp,0,vtmp);
			err=absh*v_norm_inf(vtmp);
			
			//la solucion solamente funciona si el error es menor que la tolerancia,  si es mayor se recalcula 
			//el h para la siguiente evalacion.
			if (err>rtol)
			{
				if (absh<=hmin)
					error(EF_EXIT,"ode45 (Unable to meet integration tolerances without reducing the step size)\n");
				if (nofailed==1)
				{
					nofailed=0;
					absh=max(hmin,absh*max(0.1,0.8*pow((rtol/err),pow2)));
				}
				else
				    absh=max(hmin,0.5*absh);
				
				h=tdir*absh;
				done=1;
			}
			else
				break;
		}
		//nout es el numero de resultados validos que de las matrices de salida
		oldnout=nout;
		if (outflag==3)//dos valores de tiempo de entrada inicial
		{
			nout=nout+refine;
			m=tout->m;
			if (nout>m)
			{
				//redimensionamiento de la matriz de salida
				m=tout->m;
				tout=m_resize(tout,m+chunk,1);
				
				for (i=m;i<m+chunk;i++) tout->me[i][0]=0.0;
				m=yout->m;
				
				yout=m_resize(yout,yout->m+chunk,neq);
				for(i=m;i<m+chunk;i++)
					for(j=0;j<neq;j++)
						yout->me[i][j]=0.0;
			}
			
			//calculo de los valores de salida
			ii=iv_resize(ii,nout-oldnout-1);
			for(i=oldnout+1;i<=nout-1;i++)
				ii->ive[i-oldnout-1]=i;
	
			m=ii->dim;
			iii=m_resize(iii,m,1);
			
			for(i=0;i<m;i++)
				iii->me[i][0]=((S->me[i][0])*(tnew-t))+t;
			
			for(i=0;i<m;i++)
			{    
				tmp59=ii->ive[i];
				tout->me[tmp59-1][0]=iii->me[i][0];
			}
			
			sol=m_resize(sol,neq,iii->m);
			m=y->m;
			n=iii->m;
			
mem_stat_mark(3);
			ntrp45(iii,t,y,h,f,sol);	//changed tmp for sol in return param
mem_stat_free(3);
			m=sol->n;

			for(i=0;i<m;i++)
				for(j=0;j<neq;j++)
				{
					tmp59=ii->ive[i];
					yout->me[tmp59-1][j]=sol->me[j][i];
				}
			tout->me[nout-1][0]=tnew;
			for(i=0;i<neq;i++)
				yout->me[nout-1][i]=ynew->me[i][0];
		}
		else//mas de dos valores de tiempo en la entrada inicial
		{
			while(next<=ntspan)
			{
				if(tdir*(tnew-tspan->ve[next])<0)
					break;
				else if(tnew==tspan->ve[next])
				{
					nout=nout+1;					
					tout->me[nout-1][0]=tnew;
					for(i=0;i<neq;i++)
						yout->me[nout-1][i]=ynew->me[i][0];
					next=next+1;
					break;
				}
				nout=nout+1;
				
				tout->me[nout-1][0]=tspan->ve[next];
				sol=m_resize(sol,neq,1);
				iii=m_resize(iii,1,1);
				m=y->m;
				n=iii->m;
//				
				iii->me[0][0]=tspan->ve[next];
mem_stat_mark(3);
				ntrp45(iii,t,y,h,f,sol);	//changed tmp for sol in return arg	
mem_stat_free(3);
				n=(sol->n)-1;
				for(j=0;j<neq;j++)
					yout->me[nout-1][j]=sol->me[j][n];
				next=next+1;
			}
		}
		if(nofailed==1)
		{
			temp=1.25*pow((err/rtol),pow2);
			if(temp>0.2)
				absh=absh/temp;
			else
				absh=5.0*absh;
		}

		t=tnew;
		y=m_copy(ynew,y);
		vtmp=v_resize(vtmp,f->m);
		vtmp=get_col(f,6,vtmp);
		set_col(f,0,vtmp);
	}
	//redimensionamiento de la matriz para dejar solamente los valores de salida validos

	tout=m_resize(tout,nout,1);
	yout=m_resize(yout,nout,yout->n);
	rK=m_resize(rK,nout,1+yout->n);	
	m_move(tout,0,0,nout,1,rK,0,0);
	m_move(yout,0,0,nout,yout->n,rK,0,1);
	return rK;
}
