#include "LLRL.h"
/*************************************************************************
Workspace used: 11, 5, 6, 8, 9
**************************************************************************	
	The feval function is implemented in C by passing a function in the 
	function header:
	return type (* function)(argument types)
	Ex: MAT *(*fun)()
	Esta notaci� equivale a un apuntador a un apuntador de una funcion

	FEVAL Execute the specified function
	
	FEVAL(F,x1,...,xn) evaluates the function specified by a function
	handle or function name, F, at the given arguments, x1,...,xn.
	For example, if F = @foo, FEVAL(F,9.64) is the same as foo(9.64).
	
	Copyright, Andr� Jaramillo Botero, 1999. 
**************************************************************************/
/*07-12-2002-AJB: changed assignment for block move in t_res_feval */
/*07-12-2002-AJB: eliminated mem_stat_marks, added torque argument */
/*09-23-2004-AJB: el llamado a la funci� abi ten� un argumento de m�
   para soportar la funci� de torque, se elimina temp_x por ahora */
/*09-02-2006-JCA: fixed problem with 6 df, added suport for n DF and flying base
  where Vb and dVb are base velocity and acceleration. Note: Still works with fixed base*/
/*22/02/2006, JCA, AM: Solved problem with memory in accel calls*/
/*22/02/2006, JCA, JDC:Add spatial base position: X=[roll pitch yaw x y z]'*/
/*23/02/2006, JCA, deleted redundant DH and DYN matrix.*/
/*24/02/2006 JCA Solved problems with tor reg.*/
/*28/02/2006 JCA Added torque function support*/
/*16/04/2006, JCA added support for friction*/
/*26/04/2006, JCA Added sopport for robot type, and flying base as it should be*/
/*07/05/2006 JDC, Added some additional conditions*/
/*07/05/2006 JDC, Changed interface, added DYN_base*/
/*07/05/2006 JDC, Changed prototype call for AcelBase, added DYN_base*/
/*12/05/2006* JDC, Modification for Serpentine conditions */

MAT *feval(MAT *(*fun)(), int Robot_type, MAT *tau, VEC *time_vector, double time, MAT *x, MAT *res_feval, MAT *DH, MAT *DYN,VEC *grav, VEC *fext,MAT *DYN_base)
{
	static MAT *qdd=MNULL,*temp_x=MNULL,*t_res_feval=MNULL, *xb=MNULL,*vb=MNULL,*ab=MNULL, *tor=MNULL, *Fbase=MNULL;
	int i;
	//m_resize_vars(1,1,&xb,&vb,&ab,&Fbase,NULL);
	//mem_stat_reg_vars(0,TYPE_MAT,&xb,&vb,&ab,&Fbase,NULL);
		
// validate inputs (sizes and initialization) 
	if (x == MNULL) 
		error(E_NULL,"feval (intputs)\n");
	if ((res_feval == MNULL) || (res_feval->m != 2*DH->m) || (res_feval->n != 1)) 
		res_feval=m_resize(res_feval,2*DH->m,1);
	if(grav==VNULL)
		error(E_NULL,"NULL gravity");
	else
		if (grav->dim != 6)
			error(E_SIZES,"Wrong gravity vector");
	
	if(fext==VNULL)
		error(E_NULL,"NULL fext");
	else
		if (fext->dim != 6)
			error(E_SIZES,"Wrong fext vector");
	
	temp_x = m_resize(temp_x,x->n,x->m);
	qdd = m_resize(qdd,x->n,1*DH->m);
	t_res_feval = m_resize(t_res_feval,1,2*DH->m);
	mem_stat_reg_vars(0,TYPE_MAT,&temp_x,&qdd,&t_res_feval,NULL);//&temp,, &temp_Q_Qd

	m_transp(x,temp_x);
		
	//Flying or Serpentine Robots
	if(Robot_type==2 || Robot_type==4)
	{
		m_resize_vars(6,1,&xb,&vb,&ab,NULL);
		Fbase=m_resize(Fbase,6,1);
		mem_stat_reg_vars(0,TYPE_MAT,&xb,&vb,&ab,&Fbase,NULL);
	}
		
	if (tau == MNULL)
	{
		tor=m_resize(tor,1,DH->m);
		MEM_STAT_REG(tor,TYPE_MAT);
		m_zero(tor);
	}
	else
	{
		if ((Robot_type!=2 && Robot_type!=4) && ((tau->n) != (DH->m)))
			error(E_SIZES,"feval (Torque function)");
		else
		{
			tor=m_resize(tor,1,tau->n);
			MEM_STAT_REG(tor,TYPE_MAT);

			mem_stat_mark(11);
			(*fun)(DH,DYN,tau,time_vector,time,tor);
			mem_stat_free(11);
		}
	}
		
	//Flying base robot	
	if(Robot_type==2)
	{
		
		for(i=0;i<6;i++)
		{
			xb->me[i][0]=temp_x->me[0][temp_x->n/2-6+i];
			vb->me[i][0]=temp_x->me[0][temp_x->n-6+i];
		}
			
		m_move(temp_x,0,temp_x->n/2,1,temp_x->n/2-6,temp_x,0,temp_x->n/2-6);
		temp_x=m_resize(temp_x,1,temp_x->n-12);
		
		for(i=0;i<6;i++)
			Fbase->me[i][0]=tor->me[0][tor->n-6+i];
		
		mem_stat_mark(11);
		ab=AcelBase(Robot_type,DH,DYN,DYN_base,grav,temp_x,Fbase,ab);
		mem_stat_free(11);
		
		tor=m_resize(tor,1,tor->n-6);	
	}
	
	//Serpentine robot	
	if(Robot_type==4)
	{
 		
		for(i=0;i<6;i++)
		{
			xb->me[i][0]=temp_x->me[0][temp_x->n/2-6+i];
			vb->me[i][0]=temp_x->me[0][temp_x->n-6+i];
		}
		
		m_move(temp_x,0,temp_x->n/2,1,temp_x->n/2-6,temp_x,0,temp_x->n/2-6);
		temp_x=m_resize(temp_x,1,temp_x->n-12);

		
		//Computing Friction Force of the system
		mem_stat_mark(11);
		Fbase=Friction_surface(DH,DYN,DYN_base,temp_x,xb,vb,Fbase);
		mem_stat_free(11);
		
		//Obtain Spatial Base Acceleration due to Friction Force
		mem_stat_mark(11);
		ab=AcelBase(Robot_type,DH,DYN,DYN_base,grav,temp_x,Fbase,ab);
		mem_stat_free(11);

		tor=m_resize(tor,1,tor->n-6);	
	}

	
		mem_stat_mark(11);
		qdd=accel(Robot_type,DH,DYN,temp_x,tor,qdd,xb,vb,ab,grav,fext);
		mem_stat_free(11);

	
	//Fix base or floating base
	if((Robot_type==1) || (Robot_type==3))
	{
		for(i=0;i<DH->m;i++)
			t_res_feval->me[0][i]=x->me[DH->m+i][0];
		m_move(qdd,0,0,1,DH->m,t_res_feval,0,DH->m);
		res_feval = m_resize(res_feval,2*DH->m,1);
		m_transp(t_res_feval,res_feval);
	}
	
	//Flying base -Serpentine robots
	if(Robot_type==2 || Robot_type==4)
	{
		t_res_feval = m_resize(t_res_feval,1,2*DH->m+12);

		for(i=0;i<DH->m+6;i++)
			t_res_feval->me[0][i]=x->me[DH->m+6+i][0];
		
		m_move(qdd,0,0,1,DH->m,t_res_feval,0,DH->m+6);
		res_feval = m_resize(res_feval,2*DH->m+12,1);
		m_transp(t_res_feval,res_feval);
		m_move(ab,0,0,6,1,res_feval,2*DH->m+6,0);
	}
	return res_feval;
}
