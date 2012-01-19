/*****************************************************************************************************************
	Inverse dynamics parallel Solution, log2(n)
	Thesis: Modelling and inverse dynamics solution for serpentine robots using strictly parallel algorithms
	Proyect: Time lower bound modelling and simulation of complex rigid body systems.
		Andres Jaramillo Botero ajaramil@puj.edu.co
		Juan Camilo Acosta jcacosta@puj.edu.co
		Julian David Colorado jdcolorado@puj.edu.co
		Juan Manuel Florez jmflorez@puj.edu.co
	Copyright Robotics and Automation Group, 2006
	Pontificia Universidad Javeriana, Cali, Colombia
******************************************************************************************************************/
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "../Meschach/matrix.h"
#include "../Meschach/matrix2.h"

VEC *tr2rpy(MAT *T,VEC *E);
void euler(VEC *columna, MAT *r, VEC *p);
void Rotational(float alfa,float a, float teta, float d, MAT *Rot,int k);//, VEC* p);
void skew_symetric(MAT *M, MAT *v_temp,int k);
void skew_symetric2(MAT *M, VEC *v_temp,int k);
VEC *com_step(int my_id, int num_pro,VEC *vec_send, MAT *mat_send);
VEC *com_step_inv(int my_id, int num_pro,VEC *vec_send, MAT *mat_send);
VEC *com_step2(int my_id, int num_pro,VEC *vec_send, VEC *result);
VEC *ff_ne(MAT *DH,MAT *DYN,MAT *Q_Qd_Qdd,VEC *gravity,VEC *fext,VEC *Torques,MAT *X, MAT *Vb, MAT *dVb, int my_id,int num_pro);

/*************************************************************************	
	Conversion from homogeneous transform to euler angles (y-convention)

	E=tr2rpy(T,E)
  
	T=4x4 homogeneous transformation matrix representing rotation.
	E=vector defining [roll,pitch,jaw] rotation angles (in radians)
	
	Copyright, Andr� Jaramillo Botero, 1999. 
***************************************************************************/

VEC *tr2rpy(MAT *T,VEC *E)
{
	/* validate inputs (sizes and initializations) */
	if (T == (MAT *)MNULL) 
		error(E_NULL, "tr2rpy (input matrix)\n");
	if (T->m != 4 || T->n != 4)	
		error(E_SIZES,"tr2rpy (input matrix)\n");
	if (E == (VEC *)VNULL || E->dim <= 3) E=v_resize(E,3);
	
	if (fabs(T->me[0][0]) < D_MACHEPS && fabs(T->me[1][0]) < D_MACHEPS)
	{
		E->ve[0]=0;					// roll 
		E->ve[1]=atan2(-T->me[2][0],T->me[0][0]);	// pitch 
		E->ve[2]=atan2(-T->me[1][2],T->me[1][1]);	// jaw 
	}
	else
	{
		E->ve[0]=atan2(T->me[1][0],T->me[0][0]);				// roll 
		E->ve[1]=atan2(-T->me[2][0],cos(E->ve[0])*T->me[0][0]+
				sin(E->ve[0])*T->me[1][0]);				// pitch
		E->ve[2]=atan2(sin(E->ve[0])*T->me[0][2]-cos(E->ve[0])*T->me[1][2],
			       -sin(E->ve[0])*T->me[0][1]+cos(E->ve[0])*T->me[1][1]); 	// jaw 
	}
	
	return E;
}

//**************************************************************************/
//Returns the Homogeneus Transformation from [roll pitch yaw X  Y  Z] vector'
//Obtain [r p]:Matrix Rotation, Position Vector
void euler(VEC *columna, MAT *r, VEC *p)
{
	float roll,pitch,yaw;
	roll=columna->ve[0];
	pitch=columna->ve[1];
	yaw=columna->ve[2];
	
	//Orientation Matrix
	r->me[0][0]=cos(roll)*cos(pitch);
	r->me[0][1]=cos(roll)*sin(pitch)*sin(yaw)-sin(roll)*cos(yaw);
	r->me[0][2]=cos(roll)*sin(pitch)*cos(yaw)+sin(roll)*sin(yaw);
	r->me[1][0]=sin(roll)*cos(pitch);
	r->me[1][1]=sin(roll)*sin(pitch)*sin(yaw)+cos(roll)*cos(yaw);
	r->me[1][2]=sin(roll)*sin(pitch)*cos(yaw)-cos(roll)*sin(yaw);
	r->me[2][0]=-sin(pitch);
	r->me[2][1]=cos(pitch)*sin(yaw);
	r->me[2][2]=cos(pitch)*cos(yaw);
	
	//Position Vector
	p->ve[0]=columna->ve[3];
	p->ve[1]=columna->ve[4];
	p->ve[2]=columna->ve[5];
}

/*************************************************************************
* Free Flying base Functions.
**************************************************************************/
//Returns the Rotational 3*3 Matrix
//alfa ,a teta,d are the DH parameters 
//Obtain Rot:Homogeneus Transformation Matrix
/*05/2006 JCA: fixed to avoid m_zero of ROT*/
void Rotational(float alfa,float a, float teta, float d, MAT *Rot,int k)
{
	Rot->me[0][0+3*k]=cos(teta);
	Rot->me[1][0+3*k]=sin(teta);
	Rot->me[2][0+3*k]=0;
		
	Rot->me[0][1+3*k]=-sin(teta)*cos(alfa);
	Rot->me[1][1+3*k]=cos(alfa)*cos(teta);
	Rot->me[2][1+3*k]=sin(alfa);
		
	Rot->me[0][2+3*k]=sin(alfa)*sin(teta);
	Rot->me[1][2+3*k]=-sin(alfa)*cos(teta);
	Rot->me[2][2+3*k]=cos(alfa);
}


//************************************************************
//Returns skew symetric matrix:
//Obtain values from v_temp and Obtain M skew symetric matrix
//Function form_matrix returns block matrix mtempo 
/*05/2006 JCA: fixed to avoid m_zero of M*/
/*06/2006 JCA: Block matrix support*/
void skew_symetric(MAT *M, MAT *v_temp,int k)
{
	M->me[0][0]=0;
	M->me[0][1]=-v_temp->me[2][k];
	M->me[0][2]=v_temp->me[1][k];
	
	M->me[1][0]=v_temp->me[2][k];
	M->me[1][1]=0;
	M->me[1][2]=-v_temp->me[0][k];
	
	M->me[2][0]=-v_temp->me[1][k];
	M->me[2][1]=v_temp->me[0][k];
	M->me[2][2]=0;
}
//************************************************************
//Returns skew symetric matrix:
//Obtain values from v_temp and Obtain M skew symetric matrix
//Function form_matrix returns block matrix mtempo 
/*05/2006 JCA: fixed to avoid m_zero of M*/
/*06/2006 JCA: Block matrix support*/
void skew_symetric2(MAT *M, VEC *v_temp,int k)
{
	M->me[0][0]=0;
	M->me[0][1]=-v_temp->ve[2+6*k];
	M->me[0][2]=v_temp->ve[1+6*k];
	
	M->me[1][0]=v_temp->ve[2+6*k];
	M->me[1][1]=0;
	M->me[1][2]=-v_temp->ve[0+6*k];
	
	M->me[2][0]=-v_temp->ve[1+6*k];
	M->me[2][1]=v_temp->ve[0+6*k];
	M->me[2][2]=0;
}
/************************************************************************************************************
 Complete communication step for velocities and accelerations. o(log2(n))
 Its a forward propagation
 ************************************************************************************************************
 returns velocities or accelerations depending on the call.
 %my_id: process id
 %num_pro: number of process
 %vec_send: vector to send
 %mat_send: matrix to send
 *************************************************************************************************************/
/*05/2006 JCA: Deleted redundant m_copy and v_copy and changed interface*/
/*12/06/2006 JCA: changed 3log2(n)*points-->3*log2(n), BLOCK MATRIX AND VECTORS!*/
/*20/06/2006 JCA modifications with mat_send in order to keep mat_send original values*/
/*20/06/2006 JCA modifications to avoid m_copy for matrix mult.*/
VEC *com_step(int my_id, int num_pro,VEC *vec_send, MAT *mat_send)
{
	double buffer[(int)(vec_send->dim+6*mat_send->n+1)];
	int etapa,i,j,k,np;
	MPI_Status  Status;
	MPI_Request Request;
	/*VECTOR DECLARATION*/
	static VEC *vmpi_temp=VNULL, *vec_recv=VNULL;
	
	/*MATRIX DECLARATION*/
	static MAT *mat_send3=MNULL, *mat_recv=MNULL, *mat_send2=MNULL;
	
	m_resize_vars(6,mat_send->n, &mat_send3, &mat_recv, &mat_send2, NULL);
	mem_stat_reg_vars(0,TYPE_MAT, &mat_send2, &mat_send3, &mat_recv,NULL);
	
	
	v_resize_vars(mat_send->n, &vmpi_temp, &vec_recv,NULL);
	mem_stat_reg_vars(0,TYPE_VEC, &vmpi_temp, &vec_recv,NULL);

	np=vec_send->dim/6;
	
	for(etapa=1;etapa<=(log10(num_pro-1)/log10(2)+1.001);etapa++)
	{
		//If the process must send data in this communication step
		if((int)(my_id+pow(2,etapa-1))<num_pro)
		{
			//If the process has already calculated speed or accel only sends a vector
			if(my_id<pow(2,etapa-1))
			{
				//copy vector in buffer
				for(i=0;i<vec_send->dim;i++)
					buffer[i]=vec_send->ve[i];
				// NON blocking send
				MPI_Isend(buffer,vec_send->dim,MPI_DOUBLE,(my_id+pow(2,etapa-1)),1,MPI_COMM_WORLD, &Request);
				
				//mpi send of size 1*MPI_DOUBLE to test latency of the cluster
				//MPI_Isend(buffer,1,MPI_DOUBLE,(my_id+pow(2,etapa-1)),1,MPI_COMM_WORLD, &Request);
			}
			//If the process hasnt calculated speed or accel sends a vector and a matrix
			else
			{
				//copy vector in buffer
				for(i=0;i<vec_send->dim;i++)
					buffer[i]=vec_send->ve[i];
				//copy matrix in buffer depending on comm step. 3 matrix exist to avois m_copy in calculus step
				if(etapa==1)
					for(i=0;i<6;i++)
						for(j=0;j<mat_send->n;j++)
							buffer[vec_send->dim+mat_send->n*i+j]=mat_send->me[i][j];
				else
					if(etapa%2==0)
						for(i=0;i<6;i++)
							for(j=0;j<mat_send->n;j++)
								buffer[vec_send->dim+mat_send->n*i+j]=mat_send2->me[i][j];
					else
						for(i=0;i<6;i++)
							for(j=0;j<mat_send->n;j++)
								buffer[vec_send->dim+mat_send->n*i+j]=mat_send3->me[i][j];
				// NON blocking send
				MPI_Isend(buffer,vec_send->dim+6*mat_send->n,MPI_DOUBLE,my_id+pow(2,etapa-1),1,MPI_COMM_WORLD, &Request);
				//mpi send of size 1*MPI_DOUBLE to test latency of the cluster
				//MPI_Isend(buffer,1,MPI_DOUBLE,my_id+pow(2,etapa-1),1,MPI_COMM_WORLD, &Request);
			}
		}
		//If the process must receive data in this communication step
		if(my_id>=(pow(2,etapa-1)))
		{
			//If the process is receiving a calculated speed or accel vector
			if(my_id<pow(2,etapa))
			{
				//Blocking receive
				MPI_Recv(buffer,vec_send->dim,MPI_DOUBLE,(my_id-pow(2,etapa-1)),1,MPI_COMM_WORLD, &Status);
				
				//mpi recv of size 1*MPI_DOUBLE to test latency of the cluster
				//MPI_Recv(buffer,1,MPI_DOUBLE,(my_id-pow(2,etapa-1)),1,MPI_COMM_WORLD, &Status);
				
				//copy buffer in vector
				for(i=0;i<vec_send->dim;i++)
					vec_recv->ve[i]=buffer[i];
			}
			//If the process is receiving a vector and a matrix
			else
			{
				//Blocking receive
				MPI_Recv(buffer,vec_send->dim+6*mat_send->n,MPI_DOUBLE,my_id-pow(2,etapa-1),1,MPI_COMM_WORLD, &Status);
				
				//mpi recv of size 1*MPI_DOUBLE to test latency of the cluster
				//MPI_Recv(buffer,1,MPI_DOUBLE,my_id-pow(2,etapa-1),1,MPI_COMM_WORLD, &Status);
				
				//copy buffer in vector
				for(i=0;i<vec_send->dim;i++)
					vec_recv->ve[i]=buffer[i];
				
				//copy buffer in matrix
				for(i=0;i<6;i++)
					for(j=0;j<mat_send->n;j++)
						mat_recv->me[i][j]=buffer[vec_send->dim+mat_send->n*i+j];
			}
		}
		//Calculus step		
		if(my_id>0)
		{
			//Speed or accel process
			if(etapa==(int)(((log10(my_id)/log10(2))+1.001)))
			{
				//Force=mat_mia*vec_recv+vmpi_temp2;
				//speed or accel=mat_mia*vec_recv+vmpi_temp2;
				//mv_mlt(mat_send,vec_recv,vmpi_temp);
				for(i=0;i<6;i++)
					for(j=0;j<mat_send->n/6;j++)
					{
						if(etapa==1)
    						{
							vmpi_temp->ve[6*j+i]=mat_send->me[i][6*j]*vec_recv->ve[6*j];
							for(k=1;k<6;k++)
								vmpi_temp->ve[6*j+i]=mat_send->me[i][6*j+k]*vec_recv->ve[6*j+k]+vmpi_temp->ve[6*j+i];
							}
						else
						{
							if(etapa%2==0)
							{
								vmpi_temp->ve[6*j+i]=mat_send2->me[i][6*j]*vec_recv->ve[6*j];
								for(k=1;k<6;k++)
									vmpi_temp->ve[6*j+i]=mat_send2->me[i][6*j+k]*vec_recv->ve[6*j+k]+vmpi_temp->ve[6*j+i];
							}
							else
							{
								vmpi_temp->ve[6*j+i]=mat_send3->me[i][6*j]*vec_recv->ve[6*j];
								for(k=1;k<6;k++)
									vmpi_temp->ve[6*j+i]=mat_send3->me[i][6*j+k]*vec_recv->ve[6*j+k]+vmpi_temp->ve[6*j+i];
							}
						 
						}
							//v_add(vmpi_temp,vec_send,vec_send);
						vec_send->ve[6*j+i]=vec_send->ve[6*j+i]+vmpi_temp->ve[6*j+i];
					}
					
				//v_add(vmpi_temp,vec_send,vec_send);
				
			}
			//matrix and vector calculus
			else
			{
				if(etapa<(int)(log10(my_id)/log10(2))+1.001)
				{
					//vec_send=mat_mia*vec_recv+vmpi_temp2;
					//mv_mlt(mat_send,vec_recv,vmpi_temp);
					for(i=0;i<6;i++)
						for(j=0;j<mat_send->n/6;j++)
						{
							if(etapa==1)
							{
								vmpi_temp->ve[6*j+i]=mat_send->me[i][6*j]*vec_recv->ve[6*j];
								for(k=1;k<6;k++)
									vmpi_temp->ve[6*j+i]=mat_send->me[i][6*j+k]*vec_recv->ve[6*j+k]+vmpi_temp->ve[6*j+i];
							}
							else
							{
								if(etapa%2==0)
								{
									vmpi_temp->ve[6*j+i]=mat_send2->me[i][6*j]*vec_recv->ve[6*j];
									for(k=1;k<6;k++)
										vmpi_temp->ve[6*j+i]=mat_send2->me[i][6*j+k]*vec_recv->ve[6*j+k]+vmpi_temp->ve[6*j+i];
								}
								else
								{
									vmpi_temp->ve[6*j+i]=mat_send3->me[i][6*j]*vec_recv->ve[6*j];
									for(k=1;k<6;k++)
										vmpi_temp->ve[6*j+i]=mat_send3->me[i][6*j+k]*vec_recv->ve[6*j+k]+vmpi_temp->ve[6*j+i];
								}
						 	}
							//v_add(vmpi_temp,vec_send,vec_send);
							vec_send->ve[6*j+i]=vec_send->ve[6*j+i]+vmpi_temp->ve[6*j+i];
						}
					
					//mat_send=mat_mia*mat_recv
					/*if(etapa!=1)
						m_copy(mat_send,mmpi_temp);*/
					//m_mlt(mmpi_temp,mat_recv,mat_send);
					for(i=0;i<6;i++)
						for(j=0;j<mat_send->n/6;j++)
							for(k=0;k<6;k++)
							{
								if(etapa==1)
									for(k=0;k<6;k++)
									mat_send2->me[i][6*j+k]=mat_send->me[i][6*j+0]*mat_recv->me[0][6*j+k]
												+mat_send->me[i][6*j+1]*mat_recv->me[1][6*j+k]
												+mat_send->me[i][6*j+2]*mat_recv->me[2][6*j+k]
												+mat_send->me[i][6*j+3]*mat_recv->me[3][6*j+k]
												+mat_send->me[i][6*j+4]*mat_recv->me[4][6*j+k]
												+mat_send->me[i][6*j+5]*mat_recv->me[5][6*j+k];
								else
								{
									if(etapa%2==0)
										mat_send3->me[i][6*j+k]=mat_send2->me[i][6*j+0]*mat_recv->me[0][6*j+k]
												+mat_send2->me[i][6*j+1]*mat_recv->me[1][6*j+k]
												+mat_send2->me[i][6*j+2]*mat_recv->me[2][6*j+k]
												+mat_send2->me[i][6*j+3]*mat_recv->me[3][6*j+k]
												+mat_send2->me[i][6*j+4]*mat_recv->me[4][6*j+k]
												+mat_send2->me[i][6*j+5]*mat_recv->me[5][6*j+k];
									else
										mat_send2->me[i][6*j+k]=mat_send3->me[i][6*j+0]*mat_recv->me[0][6*j+k]
												+mat_send3->me[i][6*j+1]*mat_recv->me[1][6*j+k]
												+mat_send3->me[i][6*j+2]*mat_recv->me[2][6*j+k]
												+mat_send3->me[i][6*j+3]*mat_recv->me[3][6*j+k]
												+mat_send3->me[i][6*j+4]*mat_recv->me[4][6*j+k]
												+mat_send3->me[i][6*j+5]*mat_recv->me[5][6*j+k];
								}
							}
				}
			}
		}
		//if((int)(my_id+pow(2,etapa-1))<num_pro)
			//MPI_Wait(&Request,&Status);
		MPI_Barrier(MPI_COMM_WORLD);
	}
	return(vec_send);
}
/************************************************************************************************************
 Complete communication step for Forces. o(log2(n))
Its a backward propagation
 ************************************************************************************************************
 returns velocities or accelerations depending on the call.
 %my_id: process id
 %num_pro: number of process
 %vec_send: vector to send
 %mat_send: matrix to send
 *************************************************************************************************************/
/*05/2006 JCA: Deleted redundant m_copy and v_copy and changed interface*/
/*12/06/2006 JCA: changed 3log2(n)*points-->3*log2(n)*/
/*20/06/2006 JCA modifications with mat_send in order to keep mat_send original values*/
/*20/06/2006 JCA modifications to avoid m_copy for matrix mult.*/
VEC *com_step_inv(int my_id, int num_pro,VEC *vec_send, MAT *mat_send)
{
	double buffer[(int)(vec_send->dim+6*mat_send->n+1)];
	int etapa,i,j,k,np;
	MPI_Status  Status;
	MPI_Request Request;
	
	/*VECTOR DECLARATION*/
	static VEC *vmpi_temp=VNULL, *vec_recv=VNULL;
	
	/*MATRIX DECLARATION*/
	static MAT *mat_send3=MNULL, *mat_recv=MNULL, *mat_send2=MNULL;
	
	m_resize_vars(6,mat_send->n, &mat_send3, &mat_recv, &mat_send2, NULL);
	mem_stat_reg_vars(0,TYPE_MAT, &mat_send2, &mat_send3, &mat_recv,NULL);
	
	v_resize_vars(mat_send->n, &vmpi_temp, &vec_recv,NULL);
	mem_stat_reg_vars(0,TYPE_VEC, &vmpi_temp, &vec_recv,NULL);

	for(etapa=1;etapa<=(log10(num_pro-1)/log10(2)+1.001);etapa++)
	{	
		//If the process must send data in this communication ste
		if((int)(my_id-pow(2,etapa-1))>=0)
		{	
			//If the process has already calculated forces only sends a vector
			if(my_id>=num_pro-pow(2,etapa-1))
			{
				//copy vector in buffer
				for(i=0;i<vec_send->dim;i++)
					buffer[i]=vec_send->ve[i];
				// NON blocking send
				MPI_Isend(buffer,vec_send->dim,MPI_DOUBLE,(my_id-pow(2,etapa-1)),1,MPI_COMM_WORLD, &Request);
				
				//mpi send of size 1*MPI_DOUBLE to test latency of the cluster
				//MPI_Isend(buffer,1,MPI_DOUBLE,(my_id-pow(2,etapa-1)),1,MPI_COMM_WORLD, &Request);
			}
			//If the process hasnt calculated forces sends a vector and a matrix
			else
			{
				//copy vector in buffer
				for(i=0;i<vec_send->dim;i++)
					buffer[i]=vec_send->ve[i];
				//copy matrix in buffer depending on comm step. 3 matrix exist to avois m_copy in calculus step
				if(etapa==1)
					for(i=0;i<6;i++)
						for(j=0;j<mat_send->n;j++)
							buffer[vec_send->dim+mat_send->n*i+j]=mat_send->me[i][j];
				else
					if(etapa%2==0)
						for(i=0;i<6;i++)
							for(j=0;j<mat_send->n;j++)
								buffer[vec_send->dim+mat_send->n*i+j]=mat_send2->me[i][j];
				else
					for(i=0;i<6;i++)
						for(j=0;j<mat_send->n;j++)
							buffer[vec_send->dim+mat_send->n*i+j]=mat_send3->me[i][j];
				// NON blocking send		
				MPI_Isend(buffer,vec_send->dim+6*mat_send->n,MPI_DOUBLE,my_id-pow(2,etapa-1),1,MPI_COMM_WORLD, &Request);
				
				//mpi send of size 1*MPI_DOUBLE to test latency of the cluster
				//MPI_Isend(buffer,1,MPI_DOUBLE,my_id-pow(2,etapa-1),1,MPI_COMM_WORLD, &Request);
			}
		}
		//If the process must receive data in this communication step
		if(my_id<(num_pro-pow(2,etapa-1)))
		{
			//If the process is receiving a calculated speed or accel vector
			if(my_id>=(num_pro-pow(2,etapa)))
			{
				//Blocking receive
				MPI_Recv(buffer,vec_send->dim,MPI_DOUBLE,(my_id+pow(2,etapa-1)),1,MPI_COMM_WORLD, &Status);
				
				//mpi recv of size 1*MPI_DOUBLE to test latency of the cluster
				//MPI_Recv(buffer,1,MPI_DOUBLE,(my_id+pow(2,etapa-1)),1,MPI_COMM_WORLD, &Status);
				
				//copy buffer in vector
				for(i=0;i<vec_send->dim;i++)
					vec_recv->ve[i]=buffer[i];
			}
			//If the process is receiving a vector and a matrix
			else
			{
				//Blocking receive
				MPI_Recv(buffer,vec_send->dim+6*mat_send->n,MPI_DOUBLE,my_id+pow(2,etapa-1),1,MPI_COMM_WORLD, &Status);
				
				//mpi recv of size 1*MPI_DOUBLE to test latency of the cluster
				//MPI_Recv(buffer,1,MPI_DOUBLE,my_id+pow(2,etapa-1),1,MPI_COMM_WORLD, &Status);
				
				//copy buffer in vector
				for(i=0;i<vec_send->dim;i++)
					vec_recv->ve[i]=buffer[i];
				
				//copy buffer in matrix
				for(i=0;i<6;i++)
					for(j=0;j<mat_send->n;j++)
						mat_recv->me[i][j]=buffer[vec_send->dim+mat_send->n*i+j];
			}
		}
		//Calculus step	
		if(my_id<(num_pro-1))
		{
			//Speed or accel process
			if(etapa==(int)(((log10(num_pro-1-my_id)/log10(2))+1.001)))
			{
				//Force=mat_mia*vec_recv+vmpi_temp2;
				//mv_mlt(mat_send,vec_recv,vmpi_temp);
				//v_add(vmpi_temp,vec_send,vec_send);
				for(i=0;i<6;i++)
					for(j=0;j<mat_send->n/6;j++)
					{
						if(etapa==1)
						{
							vmpi_temp->ve[6*j+i]=mat_send->me[i][6*j]*vec_recv->ve[6*j];
							for(k=1;k<6;k++)
								vmpi_temp->ve[6*j+i]=mat_send->me[i][6*j+k]*vec_recv->ve[6*j+k]+vmpi_temp->ve[6*j+i];
						}
						else
						{
							if(etapa%2==0)
							{
								vmpi_temp->ve[6*j+i]=mat_send2->me[i][6*j]*vec_recv->ve[6*j];
								for(k=1;k<6;k++)
									vmpi_temp->ve[6*j+i]=mat_send2->me[i][6*j+k]*vec_recv->ve[6*j+k]+vmpi_temp->ve[6*j+i];
							}
							else
							{
								vmpi_temp->ve[6*j+i]=mat_send3->me[i][6*j]*vec_recv->ve[6*j];
								for(k=1;k<6;k++)
									vmpi_temp->ve[6*j+i]=mat_send3->me[i][6*j+k]*vec_recv->ve[6*j+k]+vmpi_temp->ve[6*j+i];
							}
						}
						//v_add(vmpi_temp,vec_send,vec_send);
						vec_send->ve[6*j+i]=vec_send->ve[6*j+i]+vmpi_temp->ve[6*j+i];
					}
			}
			//matrix and vector calculus
			else
			{
				if(etapa<(int)(log10(num_pro-1-my_id)/log10(2))+1.001)
				{
					//vec_send=mat_mia*vec_recv+vmpi_temp2;
					//mv_mlt(mat_send,vec_recv,vmpi_temp);
					for(i=0;i<6;i++)
						for(j=0;j<mat_send->n/6;j++)
					{
						if(etapa==1)
						{
							vmpi_temp->ve[6*j+i]=mat_send->me[i][6*j]*vec_recv->ve[6*j];
							for(k=1;k<6;k++)
								vmpi_temp->ve[6*j+i]=mat_send->me[i][6*j+k]*vec_recv->ve[6*j+k]+vmpi_temp->ve[6*j+i];
						}
						else
						{
							if(etapa%2==0)
							{
								vmpi_temp->ve[6*j+i]=mat_send2->me[i][6*j]*vec_recv->ve[6*j];
								for(k=1;k<6;k++)
									vmpi_temp->ve[6*j+i]=mat_send2->me[i][6*j+k]*vec_recv->ve[6*j+k]+vmpi_temp->ve[6*j+i];
							}
							else
							{
								vmpi_temp->ve[6*j+i]=mat_send3->me[i][6*j]*vec_recv->ve[6*j];
								for(k=1;k<6;k++)
									vmpi_temp->ve[6*j+i]=mat_send3->me[i][6*j+k]*vec_recv->ve[6*j+k]+vmpi_temp->ve[6*j+i];
							}
						}
							//v_add(vmpi_temp,vec_send,vec_send);
						vec_send->ve[6*j+i]=vec_send->ve[6*j+i]+vmpi_temp->ve[6*j+i];
					}
					
					//mat_send=mat_mia*mat_recv
					/*if(etapa!=1)
					m_copy(mat_send,mmpi_temp);*/
					//m_mlt(mmpi_temp,mat_recv,mat_send);
					for(i=0;i<6;i++)
						for(j=0;j<mat_send->n/6;j++)
							for(k=0;k<6;k++)
							{
								if(etapa==1)
									for(k=0;k<6;k++)
										mat_send2->me[i][6*j+k]=mat_send->me[i][6*j+0]*mat_recv->me[0][6*j+k]
												+mat_send->me[i][6*j+1]*mat_recv->me[1][6*j+k]
												+mat_send->me[i][6*j+2]*mat_recv->me[2][6*j+k]
												+mat_send->me[i][6*j+3]*mat_recv->me[3][6*j+k]
												+mat_send->me[i][6*j+4]*mat_recv->me[4][6*j+k]
												+mat_send->me[i][6*j+5]*mat_recv->me[5][6*j+k];
								else
								{
									if(etapa%2==0)
										mat_send3->me[i][6*j+k]=mat_send2->me[i][6*j+0]*mat_recv->me[0][6*j+k]
												+mat_send2->me[i][6*j+1]*mat_recv->me[1][6*j+k]
												+mat_send2->me[i][6*j+2]*mat_recv->me[2][6*j+k]
												+mat_send2->me[i][6*j+3]*mat_recv->me[3][6*j+k]
												+mat_send2->me[i][6*j+4]*mat_recv->me[4][6*j+k]
												+mat_send2->me[i][6*j+5]*mat_recv->me[5][6*j+k];
									else
										mat_send2->me[i][6*j+k]=mat_send3->me[i][6*j+0]*mat_recv->me[0][6*j+k]
												+mat_send3->me[i][6*j+1]*mat_recv->me[1][6*j+k]
												+mat_send3->me[i][6*j+2]*mat_recv->me[2][6*j+k]
												+mat_send3->me[i][6*j+3]*mat_recv->me[3][6*j+k]
												+mat_send3->me[i][6*j+4]*mat_recv->me[4][6*j+k]
												+mat_send3->me[i][6*j+5]*mat_recv->me[5][6*j+k];
								}
							}
				}
			}
		}
		//if((int)(my_id-pow(2,etapa-1))>=0)
			//MPI_Wait(&Request,&Status);
		MPI_Barrier(MPI_COMM_WORLD);
	}
	return(vec_send);
}
/************************************************************************************************************
 Communication step for propagating a 6 dimentional vector from i process to i+1 process
 ************************************************************************************************************
 returns velocities or accelerations depending on the call.
 %my_id: process id
 %num_pro: number of process
 %vec_send: vector to send
 %result: result vector
 *************************************************************************************************************/
/*05/2006 JCA: Deleted redundant m_copy and v_copy*/
/*02/06/2006 JCA: deleted v_zero result*/
/*12/06/2006 JCA: changed 3log2(n)*points-->3*log2(n)*/
VEC *com_step2(int my_id, int num_pro,VEC *vec_send, VEC *result)
{
	double buffer[vec_send->dim];
	int i;
	MPI_Status Status;
	MPI_Request Request;
	
	if(my_id<(num_pro-1))
	{
		//copy vector in buffer
		for(i=0;i<vec_send->dim;i++)
			buffer[i]=vec_send->ve[i];
		MPI_Isend(buffer,vec_send->dim,MPI_DOUBLE,(my_id+1),1,MPI_COMM_WORLD, &Request);
		//mpi recv of size 1*MPI_DOUBLE to test latency of the cluster
		//MPI_Isend(buffer,1,MPI_DOUBLE,(my_id+1),1,MPI_COMM_WORLD, &Request);
	}
	if(my_id>0)
	{
		MPI_Recv(buffer,vec_send->dim,MPI_DOUBLE,(my_id-1),1,MPI_COMM_WORLD, &Status);
		//mpi recv of size 1*MPI_DOUBLE to test latency of the cluster
		//MPI_Recv(buffer,1,MPI_DOUBLE,(my_id-1),1,MPI_COMM_WORLD, &Status);
		//copy buffer in vector
		for(i=0;i<vec_send->dim;i++)
			result->ve[i]=buffer[i];
	}
	MPI_Barrier(MPI_COMM_WORLD);
	return result;
}
/*****************************************************************************
 Compute inverse dynamics with an o(log2(n)) parallel alghotithm. 
 
 TAU = ff_ne(DH,DYN,Q_Qd_Qdd,grav,fext, X, V, dV, Kinetic, Potential, E_Conserv)
 DH=Denavit Hartenberg parameter matrix for manipulator
 DYN=dynamic parameter matrix in the form [m,sx,sy,sz,Ixx,Iyy,Izz,Ixy,Ixz,Iyz,Jm,G,B,Tc+,Tc-]
 Q_Qd_Qdd=[Q QD QDD] (manipulator state)
 X = base spatial position
 V = Base spatial velocity.
 dV = Base spatial acceleration.
 kinetic=kinetic energy
 potential=potential energy
 E_Conserv=engergy conservation
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
*		Andr� Jaramillo Botero, ajaramil@puj.edu.co
*		Juan Camilo Acosta, jcacosta@puj.edu.co 
*		Juli� David Colorado, jdcolorado@puj.edu.co
*		Juan Manuel Florez jmflroez@puj.edu.co
******************************************************************************/
/*22/02/2006 JDC, Added Euler Function to Calcule matrix Rotation r and position vector p from Inertial frame to first body*/
/*22/02/2006 JCA  Add X entrace parameter: Spatial Base Position, Solved problem with robots with prismatic and rotational art.*/
/*24/02/2006 JDC, JCA added support for kinetic and potential energy computation.  */
/*01/03/2006 JDC, added support for Motor friction */
/*23/03/2006 JCA  CHanged interface kinetic,potential,E_Conserv*/ 
/*27/03/2006 JDC, added simple friction Model to support Friction with contact Surface - Serpentine Robots */
/*28/03/2006 JDC, Modification of Power Kinetics*/
/*03/04/2006 JCA: Modified to work in parallel mode, with velo and accel.*/
/*19/04/2006 JCA: works complete*/
/*12/06/2006 JCA: changed 3log2(n)*points-->3*log2(n)*/
VEC *ff_ne(MAT *DH,MAT *DYN,MAT *Q_Qd_Qdd,VEC *GRAV,VEC *fext,VEC *Torques,MAT *X, MAT *Vb, MAT *dVb, int my_id, int num_pro)
{
	int position=0;
	MPI_Status *Status;
	MPI_Request Request;
	
	//serial variables
	static MAT *Q=MNULL, *dQ=MNULL, *d2Q=MNULL, *r=MNULL, *p=NULL, *R=MNULL, *P=MNULL, *sk=MNULL, *Skw=MNULL, *S=MNULL, *Ticm=MNULL, *Icm=MNULL, *I=MNULL, *Itemp=MNULL, *identidad=MNULL, *rTicm=MNULL, *Ji=MNULL, *massU=MNULL, *SIcm=MNULL, *sksk=MNULL, *vf=MNULL, *C=MNULL, *MC=MNULL, *Fr1=MNULL, *T1=MNULL, *T2=MNULL;
	static VEC *F=VNULL, *H=VNULL, *ws=VNULL, *s=VNULL, *ss=VNULL, *up=VNULL, *down=VNULL, *vtempo=VNULL, *v=VNULL, *term=VNULL, *rws=VNULL, *rv1=VNULL, *qddH=VNULL, *qdH=VNULL, *temp1=VNULL, *Iv=VNULL, *sotro=VNULL, *IdV=VNULL, *vxy=VNULL, *Fr6=VNULL, *Frt=VNULL;
	
	//mpi variables.
	static VEC *speed=VNULL, *accele=VNULL, *vmpi_temp2=VNULL;
	static MAT *mat_mia=MNULL;
	
	int k,n,np,i,j,cont,cont1,etapa;
	double alfa,a,teta,d,Cin,Cin_ant,Cin_all,Pot,massneg,av,acv,aa,aca,af,acf,fin;
	
	//C matrix and vector to avoid copying into buffer
	//double matc[6][6*Q_Qd_Qdd->m+*Q_Qd_Qdd->m];
	//double speedc[6*Q_Qd_Qdd->m];
	
	n=(int) DH->m;	        /* number of joints  */
	np=(int) Q_Qd_Qdd->m;		/* number of trajectory points */
	Torques = v_resize(Torques,np);

	/* validate inputs (sizes and initializations) */
	if (DH == (MAT *)MNULL || DYN == (MAT *)MNULL)
		error(E_NULL,"ne (DH and dynamic inputs)");
	if (DH->m != DYN->m)
		error(E_SIZES,"ne (DH and dynamic inputs)\n");
	if (Q_Qd_Qdd == (MAT *)MNULL)
		error(E_NULL,"Null input (state) in newton-euler");

	/* validate size of Q_Qd_Qddd: position, veloc, accel */
	if ((int) Q_Qd_Qdd->n == 3*n && (int) Q_Qd_Qdd->m == np)
	{
		m_resize_vars(np,n, &Q, &dQ, &d2Q,NULL);
		mem_stat_reg_vars(0,TYPE_MAT, &Q, &dQ, &d2Q,NULL);
		m_move(Q_Qd_Qdd,0,0,np,n,Q,0,0);
		m_move(Q_Qd_Qdd,0,n,np,n,dQ,0,0);
		m_move(Q_Qd_Qdd,0,2*n,np,n,d2Q,0,0);
	}
	else
		error(E_SIZES,"ne (incorrect manipulator state matrix)");
	
	/* Example of input format for Q,Qd,Qdd joint
	Qd=[Qd11	Qd21	Qd31	.	Qdn1	timestep
	Qd12	Qd22	Qd32	.	Qdn2
	.		.		.		.	.
	Qd1np	Qd2np	Qd3np	.	Qdnnp]
    */
	
	/*-----------Matrix and vectors for serpentine friction-----------------*/
	
	/*vf: Tangencial and Normal Componentes of Friction Force
	C: Friction coef.
	MC: massneg*C
	Fr1: MC*vf*/
	m_resize_vars(2,2, &vf, &C, &MC, &Fr1,NULL);
	mem_stat_reg_vars(0,TYPE_MAT, &vf, &C, &MC, &Fr1,NULL);
	
	//vxy: Projecting velocity from frame i to center of gravity
	v_resize_vars(2, &vxy, NULL);
	mem_stat_reg_vars(0,TYPE_VEC, &vxy, NULL);
	
	/*---------------------------------------------------------------------*/
	
	/*H: proyection vector in the axes of motion
	term Accelerations temp vector for sk2 *v term
	qddH: d2Q * H
	qdH:  dQ * H
	temp1: force equation temp vector
	IdV: I*dV
	Fr6: 6-dimensional friction force
	Frt: Projected force at frame i, 6-dimensional*/
	v_resize_vars(6, &H, &term, &qddH, &qdH, &temp1, &IdV, &Fr6, &Frt, NULL);
	mem_stat_reg_vars(0,TYPE_VEC, &H, &term, &qddH, &qdH, &temp1, &IdV, &Fr6, &Frt,NULL);

	/*speed: vector for all traj points for i joint
	accele: vector for all traj points
	Force: vector for all traj points
	vmpi_temp2: i-1 joint speed*/
	v_resize_vars(6*np, &vmpi_temp2, &speed, &accele, &F, NULL);
	mem_stat_reg_vars(0,TYPE_VEC,  &vmpi_temp2, &speed, &accele, &F, NULL);


	/*-----------Block Matrix and vectors----------------------------------*/
	
	/*mat_mia: block matrix for comm step
	R: Rotation block matrix 6*6np
	P: traslation block matrix 6*6np*/
	m_resize_vars(6,6*np, &mat_mia, &R, &P,NULL);
	mem_stat_reg_vars(0,TYPE_MAT, &R, &P, &mat_mia,NULL);
	//r: Rotation block matrix 3*3np
	m_resize_vars(3,3*np, &r,NULL);
	mem_stat_reg_vars(0,TYPE_MAT, &r, NULL);
	//p: traslation block matrix 3*np
	m_resize_vars(3,np, &p,NULL);
	mem_stat_reg_vars(0,TYPE_MAT, &p,NULL);
	
	/*---------------------------------------------------------------------*/
	
	
	/*Skw: block 6-dimensional Skew symetric matrix
	S:  6-dimensional matrix for s.
	Icm: Inertia matrix refered to mass center frame
	I: Inertia matrix refered to body's local frame
	Sicm: temp matrix for applying parallel axes theorem to obtain the inertial operator at frame i*/
	m_resize_vars(6,6, &Skw, &S, &Icm, &I, &SIcm, NULL);
	mem_stat_reg_vars(0,TYPE_MAT, &Skw, &S, &Icm, &I, &SIcm, NULL);
	
	//sk: skew symetric temp matrix 
	//sksk: sk*sk matrix
	//Ticm: Inertial operator referred to center of gravity of body i
	//rTicm: r* Ticm
	//Itemp, 3*3 I temp matrix
	//ji: Inertial Tensor referred to center of gravity of body i
	//massU: mass * ident
	m_resize_vars(3,3, &sk, &sksk, &Ticm, &rTicm, &Itemp, &identidad, &Ji, &massU, NULL);
	mem_stat_reg_vars(0,TYPE_MAT, &sk, &sksk, &Ticm, &rTicm, &Itemp, &identidad, &Ji, &massU, NULL);
	
	//WS, up, down: temp vectors
	//ss: r* s
	v_resize_vars(3, &ws, &s, &ss, &up, &down, &vtempo, &rws, &rv1, &Iv, &sotro, &v, NULL);
	mem_stat_reg_vars(0,TYPE_VEC, &ws, &s, &ss, &up, &down, &vtempo, &rws, &rv1, &Iv, &sotro, &v,NULL);

	/**************** 4x4 **************/
	/*T1-T2-T2: Temp matrix.*/
	m_resize_vars(4,4,&T1,&T2,NULL);
	mem_stat_reg_vars(0,TYPE_MAT,&T1,&T2,NULL);
	
	
	m_ident(identidad);
	
	for(k=0;k<np;k++)
	{	
		#ifdef print_time2
		if(my_id==0 || my_id==3 || my_id==7)
		    av=MPI_Wtime();
		#endif
		
		if(DH->me[my_id][4]==0)
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
		if(my_id==0)
		{
			r->me[0][3*k]=1;r->me[0][1+3*k]=0;r->me[0][2+3*k]=0;
			r->me[1][3*k]=0;r->me[1][1+3*k]=1;r->me[1][2+3*k]=0;
			r->me[2][3*k]=0;r->me[2][1+3*k]=0;r->me[2][2+3*k]=1;
			p->me[0][k]=0;
			p->me[1][k]=0;
			p->me[2][k]=0;
			
		}
		else
		{
			alfa = DH->me[my_id-1][0];
			a=DH->me[my_id-1][1];
			if(DH->me[my_id-1][4]==0)        //Rotational
			{
				teta=Q->me[k][my_id-1];
				d=DH->me[my_id-1][3];
			}
			else                         //Prismatic
			{
				teta=DH->me[my_id-1][2];
				d=Q->me[k][my_id-1];
			}
			
			//obtain Rotational transformation
			Rotational(alfa,a,teta,d,r,k);
			p->me[0][k]=-a;					//obtain the inverse vector of position p
			p->me[1][k]=-d*sin(alfa);
			p->me[2][k]=-d*cos(alfa);
		}
		
		//Obtain the 6-dimensional R operator based on r
		m_move(r,0,3*k,3,3,R,0,6*k);	
		m_move(r,0,3*k,3,3,R,3,6*k+3);
		//Obtain the 6-dimensional P operator based on p
		
		m_move(identidad,0,0,3,3,P,0,6*k);
		m_move(identidad,0,0,3,3,P,3,6*k+3);
		
		skew_symetric(sk,p,k);
		sm_mlt(-1,sk,sk);
		
		m_move(sk,0,0,3,3,P,0,6*k+3);
//***********************************************************************************
//Computing Spatial Velocities		
//***********************************************************************************
		
		
		//m_mlt(P_trans,R_trans,PtRt);
		if(my_id!=0)
			for(i=0;i<=5;i++)
				for(j=0;j<=5;j++)	
					mat_mia->me[i][j+6*k]=P->me[0][i+6*k]*R->me[j][0+6*k] 
								+P->me[1][i+6*k]*R->me[j][1+6*k]
								+P->me[2][i+6*k]*R->me[j][2+6*k]
								+P->me[3][i+6*k]*R->me[j][3+6*k]
								+P->me[4][i+6*k]*R->me[j][4+6*k]
								+P->me[5][i+6*k]*R->me[j][5+6*k];

		
		//sv_mlt(dQ->me[k][my_id],H,dQH);
		speed->ve[0+6*k]=H->ve[0]*dQ->me[k][my_id];
		speed->ve[1+6*k]=H->ve[1]*dQ->me[k][my_id];
		speed->ve[2+6*k]=H->ve[2]*dQ->me[k][my_id];
		speed->ve[3+6*k]=H->ve[3]*dQ->me[k][my_id];
		speed->ve[4+6*k]=H->ve[4]*dQ->me[k][my_id];
		speed->ve[5+6*k]=H->ve[5]*dQ->me[k][my_id];
		/*speedc[0+6*k]=H->ve[0]*dQ->me[k][my_id];
		speedc[1+6*k]=H->ve[1]*dQ->me[k][my_id];
		speedc[2+6*k]=H->ve[2]*dQ->me[k][my_id];
		speedc[3+6*k]=H->ve[3]*dQ->me[k][my_id];
		speedc[4+6*k]=H->ve[4]*dQ->me[k][my_id];
		speedc[5+6*k]=H->ve[5]*dQ->me[k][my_id];*/
	}
	
	#ifdef print_time2
	if(my_id==0 || my_id==3 || my_id==7)
	    acv=MPI_Wtime();
	#endif
//*************************************************************************************
//SPEED PROPAGATION
//*************************************************************************************
	speed=com_step(my_id,num_pro,speed, mat_mia);//,Status, &Request);

	#ifdef print_time2
	if(my_id==0 || my_id==3 || my_id==7)
	    aa=MPI_Wtime();
	#endif
		
	vmpi_temp2=com_step2(my_id,num_pro,speed,vmpi_temp2);//,Status,Request);

	#ifdef print_speed
	printf("\nSpeed\n");
	v_output(speed);
	#endif
//*********************************************************************************
// Computing Spatial Accelerations
//*********************************************************************************
	v_zero(term);
	for(k=0;k<np;k++)
	{	
		// Obtain the last and actual angular velocity from Vi-1 and Vi
		for(cont=0;cont<=2;cont++)
		{
			ws->ve[cont]=vmpi_temp2->ve[6*k+cont];
		}
	
		skew_symetric2(sk,speed,k);	//obtain Skew symetric wi
	
		//Obtain 6-dimensional Sk
		m_move(sk,0,0,3,3,Skw,0,0);
		m_move(sk,0,0,3,3,Skw,3,3);
	
		//m_transp(r,r_trans);
		//mv_mlt(r_trans,ws,rws);
		for(cont=0;cont<=2;cont++)
		{
			rws->ve[cont]=r->me[0][cont+3*k]*ws->ve[0]
					+r->me[1][cont+3*k]*ws->ve[1]
					+r->me[2][cont+3*k]*ws->ve[2];
		}
		
		skew_symetric2(sk,rws,0); 	//obtain Skew symetric wi-1
		
		// Obtain the last and actual velocity from Vi-1 and Vi
		for(cont=0;cont<=2;cont++)
		{
			v->ve[cont]=speed->ve[6*k+cont+3];
			ws->ve[cont]=vmpi_temp2->ve[6*k+cont+3];
		}
		//mv_mlt(r_trans,ws,rv1);
		for(cont=0;cont<=2;cont++)
		{
			rv1->ve[cont]=r->me[0][cont+3*k]*ws->ve[0]
					+r->me[1][cont+3*k]*ws->ve[1]
					+r->me[2][cont+3*k]*ws->ve[2];
		}
		v_sub(v,rv1,v);
		mv_mlt(sk,v,vtempo);

		for(cont=0;cont<=2;cont++)
			term->ve[cont+3]=vtempo->ve[cont];
		
		//Obtain 6-dimensional acceleration
		sv_mlt(d2Q->me[k][my_id],H,qddH);
		sv_mlt(dQ->me[k][my_id],H,qdH);
		
		if(my_id==0)
		{
			for(cont=0;cont<=5;cont++)
				accele->ve[cont+6*k]=term->ve[cont]+qddH->ve[cont]+GRAV->ve[cont]
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
				accele->ve[cont+6*k]=term->ve[cont]+qddH->ve[cont]
							+Skw->me[cont][0]*qdH->ve[0]
							+Skw->me[cont][1]*qdH->ve[1]
							+Skw->me[cont][2]*qdH->ve[2]
							+Skw->me[cont][3]*qdH->ve[3]
							+Skw->me[cont][4]*qdH->ve[4]
							+Skw->me[cont][5]*qdH->ve[5];
		}
	}	
	#ifdef print_time2
	if(my_id==0 || my_id==3 || my_id==7)
		aca=MPI_Wtime();
	#endif
	
	accele=com_step(my_id,num_pro,accele, mat_mia);
	
	#ifdef print_time2
	if(my_id==0 || my_id==3 || my_id==7)
		af=MPI_Wtime();
	#endif

	m_zero(Itemp);
	
	#ifdef print_accele
	printf("\nAccel\n");
	v_output(accele);
	#endif
//****************************************************************************************
// Computing Spatial Forces
//****************************************************************************************
	for(k=0;k<np;k++)
	{	
		alfa = DH->me[my_id][0];
		a=DH->me[my_id][1];
		
		if(DH->me[my_id][4]==0)        //Rotational
		{
			teta=Q->me[k][my_id];
			d=DH->me[my_id][3];
		}
		else                       //Prismatic
		{
			teta=DH->me[my_id][2];
			d=Q->me[k][my_id];
		}
		Rotational(alfa,a,teta,d,r,k);
		p->me[0][k]=-a;					//obtain the inverse vector of position p
		p->me[1][k]=-d*sin(alfa);
		p->me[2][k]=-d*cos(alfa);
		
		//Obtain the 6-dimensional R operator based on r
		m_move(r,0,3*k,3,3,R,0,6*k);	
		m_move(r,0,3*k,3,3,R,3,6*k+3);
		//Obtain the 6-dimensional P operator based on p
		m_move(identidad,0,0,3,3,P,0,6*k);
		m_move(identidad,0,0,3,3,P,3,6*k+3);
		
		skew_symetric(sk,p,k);
		sm_mlt(-1,sk,sk);
		
		m_move(sk,0,0,3,3,P,0,6*k+3);
		
		//Obtain vector s from frame i+1 to center of mass of body i
		s->ve[0]=DYN->me[my_id][1];
		s->ve[1]=DYN->me[my_id][2];
		s->ve[2]=DYN->me[my_id][3];
		//v_sub(s,p,s);
		for(cont=0;cont<=2;cont++)
			s->ve[cont]=s->ve[cont]-p->me[cont][k];
	
		//mv_mlt(r,s,ss);
		for(cont=0;cont<=2;cont++)
			ss->ve[cont]=r->me[cont][0+3*k]*s->ve[0]
					+r->me[cont][1+3*k]*s->ve[1]
					+r->me[cont][2+3*k]*s->ve[2];
					
		skew_symetric2(sk,ss,0);		//Obtain Skew symetric of ss
	
		//obtain 6-dimensional operator S
		m_move(identidad,0,0,3,3,S,0,0);
		m_move(identidad,0,0,3,3,S,3,3);
		m_move(sk,0,0,3,3,S,0,3);
		
		
		skew_symetric2(sk,speed,k);	//Obtain Skew symetric of v
			
		//Obtain de Inertial Tensor referred to center of gravity of body i
		Ticm->me[0][0]=DYN->me[my_id][4];
		Ticm->me[0][1]=DYN->me[my_id][7];
		Ticm->me[0][2]=DYN->me[my_id][9];
		Ticm->me[1][0]=DYN->me[my_id][7];
		Ticm->me[1][1]=DYN->me[my_id][5];
		Ticm->me[1][2]=DYN->me[my_id][8];
		Ticm->me[2][0]=DYN->me[my_id][9];
		Ticm->me[2][1]=DYN->me[my_id][8];
		Ticm->me[2][2]=DYN->me[my_id][6];
			
		//Obtain 6-dimensional Inertial operator in center of gravity:
		//m_mlt(r,Ticm,rTicm);
		for(i=0;i<=2;i++)
			for(j=0;j<=2;j++)	
				rTicm->me[i][j]=r->me[i][0+3*k]*Ticm->me[0][j]
						+r->me[i][1+3*k]*Ticm->me[1][j]
						+r->me[i][2+3*k]*Ticm->me[2][j];
						
		//m_transp(r,r_trans);
		//m_mlt(rTicm,r_trans,Ji);
		for(i=0;i<=2;i++)
			for(j=0;j<=2;j++)	
				Ji->me[i][j]=rTicm->me[i][0]*r->me[j][0+3*k]
							+rTicm->me[i][1]*r->me[j][1+3*k]
							+rTicm->me[i][2]*r->me[j][2+3*k];
	
		m_move(Ji,0,0,3,3,Icm,0,0);
		sm_mlt(DYN->me[my_id][0],identidad,massU);
		m_move(massU,0,0,3,3,Icm,3,3);

		//applying parallel axes theorem to obtain the inertial operator at frame i
		m_mlt(S,Icm,SIcm);
	
		//m_mlt(SIcm,S_trans,I);
		for(i=0;i<=5;i++)
			for(j=0;j<=5;j++)	
				I->me[i][j]=SIcm->me[i][0]*S->me[j][0]
						+SIcm->me[i][1]*S->me[j][1]
						+SIcm->me[i][2]*S->me[j][2]
						+SIcm->me[i][3]*S->me[j][3]
						+SIcm->me[i][4]*S->me[j][4]
						+SIcm->me[i][5]*S->me[j][5];
	
		
		//*******************************************************************************
		// Computing Force equation:
		//*******************************************************************************
		m_move(I,0,0,3,3,Itemp,0,0);
		
		for(cont=0;cont<=2;cont++)
			v->ve[cont]=speed->ve[cont+6*k];
		
		mv_mlt(Itemp,v,Iv);
		mv_mlt(sk,Iv,up);
		m_mlt(sk,sk,sksk);
		mv_mlt(sksk,ss,sotro);
		sv_mlt(DYN->me[my_id][0],sotro,down);
		
		for(cont=0;cont<=2;cont++)
		{
			temp1->ve[cont]=up->ve[cont];
			temp1->ve[cont+3]=down->ve[cont];
		}
		
		//Obtain 6-dimensional force
		for(cont=0;cont<=5;cont++)
		{	
			IdV->ve[cont]=I->me[cont][0]*accele->ve[0+6*k]
					+I->me[cont][1]*accele->ve[1+6*k]
					+I->me[cont][2]*accele->ve[2+6*k]
					+I->me[cont][3]*accele->ve[3+6*k]
					+I->me[cont][4]*accele->ve[4+6*k]
					+I->me[cont][5]*accele->ve[5+6*k];
		}
		v_add(IdV,temp1,temp1);
		
		//*************************************************************************************
		//Obtain Friction Force with Contact Surface (for Serpentine Robots)
		//*************************************************************************************
			
		vf->me[0][0]=cos(X->me[k][my_id]);
		vf->me[0][1]=sin(X->me[k][my_id]);
		vf->me[1][0]=-sin(X->me[k][my_id]);
		vf->me[1][1]=cos(X->me[k][my_id]);
		
		//Matrix of Coefficients of friction
		C->me[0][0]=DYN->me[my_id][10];
		C->me[0][1]=0;
		C->me[1][0]=0;
		C->me[1][1]=DYN->me[my_id][11];
		//m_transp(vf,vf_trans);
		massneg=-DYN->me[my_id][0];
		//sm_mlt(massneg,vf_trans,massvf);
		//m_mlt(massvf,C,MC);
		for(i=0;i<=1;i++)
			for(j=0;j<=1;j++)	
				MC->me[i][j]=massneg*(vf->me[0][i]*C->me[0][j]
						+vf->me[1][i]*C->me[1][j]);
						
		m_mlt(MC,vf,Fr1);
				
		//Projecting velocity from frame i to center of gravity
		for(cont=3;cont<=4;cont++)
		{	
			vxy->ve[cont-3]=S->me[0][cont]*speed->ve[0+6*k]
				+S->me[1][cont]*speed->ve[1+6*k]
				+S->me[2][cont]*speed->ve[2+6*k]
				+S->me[3][cont]*speed->ve[3+6*k]
				+S->me[4][cont]*speed->ve[4+6*k]
				+S->me[5][cont]*speed->ve[5+6*k];
		}
		
		//6-dimensional friction force
		Fr6->ve[0]=0;
		Fr6->ve[1]=0;
		Fr6->ve[2]=-0.3333*DYN->me[my_id][11]*DYN->me[my_id][0]*0.05*0.05*dQ->me[k][my_id];
		Fr6->ve[3]=Fr1->me[0][0]*vxy->ve[0]
				+Fr1->me[0][1]*vxy->ve[1];
		Fr6->ve[4]=Fr1->me[1][0]*vxy->ve[0]
				+Fr1->me[1][1]*vxy->ve[1];
		Fr6->ve[5]=0;
		
		//Projecting force at frame i, 6-dimensional
		mv_mlt(S,Fr6,Frt);
				
		//add Friction force component to local spatial Force
		for(cont=0;cont<=5;cont++)
			F->ve[cont+6*k]=Frt->ve[cont]+temp1->ve[cont];
		
		for(i=0;i<=5;i++)
			for(j=0;j<=5;j++)	
				mat_mia->me[i][j+6*k]=R->me[i][0+6*k]*P->me[0][j+6*k] 
						+R->me[i][1+6*k]*P->me[1][j+6*k]
						+R->me[i][2+6*k]*P->me[2][j+6*k]
						+R->me[i][3+6*k]*P->me[3][j+6*k]
						+R->me[i][4+6*k]*P->me[4][j+6*k]
						+R->me[i][5+6*k]*P->me[5][j+6*k];
		
		if(my_id==num_pro-1)
		{
			for(cont=0;cont<=5;cont++)
				F->ve[cont+6*k]=mat_mia->me[cont][0+6*k]*fext->ve[0]
						+mat_mia->me[cont][1+6*k]*fext->ve[1]
						+mat_mia->me[cont][2+6*k]*fext->ve[2]
						+mat_mia->me[cont][3+6*k]*fext->ve[3]
						+mat_mia->me[cont][4+6*k]*fext->ve[4]
						+mat_mia->me[cont][5+6*k]*fext->ve[5]
						+F->ve[cont+6*k];
		}
	}
	
	#ifdef print_time2		
	if(my_id==0 || my_id==3 || my_id==7)
	    acf=MPI_Wtime();
	#endif
	
	F=com_step_inv(my_id,num_pro,F, mat_mia);//,Status, &Request);
		
	//Projecting force in the axe of motion
	for(k=0;k<np;k++)
	{
		if(DH->me[my_id][4]==0)                 //Rotational
			Torques->ve[k]=F->ve[2+6*k];
		else                               //Prismatic
			Torques->ve[k]=F->ve[5+6*k];
	}
	
	#ifdef print_time2
	if(my_id==0 || my_id==3 || my_id==7)
	{
		fin=MPI_Wtime();
		printf("avel: %f\nacvel: %f\naaccel:%f\nacaccel:%f\n",av,acv,aa,aca);
		printf("af: %f\nacf: %f\nfin: %f\n",af,acf,fin);
	}
	#endif
	
	return Torques;
}

int main( int  argc, char **argv)
{
	mem_stat_mark(1);
	int namelen, i;
	int my_id, num_pro,rec;
	double startwtime, middlewtime, finalwtime;
	
	char processor_name[MPI_MAX_PROCESSOR_NAME];
	char msg[MPI_MAX_PROCESSOR_NAME+80];
	static MAT *Q_Qd_Qdd=MNULL, *DH=MNULL, *DYN=MNULL, *DYN_base=MNULL, *xb=MNULL, *vb=MNULL, *ab=MNULL, *Ttotal=MNULL, *Xcm=MNULL;
	static VEC *gravity=VNULL, *fext=VNULL, *time=VNULL, *Torques=VNULL;
	FILE *fp;
	
	
	gravity=v_resize(gravity,6);
	MEM_STAT_REG(gravity,TYPE_VEC);
	v_zero(gravity);
	gravity->ve[5]=9.81;
	
	fext=v_resize(fext,6);
	MEM_STAT_REG(fext,TYPE_VEC);
	v_zero(fext);
	
	Q_Qd_Qdd=m_resize(Q_Qd_Qdd,1,1);
	MEM_STAT_REG(Q_Qd_Qdd,TYPE_MAT);
	
	MPI_Init( &argc, &argv );
	MPI_Comm_size( MPI_COMM_WORLD, &num_pro);
	MPI_Comm_rank( MPI_COMM_WORLD, &my_id );
	MPI_Get_processor_name(processor_name, &namelen);
	
	//if(my_id==0)
	    startwtime=MPI_Wtime();
	//df=num_pro;
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
		//Ttotal=m_finput(fp,MNULL);
		//xb=m_finput(fp,MNULL);
		//MEM_STAT_REG(xb,TYPE_MAT);
		//vb=m_finput(fp,MNULL);
		//MEM_STAT_REG(vb,TYPE_MAT);
		//ab=m_finput(fp,MNULL);
		//MEM_STAT_REG(ab,TYPE_MAT);
		
		fclose(fp);	
	}
	
	#ifdef print_id
	printf("\n---------------Inverse dynamics process %d of %d on %s -----------\n", my_id, num_pro,processor_name);
	#endif
	MPI_Barrier(MPI_COMM_WORLD);
	#ifdef print_time
	if(my_id==0 || my_id==3 || my_id==7)
		middlewtime=MPI_Wtime();
	#endif

	mem_stat_mark(2);
	Torques=ff_ne(DH,DYN,Q_Qd_Qdd,gravity,fext,Torques,Xcm,vb,ab,my_id,num_pro);
	mem_stat_free(2);
	
	#ifdef print_time
	if(my_id==0 || my_id==3 || my_id==7)
	{
		finalwtime=MPI_Wtime();
		printf("\n\nMy_id=%d",my_id);
		printf("\nstart time: %f\nMiddle time: %f\nFinal time: %f\n",startwtime,middlewtime,finalwtime);
		printf("Final-start: %f\nFinal-Middle: %f\n",finalwtime-startwtime,finalwtime-middlewtime);
	}
	#endif
	
	#ifdef print_torques
	printf("\n----------------------------------Torques:----------------------------------\n\n");
	v_output(Torques);
	#endif
	MPI_Finalize();
	mem_stat_free(1);
	return (0);
}

