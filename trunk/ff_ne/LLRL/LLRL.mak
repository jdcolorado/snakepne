#
# LLRL.mak
# 
#     This file implements the Low Level Robotics Library Makefile.
#
# Date: 2003-08-06
#
# Last modification: 2006-02-15
#
# Copyright  Robotics and Automation Group, Pontificia Universidad Javeriana - Cali.
#		             Andr� Jaramillo Botero, ajaramil@puj.edu.co
#		             Wilber Perea Castro, wpcmmx@hotmail.com
#                   Antonio Alejandro Matta G�ez, amatta@puj.edu.co
#                   Juli� David Colorado, jdcolorado@puj.edu.co
#                   Juan Camilo Acosta, jcacosta@puj.edu.co
#
# See the file "license.terms" for information on usage and redistribution of this file, and for a
# DISCLAIMER OF ALL WARRANTIES.
#
# SCCS: %Z%  %M%  %I%  %E%  %U%-------------------------------------------------------------------------------------------------------

OBJS=abi.o accel.o AcelBase.o atorque.o basetray.o base_x.o base_v.o cabi.o eul2q.o euler.o fdyn.o feval.o ff_ne.o fkine.o friction.o Friction_surface.o gravload.o \
jacob0.o jacob_base.o jacob_end.o jacobn.o jacobn_float.o Kinetic.o linkiner.o linktran.o m_fwrite.o  ne.o ntrp45.o \
ode45.o pinv.o Potential.o robots.o rotvec.o torqfun.o tr2q.o ttr2diff.o v_cross.o v_fwrite.o verlet.o vtorques.o xttr.o\
absm.o div_m.o indice.o m_fwrite.o maxm.o maxm2.o sum_mat_vec.o ikine.o serpenoid.o tr2rpy.o
 
#
# --------------------------------------------------------------------------------------------------------------------------
# The C compiler :
# --------------------------------------------------------------------------------------------------------------------------
CC=gcc
# --------------------------------------------------------------------------------------------------------------------------
# Compiler options :
# --------------------------------------------------------------------------------------------------------------------------

COPTS=-g -O2 -pg
.c.o:
	$(CC) -c $(COPTS) $<

# ------------------------------------------------ --------------------------------------------------------------------------------------------------------------------------

OBJS=abi.o accel.o --------------------------------------------------------------------------
# Libraries to link with :
# --------------------------------------------------------------------------------------------------------------------------

LIBS= -lm
RANLIB = ranlib

# --------------------------------------------------------------------------------------------------------------------------
# Objects :
# --------------------------------------------------------------------------------------------------------------------------

OBJS=abi.o accel.o AcelBase.o atorque.o basetray.o base_x.o base_v.o cabi.o eul2q.o euler.o fdyn.o feval.o ff_ne.o fkine.o friction.o Friction_surface.o gravload.o \
jacob0.o jacob_base.o jacob_end.o jacobn.o jacobn_float.o Kinetic.o linkiner.o linktran.o m_fwrite.o  ne.o ntrp45.o \
ode45.o pinv.o Potential.o robots.o rotvec.o torqfun.o tr2q.o ttr2diff.o v_cross.o v_fwrite.o verlet.o vtorques.o xttr.o\
absm.o div_m.o indice.o m_fwrite.o maxm.o maxm2.o sum_mat_vec.o ikine.o serpenoid.o tr2rpy.o 
 
all: objs
 

$(OBJS): LLRL.h
objs:$(OBJS)
	ar ru LLRL.a $(OBJS)
	$(RANLIB) LLRL.a
    
clean:
	rm -f *.o
	rm -f LLRL.a
	

# --------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------
