#
#	Inverse dynamics parallel Solution, log2(n) and precalc MAKEFILE
#	Thesis: Modelling and inverse dynamics solution for serpentine robots using strictly parallel algorithms
#	Proyect: Time lower bound modelling and simulation of complex rigid body systems.
#		Andres Jaramillo Botero ajaramil@puj.edu.co
#		Juan Camilo Acosta jcacosta@puj.edu.co
#		Julian David Colorado jdcolorado@puj.edu.co
#		Juan Manuel Florez jmflorez@puj.edu.co
#	Copyright Robotics and Automation Group, 2006
#	Pontificia Universidad Javeriana, Cali, Colombia
#
#

# The C compiler :
CC=mpicc
CC2=gcc


# Compiler options :
COPTS=-g -O7 -Wall


#LIBS
Meschach.a:
	cd ../Meschach; \
	make all; \
	
llrl.a: Meschach.a
	cd ../LLRL; \
	make -f LLRL.mak; \

#serial
serial: test-serp_flying.o 
	$(CC2) $(COPTS) -o serial test-serp_flying.o ../LLRL/LLRL.a ../Meschach/meschach.a  -lm -pg

test-serp_flying.o: llrl.a 
	$(CC2) $(COPTS) -c test-serp_flying.c -pg


#precalc comp
precalc: precalc.o 
	$(CC2) $(COPTS) -o precalc precalc.o ../LLRL/LLRL.a ../Meschach/meschach.a  -lm

precalc.o: llrl.a
	$(CC2) $(COPTS) -c precalc.c

#ff_ne2
#Precompilation routines

Solo: Meschach.a
	$(CC) -c ff_ne2.c 

print_accele: Meschach.a
	$(CC) -c ff_ne2.c -Dprint_id -Dprint_accele

print_speed: Meschach.a
	$(CC) -c ff_ne2.c -Dprint_id -Dprint_speed
		
print_torques: Meschach.a
	$(CC) -c ff_ne2.c -Dprint_id -Dprint_torques

print_id: Meschach.a
	$(CC) -c ff_ne2.c -Dprint_id

print_time: Meschach.a
	$(CC) -c ff_ne2.c -Dprint_time

print_time2: Meschach.a
	$(CC) -c ff_ne2.c -Dprint_time2

log:
	$(CC) -mpilog -o ff_ne2 ff_ne2.o  ../Meschach/meschach.a  -lm
nolog: 
	$(CC)  -o ff_ne2 ff_ne2.o  ../Meschach/meschach.a  -lm
	
clean: 
	rm -f *.o
	rm -f ../LLRL/*.o
	rm -f ../Meschach/*.o


