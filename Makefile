CC = gcc
GF = gfortran
GFFLAGS = -lgsl -lgslcblas
CFLAGS = -O3 

all: PGPEonSphere

#PGPEonSphere: funciones.mod  rungekuttamodule.mod PGPEonSphere.f90 zufall_sub.o Runge-Kutta4thorder_sub.o PGPEInteractionTerm_sub.o
#	$(GF) $(GFFLAGS) zufall_sub.o Runge-Kutta4thorder_sub.o  PGPEInteractionTerm_sub.o -o $@ PGPEonSphere.f90

PGPEonSphere: funciones.mod  rungekuttamodule.mod PGPEonSphere.f90 zufall_sub.o Runge-Kutta5thorder_sub.o PGPEInteractionTerm_sub.o 
	$(GF) $(GFFLAGS) zufall_sub.o Runge-Kutta5thorder_sub.o  PGPEInteractionTerm_sub.o -o $@ PGPEonSphere.f90


#beta_lm_dir/beta_lm.dat: beta_lm_dir/*.f90 
#	cd /home/german/Programas/SphericalGas/beta_lm_dir
#	./make

#beta_lm_dir/beta_lm_map.dat: beta_lm_dir/*.f90 
#	cd /home/german/Programas/SphericalGas/beta_lm_dir
#	make


zufall_sub.o: zufall_sub.f
	gfortran zufall_sub.f -c

#Runge-Kutta4thorder_sub.o: Runge-Kutta4thorder_sub.f90
#	gfortran Runge-Kutta4thorder_sub.f90 -c

#rungekuttamodule.mod: Runge-Kutta4thorder_sub.f90
#	gfortran Runge-Kutta4thorder_sub.f90 -c


Runge-Kutta5thorder_sub.o: Runge-Kutta5thorder_sub.f90
	gfortran Runge-Kutta5thorder_sub.f90 -c

rungekuttamodule.mod: Runge-Kutta5thorder_sub.f90
	gfortran Runge-Kutta5thorder_sub.f90 -c

PGPEInteractionTerm_sub.o: PGPEInteractionTerm_sub.f90 PGPE_modules.f90
	gfortran PGPEInteractionTerm_sub.f90 -c

funciones.mod:  PGPEInteractionTerm_sub.f90 PGPE_modules.f90
	gfortran PGPEInteractionTerm_sub.f90 -c

