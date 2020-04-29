gcc -O3 -c -o 3jsymbols.o 3jsymbols_fun.c
gfortran -c beta_lmV3.f90 -o beta_lmV3.o
gfortran -lgsl -lgslcblas -o beta_lmV3 3jsymbols.o beta_lmV3.o 
