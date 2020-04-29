// Illustration of compiler <complex> class.

#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
#include <math.h>
//using namespace std;

//typedef complex<double> dcmplx;

`

extern struct{
  double timeStep;
  int  Cmodes;
}dimensions_;


void f(double t, double complex *y,double complex *k){
  int i;


  //forced harmonic oscillator
  /*
  k[0] =  y[1];
  k[1] = -2.0*y[0]+6.0*cos(2*M_PI*t);
  */
  //

  // GPE in mode space:
  double GPE_NonLinTerm;
  
  for(){
    for(){
      for(){
	GPE_NonLinTerm += y[n]*conj(y[n_])*y[m]*beta(n,n_,m,i);
      }
    }
  }
  
  k[i] = u*GPE_NonLinTerm;

}

int rk4th_(double complex *k1,double complex *k2,double complex *k3,double complex *k4,double complex *y,double *t){
  

  double complex *y_aux;
  double dt;
  int i;
  int dim;
 
  dim = dimensions_.Cmodes;
  dt  = dimensions_.timeStep;

  y_aux = (double complex*) calloc(dim,sizeof(double complex));
 

  f(*t,y,k1);
  
  for(i=0;i<dim;i++){
    y_aux[i] = y[i] +k1[i]*dt*0.5;
  }
  f(*t+0.5*dt,y_aux,k2);

  for(i=0;i<dim;i++){
    y_aux[i] = y[i] +k2[i]*dt*0.5;
  }
  f(*t+0.5*dt,y_aux,k3);
  
  for(i=0;i<dim;i++){
    y_aux[i] = y[i] +k3[i]*dt;
  }  
  f(*t+dt,y_aux,k4);
  
  for(i=0;i<dim;i++){
    y[i] = y[i] + dt*(k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i])/6.0;
  }
  *t = *t + dt;
  
  return 0;

}


