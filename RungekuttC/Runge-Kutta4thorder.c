// Illustration of compiler <complex> class.
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <complex>
using namespace std;

typedef complex<double> dcmplx;


#define iterations 1000
#define dim 2

double dt = 1.0E-2;

void f(double t, dcmplx *y,dcmplx *k){
  int i;

  k[0] =  y[1];
  k[1] = -2.0*y[0]+6.0*cos(2*M_PI*t);

}

int RK4th(dcmplx *k1,dcmplx *k2,dcmplx *k3,dcmplx *k4,dcmplx *y,double& t){
  

  dcmplx *y_aux;
  int i;

  y_aux = (dcmplx*) calloc(dim,sizeof(dcmplx));
 

  f(t,y,k1);
  
  for(i=0;i<dim;i++){
    y_aux[i] = y[i] +k1[i]*dt*0.5;
  }
  f(t+0.5*dt,y_aux,k2);

  for(i=0;i<dim;i++){
    y_aux[i] = y[i] +k2[i]*dt*0.5;
  }
  f(t+0.5*dt,y_aux,k3);
  
  for(i=0;i<dim;i++){
    y_aux[i] = y[i] +k3[i]*dt;
  }  
  f(t+dt,y_aux,k4);
  
  for(i=0;i<dim;i++){
    y[i] = y[i] + dt*(k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i])/6.0;
  }
  t = t + dt;
  
  return 0;

}


int main(){
  
  dcmplx *k1,*k2,*k3,*k4,*y;
  int i;
  double t = 0.0;

  k1 = (dcmplx*) calloc(dim,sizeof(dcmplx));
  k2 = (dcmplx*) calloc(dim,sizeof(dcmplx));
  k3 = (dcmplx*) calloc(dim,sizeof(dcmplx));
  k4 = (dcmplx*) calloc(dim,sizeof(dcmplx));
  y  = (dcmplx*) calloc(dim,sizeof(dcmplx));


  y[0] = dcmplx(1.0,0.0);
  y[1] = dcmplx(0.0,0.0);

  printf("%15.6E %15.6E %15.6E \n",t,real(y[0]),real(y[1]));
  for(i=0;i<iterations;i++){
    RK4th(k1,k2,k3,k4,y,t);
    printf("%15.6E %15.6E %15.6E \n",t,real(y[0]),real(y[1]));
  }
  
  return 0;
   
}

