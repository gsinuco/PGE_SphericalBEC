#include <stdio.h>
#include <math.h>
#include <gsl/gsl_sf_coupling.h>

extern struct{
  int  l1,l2,l3,l4,m1,m2,m3,m4;
}jSymbols_;

void j3symbols_fun_(int *l1,int *l2,int *l3, int *l4,  int *m1, int *m2, int *m3, int *m4, double *beta,int *l_min, int *l_max){
  
  double sum;
  int l;
  
  sum = 0.0;

  //  printf("%d %d \n",*l_min,*l_max);
  
  if(*m4 == *m1-*m2+*m3){

    
    for(l=*l_min;l<=*l_max;l+=2){
      /*      
      sum += sqrt((2*(*l1)+1)*(2*(*l2)+1)*(2*(*l3)+1)*(2*(*l4)+1))*(2*l+1)*gsl_sf_coupling_3j((2*(*l1)),(2*(*l2)),(2*l),2*(*m1),-2*(*m2),-2*(*m1-*m2))*gsl_sf_coupling_3j(2*(*l3),(2*(*l4)),(2*l),2*(*m3),-2*(*m4),-2*(*m3 - *m4))*gsl_sf_coupling_3j((2*(*l1)),(2*(*l2)),(2*l),0,0,0)*gsl_sf_coupling_3j(2*(*l3),(2*(*l4)),(2*l),0,0,0); // original version
      */
      //modified begins: 16.09.2011
      
      sum += sqrt((2*(*l1)+1)*(2*(*l2)+1)*(2*(*l3)+1)*(2*(*l4)+1))*(2*l+1)*gsl_sf_coupling_3j((2*(*l1)),(2*(*l2)),(2*l),2*(*m1),-2*(*m2),-2*(*m1-*m2))*gsl_sf_coupling_3j(2*l,2*(*l3),(2*(*l4)),2*(*m1-*m2),2*(*m3),-2*(*m4))*gsl_sf_coupling_3j((2*(*l1)),(2*(*l2)),(2*l),0,0,0)*gsl_sf_coupling_3j(2*l,2*(*l3),(2*(*l4)),0,0,0); 
      
      
      // modified ends: 16.09.2011
      //printf("%f %i %i %i %f\n",sum,*l1,*l2,l,gsl_sf_coupling_3j((2*(*l1)),(2*(*l2)),(2*l),0,0,0));
    }
    //sum = pow(-1,*m1+*m4)*sum/(4*M_PI); //original version
    sum = pow(-1,*m1)*sum/(4*M_PI); //modified version 16.09.2011
  }else{
    
    printf("WARNING: m3-m4 != m2-m1");
    sum  = 0;
    
  }
  
  *beta = sum;  
}

