PROGRAM INDEXREVERT
  
  IMPLICIT NONE
  INTEGER i,j,k,l,i_,D;
  
  D = 3

  DO i_=1,81
     
     i = (i_-1)/D**3 + 1
     j = int((i_-1)/D**2 -D*(i-1)) + 1
     k = (i_ -1 - (i-1)*D**3 - (j-1)*D**2)/(D) + 1
     l = i_ -1 - (i-1)*D**3 - (j-1)*D**2 - (k-1)*D +1
     
     WRITE(*,*) i_,l,k,j,i,l + D*(k-1) + (j-1)*D**2 + (i-1)*D**3
  END DO 

     

END PROGRAM INDEXREVERT
