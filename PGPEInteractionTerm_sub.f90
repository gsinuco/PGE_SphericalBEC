! THIS PROGRAM EVALUATES THE INTERCTION TERM:
!
!  \sum_j \sum_k \sum_l c_j c_k c_l beta_ijkl
!
! OF THE GROSS-PITAEVSKII EQUATION, IN THE SPHERICAL HARMONIC BASIS.
! 
!INCLUDE "beta_lm_module.f90"


include "PGPE_modules.f90"

MODULE motionConstatn_interface
  IMPLICIT NONE
  INTERFACE
     SUBROUTINE CONSTANTS_OF_MOTION(C_modes,PSI,beta_ijkl,beta_ijkl_map, &
          &ENERGY,ANGULAR_MOMENTUM)
       
       IMPLICIT NONE  
       INTEGER,                        INTENT(IN) :: C_modes
       COMPLEX*16,       DIMENSION(:), INTENT(IN) :: PSI
       DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: beta_ijkl
       DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: beta_ijkl_map
       DOUBLE PRECISION,  INTENT(OUT):: ENERGY
       DOUBLE PRECISION, INTENT(OUT) :: ANGULAR_MOMENTUM
       
     END SUBROUTINE CONSTANTS_OF_MOTION
  END INTERFACE
END MODULE motionConstatn_interface

!!$MODULE interaction_interface
!!$  IMPLICIT NONE
!!$  INTERFACE
!!$     SUBROUTINE INTERACTIONTERM(C_modes, PSI,GPE_interaction,beta_ijkl,beta_ijkl_map)
!!$       USE funciones
!!$       IMPLICIT NONE
!!$       INTEGER,                       INTENT(IN) :: C_modes
!!$       DOUBLE PRECISION, DIMENSION(:),INTENT(IN) :: beta_ijkl
!!$       INTEGER,          DIMENSION(:),INTENT(IN) :: beta_ijkl_map
!!$       COMPLEX*16,       DIMENSION(:),INTENT(IN) :: PSI
!!$       COMPLEX*16,       DIMENSION(:),INTENT(OUT):: GPE_interaction
!!$     END SUBROUTINE INTERACTIONTERM
!!$  END INTERFACE
!!$END MODULE interaction_interface

!!$PROGRAM INTERACTION
!!$  
!!$  USE interaction_interface
!!$  IMPLICIT NONE
!!$  INTEGER C_modes,D,I
!!$  INTEGER betas,betas_calc
!!$  DOUBLE PRECISION, DIMENSION(:),  ALLOCATABLE:: beta_ijkl
!!$  INTEGER, DIMENSION(:), ALLOCATABLE :: beta_ijkl_map
!!$  DOUBLE PRECISION, DIMENSION(:),ALLOCATABLE :: PSI
!!$  DOUBLE PRECISION, DIMENSION(:),ALLOCATABLE :: GPE_interaction
!!$  
!!$  OPEN(UNIT=1,FILE='beta_lm.dat',ACTION='READ')
!!$  OPEN(UNIT=2,FILE='beta_lm_map.dat',ACTION='READ')
!!$  
!!$  C_modes = 6
!!$  D = C_modes*(C_modes-1) + C_modes
!!$  
!!$  ALLOCATE(GPE_interaction(D))
!!$  ALLOCATE(PSI(D))
!!$  
!!$  
!!$  READ(2,*) betas
!!$  READ(2,*) betas_calc
!!$  
!!$  ALLOCATE(beta_ijkl(betas_calc))
!!$  ALLOCATE(beta_ijkl_map(betas))
!!$  
!!$  DO I=1,betas_calc
!!$     READ(1,*) beta_ijkl(I)
!!$  END DO
!!$  DO I=2,betas
!!$     READ(2,*) beta_ijkl_map(I-1)
!!$  END DO
!!$  DO I=1,D
!!$     CALL RANDOM_NUMBER(PSI(i))
!!$  END DO
!!$  
!!$  CALL INTERACTIONTERM(C_modes, PSI, GPE_interaction, beta_ijkl,beta_ijkl_map)
!!$  
!!$  
!!$END PROGRAM INTERACTION

SUBROUTINE INTERACTIONTERM(C_modes, PSI,GPE_interaction,beta_ijkl,beta_ijkl_map)

  USE funciones  
  USE interactionStrength

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: C_modes
  DOUBLE PRECISION, DIMENSION(:),  INTENT(IN):: beta_ijkl
  DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: beta_ijkl_map
  COMPLEX*16, DIMENSION(:),INTENT(IN):: PSI
  COMPLEX*16, DIMENSION(:),INTENT(OUT):: GPE_interaction


  INTEGER i,j,k,l,indx
  INTEGER i_,j_,k_,l_
  INTEGER l1,l2,l3,l4
  INTEGER m1,m2,m3,m4
  INTEGER D,l_max,eltos
  INTEGER EO1,EO2,l1_min,m_,p
 
  D       = C_modes
  l_max   = INT(SQRT(1.0*C_modes))-1

!  write(*,*) l_max, C_modes

!!$  GPE_interaction(1) = -PSI(1)
!!$  GPE_interaction(2) = -PSI(2) 
!!$  GPE_interaction(3) = -PSI(3)
!!$  GPE_interaction(4) = -PSI(4) 
!!$  GPE_interaction(5) = -PSI(5)
!!$  GPE_interaction(6) = -PSI(6) 
!!$  GPE_interaction(7) = -PSI(7)
!!$  GPE_interaction(8) = -PSI(8) 
!!$  GPE_interaction(9) = -PSI(9)
  
  GPE_interaction = DCMPLX(0.0,0.0)
  eltos = 0


  DO indx=1,size(beta_ijkl_map),1

     i = int((beta_ijkl_map(indx)-1)/D**3) + 1
     j = int((beta_ijkl_map(indx)-1)/D**2 -D*(i-1)) + 1
     k = (((beta_ijkl_map(indx)-1.0) - (1.0*i-1)*D**3 - (1.0*j-1)*D**2)/(1.0*D)) + 1
     l = ((beta_ijkl_map(indx)-1.0) - (1.0*i-1.0)*D**3 - (1.0*j-1.0)*D**2 - (1.0*k-1.0)*D) + 1
          
!     WRITE(*,*) beta_ijkl_map(indx),l,k,j,i,l + D*(k-1) + (j-1)*D**2 + (i-1)*D**3
     
     l4= FLOOR(SQRT(1.0*i - 1.0))
     m4= i - 1 - l4*(l4+1)
     
     l3= FLOOR(SQRT(1.0*j - 1.0))
     m3= j - 1 - l3*(l3+1)
     
     l2= FLOOR(SQRT(1.0*k - 1.0))
     m2= k - 1 - l2*(l2+1)
     
     l1= FLOOR(SQRT(1.0*l - 1.0))
     m1= l - 1 - l1*(l1+1)   
     
     i_ = l4*(l4+1) - m4 + 1
     j_ = l3*(l3+1) - m3 + 1
     k_ = l2*(l2+1) - m2 + 1
     l_ = l1*(l1+1) - m1 + 1    
     
     
!!$     IF(i.GT.SIZE(GPE_interaction)) WRITE(*,*) "i", indx, beta_ijkl_map(indx),l,k,j,i,size(GPE_interaction),size(PSI)
!!$     IF(j.GT.SIZE(GPE_interaction)) WRITE(*,*) "j", indx, beta_ijkl_map(indx),l,k,j,i,size(GPE_interaction),size(PSI)
!!$     IF(k.GT.SIZE(GPE_interaction)) WRITE(*,*) "k", indx, beta_ijkl_map(indx),l,k,j,i,size(GPE_interaction),size(PSI)
!!$     IF(l.GT.SIZE(GPE_interaction)) WRITE(*,*) "l", indx, beta_ijkl_map(indx),l,k,j,i,size(GPE_interaction),size(PSI)
!!$     
!!$     IF(i_.GT.SIZE(GPE_interaction)) WRITE(*,*) "i_", indx, beta_ijkl_map(indx),l4,m4,size(GPE_interaction),size(PSI)
!!$     IF(j_.GT.SIZE(GPE_interaction)) WRITE(*,*) "j_", indx, beta_ijkl_map(indx),l3,m3,size(GPE_interaction),size(PSI)
!!$     IF(k_.GT.SIZE(GPE_interaction)) WRITE(*,*) "k_", indx, beta_ijkl_map(indx),l2,m2,size(GPE_interaction),size(PSI)
!!$     IF(l_.GT.SIZE(GPE_interaction)) WRITE(*,*) "l_", indx, beta_ijkl_map(indx),l1,m1,size(GPE_interaction),size(PSI)

     IF( (i.NE.j) .AND. (k.NE.l)) THEN ! .AND. (j.NE.l)) THEN

        !WRITE(*,*) "uno",indx,beta_ijkl(indx)
        GPE_interaction(i) = GPE_interaction(i)            + & 
             & PSI(j)*CONJG(PSI(k))*PSI(l)*beta_ijkl(indx) + &
             & PSI(j)*CONJG(PSI(l))*PSI(k)*beta_ijkl(indx)
        GPE_interaction(j) = GPE_interaction(j)            + &
             & PSI(i)*CONJG(PSI(k))*PSI(l)*beta_ijkl(indx) + &
             & PSI(i)*CONJG(PSI(l))*PSI(k)*beta_ijkl(indx)
        
        IF(i_.NE.i .AND. i_.NE.j) GPE_interaction(i_) = GPE_interaction(i_) + &
             & PSI(j_)*CONJG(PSI(k_))*PSI(l_)*beta_ijkl(indx)               + &
             & PSI(j_)*CONJG(PSI(l_))*PSI(k_)*beta_ijkl(indx)
        IF(j_.NE.j .AND. j_.NE.i) GPE_interaction(j_) = GPE_interaction(j_) + &
             & PSI(i_)*CONJG(PSI(k_))*PSI(l_)*beta_ijkl(indx)               + &
             & PSI(i_)*CONJG(PSI(l_))*PSI(k_)*beta_ijkl(indx)
        
!!$     ELSE IF((i.NE.j) .AND. (k.NE.l) .AND. (i.EQ.k) .AND. (j.EQ.l)) THEN
!!$        
!!$        GPE_interaction(i) = GPE_interaction(i)             + DCMPLX(0.0,-1.0) * &
!!$             & PSI(j)*CONJG(PSI(k))*PSI(l)*beta_ijkl(indx) + DCMPLX(0.0,-1.0) * &
!!$             & PSI(j)*CONJG(PSI(l))*PSI(k)*beta_ijkl(indx)
!!$
!!$        GPE_interaction(j) = GPE_interaction(j)             + DCMPLX(0.0,-1.0) * &
!!$             & PSI(i)*CONJG(PSI(k))*PSI(l)*beta_ijkl(indx) + DCMPLX(0.0,-1.0) * &
!!$             & PSI(i)*CONJG(PSI(l))*PSI(k)*beta_ijkl(indx)
!!$        
!!$        IF(i_.NE.i .AND. i_.NE.j) GPE_interaction(i_) = GPE_interaction(i_) + DCMPLX(0.0,-1.0) * &
!!$             & PSI(j_)*CONJG(PSI(k_))*PSI(l_)*beta_ijkl(indx)              + DCMPLX(0.0,-1.0) * &
!!$             & PSI(j_)*CONJG(PSI(l_))*PSI(k_)*beta_ijkl(indx)
!!$
!!$        IF(j_.NE.j .AND. j_.NE.i) GPE_interaction(j_) = GPE_interaction(j_) + DCMPLX(0.0,-1.0) * &
!!$             & PSI(i_)*CONJG(PSI(k_))*PSI(l_)*beta_ijkl(indx)              + DCMPLX(0.0,-1.0) * &
!!$             & PSI(i_)*CONJG(PSI(l_))*PSI(k_)*beta_ijkl(indx)
!!$
!!$     ELSE IF((i.NE.j) .AND. (k.NE.l) .AND. (i.EQ.l) .AND. (j.EQ.k)) THEN
!!$        
        
     ELSE IF( (i.NE.j) .AND. (k.EQ.l)) THEN
       !WRITE(*,*) "dos",indx,beta_ijkl(indx)
 
        GPE_interaction(i) = GPE_interaction(i) + PSI(j)* &
             & CONJG(PSI(k))*PSI(l)*beta_ijkl(indx) 
        GPE_interaction(j) = GPE_interaction(j) + PSI(i)* &
             & CONJG(PSI(k))*PSI(l)*beta_ijkl(indx)
        
        IF(i_.NE.i .AND. i_.NE.j) GPE_interaction(i_) = GPE_interaction(i_) + PSI(j_)* &
             & CONJG(PSI(k_))*PSI(l_)*beta_ijkl(indx) 
        IF(j_.NE.j .AND. j_.NE.i) GPE_interaction(j_) = GPE_interaction(j_) + PSI(i_)* &
             & CONJG(PSI(k_))*PSI(l_)*beta_ijkl(indx) 
     
     
     ELSE IF( (i.EQ.j) .AND. (k.NE.l)) THEN
        !WRITE(*,*) "tres",indx,beta_ijkl(indx)
        
        GPE_interaction(i) = GPE_interaction(i) +  &
             &   PSI(j)*CONJG(PSI(k))*PSI(l)*beta_ijkl(indx) + &
             &   PSI(j)*CONJG(PSI(l))*PSI(k)*beta_ijkl(indx)

        IF(i_.NE.i) GPE_interaction(i_) = GPE_interaction(i_)  + &
             &  PSI(j_)*CONJG(PSI(k_))*PSI(l_)*beta_ijkl(indx) + &
             &  PSI(j_)*CONJG(PSI(l_))*PSI(k_)*beta_ijkl(indx)       

     ELSE IF( (i.EQ.j) .AND. (k.EQ.l)) THEN
!       WRITE(*,*) "cuatro", indx,beta_ijkl(indx)
 
        GPE_interaction(i)  = GPE_interaction(i)  + PSI(j)*CONJG(PSI(k))*PSI(l)*beta_ijkl(indx) 
        
        IF(i_.NE.i) GPE_interaction(i_) = GPE_interaction(i_) +  &
             & PSI(j_)*CONJG(PSI(k_))*PSI(l_)*beta_ijkl(indx) 
     ELSE
        write(*,*) indx,i,j,k,l,beta_ijkl(indx), "no case"
     END IF

     
  END DO
  
  
!!$  
!!$  DO i = 1,D
!!$     
!!$     l4 =  FLOOR(SQRT(1.0*i - 1))
!!$     m4 =  i-1-l4*(l4+1)
!!$     !write(*,*) i,l4,m4
!!$     IF(m4.GE.0) THEN
!!$        
!!$        DO j =  i,D
!!$           
!!$           l3 =  FLOOR(SQRT(1.0*j - 1))
!!$           m3 =  j-1-l3*(l3+1)
!!$           EO1 =  MOD(l3+l4,2)
!!$           !        write(*,*) j,l3,m3
!!$           DO k =  1,D
!!$              
!!$              l2 =  FLOOR(SQRT(1.0*k - 1))
!!$              m2 =  k-1-l2*(l2+1)
!!$              
!!$              DO l = k,D,2
!!$                 
!!$                 l1 =  FLOOR(SQRT(1.0*l - 1))
!!$                 m1 =  l-1-l1*(l1+1)
!!$                 !WRITE(*,*) l4,l3,l2,l1,m4,m3,m2,m1              
!!$                 IF(EO1 .EQ.MOD(l2+l1,2)) THEN ! (L3+L4) AND (L1+L2) should have the same parity?
!!$                    
!!$                    IF( (m4-m3).EQ.(m1-m2).AND. INTERVAL_INTERSECTION(l4,l3,l2,l1)) THEN
!!$                       
!!$                       eltos = eltos + 1 
!!$                       
!!$                       i_ = l4*(l4+1) - m4 + 1
!!$                       j_ = l3*(l3+1) - m3 + 1
!!$                       k_ = l2*(l2+1) - m2 + 1
!!$                       l_ = l1*(l1+1) - m1 + 1
!!$                       
!!$                       !write(*,*) i,eltos,D
!!$                       !GPE_interaction(i) = beta_ijkl(eltos)
!!$                       
!!$                       
!!$                       IF( (i.NE.j) .AND. (k.NE.l) .AND. (j.NE.l)) THEN
!!$                          
!!$                          GPE_interaction(i) = GPE_interaction(i)             + U * DCMPLX(0.0,-1.0) * &
!!$                               & PSI(j)*CONJG(PSI(k))*PSI(l)*beta_ijkl(eltos) + U * DCMPLX(0.0,-1.0) * &
!!$                               & PSI(j)*CONJG(PSI(l))*PSI(k)*beta_ijkl(eltos)
!!$                          GPE_interaction(j) = GPE_interaction(j)             + U * DCMPLX(0.0,-1.0) * &
!!$                               & CONJG(PSI(k))*PSI(l)*beta_ijkl(eltos)        + U * DCMPLX(0.0,-1.0) * &
!!$                               & PSI(i)*CONJG(PSI(l))*PSI(k)*beta_ijkl(eltos)
!!$                          
!!$                          IF(i_.NE.i .AND. i_.NE.j) GPE_interaction(i_) = GPE_interaction(i_) + U * DCMPLX(0.0,-1.0) * &
!!$                               & PSI(j_)*CONJG(PSI(k_))*PSI(l_)*beta_ijkl(eltos)              + U * DCMPLX(0.0,-1.0) * &
!!$                               & PSI(j_)*CONJG(PSI(l_))*PSI(k_)*beta_ijkl(eltos)
!!$                          IF(j_.NE.j .AND. j_.NE.i) GPE_interaction(j_) = GPE_interaction(j_) + U * DCMPLX(0.0,-1.0) * &
!!$                               & PSI(i_)*CONJG(PSI(k_))*PSI(l_)*beta_ijkl(eltos)              + U * DCMPLX(0.0,-1.0) * &
!!$                               & PSI(i_)*CONJG(PSI(l_))*PSI(k_)*beta_ijkl(eltos)
!!$                          
!!$                       END IF
!!$                       IF( (i.NE.j) .AND. (k.EQ.l)) THEN
!!$                          GPE_interaction(i) = GPE_interaction(i) + U * DCMPLX(0.0,-1.0) * PSI(j)* &
!!$                               & CONJG(PSI(k))*PSI(l)*beta_ijkl(eltos) 
!!$                          GPE_interaction(j) = GPE_interaction(j) + U * DCMPLX(0.0,-1.0) * PSI(i)* &
!!$                               & CONJG(PSI(k))*PSI(l)*beta_ijkl(eltos)
!!$                          
!!$                          IF(i_.NE.i .AND. i_.NE.j) GPE_interaction(i_) = GPE_interaction(i_) +  U * DCMPLX(0.0,-1.0)*PSI(j_)* &
!!$                               & CONJG(PSI(k_))*PSI(l_)*beta_ijkl(eltos) 
!!$                          IF(j_.NE.j .AND. j_.NE.i) GPE_interaction(i_) = GPE_interaction(j_) +  U * DCMPLX(0.0,-1.0)*PSI(i_)* &
!!$                               & CONJG(PSI(k_))*PSI(l_)*beta_ijkl(eltos) 
!!$                       END IF
!!$                       
!!$                       IF( (i.EQ.j) .AND. (k.NE.l)) THEN
!!$                          GPE_interaction(i) = GPE_interaction(i) + U * DCMPLX(0.0,-1.0)* &
!!$                               & PSI(j)*CONJG(PSI(k))*PSI(l)*beta_ijkl(eltos) &
!!$                               & + U * DCMPLX(0.0,-1.0)*PSI(j)*CONJG(PSI(l))*PSI(k)*beta_ijkl(eltos)
!!$                          IF(i_.NE.i) GPE_interaction(i_) = GPE_interaction(i_) &
!!$                               & + U * DCMPLX(0.0,-1.0)*PSI(j_)*CONJG(PSI(k_))* &
!!$                               & PSI(l_)*beta_ijkl(eltos) + U * DCMPLX(0.0,-1.0)*PSI(j_)*CONJG(PSI(l_))*PSI(k_)*beta_ijkl(eltos)
!!$                       END IF
!!$                       
!!$                       IF( (i.EQ.j) .AND. (k.EQ.l)) THEN
!!$                          GPE_interaction(i)  = GPE_interaction(i)  +  U * DCMPLX(0.0,-1.0)*PSI(j)* &
!!$                               & CONJG(PSI(k))*PSI(l)*beta_ijkl(eltos) 
!!$                          
!!$                          IF(i_.NE.i) GPE_interaction(i_) = GPE_interaction(i_) +  U * &
!!$                               & DCMPLX(0.0,-1.0)*PSI(j_)*CONJG(PSI(k_))*PSI(l_)*beta_ijkl(eltos) 
!!$                       END IF
!!$                       
!!$                    END IF
!!$                    
!!$                 END IF
!!$              END DO
!!$           END DO
!!$        END DO
!!$     END IF
!!$     
!!$  END DO
!!$  
  
  DO i =1,D 
     l4 =  FLOOR(SQRT(1.0*i - 1))
     m4 =  i-1-l4*(l4+1)
     GPE_interaction(i) =  DCMPLX(0.0,-1.0)*(U*GPE_interaction(I) + PSI(i)*(l4*(l4+1)))
!     WRITE(*,*) l4,m4,abs(GPE_interaction(i))
  END DO


END SUBROUTINE INTERACTIONTERM


SUBROUTINE CONSTANTS_OF_MOTION(C_modes,PSI,beta_ijkl,beta_ijkl_map,ENERGY,ANGULAR_MOMENTUM)


  USE  motionConstatn_interface
  IMPLICIT NONE  
  INTEGER,                        INTENT(IN) :: C_modes
  COMPLEX*16,       DIMENSION(:), INTENT(IN) :: PSI
  DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: beta_ijkl
  DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: beta_ijkl_map
  DOUBLE PRECISION,               INTENT(OUT):: ENERGY
  DOUBLE PRECISION,               INTENT(OUT):: ANGULAR_MOMENTUM
  
  INTEGER i,j,k,l,indx
  INTEGER l1,l2,l3,l4
  INTEGER m1,m2,m3,m4
  INTEGER D,l_maX
  
  D       = C_modes
  l_max   = INT(SQRT(1.0*C_modes))-1
   
!!$  DO indx=1,size(beta_ijkl_map),1
!!$     
!!$     i = int((beta_ijkl_map(indx)-1)/D**3) + 1
!!$     j = int((beta_ijkl_map(indx)-1)/D**2 -D*(i-1)) + 1
!!$     k = (((beta_ijkl_map(indx)-1.0) - (1.0*i-1)*D**3 - (1.0*j-1)*D**2)/(1.0*D)) + 1
!!$     l = ((beta_ijkl_map(indx)-1.0) - (1.0*i-1.0)*D**3 - (1.0*j-1.0)*D**2 - (1.0*k-1.0)*D) + 1
!!$          
!!$     !     WRITE(*,*) beta_ijkl_map(indx),l,k,j,i,l + D*(k-1) + (j-1)*D**2 + (i-1)*D**3
!!$     
!!$     l4= FLOOR(SQRT(1.0*i - 1.0))
!!$     m4= i - 1 - l4*(l4+1)
!!$     
!!$     l3= FLOOR(SQRT(1.0*j - 1.0))
!!$     m3= j - 1 - l3*(l3+1)
!!$     
!!$     l2= FLOOR(SQRT(1.0*k - 1.0))
!!$     m2= k - 1 - l2*(l2+1)
!!$     
!!$     l1= FLOOR(SQRT(1.0*l - 1.0))
!!$     m1= l - 1 - l1*(l1+1)   
!!$     
!!$     ENERGY = ENERGY + CONJG(PSI(i))*PSI(j)*CONJ(PSI(k))*PSI(l)*beta_ijkl(indx)
!!$  END DO


  ENERGY = 0.0
  ANGULAR_MOMENTUM = 0.0

  DO i=1,C_modes
     l1 = FLOOR(SQRT(1.0*i-1.0))
     m1 = i - 1 -l1*(l1+1)
     ENERGY = ENERGY + l1*(l1+1)*ABS(PSI(i))**2 + ABS(PSI(i))**4
     ANGULAR_MOMENTUM = ANGULAR_MOMENTUM + l1*(l1+1)*ABS(PSI(i))**2 
  END DO

  

END SUBROUTINE CONSTANTS_OF_MOTION

