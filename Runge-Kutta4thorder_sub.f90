! gfortran Runge-Kutta4thorder_sub.f90 zufall_sub.f
! Do not forget to include:
!     USE RungeKuttaModule
! in the program using this subroutine

!include "PGPE_modules.f90"

MODULE interaction_interface
  IMPLICIT NONE
  INTERFACE
     SUBROUTINE INTERACTIONTERM(C_modes, PSI,GPE_interaction,beta_ijkl,beta_ijkl_map)
       USE funciones
       IMPLICIT NONE
       INTEGER,                       INTENT(IN) :: C_modes
       DOUBLE PRECISION, DIMENSION(:),INTENT(IN) :: beta_ijkl
       DOUBLE PRECISION, DIMENSION(:),INTENT(IN) :: beta_ijkl_map
       COMPLEX*16,       DIMENSION(:),INTENT(IN) :: PSI
       COMPLEX*16,       DIMENSION(:),INTENT(OUT):: GPE_interaction
     END SUBROUTINE INTERACTIONTERM
  END INTERFACE
END MODULE interaction_interface


MODULE GPEInteraction
  IMPLICIT NONE
  
CONTAINS
  
  FUNCTION  F(t,y, beta_ijkl,beta_ijkl_map) RESULT(K)

    USE interaction_interface
    
    IMPLICIT NONE
    COMPLEX*16,       DIMENSION(:), INTENT(IN) :: y
    DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: beta_ijkl
    DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: beta_ijkl_map
    DOUBLE PRECISION,               INTENT(IN) ::  t
    COMPLEX*16,       DIMENSION(SIZE(y))       :: K    
    
  
    !K(1) =  y(2) 
    !K(2) = -y(1)-0.25*y(2)  
    
    CALL INTERACTIONTERM(SIZE(y), y, K, beta_ijkl,beta_ijkl_map)
    
    !WRITE(*,*) K(2)
    
  END FUNCTION F
  
END MODULE GPEInteraction

MODULE RungeKuttaModule
  IMPLICIT NONE
  INTERFACE

     SUBROUTINE RKFOURTHORDER(y,t,dt,beta_ijkl,beta_ijkl_map)       
       USE GPEInteraction
       IMPLICIT NONE
       COMPLEX*16, DIMENSION(:), INTENT(IN) :: y
       DOUBLE PRECISION, INTENT(INOUT)  :: t
       DOUBLE PRECISION, INTENT(IN) :: dt
       DOUBLE PRECISION, DIMENSION(:),  INTENT(IN):: beta_ijkl
       DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: beta_ijkl_map
     END SUBROUTINE RKFOURTHORDER
     
  END INTERFACE
  
END MODULE RungeKuttaModule


!!$PROGRAM RungeKuttaProgram
!!$  USE RungeKuttaModule
!!$  IMPLICIT NONE
!!$  
!!$  COMPLEX*16, DIMENSION(:), ALLOCATABLE :: y
!!$  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: y_aux
!!$  DOUBLE PRECISION :: t,dt,pi
!!$  
!!$  INTEGER dim,i, iterations,j, seed
!!$  double precision CYCLE
!!$  data CYCLE/2.0E-8/
!!$
!!$  dim  = 2
!!$  ALLOCATE(y(dim))
!!$  ALLOCATE(y_aux(dim))
!!$  
!!$
!!$  pi = 4.0*ATAN(1.0)
!!$  iterations  = 100
!!$
!!$
!!$  dt = 2*pi/iterations
!!$  t  = 0.0  
!!$  
!!$  !  y(1) = 1.0
!!$  !  y(2) = 0.0
!!$  seed = 0
!!$  call zufalli(seed)
!!$  !CALL normalt(dim,y_aux)! 
!!$  CALL NORMALEN(size(y_aux),y_aux)
!!$  !write(*,*) size(y_aux)
!!$  WRITE(*,*) t,(y_aux(i), i=1,dim)
!!$  y(1) = y_aux(1)
!!$  y(2) = y_aux(2)
!!$  
!!$  DO j=1,4*iterations
!!$     CALL RKFOURTHORDER(y,t,dt)
!!$     WRITE(*,*) t,(REAL(y(i)), i=1,dim)
!!$  END DO
!!$  
!!$
!!$END PROGRAM RungeKuttaProgram





SUBROUTINE RKFOURTHORDER(y,t,dt, beta_ijkl,beta_ijkl_map)

  USE GPEInteraction
  
  IMPLICIT NONE
  COMPLEX*16, DIMENSION(:), INTENT(INOUT) :: y
  DOUBLE PRECISION, INTENT(INOUT)  :: t
  DOUBLE PRECISION, INTENT(IN) :: dt
  DOUBLE PRECISION, DIMENSION(:), INTENT(IN):: beta_ijkl
  DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: beta_ijkl_map
  
  COMPLEX*16, DIMENSION(:), ALLOCATABLE :: K1,K2,K3,K4,y_aux
  
  
  INTEGER i,dim
  
  ALLOCATE(K1(SIZE(Y)))
  ALLOCATE(K2(SIZE(Y)))
  ALLOCATE(K3(SIZE(Y)))
  ALLOCATE(K4(SIZE(Y)))
  ALLOCATE(y_aux(SIZE(Y)))
  
  !write(*,*) "k1"
  K1 = f(t,y,beta_ijkl,beta_ijkl_map) 
  y_aux = y +K1*dt*0.5
  
  !write(*,*) "k2"
  K2 = f(t+0.5*dt,y_aux,beta_ijkl,beta_ijkl_map)
  y_aux = y +K2*dt*0.5
  
  !write(*,*) "k3"
  K3 = f(t+0.5*dt,y_aux,beta_ijkl,beta_ijkl_map)
  y_aux = y +K3*dt
  
  !write(*,*) "k4"
  K4 = f(t+dt,y_aux,beta_ijkl,beta_ijkl_map)
  y = y + dt*(K1 + 2.0*K2 + 2.0*K3 + K4)/6.0
  
  t = t + dt;
  
END SUBROUTINE RKFOURTHORDER
  
