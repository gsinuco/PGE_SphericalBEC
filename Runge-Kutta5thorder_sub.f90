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
     
     !SUBROUTINE RUNGEKUTTA5THORDER(y,t,dt,yout,yerror)
     SUBROUTINE RUNGEKUTTA5THORDER(y,t,dt,beta_ijkl,beta_ijkl_map,yout,yerror)
       
       USE GPEInteraction
       IMPLICIT NONE
       COMPLEX*16, DIMENSION(:), INTENT(IN) :: y
       COMPLEX*16, DIMENSION(:), INTENT(OUT) :: yout
       DOUBLE PRECISION, DIMENSION(:), INTENT(OUT) :: yerror
       DOUBLE PRECISION, INTENT(IN)  :: t
       DOUBLE PRECISION, INTENT(IN) :: dt
       DOUBLE PRECISION, DIMENSION(:),  INTENT(IN):: beta_ijkl
       DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: beta_ijkl_map
     END SUBROUTINE RUNGEKUTTA5THORDER
     
  END INTERFACE
  
END MODULE RungeKuttaModule


!SUBROUTINE RUNGEKUTTA5THORDER(y,t,dt,yout,yerror)
SUBROUTINE RUNGEKUTTA5THORDER(y,t,dt,beta_ijkl,beta_ijkl_map,yout,yerror)
  
  USE GPEInteraction
  IMPLICIT NONE
  COMPLEX*16, DIMENSION(:), INTENT(IN) :: y
  COMPLEX*16, DIMENSION(:), INTENT(OUT) :: yout
  DOUBLE PRECISION, DIMENSION(:), INTENT(OUT) :: yerror
  DOUBLE PRECISION, INTENT(IN)  :: t
  DOUBLE PRECISION, INTENT(IN) :: dt
  DOUBLE PRECISION, DIMENSION(:),  INTENT(IN):: beta_ijkl
  DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: beta_ijkl_map
  COMPLEX*16, DIMENSION(:), ALLOCATABLE :: K1,K2,K3,K4,K5,K6,y_aux
  
  
  INTEGER i,dim
  
  DOUBLE PRECISION :: a2,a3,a4,a5,a6
  DOUBLE PRECISION :: b21,b31,b32,b41,b42,b43,b51,b52,b53,b54,b61,b62,b63,b64,b65
  DOUBLE PRECISION :: c1,c2,c3,c4,c5,c6

  
  DOUBLE PRECISION :: dc1,dc2,dc3,dc4,dc5,dc6
  
  a2 = 0.2
  a3 = 0.3
  a4 = 0.6
  a5 = 1.0
  a6 = 0.875
  
  b21 = 0.2
  b31 = 3.0/40.0
  b32 = 9.0/40.0
  b41 = 0.3 
  b42 = -0.9
  b43 = 1.2
  b51 = -11.0/54.0
  b52 = 2.5
  b53 = -70.0/27.0
  b54 = 35.0/27.0
  b61 = 1631.0/55296.0
  b62 = 175.0/512.0
  b63 = 575.0/13824.0
  b64 = 44275.0/110592.0
  b65 = 253.0/4096.0
  
  c1 = 37.0/378.0
  c2 = 0.0
  c3 = 250.0/621.0
  c4 = 125.0/594.0
  c5 = 0.0
  c6 = 512.0/1771.0

  dc1 = c1 - 2825.0/27648.0 
  dc3 = c3 - 18575.0/48384.0
  dc2 = 0.0
  dc4 = c4 - 13525.0/55296.0
  dc5 = -277.0/14336.0
  dc6 = c6 - 0.25

  
  ALLOCATE(K1(SIZE(Y)))
  ALLOCATE(K2(SIZE(Y)))
  ALLOCATE(K3(SIZE(Y)))
  ALLOCATE(K4(SIZE(Y)))
  ALLOCATE(K5(SIZE(Y)))
  ALLOCATE(K6(SIZE(Y)))
  ALLOCATE(y_aux(SIZE(Y)))

  
  
  !write(*,*) "k1"
  K1 = f(t,y,beta_ijkl,beta_ijkl_map) 
  !K1 = f(t,y) 
  
  y_aux = y + dt*b21*k1
  !write(*,*) y_aux
  K2 = f(t+a2*dt,y_aux,beta_ijkl,beta_ijkl_map)
  !K2 = f(t+a2*dt,y_aux)
  
  y_aux = y + dt*(b31*K1 + b32*k2)
  !write(*,*) y_aux
  K3 = f(t+a3*dt,y_aux,beta_ijkl,beta_ijkl_map)
  !K3 = f(t+a3*dt,y_aux)
  
  y_aux = y + dt*(b41*K1 + b42*K2 + b43*K3)
  K4 = f(t+a4*dt,y_aux,beta_ijkl,beta_ijkl_map)
  !write(*,*) y_aux
  !K4 = f(t+a4*dt,y_aux)
  
  y_aux = y + dt*(b51*K1 + b52*K2 + b53*K3 + b54*K4)
  K5 = f(t+a5*dt,y_aux,beta_ijkl,beta_ijkl_map)
  !write(*,*) y_aux
  !K5 = f(t+a5*dt,y_aux)
  
  y_aux = y + dt*(b61*K1 + b62*K2 + b63*K3 + b64*K4 + b65*K5)
  K6 = f(t+a6*dt,y_aux,beta_ijkl,beta_ijkl_map)
  ! write(*,*) y_aux
  !K6 = f(t+a6*dt,y_aux)
  
  
  yout = y + dt*(c1*K1  +  c3*K3  + c4*K4  + c6*K6) 
  !  write(*,*)"ks", K1,K3
  yerror =  ABS(1.0-DOT_PRODUCT(yout,yout))
  !yerror = ABS(dt*(dc1*K1  + dc2*K2  + dc3*K3  + dc4*K4  + dc5*K5  + dc6*K6))/ABS(yout)!+dt*f(t+dt,yout,beta_ijkl,beta_ijkl_map)))
  
END SUBROUTINE RUNGEKUTTA5THORDER


