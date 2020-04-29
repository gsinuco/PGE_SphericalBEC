! THIS PROGRAM EVALUATES THE INTERCTION TERM:
!
!  \sum_j \sum_k \sum_l c_j c_k c_l beta_ijkl
!
! OF THE GROSS-PITAEVSKII EQUATION, IN THE SPHERICAL HARMONIC BASIS.
! 
INCLUDE "beta_lm_module.f90"


PROGRAM INTERACTIONTERM
  USE funciones
  IMPLICIT NONE
  
  INTEGER I,J,K,L
  INTEGER l1,l2,l3,l4
  INTEGER m1,m2,m3,m4
  INTEGER C_modes, D,l_max
  INTEGER EO1,EO2,l1_min,m_,p,i_
  INTEGER betas      ! number of non-zero beta_ijkl 
  INTEGER betas_calc ! number of calculated beta_ijkl coeficientes 
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: beta_ijkl
  INTEGER, DIMENSION(:), ALLOCATABLE :: beta_ijkl_map
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE:: GPE_interaction,PSI

  OPEN(UNIT=1,FILE='beta_lm.dat',ACTION='READ')
  OPEN(UNIT=2,FILE='beta_lm_map.dat',ACTION='READ')

  C_modes = 6
  D       = C_modes*(C_modes-1) + C_modes
  l_max   = C_modes-1
  ALLOCATE(GPE_interaction(D))
  ALLOCATE(PSI(D))
  
  READ(2,*) betas
  READ(2,*) betas_calc
  
  ALLOCATE(beta_ijkl(betas_calc))
  ALLOCATE(beta_ijkl_map(betas))
   
  DO I=1,betas_calc
     READ(1,*) beta_ijkl(I)
  END DO
  DO I=2,betas
     READ(2,*) beta_ijkl_map(I-1)
  END DO
  DO I=1,D
     CALL RANDOM_NUMBER(PSI(i))
  END DO

  i_ = 0
  
  DO l4 = 0,C_modes-1,1
     DO m4 = -l4,l4
        I = stateNo(l4,m4)
        
        DO l3 = 0,C_modes-1,1
           EO1 =  MOD(l3+l4,2)
           
           DO m3 = -l3,l3                        
              
              J = stateNo(l3,m3)
              
              IF(ABS(m3-m4).LE.ABS(l3-l4)) THEN
                 m_ = ABS(l3-l4)
              ELSE
                 m_ = ABS(m3-m4)
              END IF
              
              EO2 =  MOD(m_,2)
              
              !WRITE(*,*) "#A"
              DO l2 = 0,m_,1 !  A
                 
                 l1_min = m_-l2
                 IF(EO1.NE.MOD(l2+l1_min,2)) THEN
                    l1_min = l1_min + 1                   
                 END IF
                 
                 DO l1 = l1_min,MIN0(l3+l4+l2,l_max),2
                    DO m2=MAX0(m3-m4-l1,-l2),MIN0(l2,l1+m3-m4),1
                       
                       m1 = m2-(m3-m4)
                       K  = stateNo(l2,m2)
                       L = stateNo(l1,m1)
                       
                       !p = elto(I,J,K,L,D) 
                       i_ = i_+1
                       !write(*,*) i_,beta_ijkl_map(i_),size(beta_ijkl_map),beta_ijkl(beta_ijkl_map(i_))
                       !write(*,*) size(beta_ijkl), beta_ijkl_map(p)
                       GPE_interaction(I) = GPE_interaction(I) + PSI(J)*PSI(K)*PSI(L)*&
                                                   & beta_ijkl(beta_ijkl_map(i_))
                       
                    END DO !m2
                 END DO !l1
              END DO !l1
              
              !WRITE(*,*) "#B"
              DO l2 = m_ + 1, MIN0(l3+l4,l_max),1 !  B
                 l1_min = 0
                 IF(EO1.NE.MOD(l2+l1_min,2)) THEN
                    l1_min = l1_min + 1                   
                 END IF
                 
                 DO l1 = l1_min,MIN0(l3+l4+l2,l_max),2
                    DO m2=MAX0(m3-m4-l1,-l2),MIN0(l2,l1+m3-m4),1
                       m1 = m2-(m3-m4)
                       
                       K  = stateNo(l2,m2)
                       L = stateNo(l1,m1)
                       
                       
                       !p = elto(I,J,K,L,D) 
                       i_ = i_+1
                       
                       !write(*,*) i_,beta_ijkl_map(i_),size(beta_ijkl_map),beta_ijkl(beta_ijkl_map(i_))
                       GPE_interaction(I) = GPE_interaction(I) + PSI(J)*PSI(K)*PSI(L)*&
                            & beta_ijkl(beta_ijkl_map(i_))
                       
                    END DO !m2
                 END DO !l1
              END DO  !l2                   
              
              !WRITE(*,*) "#C"
              DO l2 = l3+l4+1,l_max,1 !  C
                 
                 l1_min = l2 - (l3+l4) 
                 IF(EO1.NE.MOD(l2+l1_min,2)) THEN
                    l1_min = l1_min + 1                   
                 END IF
                 
                 DO l1 = l1_min,MIN0(l3+l4+l2,l_max),2
                    DO m2=MAX0(m3-m4-l1,-l2),MIN0(l2,l1+m3-m4),1
                       m1 = m2-(m3-m4)
                       K  = stateNo(l2,m2)
                       L = stateNo(l1,m1)                       
                       
                       !p = elto(I,J,K,L,D) 
                       i_ = i_+1
                       
                       !write(*,*) i_,beta_ijkl_map(i_),size(beta_ijkl_map),beta_ijkl(beta_ijkl_map(i_))
                       GPE_interaction(I) = GPE_interaction(I) + PSI(J)*PSI(K)*PSI(L)*&
                            & beta_ijkl(beta_ijkl_map(i_))
                    END DO !m2
                 END DO !l1
              END DO  !l2              



           END DO!m3
        END DO !l3
        WRITE(*,*) I, GPE_interaction(I)
     END DO!m4
  END DO!l4

END PROGRAM INTERACTIONTERM
