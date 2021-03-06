INCLUDE "beta_lm_module.f90"

!-------------------


!!$ THIS PROGRAM IS PART OF A PROJECT TO EVALUATE TIME EVOLUTION OF A GAS
!!$ CONFINED TO THE SURFACE OF AN SPHERE ACCORDING TO THE PROJECTED GROSS-PITAEVSKII EQUATION.
!!$
!!$ IN THE BASIS OF SPHERICAL HARMONICS, THE WAVE FUNCTIONS READS AS:
!!$
!!$ |PSI> = SUM_L^M C_|l,m>  Y_l^m
!!$
!!$ AND THE PGPE EQUATION TAKES THE FORM:
!!$
!!$ i d C_|l,m>/dt  =  - C_|l,m> (l*(l+1)) + U_0 * F_|l,m>
!!$
!!$ with 

!!$ F_|l,m> =  SUM_l3,m3 SUM_l2,m2 SUM_l1,m1 C_|l3,m3> C_|l2,m2>* C_|l1,m1> \int Y_l^m* Y_l3^m3 Y_l2^m2* Y_l1^m1 d \omega
!!$
!!$
!!$ WITH THIS PROGRAM, THE INTEGRAL OF THE PRODUCT OF FOUR SPHERICAL HARMONICS IS CALCULATED.

!!$ SUCH AN INTEGRAL IS NAMED BETA_IJKL, WHERE EACH INDEX CORRESPONDS TO A STATE |l,m>. THE VALUE IS GIVEN BY:
!!$
!!$ BETA_IJKL = BETA_P = \Sum_{l = \ABS(l1-l2)}^{(l1+l2)} \sqrt{(2l1+1)(2l2+1)(2l3+1)(2l4+1)} (2l+1) (-1)^{-m_3} 
!!$ \left( l1 l2 l \n 0 0 0 \right) \left( l3 l4 l \n 0 0 0\right) \left( l1 l2 l \n m1 -m2 (m2-m1) \right) \left( l3 l4 l \n -m3 m4 m3-m4 \right) 
!!$
!!$ Because the number of coeficientes to evaluate growths as (C_modes*(C_mods+1))**4, we use several to reduce the number of calculations.
!!$
!!$ Given l4 m4, l3,m3, determine the combinations of l2,m2,l1,m1 with non-zero coeficients in the sum for beta_ijkl. This happen if any of the following happen:
!!$
!!$ - abs(l2-l1) and l2+l1 \in [abs(l3-l4),l3+l4]
!!$ - abs(l2-l1)  \in [abs(l3-l4),l3+l4]
!!$ - l2+l1  \in [abs(l3-l4),l3+l4]
!!$ - abs(l2-l1) <  abs(l3-l4) and l1+l2 > l3+l4 (\eq  abs(l3-l4) and l3+l4 \in [abs(l2-l1),l1+l2])
!!$
!!$
!!$ In a l1-l2 plane, this conditions define the region:
!!$
!!$
!!$ l2
!!$ |   /                   /
!!$ |  /       C           /
!!$ | /                   / 
!!$ |/_________________  /
!!$ |         B         /
!!$ |_________________ / 
!!$ |\                / 
!!$ |  \      A      / 
!!$ |    \          /
!!$ |      \       /
!!$ |        \    /
!!$ --------------------------------- l1
!!$
!!$ Also, l1+l2 and l3+l4 should have the same parity, because 
!!$
!!$
!!$ \left( l1 l2 l \n 0 0 0 \right)  = 0 is l1+l2+l is odd
!!$
!!$
!!$
!!$ In addition, given m3 and m4, then m1 and m2 should satisfy:
!!$
!!$ m2-m1 = m3-m4
!!$
!!$ Also, m3-m4 determines the minimun value of l1+l2:
!!$
!!$
!!$      -l1         m2           l1
!!$       .      /   |          .
!!$ ............/...........................l2                 
!!$       .    /     |          .
!!$       .   /      |          .
!!$       .  /       |          .
!!$ --------/----------------------------m1
!!$       ./         |          .
!!$       /          |          . 
!!$      / .         |          .
!!$ ................................-l2
!!$       .          |
!!$       .          |
!!$
!!$ The diagonal line represents m2-m1 = cte. The, box defined by l1 and l2 should be large enough to be crossed by the diagonal.
!!$
!!$ This conditions is satisfied if abs(l1+l2) >= cte

!!$ All this is taken into account in the first loop which builds an array containing the index corresponding to ijkl for which beta !=0.
!!$
!!$ Then, a second series of loops calculate the functions, taking into account the parity condition only for one of the beta's related by symmetry operator.
!!$
!!$ The output of this program are two vectors (written into a file):
!!$
!!$ beta(i) containing as many elements as needed to be calculated
!!$ beta_done(i) indicating the value of beta_ijkl in beta(i)

!!$  This arrays should be used into a program to evaluate the non-linear term of the GPE  
!!$ the first element of beta_done(1) = -1, does not correspond to any beta_ijkl.

!-----------------

PROGRAM BETA_
  USE funciones
  IMPLICIT NONE
  INTEGER :: l1,l2,l3,l4,m1,m2,m3,m4
  INTEGER*8 :: C_modes,D,l_max
!  DOUBLE PRECISION, DIMENSION(1:1) :: PAD1
!  DOUBLE PRECISION, DIMENSION(1:1) :: PAD1_
  INTEGER :: EO1,lmin,lmax
  INTEGER :: i,j,k,l
  INTEGER*8 :: eltos
  DOUBLE PRECISION :: beta

  OPEN(unit=1,file="beta_lm_cmodes3_lmax.dat",action="write")
  OPEN(unit=2,file="beta_lm_cmodes3_calc_lmax.dat",action="write")

  C_modes = 3
  D       = C_modes*C_modes
  l_max   = C_modes-1
 
  eltos = 0
  
  DO i = 1,D
     
     l4 =  FLOOR(SQRT(1.0*i - 1))
     m4 =  i-1-l4*(l4+1)
     !write(*,*) i,l4,m4
     IF(m4.GE.0) THEN
        
        DO j =  i,D
           
           l3 =  FLOOR(SQRT(1.0*j - 1))
           m3 =  j-1-l3*(l3+1)
           EO1 =  MOD(l3+l4,2)
           !        write(*,*) j,l3,m3
           DO k =  1,D
              
              l2 =  FLOOR(SQRT(1.0*k - 1))
              m2 =  k-1-l2*(l2+1)
              
              DO l = k,D,2
                 
                 l1 =  FLOOR(SQRT(1.0*l - 1))
                 m1 =  l-1-l1*(l1+1)
                 !WRITE(*,*) l4,l3,l2,l1,m4,m3,m2,m1              
                 IF(EO1 .EQ.MOD(l2+l1,2)) THEN ! (L3+L4) AND (L1+L2) should have the same parity?
                    
                    IF( (m4-m3).EQ.(m1-m2).AND. INTERVAL_INTERSECTION(l4,l3,l2,l1)) THEN
                       
                       eltos = eltos + 1 
                       
                       lmin = MOD(MAX0(ABS(l1-l2),ABS(l3-l4)),2)  ! TO ENSURE L3+L4+LMIN EVEN
                       IF(lmin.EQ.EO1) THEN
                          lmin = MAX0(ABS(l1-l2),ABS(l3-l4))
                       ELSE
                          lmin = MAX0(ABS(l1-l2),ABS(l3-l4)) + 1
                       END IF
                       !lmax = MIN0(l1+l2,l3+l4)
                       lmax = l_max
                       
                       CALL j3symbols_fun(l1,l2,l3,l4,m1,m2,m3,m4,beta,lmin,lmax)
                       
                       !PAD1   = (/ elto(i,j,k,l,D) /)
                       
                       !PAD1_   = (/ beta /)
                       
                       !IF(PAD1(1).LT.0) WRITE(*,*) i,j,k,l,D,intPAD1!*(i-1)! l + D*(k-1) + D*D*(j-1) + D*D*D*(i-1)
                       WRITE(1,*)  elto( stateNo(l4,m4), stateNo(l3,m3), stateNo(l2,m2), stateNo(l1,m1),D), beta
                    END IF
                 END IF
              END DO
           END DO
        END DO
     END IF
  END DO

  WRITE(2,*) eltos
  WRITE(*,*) "# Number of ijkl combinations evaluated: ", eltos
!!$  WRITE(*,*) "# number of non-zero ijkl combinations: ", SIZE(beta_done)-1
  WRITE(*,*) "# Total number of ijkl combinations: ", D**4
!!$  !WRITE(*,*) "Number of ijkl combinations:", SIZE(beta_done)-1,beta_done(2:SIZE(beta_done))
!!$  DEALLOCATE(beta_table) 
!!$  DEALLOCATE(beta_list) 
!!$  DEALLOCATE(beta_done)

  201 FORMAT(I10," ",1E25.17)
END PROGRAM BETA_


