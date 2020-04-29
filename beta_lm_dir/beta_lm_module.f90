MODULE funciones
  IMPLICIT NONE
  INTERFACE
     
     SUBROUTINE I_mrgrnk(XDONT, IRNGT)
       Integer, Dimension (:), Intent (In)  :: XDONT
       Integer, Dimension (:), Intent (Out) :: IRNGT    
     END SUBROUTINE I_mrgrnk

  END INTERFACE

CONTAINS

  FUNCTION stateNo(l,m)
    
    IMPLICIT NONE
    INTEGER  ::  l,m 
    DOUBLE PRECISION  :: stateNo
    
    stateNo = l*(l+1) + m  + 1 ! |l,m> = |stateNo>
    
    RETURN
    
  END FUNCTION stateNo

  
  FUNCTION elto(i,j,k,l,D)
    
    IMPLICIT NONE
    DOUBLE PRECISION   :: i,j,k,l
    INTEGER*8 :: D
    DOUBLE PRECISION  :: elto
    
    elto = l + D*(k-1) + D*D*(j-1) + D*D*D*(1.0*i-1) ! beta_ijkl = beta_elto
    
    RETURN
    
  END FUNCTION elto
  
  FUNCTION done_index(p,bt,bl) 

    ! This function evaluate if the value beta_p has been evaluated:
    !
    ! p: INTEGER indicates the beta elto to be calculated
    ! bt: ARRAY is the list of beta_eltos
    ! bi: ARRAY is the list of beta_indices to order beta_eltos
    !
    !the function implements a bisection to find the index such that the bl(index)-th element of the list of the table corresponds to p.
    
    ! it also takes care in case the root is next to the interval boundaries
    
    IMPLICIT NONE
    INTEGER :: p,i
    INTEGER, DIMENSION(:),INTENT(IN) :: bt,bl
    INTEGER :: index_l,index_r,index_aux,done_index

    index_l   = 1
    index_r   = SIZE(bt)
    
    IF(bt(bl(index_l)).EQ.p) THEN
       done_index = index_l 
       index_aux = index_l 
       RETURN
    END IF
    
    IF(bt(bl(index_r)).EQ.p) THEN    
       done_index = index_r  
       index_aux = index_l 
       RETURN
    END IF

    index_l   = 2

    index_r   = SIZE(bt)-1
    
    
    DO I=0,CEILING(LOG(1.0*SIZE(bt))/LOG(2.0))+1,1  ! How many times can an interval be splitted in two?
       IF(bt(bl(index_r-1)).EQ.p) THEN 
          !write(*,*) "done!", index_r - 1,p,bt(bl(index_r-1))
          done_index = index_r - 1
          RETURN
       END IF
       IF(bt(bl(index_r+1)).EQ.p) THEN 
          !write(*,*) "done!", index_r + 1,p,bt(bl(index_r+1))
          done_index = index_r + 1
          RETURN
       END IF
       IF(bt(bl(index_l-1)).EQ.p) THEN 
          !write(*,*) "done!", index_l - 1,p,bt(bl(index_l-1))
          done_index = index_l + 1
          RETURN
       END IF
       IF(bt(bl(index_l+1)).EQ.p) THEN 
          !write(*,*) "done!", index_l + 1,p,bt(bl(index_l+1))
          done_index = index_l + 1
          RETURN
       END IF
       index_aux = CEILING(0.5*(index_r+index_l))
       !write(*,*) p,index_r,index_l,index_aux,bt(index_aux)
       IF(bt(bl(index_aux)).GT.p) THEN
          index_r = index_aux
       END IF
       
       IF(bt(bl(index_aux)).LT.p) THEN 
          index_l =  index_aux
       END IF
    
       IF(bt(bl(index_aux)).EQ.p) THEN 
          !write(*,*) "done!", index_aux,p,bt(bl(index_aux))
          done_index = index_aux
          RETURN
       END IF      

    END DO

    done_index = index_aux

    RETURN

  END FUNCTION done_index

  FUNCTION INTERVAL_INTERSECTION(l4,l3,l2,l1)
    
    IMPLICIT NONE
    INTEGER  ::  l4,l3,l2,l1
    LOGICAL  :: INTERVAL_INTERSECTION
    
    INTERVAL_INTERSECTION = .FALSE.

    IF(ABS(L1-L2).GE.ABS(L3-L4) .AND. ABS(L1-L2).LE.(L3+L4) ) INTERVAL_INTERSECTION = .TRUE. 
    IF(ABS(L3-L4).GE.ABS(L1-L2) .AND. ABS(L3-L3).LE.(L1+L2) ) INTERVAL_INTERSECTION = .TRUE. 

  

  END FUNCTION INTERVAL_INTERSECTION
  
  SUBROUTINE WRITE_VECTOR(A)
    IMPLICIT NONE
    
    INTEGER, DIMENSION(:), INTENT(IN) :: A
    INTEGER I
    
    DO I = 1,SIZE(A)
       WRITE(*,*) A(I)
    END DO
    RETURN
  END SUBROUTINE WRITE_VECTOR
  
  
END MODULE funciones
