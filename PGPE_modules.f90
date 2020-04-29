MODULE funciones
  
  IMPLICIT NONE
CONTAINS
  
  FUNCTION stateNo(l,m)
    
    IMPLICIT NONE
    INTEGER  ::  l,m 
    INTEGER  :: stateNo
    
    stateNo = l*(l+1) + m  + 1 ! |l,m> = |stateNo>
    
    RETURN
    
  END FUNCTION stateNo

  FUNCTION INTERVAL_INTERSECTION(l4,l3,l2,l1)
    
    IMPLICIT NONE
    INTEGER  ::  l4,l3,l2,l1
    LOGICAL  :: INTERVAL_INTERSECTION
    
    INTERVAL_INTERSECTION = .FALSE.
    
    IF(ABS(L1-L2).GE.ABS(L3-L4) .AND. ABS(L1-L2).LE.(L3+L4) ) INTERVAL_INTERSECTION = .TRUE. 
    IF(ABS(L3-L4).GE.ABS(L1-L2) .AND. ABS(L3-L3).LE.(L1+L2) ) INTERVAL_INTERSECTION = .TRUE. 
    
    
    
  END FUNCTION INTERVAL_INTERSECTION

END MODULE funciones



MODULE interactionStrength
  DOUBLE PRECISION :: U = 5.0
END MODULE interactionStrength
