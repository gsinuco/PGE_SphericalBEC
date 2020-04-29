MODULE RKadaptiveCts
  DOUBLE PRECISION :: SAFETY = 0.9
  DOUBLE PRECISION :: ERRCON = 1.84E-9
  DOUBLE PRECISION :: tolerance =   1.0E-5
  DOUBLE PRECISION :: STEP_ERROR
  LOGICAL          :: DONE = .FALSE.
END MODULE RKadaptiveCts

PROGRAM PGPESHPERE
  
  USE RungeKuttaModule
  USE interactionStrength
  USE RKadaptiveCts
  USE motionConstatn_interface
 
  IMPLICIT NONE
  
  TYPE :: C_field      ! Coherent Field in real space
     COMPLEX*16, DIMENSION(:,:), allocatable :: phi 
     INTEGER n_c     
  END TYPE C_field
  
  TYPE :: I_field      ! Incoherent field in real space
     COMPLEX*16, DIMENSION(:,:), allocatable :: phi
     INTEGER n_i     
  END TYPE I_field
  
  TYPE :: PGPE_field   ! Field in mode space   
     COMPLEX*16, DIMENSION(:), allocatable :: psi
     INTEGER n_c     
  END TYPE PGPE_field

  TYPE :: sphere       ! Sphere Characteristics
     DOUBLE PRECISION R;
     INTEGER n_theta
     INTEGER n_phi     
  END TYPE sphere
  
  TYPE :: PGPE         ! Projection characterisitics
     DOUBLE PRECISION  E_cut     
     INTEGER C_modes
     INTEGER I_modes
  END TYPE PGPE    

  TYPE(sphere) s
  TYPE(PGPE) P
  TYPE(C_field) C
  !TYPE(I_field) I
  TYPE(PGPE_field) Projection
  
  DOUBLE PRECISION time,timeStep,PI,AtomNumber,r,NORM
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: RandomAmplitudes, RandomPhases
  INTEGER       :: ntheta,nphi, seed, iterations, k,l,Cmodes,i,l_max,C_modes 

  DOUBLE PRECISION, DIMENSION(:),  ALLOCATABLE:: beta_ijkl
  DOUBLE PRECISION, DIMENSION(:),  ALLOCATABLE :: beta_ijkl_map
  INTEGER :: betas, betas_calc
  

  COMPLEX*16, DIMENSION(:), ALLOCATABLE :: yout
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: yerror


  DOUBLE PRECISION :: ENERGY,ANGULAR_MOMENTUM

  character*32 cc


  common/geometry/r,ntheta,nphi
  common/dimensions/timeStep,Cmodes
  
 
  PI = 4.0*ATAN(1.0)
  AtomNumber = 1.0 

  s%R       = 1
  s%n_theta = 3
  s%n_phi   = 3

  C_modes =  10 ! as in beta_lm_dir/beta_lm.f90
  l_max = C_modes - 1
  
  P%C_modes = (C_modes)**2
  !P%I_modes = 2
  P%E_cut   = 10

  r      = s%R
  ntheta = s%n_theta
  nphi   = s%n_phi
  Cmodes = P%C_modes
  


  
  cc = '#Allocating memory:'
!  WRITE(*,*) cc

  ALLOCATE(C%phi(s%n_theta,s%n_phi))    ! Field in real space
  !ALLOCATE(I%phi(s%n_theta,s%n_phi)) 
  ALLOCATE(Projection%psi(P%C_modes))   ! Field in mode space
  ALLOCATE(RandomAmplitudes(P%C_modes)) ! needed to set initial field  
  ALLOCATE(RandomPhases(P%C_modes))     ! needed to set initial field

  ALLOCATE(yout(P%C_modes))   ! Advanced field
  ALLOCATE(yerror(P%C_modes)) ! estimated error

  C%phi  = DCMPLX(0.0,0.0)
  C%n_c  = P%C_modes

  Projection%psi = DCMPLX(0,0)
  Projection%n_c    = P%C_modes

  cc = '#Populating field modes:'
!  WRITE(*,*) cc
  !----------------------------------------------------------------
  !-------------------- POPULATE THE FIELD MODES ------------------
  !----------------------------------------------------------------
  
  
  ! GENERATION OF THE  THE INITIAL RANDOM STATE:
  !
  ! zufall random number generator package is distributed by NETLIB:
  !
  ! http://www.netlib.org/random/index.html
  ! 
  ! Phases are distributed uniformily on the interval [0,2*pi]
  ! Amplitudes have a Gaussian distribution of sigma = 1 and mean = 0.
  !
  ! 
  ! First: Initialize the seed for the random generators:
  
  seed = 0
  call ZUFALLI(seed)
  
  ! Then, call the random generators:
  
  !CALL NORMALEN(size(RandomAmplitudes),RandomAmplitudes) ! for the amplitudes
  CALL ZUFALL(size(RandomAmplitudes),RandomAmplitudes)           ! for the amplitude
  !RandomAmplitudes(P%C_modes-10:P%C_modes) = 0.0        ! CAUTION: cutting population.
  CALL ZUFALL(size(RandomPhases),RandomPhases)           ! for the phases
  RandomPhases = 2*pi*RandomPhases                       ! rescale the interval of random phases 
  
  ! With phases and amplitudes, set the initial values of the mode coeficientes: 

  RandomAmplitudes = 0.5*RandomAmplitudes
  RandomAmplitudes(1) = 0.5  + 0.5*RandomAmplitudes(1) 

  Projection%psi = DCMPLX(0.0,0.0)
  DO k =1,4!P%C_modes
     Projection%psi(k) = (DCMPLX(RandomAmplitudes(k),0.0))*DCMPLX(COS(RandomPhases(k)) &
          & ,SIN(RandomPhases(k)))
     !WRITE(*,*) k,REAL(projection%psi(k)),AIMAG(projection%psi(k))
  END DO

  !Normalization of the wavefunction

  NORM = DOT_PRODUCT(Projection%psi,Projection%psi)

  Projection%psi = SQRT(AtomNumber)*Projection%psi/SQRT(NORM)

  NORM = DOT_PRODUCT(Projection%psi,Projection%psi)
  !WRITE(*,*) NORM
 
 
  DEALLOCATE(RandomAmplitudes)  ! Not needed any longer
  DEALLOCATE(RandomPhases)      ! Not needed any longer

  cc = '#time evolution:'
!  WRITE(*,*) cc
  
  !----------------------------------------------------------------
  !--------------------- EVOLUTION OF THE FIELD -------------------
  !----------------------------------------------------------------
  !
  ! RKFOURTHORDER: Time evolution of the differential equation df_dt = G(t), with a Runge-Kutta fourth order.
  !                This is my version of Runge-Kutta. 
  !                It is located in file:
  !                   Runge-Kutta4thorder_sub.f90
  !                do not forget to include the statement:
  !                   USE RungeKuttaModule
  !                at the begining of this program 
  !

  ! To evaluate the non-linear term of the GPE, first read and store the needed coeficients

  OPEN(UNIT=1,FILE='beta_lm_dir/beta_lm_cmodes10_lmax.dat',ACTION='READ')
  OPEN(UNIT=2,FILE='beta_lm_dir/beta_lm_cmodes10_calc_lmax.dat',ACTION='READ')
  OPEN(UNIT=3,FILE='PGPEonSphere_cmodes10_re_HT.dat',ACTION='WRITE')
  OPEN(UNIT=4,FILE='PGPEonSphere_cmodes10_im_HT.dat',ACTION='WRITE')
  OPEN(UNIT=5,FILE='PGPEonSphere_cmodes10_psi_HT.dat',ACTION='WRITE')
  
  !READ(2,*) betas 
  READ(2,*) betas_calc
  
  ALLOCATE(beta_ijkl(betas_calc))
  ALLOCATE(beta_ijkl_map(betas_calc))
  
  DO I=1,betas_calc
     READ(1,*) beta_ijkl_map(I),beta_ijkl(I)
!     WRITE(*,*) betas_calc,I, beta_ijkl_map(I),beta_ijkl(I)
!     write(*,*) beta_ijkl_map(I),beta_ijkl(I)
  END DO
 
  ! then, time evolve the field:

  iterations =  500
  timeStep = 2*pi/iterations
  time  = 0.0  
  !WRITE(*,*) time,(REAL(Projection%psi(l)), l=1,P%C_modes)
  DO k=1,20*iterations
!     CALL RKFOURTHORDER(Projection%psi,time,timeStep,beta_ijkl,beta_ijkl_map)
     DONE = .FALSE.
     
     DO WHILE (DONE .EQV. .FALSE.)
        CALL RUNGEKUTTA5THORDER(Projection%psi,time,timeStep,beta_ijkl,beta_ijkl_map,yout,yerror)
        !write(*,*) yerror
        yerror = yerror/tolerance
        STEP_ERROR = MAXVAL(yerror,1)
        !write(*,*) step_error
        IF(STEP_ERROR .LE. 1.0 ) THEN
           DONE = .TRUE.
        ELSE
           !write(*,*) yerror(1)
           IF(timeStep.GT.0.0) THEN
              timeStep = DMAX1(SAFETY*timeStep*(STEP_ERROR**(-0.25)),0.1*timeStep)
              !write(*,*) "new dt:", dt
           ELSE
              timeStep = DMIN1(SAFETY*timeStep*(STEP_ERROR**(-0.25)),0.1*timeStep)
           END IF
           IF(timeStep.LT.5E-3) THEN
              timeStep = 5E-3
              DONE = .TRUE.
           END IF
        END IF
        
     END DO
     Projection%psi = yout
     time = time + timeStep
     IF(STEP_ERROR.GT.ERRCON) THEN
        timeStep =  SAFETY*timeStep*(STEP_ERROR**(-0.2))
     ELSE
        timeStep  = 2.0*timeStep
     END IF

     NORM = DOT_PRODUCT(Projection%psi,Projection%psi)
     CALL CONSTANTS_OF_MOTION(C_modes,Projection%psi,beta_ijkl,beta_ijkl_map,ENERGY,ANGULAR_MOMENTUM)
     WRITE(*,*) time, timeStep,NORM,ENERGY,ANGULAR_MOMENTUM
     WRITE(3,*) time,(REAL(Projection%psi(l)), l=1,P%C_modes)
     WRITE(4,*) time,(AIMAG(Projection%psi(l)), l=1,P%C_modes)
     WRITE(5,*) time,(ABS(Projection%psi(l)),l=1,P%C_modes)
  END DO
  
!  WRITE(*,*) '#End of PGPEonSphere'
  
END PROGRAM PGPESHPERE

