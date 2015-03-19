! a:                
! parameter in the proposed wave function phi = sqrt(a/PI)*exp(-a*x^2)
!
! local_energy:     
! the local energy calculated with parameter a using metrop-hast.
SUBROUTINE metropolis_hastings(a, local_energy)

    implicit none
    integer, intent (out) :: local_energy
    integer :: i
    real(8), intent (in) :: a
    real(8) :: x, xn, phi_old, phi_new, p, PI
    real(8) :: rnd(1:2)
    
    CALL init_rnd()
    PI = 4*atan(1d+0)
    
    ! loop big amount of times to obtain equilibrium
    DO i=1,5000000
        CALL random_number(rnd)
        
        ! update x with a random step
        x = -5
        xn = x + 0.5*rnd(1)
        
        ! calculate the wave function squared for the old and new value for x
        phi_old = sqrt(a/PI)*exp(-a*x**2)
        phi_new = sqrt(a/PI)*exp(-a*xn**2)
        
        ! calculate stepping probability
        p = (phi_new/phi_old)**2
        
        ! make the step if p > 1 or p > random number
        IF (p > rnd(2)) THEN
            x = xn
        END IF
        
    END DO
    
    ! loop 1000 times to average local energy
    DO i=1,1000
        CALL random_number(rnd)
        
        ! update x with a random step
        x = -5
        xn = x + 0.5*rnd(1)
        
        ! calculate the wave function squared for the old and new value for x
        phi_old = sqrt(a/PI)*exp(-a*x**2)
        phi_new = sqrt(a/PI)*exp(-a*xn**2)
        
        ! calculate stepping probability
        p = (phi_new/phi_old)**2
        
        ! make the step if p > 1 or p > random number
        IF (p > rnd(2)) THEN
            x = xn
        END IF
        ! ???????????????????????????????????????
        local_energy = local_energy + 3
        ! ???????????????????????????????????????
    END DO

END SUBROUTINE metropolis_hastings

! subroutine for initializing the random value generator
SUBROUTINE init_rnd()
    implicit none
    integer, dimension(:), allocatable :: seed
    integer, dimension(8) :: dtVals
    integer :: seedSize

    ! get a seed from system time
    CALL date_and_time(VALUES=dtVals)

    ! get size needed for random_seed
    CALL random_seed(SIZE=seedSize)

    ! intialize random generator with seed
    allocate(seed(seedSize))
    CALL random_seed(PUT=dtVals((9-seedSize):8))

END SUBROUTINE init_rnd