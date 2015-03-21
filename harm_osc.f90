! a:                
! parameter in the proposed wave function phi = sqrt(a/PI)*exp(-a*x^2)
!
! energy:     
! the local energy calculated with parameter a using metrop-hast.
SUBROUTINE metropolis_hastings(a, energy)

    implicit none
    real(8), intent (out) :: energy
    real(8), intent (in) :: a
    
    integer :: i
    real(8) :: x, xn, phi_old, phi_new, p, PI
    real(8) :: rnd(1:3)
    
    CALL init_rnd()
    PI = 4*atan(1d+0)
    ! loop to obtain equilibrium
    x = -10
    energy = 0
    DO i=1,5000
        CALL random_number(rnd)
        
        ! update x with a normally distributed random step
        xn = x + sqrt(-2*log(rnd(1)))*sin(2*PI*rnd(2))
        
        ! calculate the wave function squared for the old and new value for x
        phi_old = sqrt(a/PI)*exp(-a*x**2)
        phi_new = sqrt(a/PI)*exp(-a*xn**2)
        
        ! calculate stepping probability
        p = (phi_new/phi_old)**2
        
        ! make the step if p > 1 or p > random number
        IF (p > rnd(3)) THEN
            x = xn
        END IF
        
    END DO
    
    ! loop 5000 times to calculate the integral
    DO i=1,100000
        CALL random_number(rnd)
        
        ! update x with a normally distributed random step
        xn = x + sqrt(-2*log(rnd(1)))*sin(2*PI*rnd(2))
        
        ! calculate the wave function squared for the old and new value for x
        phi_old = sqrt(a/PI)*exp(-a*x**2)
        phi_new = sqrt(a/PI)*exp(-a*xn**2)
        
        ! calculate stepping probability
        p = (phi_new/phi_old)**2
        
        ! make the step if p > 1 or p > random number
        IF (p > rnd(3)) THEN
            x = xn
        END IF
        
        ! calculate the energy using the metropolis-hastings generated points
        ! and the local energy
        energy = energy + (a + (0.5 - 2*a**2)*x**2)/100000
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