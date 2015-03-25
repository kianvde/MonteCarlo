!--metropolis-hastings loop--
!
! NB run subroutine init_rnd() once before using the loop
!
! INPUT
! r_in :                [real(6)] => vector containing r1 and r2 [r1, r2]
! problem parameters:   [real(4)] => vector containing problem parameters [s, a, b]
!       s: separation between protons
!       a: trial wave function parameter
!       b: trial wave function parameter
! sigma:                [real]    => standard deviation of update step
! N:                    [integer] => number of updates
! compute:              [boolean] => specify if the energy and the descent factor should be computed
!
!
! OUTPUT
! energy:               [real]    => the averaged local energy
! descent:              [real]    => beta update parameter
! r_out:                [real(6)] => vector containing r1 and r2 [r1, r2]         

                            !  |------------------input-------------------|  |-------output-------|
SUBROUTINE metropolis_hastings(r_in, problem_parameters, sigma, N, compute, energy, descent, r_out)
    implicit none

    ! input/output
    real(8), intent(in)  :: r_in(6)
    real(8), intent(in)  :: problem_parameters(3)
    real(8), intent(in)  :: sigma
    integer, intent(in)  :: N
    logical, intent(in)  :: compute
    real(8), intent(out) :: energy, descent
    real(8), intent(out) :: r_out(6)
    
    ! computation variables
    integer :: i
    real(8) :: rnd(7)
    real(8) :: PI, p
    real(8) :: s, a, b
    
    real(8) :: r1(3), r2(3), r1n(3), r2n(3)
    real(8) :: r1Lu(3), r1Ru(3), r2Lu(3), r2Ru(3), r12u(3), t(3), temp(3)
    real(8) :: r1L(2), r1R(2), r2L(2), r2R(2), r12(2), p1L(2), p1R(2), p2L(2), p2R(2), p12(2)
    real(8) :: sd_A, sd_A_avg, sd_AE_avg, local_energy
    
    ! extract r1 and r2 and the problem parameters
    r1 = r_in(1:3)
    r2 = r_in(4:6)
    s = problem_parameters(1)
    a = problem_parameters(2)
    b = problem_parameters(3)
    
    ! calculate useful constants
    t  = (/ s/2, 0d+0, 0d+0 /)
    PI = 4*atan(1d+0)
    
    ! for loop metropolis-hastings
    energy = 0
    sd_A_avg = 0
    sd_AE_avg = 0
    DO i=1,N
        ! generate an array of random numbers for proposing update for r1 and r2
        rnd = 0
        CALL random_number(rnd)

        ! propose update for r1 and r2 according to a step with uniform 
        ! random angles (phi, theta) and normal distributed radial 
        ! step size (R) with standard deviation s/4
        r1n = r1 + sigma*sqrt(-2*log(rnd(1)))*cos(2*PI*rnd(2))* & 	! 1+R*
            (/  sin(PI*rnd(3))*cos(2*PI*rnd(4)), &                  ! (sin(theta)cos(phi),
                sin(PI*rnd(3))*sin(2*PI*rnd(4)), &                  !  sin(theta)sin(phi),
                cos(PI*rnd(3)) /)                                   !  cos(theta))

        r2n = r2 + sigma*sqrt(-2*log(rnd(1)))*sin(2*PI*rnd(2))* &
            (/  sin(PI*rnd(5))*cos(2*PI*rnd(6)), &
                sin(PI*rnd(5))*sin(2*PI*rnd(6)), &
                cos(PI*rnd(5)) /)

        ! calculate useful quantities in the vector constructors the second 
        ! value is the value for updated positions            
        r1L = (/ sum((r1+t)**2)**0.5, sum((r1n+t)**2)**0.5 /)
        r1R = (/ sum((r1-t)**2)**0.5, sum((r1n-t)**2)**0.5 /)
        r2L = (/ sum((r2+t)**2)**0.5, sum((r2n+t)**2)**0.5 /)
        r2R = (/ sum((r2-t)**2)**0.5, sum((r2n-t)**2)**0.5 /)
        r12 = (/ sum((r1-r2)**2)**0.5, sum((r1n-r2n)**2)**0.5 /)
        r1Lu = (r1+t)/r1L(1)
        r1Ru = (r1-t)/r1R(1)
        r2Lu = (r2+t)/r2L(1)
        r2Ru = (r2-t)/r2R(1)
        r12u = (r1-r2)/r12(1)
        p1L = (/ exp(-r1L(1)/a), exp(-r1L(2)/a) /)
        p1R = (/ exp(-r1R(1)/a), exp(-r1R(2)/a) /)
        p2L = (/ exp(-r2L(1)/a), exp(-r2L(2)/a) /)
        p2R = (/ exp(-r2R(1)/a), exp(-r2R(2)/a) /)
        p12 = (/ exp(r12(1)/(2+2*b*r12(1))), exp(r12(2)/(2+2*b*r12(2))) /)
        
        ! calculate the fraction of the new wave function divided by the old
        ! one
        p = ((exp(-sum((r1n+t)**2)**0.5/a) + exp(-sum((r1n-t)**2)**0.5/a))/ &   ! factor due to electron 1
            (p1L(1) + p1R(1)))* &
            ((exp(-sum((r2n+t)**2)**0.5/a) + exp(-sum((r2n-t)**2)**0.5/a))/ &   ! factor due to electron 2
            (p2L(1) + p2R(1)))* &
            (exp(r12(2)/(2+2*b*r12(2))) / p12(1))                                      ! interaction factor
        
        IF (compute) THEN

            ! calculate the local energy
            temp = (p1L(1)*r1Lu+p1R(1)*r1Ru)/(p1L(1)+p1R(1)) - (p2L(1)*r2Lu+p2R(1)*r2Ru)/(p2L(1)+p2R(1))
            
            local_energy = - 1/(a**2) + (p1L(1)/r1L(1) + p1R(1)/r1R(1))/(a*(p1L(1)+p1R(1))) + &
                (p2L(1)/r2L(1) + p2R(1)/r2R(1))/(a*(p2L(1)+p2R(1))) & 
                - 1/r1L(1) - 1/r1R(1) - 1/r2L(1) - 1/r2R(1) + 1/r12(1) + &
                dot_product(temp, r12u/(2*a*(1+b*r12(1))**2)) &
                - ((4*b+1)*r12(1)+4)/((4*(1+b*r12(1))**4)*r12(1))
            energy = energy + local_energy
  
            ! calculate terms for steepest descent method
            sd_A = -r12(1)**2/(2*(1+b*r12(1))**2)
            sd_A_avg = sd_A_avg + sd_A
            sd_AE_avg = sd_AE_avg + sd_A*local_energy
        END IF

        ! update r1 and r2 according to metropolis hastings algorithm
        IF (p**2 > rnd(7)) THEN
            r1 = r1n
            r2 = r2n
        END IF   
    END DO
    
    ! compute local energy and steepest descent method factor
    IF (compute) THEN
        energy = energy/N
        
        ! calculate descent for beta update
        descent = (2d+0/N)*(sd_AE_avg - energy*sd_A_avg)
    END IF
    
    ! calculate r_out
    r_out = (/ r1, r2 /)
    
END SUBROUTINE metropolis_hastings


!---subroutine for initializing the random value generator---!
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