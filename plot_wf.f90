
SUBROUTINE metropolis_hastings(bucket)
    implicit none

    integer, intent (out) :: bucket(0:24, 0:24, 0:24)
    integer :: i, bi_x, bi_y, bi_z
    real(8) :: r1(0:2), r2(0:2), r1n(0:2), r2n(0:2), s(0:2)
	
    real(8) :: rnd(1:9)
    real(8) :: PI, r12, r12n, p, p_r1, p_r2, p_r12

    s  = (/ 2, 0, 0 /)
    r1 = (/ -1, 0, 0 /)
    r2 = (/ 1, 0, 0 /)

    PI = 4*atan(1d+0)
    CALL init_rnd()
    DO i=1,50000000

        ! generate an array of random numbers for proposing update for r1 and r2
        rnd = 0
        CALL random_number(rnd)

        ! propose update for r1 and r2 according to a step with uniform 
        ! random angles (phi, theta) and normal distributed radial 
        ! step size (R) ?WITH STANDARD DEVIATION...?
        r1n = r1 + 0.1*sqrt(-2*log(rnd(1)))*cos(2*PI*rnd(2))* & 	! 1+R*
            (/  sin(PI*rnd(3))*cos(2*PI*rnd(4)), &                  ! (sin(theta)cos(phi),
                sin(PI*rnd(3))*sin(2*PI*rnd(4)), &                  !  sin(theta)sin(phi),
                cos(PI*rnd(3)) /)                                   !  cos(theta))

        r2n = r2 + 0.1*sqrt(-.02*log(rnd(5)))*cos(2*PI*rnd(6))* &
            (/  sin(PI*rnd(7))*cos(2*PI*rnd(8)), &
                sin(PI*rnd(7))*sin(2*PI*rnd(8)), &
                cos(PI*rnd(7)) /)

        ! calculate the fraction of the new wavefunction divided by the old
        ! one
        r12n = sum((r1n-r2n)**2)*0.5
        r12 = sum((r1-r2)**2)*0.5
        p_r1 =  exp(-sum((r1n+s)**2)**0.5)+exp(-sum((r1n-s)**2)**0.5)/ &	! factor due to electron 1
                exp(-sum((r1+s)**2)**0.5)+exp(-sum((r1-s)**2)**0.5)
        p_r2 =  exp(-sum((r2n+s)**2)**0.5)+exp(-sum((r2n-s)**2)**0.5)/ &	! factor due to electron 2
                exp(-sum((r2+s)**2)**0.5)+exp(-sum((r2-s)**2)**0.5)
        p_r12 = exp(r12n/(2+2*r12n)) / exp(r12/(2+2*r12))                   ! interaction factor
        p = p_r1*p_r2*p_r12

        ! update r1 and r2 according to metropolis hastings algorithm
        IF (p > rnd(9)) THEN
            r1 = r1n
            r2 = r2n
        END IF

        bi_x = int(25*(r1(0)+4)/8)
        bi_y = int(25*(r1(1)+4)/8)
        bi_z = int(25*(r1(2)+4)/8)
        IF (((bi_x >=0 ) .and. (bi_x < 25)) .and. &
            ((bi_y >=0 ) .and. (bi_y < 25)) .and. &
            ((bi_z >=0 ) .and. (bi_z < 25))) THEN
            bucket(bi_x,bi_y,bi_z) = bucket(bi_x,bi_y,bi_z) + 1
        END IF

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
    CALL RANDOM_SEED(PUT=dtVals((9-seedSize):8))

END SUBROUTINE init_rnd