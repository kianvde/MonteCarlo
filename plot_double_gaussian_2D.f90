SUBROUTINE metropolis_hastings(bucket)
    implicit none

    integer, intent (out) :: bucket(0:99, 0:99)
    integer :: i, bi_x, bi_y
    real(8) :: p, p1, p2, PI
    real(8) :: rnd(1:7), r1(0:1), r1n(0:1), r2(0:1), r2n(0:1), s(0:1)
    
    s = (/ 1, 0 /)
    
    PI = 4*atan(1d+0)
    DO i=1,5000000
        CALL random_number(rnd)
        r1n = r1 + 0.1*sqrt(-2*log(rnd(1)))*cos(2*PI*rnd(2))* & 	
            (/ cos(2*PI*rnd(3)), sin(2*PI*rnd(3)) /)
        r2n = r2 + 0.1*sqrt(-2*log(rnd(4)))*cos(2*PI*rnd(5))* & 	
            (/ cos(2*PI*rnd(6)), sin(2*PI*rnd(6)) /)            
        
        p1 = (exp(-sum((r1n-s)**2))/exp(-sum((r1-s)**2)))**2
        p2 = (exp(-sum((r2n+s)**2))/exp(-sum((r2+s)**2)))**2
        p = p1*p2
        
        IF (p > rnd(7)) THEN
            r1 = r1n
            r2 = r2n
        END IF
        
        bi_x = int(100*(r2(0)+2)/4)
        bi_y = int(100*(r2(1)+2)/4)
        IF (((bi_x >=0 ) .and. (bi_x < 100)) .and. &
            ((bi_y >=0 ) .and. (bi_y < 100)))  THEN
            bucket(bi_x,bi_y) = bucket(bi_x,bi_y) + 1
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
    CALL random_seed(PUT=dtVals((9-seedSize):8))

END SUBROUTINE init_rnd