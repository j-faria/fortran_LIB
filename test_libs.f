program test_libs

    use lib_assert
    use lib_periodogram, only: bls
    use lib_array, only: linspace
    use lib_plot
    use lib_statistics, only: mean
    use lib_matrix
    use lib_utils
    
    implicit none

    integer,parameter   :: dp = selected_real_kind(p=15,r=307)
    character(len=20)   :: which_test
    logical             :: timing
    integer             :: argcount
    
    ! for timing measurements
    integer :: t1,t2
    integer :: clock_rate, clock_max
    
    
    ! to test lib_matrix
    real(dp), dimension(3,3) :: m, new_m
    real(dp), dimension(5,5) :: other_m, new_other_m
    real(dp), dimension(3,3) :: vec
    real(dp), dimension(3)   :: val
    real(dp), dimension(3,3) :: l, u, p
    real(dp) :: det

    ! to test lib_utils
    real(dp) :: a    
    
    
    real(dp), dimension(100)    :: x, y
    real(dp), dimension(1000)   :: f, pdgm, period
    real(dp)                    :: ofac, hifac, prob
    integer                     :: jmax, nf

    integer :: i, j
    
    
    !--- process command line argument to determine which test to do
    argcount = iargc()
    if (argcount < 1) then
        write(*,'(a,x,/)') "ERROR! usage: test which_library"
        stop
    endif
    call getarg(1, which_test)
    timing = .false.
    if (argcount == 2) timing = .true.

    ! do tests
    select case (which_test)
    
        case('matrix')
                
            m(:,1) = (/1,4,7/)
            m(:,2) = (/2,5,8/)
            m(:,3) = (/3,6,9/)
            !write(*,'(3f5.1)') m
            call centered(m, new_m)
            call assert(new_m(1,1)==-1.0, 'error in CENTERED')
            call assert(new_m(3,3)==1.0, 'error in CENTERED')
            
            m(:,1) = (/1,3,3/)
            m(:,2) = (/1,4,3/)
            m(:,3) = (/1,3,4/)
            !write(*,'(3f5.1)') m
            call inverse(m, new_m)
            !write(*,'(3f23.17)') new_m
            call assert_almost_equal(new_m(1,3), -1.0_dp, 'error in INVERSE1', eps=2*epsilon(1.0_dp))
            call assert_almost_equal(new_m(3,3), 1.0_dp, 'error in INVERSE2')
            call assert_almost_equal(new_m(1,1), 7.0_dp, 'error in INVERSE3', eps=2*epsilon(1.0_dp))
            
            m(:,1) = (/6,1,1/)
            m(:,2) = (/4,-2,5/)
            m(:,3) = (/2,8,7/)
            det = determinant(m)
            call assert_almost_equal(det, -306.0_dp, 'error in DETERMINANT')
                
                
            m(:,1) = (/1,1,0/)
            m(:,2) = (/0,2,0/)
            m(:,3) = (/0,-1,4/)
            !write(*,'(3f6.2)') m
            call eigen(m, val, vec)
            !write(*,*) val
            call assert_almost_equal(val(3), 4.0_dp, 'error in EIGEN')

                
            !write(*,'(3f6.2)') m
            call sqrt_matrix(m, new_m)
            !write(*,'(3f15.8)') new_m
            call assert_almost_equal(new_m(1,1), sqrt(2.0_dp), 'error in SQRT_MATRIX')
            
            m = identity(3)
            !write(*,'(3f6.2)') m
            call assert_almost_equal(m(2,2), 1.0_dp, 'error in IDENTITY')

                
                
                
                
                
                
!!                m(:,1) = (/1,3,5/)
!!                m(:,2) = (/2,4,7/)
!!                m(:,3) = (/1,1,0/)
!!                write(*,'(3f5.1)') m
!!                call lu_decomp(m, new_m)
!!                write(*,*)
!!                write(*,'(3f5.1)') new_m
!!                write(*,*)
!!                
!!                l = 0
!!                u = 0
!!                do i=1,3
!!                    do j=1,3
!!                        if (i>=j) then
!!                            u(i,j) = new_m(i,j)
!!                            if(i==j) then
!!                                l(i,j) = 1
!!                                p(i,j) = 3
!!                            end if
!!                        else 
!!                            l(i,j) = new_m(i,j)
!!                        end if
!!                    end do
!!                end do
!!!                write(*,'(3f5.1)') u
!!!                write(*,*)
!!!                write(*,'(3f5.1)') l
!!                
!!!                other_m(:,1) =  (/1.0 , 1.2 , 1.4 , 1.6 , 1.8 , 2.0 , 2.2 , 2.4 , 2.6 /)
!!!                other_m(:,2) =  (/ 1.2 , 1.0 , 1.2 , 1.4 , 1.6 , 1.8 , 2.0 , 2.2 , 2.4 /)
!!!                other_m(:,3) =  (/ 1.4 , 1.2 , 1.0 , 1.2 , 1.4 , 1.6 , 1.8 , 2.0 , 2.2 /)
!!!                other_m(:,4) =  (/ 1.6 , 1.4 , 1.2 , 1.0 , 1.2 , 1.4 , 1.6 , 1.8 , 2.0 /)
!!!                other_m(:,5) =  (/ 1.8 , 1.6 , 1.4 , 1.2 , 1.0 , 1.2 , 1.4 , 1.6 , 1.8 /)
!!!                other_m(:,6) =  (/ 2.0 , 1.8 , 1.6 , 1.4 , 1.2 , 1.0 , 1.2 , 1.4 , 1.6 /)
!!!                other_m(:,7) =  (/ 2.2 , 2.0 , 1.8 , 1.6 , 1.4 , 1.2 , 1.0 , 1.2 , 1.4 /)
!!!                other_m(:,8) =  (/ 2.4 , 2.2 , 2.0 , 1.8 , 1.6 , 1.4 , 1.2 , 1.0 , 1.2 /)
!!!                other_m(:,9) =  (/ 2.6 , 2.4 , 2.2 , 2.0 , 1.8 , 1.6 , 1.4 , 1.2 , 1.0 /)
!!!                write(*,'(9f5.1)') other_m
!!!                call lu_decomp(other_m, new_other_m)
!!!                write(*,*)
!!!                write(*,'(9f5.1)') new_other_m
!!                
!!                write(*,'(3f5.1)') matmul(p, matmul(l, u))
!!                
                write(*,*) 'tested lib_matrix without errors'
        
        
        case('utils')
        
            a = create_nan()
            call assert(isnan(a), 'error in CREATE_NAN')
            a = create_inf(-1.0_dp)
            call assert(isinf(a), 'error in CREATE_INF')
            a = create_inf(1.0_dp)
            call assert(isinf(a), 'error in CREATE_INF')
            
            write(*,*) 'tested lib_utils without errors'
            
            
                
    end select
    
    
    
!    ! build signal
!    call linspace(1.0_dp, 10.0_dp, x)
!    y = sin(x)
!    
!    ! calculate periodogram
!    ofac = 6.0
!    call bls(x, y, ofac, hifac, f, pdgm, nf, jmax, prob)
!    
!    period = 1.0/f
!    call plot(period(1:nf), pdgm(1:nf), ' 5-')
!!    call plot(x, y, ' 5-')


    

    
    
    

end program test_libs
