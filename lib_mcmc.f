! Markov-chain Monte Carlo related routines
! João Faria (c) 2013

module lib_mcmc
!  Defines the following routines:
!       gengaussvect
!       creategrid
!       fdremp
!       inverserandlin
!       normlogdens
!       gausskernel
!       gauss

  use lib_random, only: seed_random_number, random_normal
  use lib_statistics, only: mean, stdev, variance, cumsum
  use lib_matrix
  use lib_utils
    
  implicit none
   
  integer,parameter :: sp = selected_real_kind(p=6, r=37)
  integer,parameter :: dp = selected_real_kind(p=15, r=307)
  
  real(dp), parameter :: pi = 3.1415926535897932384626433832795029_dp
  
  ! the most negative argument to the exp() function is log(tiny(1.0_dp)) which
  ! is around -708. We set the limit at 10% of this value to get reasonable results
  real(dp), parameter :: smallest_exp_arg = 0.1_dp * log(tiny(1.0_dp))
  
contains


!*******************************************************************************
! generation of a random vector X from a multivariate normal distribution with 
! mean zero and covariance matrix SIGMA.
! note: in order to improve the speed of the algorithms, the subroutine uses the
!       square root matrix of sigma
!*******************************************************************************
  subroutine gengaussvect(sigma, x)

    real(dp), intent(in)    :: sigma(:,:)
    real(dp), intent(out)   :: x(:)
    integer :: n, i
    
    n = size(sigma, dim=1)
    if (size(sigma, dim=2)/=n) then
     print *, 'GENGAUSSVECT : gam is not square'
    end if

    do i=1,n
     x(i) = random_normal()
    end do

    x = matmul(sigma, x)

  end subroutine gengaussvect

!*******************************************************************************
! square root of a matrix
! calculates matrix B such that B*B = A
!*******************************************************************************
! lib_matrix -> sqrt_matrix


!*******************************************************************************
! extraction of the eigenvalues and the eigenvectors of a matrix
!*******************************************************************************
! lib_matrix -> eigen

!*******************************************************************************
! variance matrix of a matrix of data
!*******************************************************************************
! lib_matrix -> variance_matrix

!*******************************************************************************
! mean of a vector
!*******************************************************************************
! lib_statistics -> mean


!*******************************************************************************
! variance of a vector
!*******************************************************************************
! lib_statistics -> variance


!*******************************************************************************
! standard deviation
!*******************************************************************************
! lib_statistics -> stdev


!*******************************************************************************
! inversion of a matrix (used LAPACK dgetri routine)
!*******************************************************************************
! lib_matrix -> inverse

!*******************************************************************************
! determinant of a matrix (numerical recipes)
!*******************************************************************************
! lib_matrix -> determinant


!*******************************************************************************
! scaled empirical cumulative distribution on a grid
!*******************************************************************************
  subroutine fdremp(vect,grille)

    real(dp), intent(inout) :: vect(:)
    real(dp), intent(in)    :: grille(:)
    real(dp), allocatable   :: vect2(:)
    integer :: n, i

    n=size(vect)
    allocate(vect2(n))
    vect2(1)=0

    do i=2,n
     vect2(i)=0.5*(vect(i-1)+vect(i))*(grille(i)-grille(i-1))
    end do

    if (sum(vect2)==0) then
     print*, 'warning!! the function is null on the whole grid'
     stop
    end if

    vect=vect2
    call cumsum(vect)
    vect=vect/vect(n)
    deallocate(vect2)

  end subroutine fdremp

!*******************************************************************************
! inversion of the empirical cdf
!*******************************************************************************
  function inverserandlin(grille, fdr)

    real, intent(in) :: grille(:), fdr(:)
    real :: inverserandlin, alea
    integer :: i,n

    n=size(fdr)
    call random_number(alea)
    i=2
    do while ((alea>fdr(i)) .and. (i<n)) 
     i=i+1
    end do

    inverserandlin = grille(i-1) + (alea-fdr(i-1))*((grille(i)-grille(i-1))/(fdr(i)-fdr(i-1)))

  end function inverserandlin

!*******************************************************************************
! log-density of the uniform distribution on [a,b]
!*******************************************************************************
  real(dp) function loguniform(x, a, b)
    
    real(dp), intent(in) :: x, a, b

    if(x<a) then
        loguniform = create_inf(-1.0_dp)
        return
    else if (x>b) then
        loguniform = create_inf(-1.0_dp)
        return
    else
        loguniform=-log(b-a)
    endif

  end function loguniform
  
!*******************************************************************************
! log-density of the jeffreys distribution on [a,b]
!*******************************************************************************
  real(dp) function logjeffreys(x, a, b)
    
    real(dp), intent(in) :: x, a, b

    if(x<a) then
        logjeffreys = create_inf(-1.0_dp)
        return
    else if (x>b) then
        logjeffreys = create_inf(-1.0_dp)
        return
    else
        logjeffreys = -log(x*log(b/a))
    endif

  end function logjeffreys

!*******************************************************************************
! log-density of the multivariate normal distribution
!*******************************************************************************
  real(dp) function normlogdens(x, mu, sigma)

    real(dp), intent(in)    :: x(:), mu(:), sigma(:, :)
    real(dp), allocatable   :: teta(:, :)

    integer :: p, n, m

    if(determinant(sigma)==0.) then
        print*, 'normlogdens: variance matrix has zero determinant. log-density set to minus infinity'
        normlogdens=-1.*(10.**30.)
        return
    endif

    n = size(sigma, dim=1)
    m = size(sigma, dim=2)
    allocate(teta(n,m))
    call inverse(sigma, teta)
    
    p = size(mu)
    normlogdens = -0.5*log(determinant(sigma))-0.5*p*log(2.*pi)-0.5*dot_product(x-mu, matmul(teta, (x-mu)))
    
    deallocate(teta)

  end function normlogdens

!*******************************************************************************
!! density estimation by a gaussian kernel
!*******************************************************************************
  real(dp) function gausskernel(x, serie)

    real(dp), intent(in) :: x, serie(:)

    real(dp) :: std, h
    integer :: n, i

    n = size(serie)

    std = stdev(serie)

    h = std/((1.0_dp*n)**0.2_dp)

    gausskernel = 0

    do i=1,n
        gausskernel = gausskernel + gauss((x-serie(i))/h, 0.0_dp, 1.0_dp)
    end do

    gausskernel = gausskernel/(n*h)

  end function gausskernel

!*******************************************************************************
! density distribution of an univariate normal distribution
!   x       : variable
!   mean    : mean of the distribution
!   std     : standard deviation of the distribution
!*******************************************************************************
  real(dp) function gauss(x, mean, std)

    real(dp), intent(in) :: x, mean, std

    if(std==0.) then
     print*, 'gauss: variance is zero. density set to zero'
     gauss = 0.
     return
    endif

    gauss = (1.0_dp/ (sqrt(2.0_dp*pi)*std) )*exp( -0.5_dp * (((x-mean)/std)**2) )

  end function gauss


!*******************************************************************************
! execute niter iterations of the metropolis algorithm (n-dimensional). The jump
! function is multivariate normal with covariance matrix c*sigma
!
!   f:      natural log of the target function
!   init:   starting vector (in)
!   sigma:  starting covariance matrix of the gaussian jump function (in)
!           final covariance matrix of the jump function (out)
!   niter:  number of iterations
!   nc:     the multiplicative constant c will be adjusted every nc iterations,
!           to preserve an adequate acceptation rate
!   ncov:   the covariance matrix sigma will be updated every ncov iterations
!   nsave:  output current sample every nsave iterations
!   out:    sample of simulated values
!   file_unit: optionally, write (all!) samples to this file unit
!*******************************************************************************
  subroutine metropolis(f, init, profil, sigma, niter, &
                        nc, ncov, nsave, out, file_unit)

    interface 
        real(kind=8) function f(teta)
            real(kind=8), intent(in):: teta(:)
        end function
    end interface

    real(dp), intent(in)    :: init(:)
    real(dp), intent(inout) :: sigma(:,:)
    real(dp), intent(out)   :: out(:,:)
    
    integer, intent(in)     :: niter, nc, ncov, nsave, profil(:)
    integer, intent(in), optional :: file_unit


    logical :: output = .false.
    
    real(dp), allocatable   :: gen(:), sigma_new(:,:)
    real(dp), allocatable   :: parnew(:), genvect(:)
    integer, allocatable    :: varpart(:), ctpart(:)
    real(dp)    :: fact_rejet, r, random, ratio
    integer     :: n, dim, count1, count2, rejet, i

    ! number of updates
    print*, 'metropolis: multiplicative constant will be updated ', floor((1.*niter)/(1.*nc)), ' times'
    print*, '            covariance matrix will be updated ', floor((1.*niter)/(1.*ncov)), ' times'

    ! should output to file?
    if (present(file_unit)) output = .true.

    ! determine variable and fixed components of the parameters vector
    n = size(init)
    dim = count(profil==1)
    allocate(varpart(dim))
    
    if(dim /= n) then
        allocate(ctpart(n-dim))
    else
        allocate(ctpart(n))
    endif
    
    count1=0
    count2=0
    
    do i=1, size(profil)
        if(profil(i)==1) then
            count1 = count1 + 1
            varpart(count1) = i
        else
            count2 = count2 + 1
            ctpart(count2) = i
        endif
    enddo

    ! initialization 
    call seed_random_number()   ! seed RNG
    
    allocate(sigma_new(n, n))
    call sqrt_matrix(sigma, sigma_new)
    sigma = sigma_new
    
    allocate(gen(n))
    gen = 0.0_dp

    allocate(genvect(n))
    allocate(parnew(n))
    parnew = init
    
    out(1,:) = init
    fact_rejet = 2.4/sqrt(real(dim))
    rejet = 0
    
    ! information
    do i=1,dim
        write(*,*) sigma(i,i)
    end do
    write(6, '(8f15.6)') init
    write(*,'(a)',advance='no') 'press Enter to start: '
    read(*,*)
    
    
    ! start the iterations
    main: do i=1,niter
  
        ! generation of the jump on variable components
        call gengaussvect(fact_rejet*sigma, gen)
        if(dim/=n) genvect(ctpart) = profil(ctpart)
        genvect(varpart) = gen


        ! computation of the ratio of target functions
        ! be careful!! f is the natural logarithm of the target function
        ratio = f(parnew+genvect) - f(parnew)
        if(ratio < smallest_exp_arg) then
            r = 0.0_dp
        else if(ratio > 0) then
            r = 1.0_dp
        else
            r = exp(ratio)
        endif
  
        ! accept or reject the new value  
        call random_number(random)
        if (random < r ) then
            parnew = parnew + genvect
            !write(*,*) 'accepting'
        else
            rejet = rejet + 1
        endif
        out(i+1,:) = parnew
        
        
!        write(6, '(8f15.6)') parnew
!        write(6, '(8f15.6)') genvect
!        write(*,'(" 1 ")') !BP
!        read(*,*)
        
        !output to file
        if (output) write(file_unit, '(8f15.6)') parnew
      
        write(*,*) i, rejet, (i-rejet)/real(i)

        ! updates of the multiplicative constant
        if (mod(i,nc)==0) then
            if ((1. - real(rejet)/real(nc)) < 0.23) then
                fact_rejet = fact_rejet *0.9
            else if ((1. - real(rejet)/real(nc)) > 0.44) then
                fact_rejet = fact_rejet * 1.1
            end if
            !write(*,*) fact_rejet, 1. - real(rejet)/real(nc)
            !rejet=0
        end if

        ! updates of the covariance matrix
        if (mod(i,ncov)==0) then
            call variance_matrix(out(i-ncov+1:i, varpart), sigma_new)
            if (.not.(all(sigma_new==0.))) then
                sigma=sigma_new
                call sqrt_matrix(sigma, sigma_new)
                sigma=sigma_new
            endif
        end if

    end do main
    
    deallocate(sigma_new)
    deallocate(gen)
    deallocate(genvect)
    deallocate(parnew)
    deallocate(varpart, ctpart)

  end subroutine metropolis


end module lib_mcmc
