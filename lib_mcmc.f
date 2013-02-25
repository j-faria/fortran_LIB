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
  
contains


!*******************************************************************************
! generation of a random vector X from a multivariate normal distribution with 
! mean zero and covariance matrix GAM.
! note: in order to improve the speed of the algorithms, the subroutine uses the
!       square root matrix of gamma
!*******************************************************************************
subroutine gengaussvect(gam, x)

    real(dp), intent(in)    :: gam(:,:)
    real(dp), intent(out)   :: x(:)
    integer :: n, i
    
    n = size(gam, dim=1)
    if (size(gam, dim=2)/=n) then
     print *, 'GENGAUSSVECT : gam is not square'
    end if

    do i=1,n
     x(i) = random_normal()
    end do

    x = matmul(gam, x)

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
! function is multivariate normal with covariance matrix c*gamma
!
!   f:      natural log of the target function
!   init:   starting vector (in)
!   gamma:  starting covariance matrix of the gaussian jump function (in)
!           final covariance matrix of the jump function (out)
!   niter:  number of iterations
!   nc:     the multiplicative constant c will be adjusted every nc iterations,
!           to preserve an adequate acceptation rate
!   ncov:   the covariance matrix gamma will be updated every ncov iterations
!   out:    sample of simulated values
!*******************************************************************************
subroutine metropolis(f, x, init, profil, gamma, niter, nc, ncov, out)

    interface 
      real(kind=8) function f(x, teta)
            real(kind=8), intent(in):: x(:,:), teta(:)
      end function
    end interface

    real(dp), intent(in) :: x(:,:)
    real(dp), intent(in)    :: init(:)
    real(dp), intent(inout) :: gamma(:,:)
    real(dp), intent(out)   :: out(:,:)
    
    integer, intent(in)     :: niter, nc, ncov, profil(:)
    
    real(dp), allocatable   :: gen(:), gamnew(:,:)

    real(dp), allocatable   :: parnew(:), genvect(:)
    integer, allocatable    :: varpart(:), ctpart(:)
    real(dp)    :: fact_rejet, r, random
    integer     :: n, dim, count1, count2, rejet, i

    ! number of updates
    print*, 'metropolis: multiplicative constant will be updated ', floor((1.*niter)/(1.*nc)), ' times'
    print*, '            covariance matrix will be updated ', floor((1.*niter)/(1.*ncov)), ' times'

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
    
    allocate(gamnew(n, n))
    call sqrt_matrix(gamma, gamnew)
    gamma = gamnew
    deallocate(gamnew)

    allocate(gen(n))
    gen = 0.0_dp

    allocate(genvect(n))
    allocate(parnew(n))
    parnew = init
    
    out(1,:) = init
    fact_rejet = 2.4/sqrt(real(dim))
    rejet = 0

    ! start the iterations
    main: do i=1,niter
  
      ! generation of the jump on variable components
      call gengaussvect(fact_rejet*gamma, gen)
      if(dim/=n) genvect(ctpart) = profil(ctpart)
      genvect(varpart) = gen


      ! computation of the ratio of target functions
      ! be careful!! f is the natural logarithm of the target function
      if(f(x, parnew+genvect)-f(x, parnew)<=-30.*log(10.)) then
        r = 0.0_dp
      else
        if(f(x, parnew+genvect)-f(x, parnew)>0) then
          r = 1.0_dp
        else
          r =	exp(f(x,parnew+genvect)-f(x,parnew))
        endif
      endif
  
      ! accept or reject the new value  
      call random_number(random)
      if (random < r ) then
        parnew = parnew + genvect
      else
        rejet = rejet + 1
      endif
      out(i+1,:) = parnew

      ! updates of the multiplicative constant
      if (mod(i,nc)==0) then
        if ((1. - real(rejet)/real(nc)) < 0.23) then
          fact_rejet = fact_rejet *0.9
        elseif ((1. - real(rejet)/real(nc)) > 0.44) then
          fact_rejet = fact_rejet * 1.1
        end if
        rejet=0
      end if

      ! updates of the covariance matrix
      if (mod(i,ncov)==0) then
        call variance_matrix(out(i-ncov+1:i, varpart), gamnew)
        if (.not.(all(gamnew==0.))) then
          gamma=gamnew
          call sqrt_matrix(gamma, gamnew)
          gamma=gamnew
        endif
      end if

    end do main
    
    deallocate(gen)
    deallocate(genvect)
    deallocate(parnew)
    deallocate(varpart, ctpart)

end subroutine metropolis


end module lib_mcmc
