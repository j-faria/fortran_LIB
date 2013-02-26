! Matrix related routines
! João Faria 2013

module lib_matrix
!  Defines the following routines:
!   - centered matrix
!   - inverse matrix

  use lib_statistics, only: mean
  
  implicit none
  save

  private

  integer,parameter :: idp = selected_int_kind(13)
  integer,parameter :: sp = selected_real_kind(p=6,r=37)
  integer,parameter :: dp = selected_real_kind(p=15,r=307)
  
  public :: centered
  public :: inverse
  public :: determinant
  public :: lu_decomp
  public :: eigen
  public :: sqrt_matrix
  public :: variance_matrix
  
contains

!*******************************************************************************
! compute the centered matrix, ie. the matrix elements minus the column means
!*******************************************************************************
  subroutine centered(matrix, centered_matrix)

    real(dp), intent(in)    :: matrix(:,:)
    real(dp), intent(out)   :: centered_matrix(:, :)
    real(dp), allocatable   :: means(:)
    integer :: m, n, i, j

    m = size(matrix, dim=1)
    n = size(matrix, dim=2)
    allocate(means(m))

    do i=1,m
        means(i) = mean(matrix(i,:))
    end do
    
    do i=1,n
        centered_matrix(:, i) = matrix(:, i) - means(:)
    end do

    deallocate(means)
    return
    
  end subroutine centered
  
  
!*******************************************************************************
! inversion of a matrix (uses LAPACK dgetri and dgetrf routines)
! doesn't do proper error checking
!*******************************************************************************
  subroutine inverse(matrix, inverse_matrix)

    real(dp), intent(in)    :: matrix(:, :)
    real(dp), intent(out)   :: inverse_matrix(:, :)
    integer, allocatable    :: ipiv(:)
    real(dp), allocatable   :: work(:)
    integer :: n, lwork, lda, info
    integer :: i, j

    n = size(matrix, dim=1)
    if (size(matrix, dim=2) /= n) then
        print *, 'inverse: matrix must be square'
        stop
    end if
    
    lda = n
    lwork = n*n
    allocate(ipiv(n))
    allocate(work(lwork))

    ! store MATRIX in INVERSE_MATRIX to prevent it from being overwritten by LAPACK
    inverse_matrix = matrix
    ! DGETRF computes an LU factorization of a general M-by-N matrix, using 
    ! partial pivoting with row interchanges
    call dgetrf(n, n, inverse_matrix, lda, ipiv, info )
    ! DGETRI computes the inverse of a matrix using the LU factorization 
    ! computed by DGETRF
    call dgetri(n, inverse_matrix, n, ipiv, work, lwork, info)

    deallocate (ipiv)
    deallocate (work)
    
    return
    
  end subroutine inverse
  
  
!*******************************************************************************
! determinant of a square matrix (code by Ashwith J. Rego)
!  - the matrix is converted to upper traingular form
!  - the determinant of a triangular matrix is obtained by finding the product 
!    of the diagonal elements
!*******************************************************************************
  real(dp) function determinant(matrix)

    real(dp) :: matrix(:,:)
    real(dp) :: m, temp
    integer :: i, j, k, l, n
    logical :: detexists = .true.
    
    n = size(matrix, dim=1)
    l = 1
    !convert to upper triangular form
    do k = 1, n-1
        if (matrix(k,k) == 0) then
            detexists = .false.
            do i = k+1, n
                if (matrix(i,k) /= 0) then
                    do j = 1, n
                        temp = matrix(i,j)
                        matrix(i,j)= matrix(k,j)
                        matrix(k,j) = temp
                    end do
                    detexists = .true.
                    l=-l
                    exit
                endif
            end do
            if (.not. detexists) then
                determinant = 0
                write(*,*) "DETERMINANT: the determinant doesn't exist, set to 0"
                return
            end if
        endif
        do j = k+1, n
            m = matrix(j,k)/matrix(k,k)
            do i = k+1, n
                matrix(j,i) = matrix(j,i) - m*matrix(k,i)
            end do
        end do
    end do
   
    !calculate determinant by finding product of diagonal elements
    determinant = l
    do i = 1, n
        determinant = determinant * matrix(i,i)
    end do
   
  end function determinant
  
!*******************************************************************************
! LU decomposition of a matrix (uses LAPACK dgetrf routine)
!*******************************************************************************
  subroutine lu_decomp(matrix, lu_matrix)

    real(dp), intent(in)    :: matrix(:, :)
    real(dp), intent(out)   :: lu_matrix(:, :)
    integer, allocatable    :: ipiv(:)
    integer :: n, lda, info
    integer :: i, j

    n = size(matrix, dim=1)
    if (size(matrix, dim=2) /= n) then
        print *, 'lu_decomp: matrix must be square'
        stop
    end if
    
    lda = n
    allocate(ipiv(n))

    ! store MATRIX in LU_MATRIX to prevent it from being overwritten by LAPACK
    lu_matrix = matrix
    ! DGETRF computes an LU factorization of a general M-by-N matrix, using 
    ! partial pivoting with row interchanges
    call dgetrf(n, n, lu_matrix, lda, ipiv, info )
!    lu_matrix = transpose(lu_matrix)
    
    deallocate (ipiv)
    
    return
    
  end subroutine lu_decomp
  
!*******************************************************************************
! calculate eigenvalues and eigenvectors of a matrix
!   if any eigenvalue is complex, this routine only returns the real part and
!   issues a warning. Use EIGEN_COMPLEX instead.
!*******************************************************************************
  subroutine eigen(matrix, val, vec)
    real(dp), intent(in)    :: matrix(:,:)
    real(dp), intent(out)   :: val(:), vec(:,:)
    real(dp), allocatable   :: wr(:), wi(:), work(:), vl(:,:), vr(:,:)
    integer :: n, i, j, imag
    integer :: lda, ldvl, ldvr
    integer :: info, lwork
    character :: jobvl, jobvr
  
    n = size(matrix, dim=1)
    if (size(matrix, dim=2) /= n) then
        print *, 'EIGEN: matrix must be square'
        stop
    end if
    
    allocate(vl(n,n))
    allocate(vr(n,n))
    allocate(wi(n))
    allocate(wr(n))
    
    lda = n
    ldvl = n
    ldvr = n
    lwork = 4*n
    allocate(work(lwork))
    
    ! compute eigen(values,vectors)
    call dgeev('v', 'v', n, matrix, lda, wr, wi, vl, ldvl, vr, ldvr, work, lwork, info)

    ! check for convergence
    if( info.gt.0 ) then
        write(*,*) 'EIGEN: failed to compute eigenvalues'
        stop
    end if
    
    imag = 0
    do i = 1, n
        if( wi(i) == 0 ) then
            val(i) = wr(i)
        else
            val(i) = wr(i)
            imag = imag + 1
        end if
    end do
    if (imag /= 0) write(6,'(a,i2,x,2a)') 'EIGEN:', imag, 'eigenvalues are complex, ', &
                                           'only returning real part. Consider using EIGEN_COMPLEX'
 
!      DO I = 1, N
!         J = 1
!         DO WHILE( J.LE.N )
!            IF( WI( J ).EQ. 0 ) THEN
!               WRITE(*,9998,ADVANCE='NO') VL( I, J )
!               vec(i,j) = vl(i,j)
!               J = J + 1
!            END IF
!         END DO
!         write(*,*)
!      END DO

! 9998 FORMAT( 11(:,1X,F6.2) )

  
  end subroutine eigen


!*******************************************************************************
! Compute the principal square root X of an NxN real matrix A, such that X*X=A
!   algorithm from: Bjorck & Hammarling, A Schur method for the square root of a 
!   matrix, Lin Algebra Appl, 52-53:127--140 (1983)
!*******************************************************************************
  subroutine sqrt_matrix(matrix_A, matrix_B)

    real(dp), intent(in)    :: matrix_a(:,:)
    real(dp), intent(out)   :: matrix_b(:,:)

    real(dp), allocatable :: q(:,:), s(:,:), u(:,:), wr(:), wi(:), work(:)
    real(dp), allocatable :: quqh(:,:), qu(:,:), rwork(:)
    logical, allocatable  :: bwork(:)
    integer n, lwork, lda
    integer i, j, k, sdiag
    integer sdim, info

    n = size(matrix_a, dim=1)
    if (size(matrix_a, dim=2) /= n) then
      print *, 'sqrt_matrix: matrix must be square'
      stop
    end if

    lwork = 16*n
    lda = n
    sdim = 0
    
    ! allocate working arrays
    allocate(q(n,n))
    allocate(s(n,n))
    allocate(u(n,n))
    allocate(wr(n))
    allocate(wi(n))
    allocate(work(lwork))
    allocate(rwork(n))
    allocate(bwork(n))

    ! store A in S to prevent it from being overwritten by LAPACK
    s = matrix_a

    ! compute the schur form Q*S*(Q**h)
    call dgees('v', 'n', dummy_select, n, s, lda, sdim, wr, wi, q, n, work, lwork, rwork, bwork, info)    
    
    ! compute the square root U of S one super-diagonal at a time.
    do j = 1,n-1 ! set the lower triangle to zero
        do i = j+1,n
            u(i,j) = 0
        end do
    end do
    
    do i = 1,n ! set the diagonal elements
        u(i,i) = sqrt(s(i,i))
    enddo
    
    do sdiag = 1,n-1 ! loop over the n-1 super-diagonals
        do i = 1,n-sdiag
            j = i+sdiag
            u(i,j) = s(i,j)
            do k = i+1,j-1
                u(i,j) = u(i,j) - u(i,k)*u(k,j)
            end do
            u(i,j) = u(i,j) / (u(i,i) + u(j,j))
        end do
    end do
    
    !  compute B = Q*U*(Q**h).
    qu = matmul(q, u)
    quqh = matmul(qu, transpose(q))
    matrix_b = quqh


end subroutine sqrt_matrix
!**********************************************
    ! dummy select function
    logical function dummy_select(arg1, arg2)
        real(dp), intent(in) :: arg1, arg2
        ! the return value is always .true.
        dummy_select = .true.
    end function dummy_select
!**********************************************




!*******************************************************************************
! variance matrix of a matrix of data
!*******************************************************************************
  subroutine variance_matrix(matrix, var_matrix)

    real(dp), intent(in)    :: matrix(:,:)
    real(dp), intent(out)   :: var_matrix(:,:)
    real(dp), allocatable   :: m(:,:)

    integer :: n, nobs

    n = size(matrix, dim=2)
    nobs = size(matrix, dim=1)
    allocate(m(nobs, n))

    ! mean center the matrix
    call centered(matrix, m)

    ! calculate variance matrix
    var_matrix = matmul(transpose(m), m) * (1./real(nobs))

    deallocate(m)

  end subroutine variance_matrix


end module lib_matrix




!!*******************************************************************************
!! determinant of a matrix (numerical recipes)
!!*******************************************************************************
!real(dp) function det(matrix)
!    
!    real(dp), intent(in)    :: matrix(:,:)

!    real(dp), allocatable   :: tabbis(:,:)
!    real(dp) :: d
!    integer, allocatable    :: indx(:)
!    integer :: n, j

!    n = size(matrix(1,:))
!    
!    allocate(indx(n),tabbis(n,n))
!    
!    tabbis = matrix
!    
!    call lu_decomposition(tabbis)
!    
!    d = 1.0_dp
!    do j=1,n
!        d = d * tabbis(j,j)
!    enddo
!    
!    det = d
!    
!    deallocate(tabbis, indx)
!    
!end function
