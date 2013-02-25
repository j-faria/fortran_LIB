! NaN, Inf, etc routines
! João Faria 2013

module lib_utils

  implicit none
  save

  private

  integer,parameter :: sp = selected_real_kind(p=6,r=37)
  integer,parameter :: dp = selected_real_kind(p=15,r=307)
  
  public :: create_nan
  public :: create_inf
  public :: isinf
  
contains
  
  real(dp) function create_nan()
  ! create an IEEE Not A Number
    real(dp) :: x = -1.0d0
    create_nan = sqrt(x)
  end function create_nan
  
  real(dp) function create_inf(x)
  ! create +/-Infinity
    real(dp), intent(in) :: x
    if (x > 0) create_inf = huge(x) + huge(x)
    if (x < 0) create_inf = - huge(x) - huge(x)
  end function create_inf
  
  logical function isinf(x)
    real(dp), intent(in) :: x
    isinf = (2*x==x .and. x /= 0)
  end function isinf

end module lib_utils
