module lib_simplex

  implicit none
  save

  private
  
  integer,parameter :: sp = selected_real_kind(p=6,r=37)
  integer,parameter :: dp = selected_real_kind(p=15,r=307)

  public :: minim

contains

  subroutine minim(p, nop, func, maxfn, iprint, stopcr, nloop, iquad,  &
                 simp, var, functn, ifault)

!     a program for function minimization using the simplex (Nelder-Mead) method
!     for details, see nelder & mead, the computer journal, january 1965

!     !-- change history at bottom of file --!

!     arguments:-
!     p()     = input, starting values of parameters
!               output, final values of parameters
!     step()  = input, initial step sizes
!     nop     = input, no. of parameters, incl. any to be held fixed
!     func    = output, the function value corresponding to the final
!                 parameter values.
!     maxfn     = input, the maximum no. of function evaluations allowed.
!               say, 20 times the number of parameters, nop.
!     iprint  = input, print control parameter
!                 < 0 no printing
!                 = 0 printing of parameter values and the function
!                     value after initial evidence of convergence.
!                 > 0 as for iprint = 0 plus progress reports after every
!                     iprint evaluations, plus printing for the initial simplex.
!     stopcr  = input, stopping criterion.
!               the criterion is applied to the standard deviation of
!               the values of func at the points of the simplex.
!     nloop   = input, the stopping rule is applied after every nloop
!               function evaluations.   normally nloop should be slightly
!               greater than nop, say nloop = 2*nop.
!     iquad   = input, = 1 if fitting of a quadratic surface is required
!                      = 0 if not
!               n.b. the fitting of a quadratic surface is strongly
!               recommended, provided that the fitted function is
!               continuous in the vicinity of the minimum.   it is often
!               a good indicator of whether a premature termination of
!               the search has occurred.
!     simp    = input, criterion for expanding the simplex to overcome
!               rounding errors before fitting the quadratic surface.
!               the simplex is expanded so that the function values at
!               the points of the simplex exceed those at the supposed
!               minimum by at least an amount simp.
!     var()   = output, contains the diagonal elements of the inverse of
!               the information matrix.
!     functn  = input, name of the user's subroutine - arguments (p,func)
!               which returns the function value for a given set of
!               parameter values in array p.
!               NOTE: functn must be declared external in the calling program
!     ifault  = output, = 0 for successful termination
!                 = 1 if maximum no. of function evaluations exceeded
!                 = 2 if information matrix is not +ve semi-definite
!                 = 3 if nop < 1
!                 = 4 if nloop < 1
!     step()  = input, optional, step sizes for building initial simplex
!               by default, the code uses a step of 0.1*p
!     varypar = input, optional, array of logicals that determine which 
!               parameters are varied. By default, varypar=.true.
!
!     n.b. p, step and var (if iquad = 1) must have dimension at least nop
!          in the calling program.

!*****************************************************************************

    integer, intent(in)                 :: nop, maxfn, iprint, nloop, iquad
    integer, intent(out)                :: ifault
    real(dp), intent(in)                :: stopcr, simp
    real(dp), intent(inout)             :: p(:)
    real(dp), intent(out)               :: var(:), func
    real(dp)   :: step(nop)
    logical    :: varypar(nop) 

    interface
      subroutine functn(p, func)
        implicit none
        integer,parameter      :: dp = selected_real_kind(p=15,r=307)
        real (dp), intent(in)  :: p(:)
        real (dp), intent(out) :: func
      end subroutine functn
    end interface

!     local variables

    real (dp)   :: g(nop+1,nop), h(nop+1), pbar(nop), pstar(nop), pstst(nop), &
               aval(nop), pmin(nop), temp(nop), bmat(nop*(nop+1)/2),  &
               vc(nop*(nop+1)/2), ymin, rmax, hstst, a0, hmin, test, hmean, &
               hstd, hstar, hmax, savemn

    real (dp), parameter :: zero = 0._dp, half = 0.5_dp, one = 1._dp, two = 2._dp
    integer     :: i, i1, i2, iflag, ii, ij, imax, imin, irank, irow, j, j1, jj, &
               k, l, loop, nap, neval, nmore, np1, nullty
    integer :: error

!     a = reflection coefficient, b = contraction coefficient, and
!     c = expansion coefficient.

    real (dp), parameter :: a = 1._dp, b = 0.5_dp, c = 2._dp

!     set lout = logical unit no. for output

    integer, parameter :: lout = 6

!     if progress reports have been requested, print heading

    if (iprint > 0) write (lout,5000) iprint

!     check input arguments

    ifault = 0
    if (nop <= 0) ifault = 3
    if (nloop <= 0) ifault = 4
    if (ifault /= 0) return

!     handle optional inputs
    step = 1._dp + 0.1_dp * p
    varypar = .true.

!     set nap = no. of parameters to be varied

    nap = count(varypar)
    neval = 0
    loop = 0
    iflag = 0

!     if nap = 0 evaluate function at the starting point and return

    if (nap <= 0) then
      call functn(p,func)
      return
    end if

!     set up the initial simplex

20  g(1,:) = p
    irow = 2
    do i = 1, nop
      if (varypar(i)) then
        g(irow,:) = p
        g(irow,i) = p(i) + step(i)
        irow = irow + 1
      end if
    end do

    np1 = nap + 1
    do i = 1, np1
      p = g(i,:)
      call functn(p,h(i))
      neval = neval + 1
      if (iprint > 0) then
        write (lout,5100) neval, h(i), p
      end if
    end do

!     start of main cycle.

!     find max. & min. values for current simplex (hmax & hmin).

    main_loop: do
      loop = loop + 1
      imax = 1
      imin = 1
      hmax = h(1)
      hmin = h(1)
      do i = 2, np1
        if (h(i) > hmax) then
          imax = i
          hmax = h(i)
        else
          if (h(i) < hmin) then
            imin = i
            hmin = h(i)
          end if
        end if
      end do

!     find the centroid of the vertices other than p(imax)

      pbar = zero
      do i = 1, np1
        if (i /= imax) then
          pbar = pbar + g(i,:)
        end if
      end do
      pbar = pbar / nap

!     reflect maximum through pbar to pstar,
!     hstar = function value at pstar.

      pstar = a * (pbar - g(imax,:)) + pbar
      call functn(pstar,hstar)
      neval = neval + 1
      if (iprint > 0) then
        if (mod(neval,iprint) == 0) write (lout,5100) neval, hstar, pstar
      end if

!     if hstar < hmin, reflect pbar through pstar,
!     hstst = function value at pstst.

      if (hstar < hmin) then
        pstst = c * (pstar - pbar) + pbar
        call functn(pstst,hstst)
        neval = neval + 1
        if (iprint > 0) then
          if (mod(neval,iprint) == 0) write (lout,5100) neval, hstst, pstst
        end if

!     if hstst < hmin replace current maximum point by pstst and
!     hmax by hstst, then test for convergence.

        if (hstst >= hmin) then   ! replace maximum point by pstar & h(imax) by hstar.
          g(imax,:) = pstar
          h(imax) = hstar
        else
          g(imax,:) = pstst
          h(imax) = hstst
        end if
        go to 250
      end if

!     hstar is not < hmin.
!     test whether it is < function value at some point other than
!     p(imax).   if it is replace p(imax) by pstar & hmax by hstar.

      do i = 1, np1
        if (i /= imax) then
          if (hstar < h(i)) then  ! replace maximum point by pstar & h(imax) by hstar.
            g(imax,:) = pstar
            h(imax) = hstar
            go to 250
          end if
        end if
      end do

!     hstar > all function values except possibly hmax.
!     if hstar <= hmax, replace p(imax) by pstar & hmax by hstar.

      if (hstar <= hmax) then
        g(imax,:) = pstar
        hmax = hstar
        h(imax) = hstar
      end if

!     contracted step to the point pstst,
!     hstst = function value at pstst.

      pstst = b * g(imax,:) + (one-b) * pbar
      call functn(pstst,hstst)
      neval = neval + 1
      if (iprint > 0) then
        if (mod(neval,iprint) == 0) write (lout,5100) neval, hstst, pstst
      end if

!     if hstst < hmax replace p(imax) by pstst & hmax by hstst.

      if (hstst <= hmax) then
        g(imax,:) = pstst
        h(imax) = hstst
        go to 250
      end if

!     hstst > hmax.
!     shrink the simplex by replacing each point, other than the current
!     minimum, by a point mid-way between its current position and the
!     minimum.

      do i = 1, np1
        if (i /= imin) then
          do j = 1, nop
            if (varypar(j)) g(i,j) = (g(i,j) + g(imin,j)) * half
            p(j) = g(i,j)
          end do
          call functn(p,h(i))
          neval = neval + 1
          if (iprint > 0) then
            if (mod(neval,iprint) == 0) write (lout,5100) neval, h(i), p
          end if
        end if
      end do

!     if loop = nloop test for convergence, otherwise repeat main cycle.

  250 if (loop < nloop) cycle main_loop

!     calculate mean & standard deviation of function values for the
!     current simplex.

      hmean = sum( h(1:np1) ) / np1
      hstd = sum( (h(1:np1) - hmean) ** 2 )
      hstd = sqrt(hstd / np1)

!     if the rms > stopcr, set iflag & loop to zero and go to the
!     start of the main cycle again.

      if (hstd > stopcr .and. neval <= maxfn) then
        iflag = 0
        loop = 0
        cycle main_loop
      end if

!     find the centroid of the current simplex and the function value there.

      do i = 1, nop
        if (varypar(i)) then
          p(i) = sum( g(1:np1,i) ) / np1
        end if
      end do
      call functn(p,func)
      neval = neval + 1
      if (iprint > 0) then
        if (mod(neval,iprint) == 0) write (lout,5100) neval, func, p
      end if

!     test whether the no. of function values allowed, maxfn, has been
!     overrun; if so, exit with ifault = 1.

      if (neval > maxfn) then
        ifault = 1
        if (iprint < 0) return
        write (lout,5200) maxfn
        write (lout,5300) hstd
        write (lout,5400) p
        write (lout,5500) func
        return
      end if

!     convergence criterion satisfied.
!     if iflag = 0, set iflag & save hmean.
!     if iflag = 1 & change in hmean <= stopcr then search is complete.

      if (iprint >= 0) then
        write (lout,5600)
        write (lout,5400) p
        write (lout,5500) func
      end if

      if (iflag == 0 .or. abs(savemn-hmean) >= stopcr) then
        iflag = 1
        savemn = hmean
        loop = 0
      else
        exit main_loop
      end if

    end do main_loop

    if (iprint >= 0) then
      write (lout,5700) neval
      write (lout,5800) p
      write (lout,5900) func
    end if
    if (iquad <= 0) return

!------------------------------------------------------------------

!     quadratic surface fitting

    if (iprint >= 0) write (lout,6000)

!     expand the final simplex, if necessary, to overcome rounding
!     errors.

    hmin = func
    nmore = 0
    do i = 1, np1
      do
        test = abs(h(i)-func)
        if (test < simp) then
          do j = 1, nop
            if (varypar(j)) g(i,j) = (g(i,j)-p(j)) + g(i,j)
            pstst(j) = g(i,j)
          end do
          call functn(pstst,h(i))
          nmore = nmore + 1
          neval = neval + 1
          if (h(i) >= hmin) cycle
          hmin = h(i)
          if (iprint >= 0) write (lout,5100) neval, hmin, pstst
        else
          exit
        end if
      end do
    end do

!     function values are calculated at an additional nap points.

    do i = 1, nap
      i1 = i + 1
      pstar = (g(1,:) + g(i1,:)) * half
      call functn(pstar,aval(i))
      nmore = nmore + 1
      neval = neval + 1
    end do

!     the matrix of estimated second derivatives is calculated and its
!     lower triangle stored in bmat.

    a0 = h(1)
    do i = 1, nap
      i1 = i - 1
      i2 = i + 1
      do j = 1, i1
        j1 = j + 1
        pstst = (g(i2,:) + g(j1,:)) * half
        call functn(pstst,hstst)
        nmore = nmore + 1
        neval = neval + 1
        l = i * (i-1) / 2 + j
        bmat(l) = two * (hstst + a0 - aval(i) - aval(j))
      end do
    end do

    l = 0
    do i = 1, nap
      i1 = i + 1
      l = l + i
      bmat(l) = two * (h(i1) + a0 - two*aval(i))
    end do

!     the vector of estimated first derivatives is calculated and
!     stored in aval.

    do i = 1, nap
      i1 = i + 1
      aval(i) = two * aval(i) - (h(i1) + 3._dp*a0) * half
    end do

!     the matrix q of nelder & mead is calculated and stored in g.

    pmin = g(1,:)
    do i = 1, nap
      i1 = i + 1
      g(i1,:) = g(i1,:) - g(1,:)
    end do

    do i = 1, nap
      i1 = i + 1
      g(i,:) = g(i1,:)
    end do

!     invert bmat

    call syminv(bmat, nap, bmat, temp, nullty, ifault, rmax)
    if (ifault == 0) then
      irank = nap - nullty
    else                                 ! bmat not +ve definite
                                         ! resume search for the minimum
      if (iprint >= 0) write (lout,6100)
      ifault = 2
      if (neval > maxfn) return
      if (iprint >= 0) write (lout,6200)
      step = half * step
      go to 20
    end if

!     bmat*a/2 is calculated and stored in h.

    do i = 1, nap
      h(i) = zero
      do j = 1, nap
        if (j <= i) then
          l = i * (i-1) / 2 + j
        else
          l = j * (j-1) / 2 + i
        end if
        h(i) = h(i) + bmat(l) * aval(j)
      end do
    end do

!     find the position, pmin, & value, ymin, of the minimum of the
!     quadratic.

    ymin = dot_product( h(1:nap), aval(1:nap) )
    ymin = a0 - ymin
    do i = 1, nop
      pstst(i) = dot_product( h(1:nap), g(1:nap,i) )
    end do
    pmin = pmin - pstst
    if (iprint >= 0) then
      write (lout,6300) ymin, pmin
      write (lout,6400)
    end if

!     q*bmat*q'/2 is calculated & its lower triangle stored in vc

    do i = 1, nop
      do j = 1, nap
        h(j) = zero
        do k = 1, nap
          if (k <= j) then
            l = j * (j-1) / 2 + k
          else
            l = k * (k-1) / 2 + j
          end if
          h(j) = h(j) + bmat(l) * g(k,i) * half
        end do
      end do

      do j = i, nop
        l = j * (j-1) / 2 + i
        vc(l) = dot_product( h(1:nap), g(1:nap,j) )
      end do
    end do

!     the diagonal elements of vc are copied into var.

    j = 0
    do i = 1, nop
      j = j + i
      var(i) = vc(j)
    end do
    if (iprint < 0) return
    write (lout,6500) irank
    call print_tri_matrix(vc, nop, lout)

    write (lout,6600)
    call syminv(vc, nap, bmat, temp, nullty, ifault, rmax)

!     bmat now contains the information matrix

    write (lout,6700)
    call print_tri_matrix(bmat, nop, lout)

    ii = 0
    ij = 0
    do i = 1, nop
      ii = ii + i
      if (vc(ii) > zero) then
        vc(ii) = one / sqrt(vc(ii))
      else
        vc(ii) = zero
      end if
      jj = 0
      do j = 1, i - 1
        jj = jj + j
        ij = ij + 1
        vc(ij) = vc(ij) * vc(ii) * vc(jj)
      end do
      ij = ij + 1
    end do

    write (lout,6800)
    ii = 0
    do i = 1, nop
      ii = ii + i
      if (vc(ii) /= zero) vc(ii) = one
    end do
    call print_tri_matrix(vc, nop, lout)

!     exit, on successful termination.

    write (lout,6900) nmore
    return

5000 format (' progress report every',i4,' function evaluations'/  &
             ' eval.   func.value.          parameter values')
5100 format (/' ', i4, '  ', g12.5, '  ', 5g11.4, 3(/t22, 5g11.4))
5200 format (' no. of function evaluations > ',i5)
5300 format (' rms of function values of last simplex =', g14.6)
5400 format (' centroid of last simplex =',4(/' ', 6g13.5))
5500 format (' function value at centroid =', g14.6)
5600 format (/' evidence of convergence')
5700 format (/' minimum found after',i5,' function evaluations')
5800 format (' minimum at',4(/' ', 6g13.6))
5900 format (' function value at minimum =', g14.6)
6000 format (/' fitting quadratic surface about supposed minimum'/)
6100 format (/' matrix of estimated second derivatives not +ve defn.'/ &
             ' minimum probably not found'/)
6200 format (/t11, 'search restarting'/)
6300 format (' minimum of quadratic surface =',g14.6,' at',4(/' ', 6g13.5))
6400 format (' if this differs by much from the minimum estimated ',      &
             'from the minimization,'/' the minimum may be false &/or the ',  &
             'information matrix may be inaccurate'/)
6500 format (' rank of information matrix =',i3/ &
             ' inverse of information matrix:-')
6600 format (/' if the function minimized was -log(likelihood),'/         &
             ' this is the covariance matrix of the parameters.'/         &
             ' if the function was a sum of squares of residuals,'/       &
             ' this matrix must be multiplied by twice the estimated ',   &
             'residual variance'/' to obtain the covariance matrix.'/)
6700 format (' information matrix:-'/)
6800 format (/' correlation matrix:-')
6900 format (/' a further',i4,' function evaluations have been used'/)

  end subroutine minim




subroutine syminv(a, n, c, w, nullty, ifault, rmax)

!     algorithm as7, applied statistics, vol.17, 1968.

!     arguments:-
!     a()    = input, the symmetric matrix to be inverted, stored in
!                lower triangular form
!     n      = input, order of the matrix
!     c()    = output, the inverse of a (a generalized inverse if c is
!                singular), also stored in lower triangular form.
!                c and a may occupy the same locations.
!     w()    = workspace, dimension at least n.
!     nullty = output, the rank deficiency of a.
!     ifault = output, error indicator
!                 = 1 if n < 1
!                 = 2 if a is not +ve semi-definite
!                 = 0 otherwise
!     rmax   = output, approximate bound on the accuracy of the diagonal
!                elements of c.  e.g. if rmax = 1.e-04 then the diagonal
!                elements of c will be accurate to about 4 dec. digits.

!     latest revision - 1 april 1985

!***************************************************************************

real (dp), intent(in out) :: a(:), c(:), w(:)
integer, intent(in)       :: n
integer, intent(out)      :: nullty, ifault
real (dp), intent(out)    :: rmax

real (dp), parameter :: zero = 0._dp, one = 1._dp
integer              :: i, icol, irow, j, jcol, k, l, mdiag, ndiag, nn, nrow
real (dp)            :: x

nrow = n
ifault = 1
if (nrow > 0) then
  ifault = 0

!     cholesky factorization of a, result in c

  call chola(a, nrow, c, nullty, ifault, rmax, w)
  if (ifault == 0) then

!     invert c & form the product (cinv)'*cinv, where cinv is the inverse
!     of c, row by row starting with the last row.
!     irow = the row number, ndiag = location of last element in the row.

    nn = nrow * (nrow+1) / 2
    irow = nrow
    ndiag = nn
    10 if (c(ndiag) /= zero) then
      l = ndiag
      do i = irow, nrow
        w(i) = c(l)
        l = l + i
      end do
      icol = nrow
      jcol = nn
      mdiag = nn

      30 l = jcol
      x = zero
      if (icol == irow) x = one / w(irow)
      k = nrow
      40 if (k /= irow) then
        x = x - w(k) * c(l)
        k = k - 1
        l = l - 1
        if (l > mdiag) l = l - k + 1
        go to 40
      end if

      c(l) = x / w(irow)
      if (icol == irow) go to 60
      mdiag = mdiag - icol
      icol = icol - 1
      jcol = jcol - 1
      go to 30
    end if ! (c(ndiag) /= zero)

    l = ndiag
    do j = irow, nrow
      c(l) = zero
      l = l + j
    end do

    60 ndiag = ndiag - irow
    irow = irow - 1
    if (irow /= 0) go to 10
  end if
end if
return

end subroutine syminv



subroutine chola(a, n, u, nullty, ifault, rmax, r)

!     algorithm as6, applied statistics, vol.17, 1968, with
!     modifications by a.j.miller

!     arguments:-
!     a()    = input, a +ve definite matrix stored in lower-triangular
!                form.
!     n      = input, the order of a
!     u()    = output, a lower triangular matrix such that u*u' = a.
!                a & u may occupy the same locations.
!     nullty = output, the rank deficiency of a.
!     ifault = output, error indicator
!                 = 1 if n < 1
!                 = 2 if a is not +ve semi-definite
!                 = 0 otherwise
!     rmax   = output, an estimate of the relative accuracy of the
!                diagonal elements of u.
!     r()    = output, array containing bounds on the relative accuracy
!                of each diagonal element of u.

!     latest revision - 1 april 1985

!***************************************************************************

real (dp), intent(in)   :: a(:)
integer, intent(in)     :: n
real (dp), intent(out)  :: u(:)
integer, intent(out)    :: nullty, ifault
real (dp), intent(out)  :: rmax, r(:)

!     eta should be set equal to the smallest +ve value such that
!     1._dp + eta is calculated as being greater than 1._dp in the accuracy
!     being used.

real (dp), parameter :: eta = epsilon(1.0_dp), zero = 0._dp
integer              :: i, icol, irow, j, k, l, m
real (dp)            :: rsq, w

ifault = 1
if (n > 0) then
  ifault = 2
  nullty = 0
  rmax = eta
  r(1) = eta
  j = 1
  k = 0

!     factorize column by column, icol = column no.

  do  icol = 1, n
    l = 0

!     irow = row number within column icol

    do  irow = 1, icol
      k = k + 1
      w = a(k)
      if (irow == icol) rsq = (w*eta) ** 2
      m = j
      do  i = 1, irow
        l = l + 1
        if (i == irow) exit
        w = w - u(l) * u(m)
        if (irow == icol) rsq = rsq + (u(l)**2*r(i)) ** 2
        m = m + 1
      end do

      if (irow == icol) exit
      if (u(l) /= zero) then
        u(k) = w / u(l)
      else
        u(k) = zero
        if (abs(w) > abs(rmax*a(k))) go to 60
      end if
    end do

!     end of row, estimate relative accuracy of diagonal element.

    rsq = sqrt(rsq)
    if (abs(w) > 5.*rsq) then
      if (w < zero) return
      u(k) = sqrt(w)
      r(i) = rsq / w
      if (r(i) > rmax) rmax = r(i)
    else
      u(k) = zero
      nullty = nullty + 1
    end if

    j = j + icol
  end do
  ifault = zero
end if
60 return

end subroutine chola



subroutine print_tri_matrix(a, n, lout)

integer, intent(in)    :: n, lout
real (dp), intent(in)  :: a(:)

!     local variables
integer  :: i, ii, i1, i2, l

l = 1
do l = 1, n, 6
  ii = l * (l-1) / 2
  do i = l, n
    i1 = ii + l
    ii = ii + i
    i2 = min(ii,i1+5)
    write (lout,'(1x,6g13.5)') a(i1:i2)
  end do
end do
return

end subroutine print_tri_matrix

end module lib_simplex

! CHANGE HISTORY
!     programmed by d.e.shaw,
!     csiro, division of mathematics & statistics
!     p.o. box 218, lindfield, n.s.w. 2070

!     with amendments by r.w.m.wedderburn
!     rothamsted experimental station
!     harpenden, hertfordshire, england

!     further amended by alan miller
!     csiro division of mathematical & information sciences
!     private bag 10, clayton, vic. 3169

!     fortran 90 conversion by alan miller, june 1995
!     alan.miller @ vic.cmis.csiro.au
!     latest revision - 5 december 1999

!     some amendments by joão faria, jan 2013
!     joao.faria @ astro.up.pt
