! periodogram related routines
! joão faria 2013

module lib_periodogram

    use lib_assert
    use lib_statistics, only: mean
    use lib_array, only: linspace
    
    implicit none
    save


    integer,parameter     :: sp = selected_real_kind(p=6,r=37)
    integer,parameter     :: dp = selected_real_kind(p=15,r=307)
    real(dp), parameter   :: pi = 3.1415926535897932384626433832795029
    real(dp), parameter   :: twopi = 2.0_dp * pi
    real(dp), parameter   :: fourpi = 4.0_dp * pi  

contains

    subroutine ls(x, y, n, ofac, hifac, px, py, np, nout, jmax, prob)
! NOTE: this is Numerical Recipes "period" routine (with some f95 edits)

!   given n data points with abscissas x(1:n) (which need not be equally spaced)
!   and ordinates y(1:n), and given a desired oversampling factor ofac (a typical
!   value being 4 or larger), this routine fills array px with an increasing 
!   sequence of frequencies (not angular frequencies) up to hifac times the 
!   “average” nyquist frequency, and fills array py with the values of the 
!   lomb-scargle normalized periodogram at those frequencies. the arrays x and y
!   are not altered. np, the dimension of px and py, must be large enough to 
!   contain the output, or an error (pause) results. the routine also returns 
!   jmax such that py(jmax) is the maximum element in py, and prob, an estimate 
!   of the significance of that maximum against the hypothesis of random noise. 
!   a small value of prob indicates that a significant periodic signal is present.

	    integer i,j,jmax,n,nout,np,nmax
	    double precision hifac,ofac,prob,px(np),py(np),x(n),y(n)
	    parameter(nmax=2000)
	    double precision ave,c,cc,cwtau,effm,expy,pnow,pymax,s,ss
	    double precision sumc,sumcy,sums,sumsh,sumsy,swtau,var
	    double precision wtau,xave,xdif,xmax,xmin,yy
	    double precision arg,wtemp,wi(nmax),wpi(nmax)
	    double precision wpr(nmax),wr(nmax),twopid
	    twopid = twopi
	    nout=0.5d0*ofac*hifac*n
	    if(nout.gt.np) then
	        write(*,*) 'output arrays too short : in period'
	        stop
        end if
        
	    call avevar(y,n,ave,var)
	    xmax=x(1)
	    xmin=x(1)
	    do j=1,n
	       if(x(j).gt.xmax)xmax=x(j)
	       if(x(j).lt.xmin)xmin=x(j)
	    enddo
	    xdif=xmax-xmin
	    xave=0.5d0*(xmax+xmin)
	    pymax=0.d0
	    pnow=1d0/(xdif*ofac)
	    do j=1,n
	       arg=twopid*((x(j)-xave)*pnow)
	       wpr(j)=-2.d0*sin(0.5d0*arg)**2
	       wpi(j)=sin(arg)
	       wr(j)=cos(arg)
	       wi(j)=wpi(j)
	    enddo
	    do i=1,nout
	       px(i)=pnow
	       sumsh=0d0
	       sumc=0d0
	       do j=1,n
	          c=wr(j)
	          s=wi(j)
	          sumsh=sumsh+s*c
	          sumc=sumc+(c-s)*(c+s)
	       enddo
	       wtau=0.5d0*atan2(2d0*sumsh,sumc)
	       swtau=sin(wtau)
	       cwtau=cos(wtau)
	       sums=0d0
	       sumc=0d0
	       sumsy=0d0
	       sumcy=0d0
	       do j=1,n
	          s=wi(j)
	          c=wr(j)
	          ss=s*cwtau-c*swtau
	          cc=c*cwtau+s*swtau
	          sums=sums+ss**2
	          sumc=sumc+cc**2
	          yy=y(j)-ave
	          sumsy=sumsy+yy*ss
	          sumcy=sumcy+yy*cc
	          wtemp=wr(j)
	          wr(j)=(wr(j)*wpr(j)-wi(j)*wpi(j))+wr(j)
	          wi(j)=(wi(j)*wpr(j)+wtemp*wpi(j))+wi(j)
	       enddo
	       py(i)=0.5*(sumcy**2/sumc+sumsy**2/sums)/var
	       if(py(i).ge.pymax)then
	          pymax=py(i)
	          jmax=i
	       endif
	       pnow=pnow+1./(ofac*xdif)
	    enddo
	    expy=exp(-pymax)
	    effm=2d0*nout/ofac
	    prob=effm*expy
	    if(prob.gt.0.01d0)prob=1d0-(1d0-expy)**effm
	    return
    end subroutine ls


! AVEVAR from Numerical Recipes
    subroutine avevar(data,n,ave,var)
! given array data(1:n), returns its mean as ave and its variance as var.
	    integer j,n
	    double precision ave,var,data(n),s,ep
	    ave=0d0
	    do j=1,n
	       ave=ave+data(j)
	    enddo
	    ave=ave/n
	    var=0d0
	    ep=0d0
	    do j=1,n
	       s=data(j)-ave
	       ep=ep+s
	       var=var+s*s
	    enddo
	    var=(var-ep**2/n)/(n-1)
	    return
	end subroutine avevar





!!    subroutine ls_fast(x,y,n,ofac,hifac,wk1,wk2,nwk,nout,jmax,prob)
!!! NOTE: this is Numerical Recipes "fasper" routine

!!!   Same as ls but calculated with the fast algorithm of Press & Rybicki (1989)
!!!   with an operation count of order N*log(N)
!!!   Acording to Press et al. (1992) it is already faster than ls for data sets 
!!!   as small as 60 or 100 points.

!!        integer jmax,n,nout,nwk,macc
!!        real hifac,ofac,prob,wk1(nwk),wk2(nwk),x(n),y(n)
!!        parameter (macc=4)
!!!u    uses avevar,realft,spread
!!        integer j,k,ndim,nfreq,nfreqt
!!        real ave,ck,ckk,cterm,cwt,den,df,effm,expy,fac,fndim,hc2wt,hs2wt, &
!!             hypo,pmax,sterm,swt,var,xdif,xmax,xmin
!!        nout=0.5*ofac*hifac*n
!!        nfreqt=ofac*hifac*n*macc
!!        nfreq=64
!!1       if (nfreq.lt.nfreqt) then
!!            nfreq=nfreq*2
!!            goto 1
!!        endif
!!        ndim=2*nfreq
!!        if(ndim.gt.nwk) then
!!            write(*,*) 'workspaces too small : in fasper'
!!            stop
!!        end if
!!        call avevar(y,n,ave,var)
!!        xmin=x(1)
!!        xmax=xmin
!!        do j=2,n
!!            if(x(j).lt.xmin)xmin=x(j)
!!            if(x(j).gt.xmax)xmax=x(j)
!!        end do    
!!      xdif=xmax-xmin
!!        do j=1,ndim
!!            wk1(j)=0.
!!            wk2(j)=0.
!!        end do
!!        fac=ndim/(xdif*ofac)
!!        fndim=ndim
!!        do j=1,n
!!            ck=1.+mod((x(j)-xmin)*fac,fndim)
!!            ckk=1.+mod(2.*(ck-1.),fndim)
!!            call spread(y(j)-ave,wk1,ndim,ck,macc)
!!            call spread(1.,wk2,ndim,ckk,macc)
!!        end do
!!        call realft(wk1,ndim,1)
!!        call realft(wk2,ndim,1)
!!        df=1./(xdif*ofac)
!!        k=3
!!        pmax=-1.
!!        do j=1,nout
!!            hypo=sqrt(wk2(k)**2+wk2(k+1)**2)
!!            hc2wt=0.5*wk2(k)/hypo
!!            hs2wt=0.5*wk2(k+1)/hypo
!!            cwt=sqrt(0.5+hc2wt)
!!            swt=sign(sqrt(0.5-hc2wt),hs2wt)
!!            den=0.5*n+hc2wt*wk2(k)+hs2wt*wk2(k+1)
!!            cterm=(cwt*wk1(k)+swt*wk1(k+1))**2/den
!!            sterm=(cwt*wk1(k+1)-swt*wk1(k))**2/(n-den)
!!            wk1(j)=j*df
!!            wk2(j)=(cterm+sterm)/(2.*var)
!!            if (wk2(j).gt.pmax) then
!!                pmax=wk2(j)
!!                jmax=j
!!            endif
!!            k=k+2
!!        end do
!!        expy=exp(-pmax)
!!        effm=2.*nout/ofac
!!        prob=effm*expy
!!        if(prob.gt.0.01)prob=1.-(1.-expy)**effm
!!        return
!!        
!!    end subroutine ls_fast


! SPREAD from Numerical Recipes
    subroutine spread(y,yy,n,x,m)
!   given an array yy of length n , extirpolate (spread) a value y into m actual array elements
!   that best approximate the “fictional” (i.e., possibly noninteger) array element number x.
!   the weights used are coefficients of the lagrange interpolating polynomial.

        integer m,n
        real x,y,yy(n)
        integer ihi,ilo,ix,j,nden,nfac(10)
        real fac
        save nfac
        data nfac /1,1,2,6,24,120,720,5040,40320,362880/
        if(m.gt.10) then
            write(*,*) 'factorial table too small : in spread'
            stop
        end if
        ix=x
        if(x.eq.float(ix)) then
            yy(ix)=yy(ix)+y
        else
            ilo=min(max(int(x-0.5*m+1.0),1),n-m+1)
            ihi=ilo+m-1
            nden=nfac(m)
            fac=x-ilo
            do j=ilo+1,ihi
              fac=fac*(x-j)
            end do
            yy(ihi)=yy(ihi)+y*fac/(nden*(x-ihi))
            do j=ihi-1,ilo,-1
              nden=(nden/(j+1-ilo))*(j-ihi)
              yy(j)=yy(j)+y*fac/(nden*(x-j))
            end do
        endif
    
        return
    end subroutine spread
      

    subroutine bls(x, y, ofac, hifac, px, py, nf, jmax, prob)
!   this computes the Bayesian Lomb-Scargle periodogram as defined in Gregory (2005)
!   inputs are arrays of dimension nd corresponding to abscissas x(1:nd) (which 
!   need not be equally spaced) and ordinates y(1:nd). Additionally, an oversampling
!   factor ofac sets the frequency resolution as df=1/(ofac*T) where T=x(nd)-x(1)

        implicit none
        
        real(dp), dimension(:), intent(in)      :: x, y
        real(dp), dimension(:), intent(inout)   :: px, py
        real(dp), intent(in)    :: ofac, hifac
        real(dp), intent(out)   :: prob
        integer, intent(out)    :: nf, jmax       
        
        integer     :: ii, jj
        integer     :: nd
        real(dp)    :: y_mean, nd2bar
        real(dp)    :: ti, tf, T, dt_min
        real(dp)    :: f_min, f_max, df, SNR
        real(dp)    :: theta, twopif, fourpif, cos2pift, sin2pift, cos4pift, sin4pift
        real(dp)    :: r, i, c, s, diff, z
        real(dp), allocatable   :: ywork(:), f(:), h2bar(:), p(:)
        
        
        nd = size(x)
        call assert(nd == size(y), 'in bls, arrays x and y must have same dimensions')
        
        y_mean = mean(y)
        allocate(ywork(nd))
        
        ! x need not be sorted:
        ti = minval(x)
        tf = maxval(x)
        T = tf - ti
        
        ! minimum separation in abcissas
        dt_min = T
        do ii=1,nd
            if (x(ii+1) - x(ii) < dt_min) dt_min = x(ii+1) - x(ii)
        end do
        
        ! the bayesian version of the periodogram can produce arbitrarily sharp 
        ! spectral features. In general it's necessary to set the frequency 
        ! resolution smaller than for the LS. The guide is eq 13.4 of Gregory (2005)
        !   df ~ 1/(1.6*SNR*duration*sqrt[nd])
        SNR = (0.5_dp * (maxval(y) - minval(y)) ) / sqrt(2.0_dp)	
        df = 1. / (1.6 * SNR * T * sqrt(dble(nd)))

        f_min = 1.0_dp / T
        f_max = 1.0_dp/(2.0_dp*dt_min)
        nf = ofac * int(abs( (f_max - f_min) / df ))
        call assert(size(px) >= nf, 'bls, px and py must be big enough to hold every frequency')
        allocate(p(nf), h2bar(nf))

        ywork = y - y_mean
        nd2bar = sum(y*y) 
                
!		for(i=0; i<nd[0]; i++) {
!			d[i] -= mean_data;
!			nd2bar += d[i]*d[i];
!		}

		! --- initialize the frequencies
		allocate(f(nf))
		call linspace(f_min, f_max, f)
!		for (i=0; i<numberf; i++)
!			f[i] = f_min + i*df;

		
		! --- compute the lombscargle periodogram
        z = 1.0_dp

		! step through each frequency
		do jj=1,nf
!		for(j=0; j<numberf; j++) {
		    twopif = twopi * f(jj)
		    fourpif = fourpi * f(jj)		
!			twopif  = m_2pi * f[j];
!			fourpif = m_4pi * f[j];

			! calculate theta
			cos4pift = 0.0_dp
			sin4pift = 0.0_dp
			
			do ii=1,nd
!			for(i=0; i<nd[0]; i++) {
                cos4pift = cos4pift + dcos(fourpif * x(ii)) * z*z
				sin4pift = sin4pift + dsin(fourpif * x(ii)) * z*z
			end do
			
			theta = 0.5_dp * datan2(sin4pift, cos4pift)

			! lomb-scargle and p(f|d,i) at this frequency
			r = 0.0_dp
			i = 0.0_dp
			c = 0.0_dp
			s = 0.0_dp
			
			do ii=1,nd
!			for(i=0; i<nd[0]; i++) {
			    cos2pift = dcos(twopif * x(ii) - theta)
				sin2pift = dsin(twopif * x(ii) - theta)
				r = r + ywork(ii) * cos2pift * z
				i = i + ywork(ii) * sin2pift * z
				c = c + cos2pift**2 * z**2
				s = s + sin2pift**2 * z**2
			end do

			! this is the lomb-scargle periodogram:
			h2bar(jj) = (r*r)/c + (i*i)/s

			diff = nd2bar - h2bar(jj)
			
			! this is the bayesian generalization 
			! see bretthorst (2000, 2001) or gregory (2005)
			p(jj) = (1.0_dp/sqrt(c*s)) * diff**((2-nd)/2.0_dp)
		
		end do

        px = f
		py = p
		
		prob = 1.0
		jmax = 1
		
		return
		
	end subroutine bls

end module lib_periodogram
		
		
		
		
		
