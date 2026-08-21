
      subroutine ITERATIONS(kount,iterend,newstep, &
      time,deltime,finaltime,nzit,dencell,tcell,pressure,cellMW,ktime,ktimemax, &
      iterperstep,unsteadytol,resave,resavenorm,resmaxnorm,nqresmaxnorm,isend,ireducestep)

!=============================================================
      use dimensmod
      use nitermod
      use monitormod
      use numericmod
      use phiminmaxmod
      use logicmod
                        
      logical :: isend,ireducestep

! Unsteady
! deltime = finaltime / ntimestep is time step
! iterperstep is iterations/time step
! ktime cycles 1 -> iterperstep -> 1 ... for each deltime

! finaltime                                     24.0s
! ntimestep                                      12
! iterperstep              4
! deltime              2.0s
! ktime        1  2  3  4  1  2  3  4  1  2  3   4
! newstep      0           1           2                 
! newzit       0  1  2  3  0  1  2  3  0  ...

DO
      ! Termination
      if(isend) exit
      
      ireducestep = .false.

      ! Iteration number, ktime, newstep updated here.
      niter = niter + 1
      ktime = ktime + 1               
      if(ktime > iterperstep) ktime = 1
      if(ktime == 1) niter = 1
      if(ktime  ==  1) newstep = newstep + 1

      !  completed time at end of time step deltime
      if(ktime  ==  iterperstep)  time  = time + deltime

! Output intervals    (must be for ktime = iterperstep)
      newzit    = mod(ktime,iterperstep)
      nzit    = max( mod(newstep,10), newzit)


      ! ktime = 1 for 1st iteration

      ! Notice for Last iteration
      if(niter >= iterend)  isend = .true.
      if(time  >= finaltime) isend = .true.

!      call CycleTIMER(niter)

      if(ktime == 1) call UNSTEADYPHIUPDATE(iterperstep,nqLast)

      residmax = 0.0   ! for scalars
      
      call PROPERTIES(dencell,tcell,cellMW)

      ! SPECIES REACTIONS 
      call toPhiMaxNQ

       !  Update Phimag at iupdatephimag niter or newstep intervals (initial update for new steady at iter = 2)
       if(iupdatephimag > 0) then
       if(ktime==1 .and. mod(newstep,iupdatephimag) == 0) call PHIMAGUPDATE
       endif

      ! Convergence test; skips to ktime = iterperstep - 1 if criterion met
      if(unsteadytol > 0.0) call SKIP(ktime,ktimemax,iterperstep,unsteadytol,delphi,nqmaxerr,phimax)

       ! Jacobian solver  
       call CFDJACOBSOLVE(deltime,isend)

      if(ktime == iterperstep) call SPECIESOUTPUTS(resave,resavenorm,resmaxnorm,nqresmaxnorm) 
      if(itimeplot /= 0 .and. kount==1 .and. ktime == iterperstep) call TIMEOUTPUT(newstep,time,Tcell,Pressure,dencell,cellMW)

      ! ITERATION OUTPUT
      if(ktime == iterperstep) then
        if(delphi > unsteadytol) then
          !call ITEROUT(niter,nzit,Tcell,resavenorm,resmaxnorm,nqresmaxnorm,time)
          call ITEROUT(niter,nzit,Tcell,delphi,nqmaxerr,phimax,time)
          
          ireducestep = .true.
          isend = .true.

        endif
      endif
ENDDO

      return
      end subroutine ITERATIONS


!===============================================================
                

!=======================================

      subroutine PRELIMALLOC

!=======================================
      use NUMERICMOD
      use TMOD
      use FOOTMOD
      use filesmod
      use dimensmod      
      use sizemod
      use logicmod
      use PHIMINMAXMOD
      use FUELMOD
            
!--------------------------------------
      allocate (resid(nqlast), residmax(nqlast), urfs(nqlastm), SorsURF(nqSp1:nqlast), &
                 phimagnitude(nqlast), reserveList(nqLastm), stat=istat) 
      call STATCHECK(istat,'          ')      

      allocate (phiresmax(nqlast), stat=istat) 
      call STATCHECK(istat,'          ')     

      resid     = 0.0
      residmax  = 0.0
      phiresmax = 0.0


      return
999   call errorline(filein,'13    ')
      end

!===============================================================


      subroutine STATCHECK(istat,location)
      use filesmod

      character (LEN=*) location
      if(istat.ne.0) then
      write(file_err,*)' Not enough memory for this grid '
      write(file_err,*) location
!      call WINDOWSTRING(' Not enough memory for this grid ')
      stop
      endif
      
      return
      end 
 
!===============================================================                    

      subroutine ALLOC1
      
!=======================================
      use dimensmod      
      use LOGICMOD
      use DENMOD
      use TMOD
      use SIMPLERMOD
      use Enthsourcemod
      use PHIMOD
      use PHIOMOD
      use FUELMOD
      use sizemod
      use filesmod
      use MOLWTMOD
      use TWODIMMOD
            
!--------------------------------------      

      allocate (phi(5:nqlastm), stat=istat) 
      call STATCHECK(istat,'phi       ')     

      phi = 0.0

      allocate (phio(5:nqlast), stat=istat) 
      call STATCHECK(istat,'phio      ')     

      phio = 0.0

      allocate (phihold(nqsp1:nqlast), stat=istat) 
      call STATCHECK(istat,'phihold   ')     

      phihold = 0.0
      dencell = 0.0
      cellMW = 0.0      
      Tcell = 0.0
       
      return
      end 

!=================================================================


      subroutine SCALARCOEFFS(ierr,nq,phi,urf,resmean,fnphi,ap1,dencell,denocell,phio,ap,su,smp,deltime,resmax)

! resmean is unsigned mean(absolute)
! resmax is signed max(absolute) 
!==================================================================

      use filesmod
      use dimensmod
      use paramsmod
      use iterperstepmod
      use twodimmod
      use speciesmod
      use logicmod
             
      real(8) :: fnphi, ap1, apsave
      real(8) :: su, ap

      real :: smp

      data nqprev / 0 /, resmaxall / 0.0 / 
!-----------------------------------------------------------------
      
! assembly of coefficients
      ierr   = 0
      n0     = 0
      resmean = 0.0
      resmax  = 0.0


      if(nq /= nqprev) then
      resmaxall = 0.0
      nqprev = nq
      endif

      n0 = n0 + 1
      apsave = ap
      vol = 1.0
      sp   = 0.0

! transfer from BOUNDS  (sp comes in as -ap;  su is transferred)
!-----------------------------------
! Unsteady

      if(unsteady) then
      aop = denocell * vol / deltime
      anp = dencell  * vol / deltime
      phiold = phi
      if(iterperstep > 1) phiold = phio
      su  = su  +  aop * phiold
      sp  = sp - aop
      smp = smp +  anp - aop
      endif

!-----------------------------------
!  false source term
      cp = max(smp, 0.0)
      su = su + cp * phi
      sp = sp - cp

!-----------------------------------
! AP
      ap = ap - sp

       if(ap  <=  rtiny) then
       write(file_err,'(" STOP: EQN NOT SOLVED: ap <= epsilon for nq ",i5)') nq
       write(file_err,*) ' ap sp apsave ', ap,sp, apsave
       stop
!       ierr = 1
!       return 
       endif
!-------------------------------
! Residual
      resor = 0.0
      resor =  resor - ap * phi  +  su     

      ! fnphi =  S - Ap*Q + sum(AjQj) = 0  for MNEWT  (single i,j,k)
      fnphi = resor    
      ap1   = ap


      resphi  = resor/ap  
      resmean = resmean + abs(resphi)                ! without sign
      if(abs(resphi) > abs(resmax)) then
      resmax = resphi                           ! with sign
       if(abs(resmax) > abs(resmaxall)) then
       resmaxall = resmax
       endif
      endif


      resmean = resmean / max(float(n0),1.0)


!  relaxation
      ap = ap / urf
      su = su + (1. - urf) *  ap * phi
      
      return
      end
!==========================================================

      subroutine UNSTEADYPHIUPDATE(iterperstep,nqLast)

!======================================================  
      use DENMOD
      use phimod
      use phiomod
            
     ! Density update:  deno is density used for previous time step
      denocell = dencell

      if(iterperstep  >  1) phio(5:nqlast) = phi(5:nqlast)

      return
      end

!========================================================

      subroutine CFDJACOBSOLVE(deltime,isend)

!======================================================
      use dimensmod
      use paramsmod
      use filesmod
      use LOGICMOD
      use PHIMOD
      use DENMOD
      use TMOD
      use PHIOMOD
      use JACMOD
      use FNQMOD
      use SIMPLERMOD
      use NUMERICMOD
      use SPECIESMOD
      use MONITORMOD
      use PHIMINMAXMOD
      use NEWLOOPMOD
      use focusmod
      use PHIMINRESMOD
 
                         
      real(8) :: tsu(nqSp1:nqlast), tsp(nqSp1:nqlast)
      real(8) :: ap1
      real :: Yprev(nqsp1:nqlastm)
      real :: Yspecies(nqsp1:nqlastm)

      logical :: isend
      real :: ressign(nqsp1:nqlast) 
!-----------------------------------------------------

      resid(nqsp1:nqlast)    = 0.0 
      residmax(nqsp1:nqlast) = 0.0 

      sumsource = 0.0       
      resmon(nqsp1:nqlast)  = 0.0
       
  
! resnormijk = residual(nq)/Phimax(nq) averaged for all nq at specific i,j,k
! Criterion for skipping cells or Focus is ratio of resnormijk / resnormean 
! resnormean is average of residual(nq)/Phimax(nq) for all i,j,k  
  
      su = 0.0
      ap = 0.0
      smp = 0.0


        ! N-R method used for focus
        resnormijk = rbig    ! to start cycle for each i j k

         ! Reactions
         vol = 1.0
         Yspecies(nqsp1:nqlastm) = Phi(nqsp1:nqlastm)
         call ALLREACTSOURCE(1,-1,Tcell, Yspecies, dencell, tsu,tsp,vol)

         dfdq = dfdq * vol 

         ! fnq(nq) and ap1      fnq =   S - Ap*Q + sum(AjQj) = 0
         fnq = 0.0
          
         resnormijk = 0.0

          ! Species Loop -------------------------------
          do nq = nqsp1, nqlast
          if(.not. solvenq(nq)) cycle
          
          su =    tsu(nq) * vol    ! kg/s
          ap =  - tsp(nq) * vol
          
          if(ap < 0.0 .or. su < 0.0)then
          write(file_err,*)'STOP: ap or su -ve before SCALARCOEFFS in CFDJACOBSOLVE'
          write(file_err,*)'nq,ap,su ', nq,ap,su      
          write(file_err,*)'tsu(nq),-tsp(nq) ',tsu(nq),-tsp(nq) 
          stop
          endif
          
          ! neighbouring coeffs and residual at one location
          ! residual = |ressign(nq)| residual is +ve; ressign is +ve or -ve 
          call SCALARCOEFFS(ierr,nq,phi(nq), &            
                urfs(nq),residual,fnq(nq),ap1,dencell,denocell,phio(nq),ap,su,smp,deltime,ressign(nq))

          ! Self differential is -Ap; overwrite chemical source already in Ap
          dfdq(nq,nq)  = - Ap1

          ! for Focus
          if(phimaxnq(nq) > phiminres) then
          resnormijk = resnormijk + residual/Phimaxnq(nq)  ! sum normalised resids for all nq at one i,j,k  
          endif
          enddo  ! End nq Loop ----------------------------

          ! Sum source output (old sources with old set of phis)
          do nq = nqsp1, nqlast
          if(.not. solvenq(nq)) cycle
          sumsource(nq)   = sumsource(nq) + (tsu(nq) + tsp(nq) * phi(nq)) * vol
          enddo

          if(resnormijk > rbig) then   ! solution has diverged
          write(file_err,*)' Huge residual. Solution has diverged. Stop '
          ! stop
          endif

          do nq = nqsp1, nqlast
          if(.not. solvenq(nq)) cycle
          resid(nq) = resid(nq) + abs(ressign(nq))    ! summing over all i,j,k for particular nq,iregn  
          enddo 

          ! Solve all variables at i,j,k
          Yprev(nqsp1:nqlast)    = Yspecies(nqsp1:nqlast)
          call MNEWT(Yspecies,fnq,dfdq,maxspecies, URFspeciesdflt)

          ! Correct Over-Under shoots. (Under-relax = 0.1) 
          do nq = nqsp1, nqLast        
          if(.not. solvenq(nq)) Yspecies(nq) = 0.0   ! removes mnewt junk
          call PHIMINMAXSUB(nq,Yprev(nq),Yspecies(nq),nmovenq(nq),novernq(nq),nundernq(nq))
          enddo
          Phi(nqsp1:nqlast) = Yspecies(nqsp1:nqlast)       ! all variables
          phi(nqBalance)    = Ybalance(Yspecies)   

      resmon(nqsp1:nqlast) = resmon(nqsp1:nqlast) + resid(nqsp1:nqlast)

!------------------------------------------------------------------------------------
! End reports
      if(isend) then

      call ELEMENTSOURCEBAL

      do nq = nqsp1, nqlast 
      !if(novernq(nq)  >  0) write(fileRPT,'(i6," cells  Allowed overshoot for ",  a24)') novernq(nq),  spname(nq)
      !if(nmovenq(nq)  > 0)  write(fileRPT,'(i6," cells  MOVEMENT LIMITED for ",   a24)') nmovenq(nq),  spname(nq)
      !if(nundernq(nq) > 0)  write(fileRPT,'(i6," cells  UNDERSHOOT LIMITED for ", a24)') nundernq(nq), spname(nq)
      enddo      

      endif
!--------------------

      return
      end
!======================================================

      SUBROUTINE mnewt(x,fvec,fjac,np,urf)

!======================================================
! from Newton-Raphson method for non-linear systems of equations 
! #9.6 Numerical Recipes 2nd ed. Press et. al. Cambridge 1992

      use filesmod  
      use DIMENSMOD
      implicit  integer(i - n), real(8) (a - h, o - z)

      real    :: x(np), urf ! (np)
      integer :: indx(np),inot(np)
      real(8) :: fvec(np),p(np) 
      real(8) :: fjac(NP,NP)

!---------------------------------------------

        p(1:np) = -fvec(1:np)

        call ludcmp(fjac,np,indx,d,inot)
        call lubksb(fjac,np,indx,p)

! urf correction
        x(1:np) = x(1:np) + float(1 - inot(1:np)) * p(1:np) * urf  ! (i)

      return
      END
!==================================

      SUBROUTINE ludcmp(a,n,indx,d,inot)

!===================================
      use filesmod
      
      implicit  integer(i - n), real(8) (a - h, o - z)
       
      real(8), PARAMETER :: TINY8 = 1.0d-40
      INTEGER :: indx(n),inot(n)
      real(8) :: a(N,N)
      real(8) :: suma, aamax, dum       
      real(8) :: vv(N)
!----------------------------------

      inot = 0

      d = 1.0
      
      do i = 1,n
       aamax = 0.0
       do j = 1,n
       aamax = max( aamax, abs(a(i,j)) )
       enddo

        
       if (aamax < TINY8) then
       inot(i) = 1  
       aamax   =  TINY8              !   a(i,i)  = -TINY      ! small self-coeff
       endif

      vv(i) = 1./aamax
      enddo

      do j = 1,n
       do i = 1,j-1
       suma = a(i,j)
        do k = 1,i-1
        suma = suma - a(i,k)*a(k,j)
        enddo
       a(i,j) = suma
       enddo
       aamax = 0.0
        
       do i = j,n
       suma = a(i,j)
        do k = 1,j-1
        suma = suma - a(i,k)*a(k,j)
        enddo
       a(i,j) = suma
       dum = vv(i)*abs(suma)
        if (dum >= aamax) then
        imax = i                     ! must solve all species to avoid bad imax value here.
        aamax = dum
        endif
       enddo
        
       if(j /= imax)then
       do k = 1,n
       dum = a(imax,k)
       a(imax,k) = a(j,k)
       a(j,k) = dum
       enddo
       d = -d
       vv(imax) = vv(j)
       endif
        
       indx(j) = imax                        
       if(a(j,j) == 0.0) a(j,j) = TINY8
        
       if(j /= n)then
       dum = 1.0/a(j,j)
        do i = j+1,n
        a(i,j) = a(i,j)*dum
        enddo
       endif
      enddo
      
      return
      END
!===================================

      SUBROUTINE lubksb(a,n,indx,b)

!===================================

      implicit  integer(i - n), real(8) (a - h, o - z)

      INTEGER indx(n)
      REAL(8) :: b(n), a(n,n)
      real(8) :: sum
      
      ii=0

      do i=1,n
        LL=indx(i)
        sum=b(LL)
        b(LL)=b(i)
        if (ii.ne.0)then
          do j=ii,i-1
          sum=sum-a(i,j)*b(j)
          enddo
        else if (sum /= 0.0) then
          ii=i
        endif
        b(i)=sum
      enddo

      do i=n,1,-1
        sum=b(i)
        do j=i+1,n
        sum=sum-a(i,j)*b(j)
        enddo
        b(i)=sum/a(i,i)
      enddo
      return
      END
!==========================================================================


      subroutine PHIMINMAXSUB(nq,phiprev,phi,nmovement,nover,nunder)

!===================================================
      use DIMENSMOD
      use PHIMINMAXMOD
      
! Movement limitation if fracmove's > 0.0
      if(fracmoveup   > 0.0) then
      phiNewMax = min((phiprev + fracmoveup   * phimagnitude(nq)), 1.0)
       if(phi > phiNewMax ) then
       phi = phiNewMax
       if(phi > criterionmove * phimagnitude(nq)) nmovement = nmovement + 1   ! only record if phi is large enough
       endif 
      endif

      if(fracmovedown > 0.0) then
      phiNewMin = max((phiprev - fracmovedown * phiprev), 0.0)
      if(phi < phiNewMin ) then
       phi = phiNewMin
       if(phi > criterionmove * phimagnitude(nq)) nmovement = nmovement + 1
       endif 
      endif


! Max limit
      if(phi > 1.0) then
!       if(nq <= nqLast .and. phi > 1.0) then
       nover = nover + 1
       phi  = 1.0      
!       endif

! Min limit
      elseif(phi  <  0.0) then
      nunder = nunder + 1    
      phi = 0.0
      endif

      return
      end
!===================================================

      subroutine toPhiMaxNQ

! Puts maximum phi into Phimaxnq for Laminar
! CMC uses conditioned values

!======================================================================
      use dimensmod
      use filesmod
      use phimod
      use MONITORMOD

      do nqspec =  nqsp1,nqLast
      Phimaxnq(nqspec) = phi(nqspec)
      enddo

      return
      end
!===================================================================                            

      subroutine PHIMAGUPDATE

! phimagnitude(nq) from Phimaxnq
!=======================================================================

      use dimensmod
      use filesmod
      use PHIMINMAXMOD
      use monitormod

      do nq = nqsp1, nqlast
      phimagnitude(nq) = PHIMAGFN(Phimaxnq(nq))
      enddo

      return
      end

!===================================================================                            

      function PHIMAGFN(phimax)

! sets effective phimagnitude from phimax 
!===================================================

      use Phiminmaxmod
            
      PHIMAGFN = min(phimax, 1.0)   ! phimagitude cannot be > 1.0
      PHIMAGFN = max(phimagfn, phimagnitudemin)      ! has to be >= phimagnitudemin

      return
      end
!==================================================

      subroutine ALLOCDFDQ(nqsp1a,nqlasta)

!=======================================================================
      use dimensmod
      use JACMOD
      use FNQMOD
      use filesmod
      use sizemod

      write(fileRPT,'(/" dfdQ Allocations ")')
      write(fileRPT,'(" nqlasta      ",i6)') nqlasta
      write(fileRPT,'(" nqsp1a       ",i6)') nqsp1a

      allocate (dfdQ(nqsp1a:nqlasta,nqsp1a:nqlasta), fnq(nqsp1a:nqlasta), stat=istat)
      call STATCHECK(istat,'dfdQ      ')

      dfdq = 0.0

      return
      end
!=======================================================================

      subroutine SKIP(ktime,ktimemax,iterperstep,unsteadytol,delphi,nqmaxerr,phimax)

!==============================================================  
      use dimensmod
      use phimod
      use PHIMINRESMOD

      
         if(ktime < iterperstep-1 .and. ktime>1) then

         delphi = 0.0
         nqmaxerr = 0
         do nq = nqsp1,nqlast 
           if(phi(nq) < phiminres) cycle
           delphinq = abs(phi(nq)-phihold(nq))/phi(nq)
           if(delphinq > delphi) then
             delphi = delphinq
             nqmaxerr = nq
             phimax = phi(nq)
           endif      
         enddo  
          if(delphi < unsteadytol) then
          ktimemax = max(ktimemax,ktime)
          ktime = iterperstep - 1
          endif
         endif

     phihold = phi(nqsp1:nqlast)   

     return
     end 

!==============================================================  
        logical function EQAL(x,y)

!=========================================================

        parameter( epsilon = 1.e-6)

        if ( abs(x-y)  <=  abs( x * epsilon) ) then
            eqal = .true.
        else
            eqal = .false.
        endif

        return
        end

!====================================================================

        logical function NEQAL(x,y)

!=========================================================

        parameter( epsilon = 1.e-6)

        if ( abs(x-y)  >  abs( x * epsilon) ) then
           neqal = .true.
        else
           neqal = .false.
        endif

        return
        end

!====================================================================

        logical function NEQAL8(x,y)

!=========================================================

        real(8), parameter :: epsilon = 1.e-6
        real(8) :: x,y

        if ( abs(x-y)  >  abs( x * epsilon) ) then
           neqal8 = .true.
        else
           neqal8 = .false.
        endif

        return
        end

!====================================================================
