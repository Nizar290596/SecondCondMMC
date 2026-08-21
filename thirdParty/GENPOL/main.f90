
!   GENPOLsingle     
!   Copyright:  J. H. Kent.
       
      Program GENPOLsingle

      use dimensmod
      use filesmod
      use phimod
      use speciesmod

!------------------------------------------------

      call STARTUP
      
      ! READ DATA
      open(unit = filedataset,   file = 'soot/dataset.in',   status='old')
      read(filedataset,*) iparticle
      read(filedataset,*) T     ! T(K)
      read(filedataset,*) P      ! P(kPa)
      read(filedataset,*) dt 
      do nq = nqsp1,nqlastm
        read(filedataset,*) phi(nq)
      enddo 
      close(filedataset)
      
      call PARTICLEDATAINPUT(phi(nqsp1:nqlastm),T,P)

      call EACHPARTICLE(dt,phi(nqsp1:nqlastm),T,P)

      call PARTICLEDATAOUTPUT(phi(nqsp1:nqlastm))

!      ! Phi mass fractions output
      open(unit = filephiout,  file=  'soot/phiout.txt',   status='replace')
      do nq = nqsp1, nqlastm
      write(filephiout,'((1pg12.3,t20,a24))') phi(nq),spName(nq) 
      enddo
      close(filephiout)
      
      !call PLOTIME    ! optional timed output

      call SHUTDOWN
      end program 

!===============================================================

      subroutine STARTUP

!======================================================================
      use bin1mod
      use bin2mod
      use dimensmod
      use elementsmod
      use filesmod
      use focusmod
      use isendmod
      use iterperstepmod
      use logicmod
      use molwtmod
      use monitormod
      use monitorposmod
      use nitermod
      use numericmod
      use paramsmod
      use phiminmaxmod 
      use phiminresmod
      use phimod
      use phiomod
      use phitimemod
      use pressmod 
      use reactmod
      use TAGMOD
      use simplermod
      use sizemod
      use sootmod
      use sectionalmod
      use speciesmod
      use thermomod
      use tmod
      use twodimmod
      use DENMOD

     
      
!=================================================================
! Start-up
      time = 0.0
      maxit=0 
      ierrsolve = 0

      call OPENFILES()

      call SPECIESCOUNT    ! defines nqlast, nqlastm
      call REACTIONCOUNT
      call ALLOCSPEC

      call PRELIMALLOC

      call INPUT(iterperstep,unsteadytol,deltimenom)

      if(maxit==0) isend = .true.
      call OPENFILES2()
      call ALLOC1

      call READSPECIESNUMERICAL(iterperstep)
      call SPECIESIN

      call THERMOIN            
      call Hquadratic

! REACTIONS
      call REACTIN(1,idum) ! needed for gasification reactions even if Flamelet
      call ALLOCREACT
      call REACTIN(2,ierrsolve)

      ! Put zero iterations if Species Not-solved error
       if(ierrsolve /= 0) then
       write(file_err,'(/" STOP: Species-Needed error. See above. ")')
       maxit = 0
       isend = .true.
       endif

      call ALLOCDFDQ(nqsp1,nqlast)

      ktime = iterperstep

!  INITIAL OUTPUT
      call INITOUT(deltime,maxit)

      close(fileIN)
      close(filescheme)
      close(filethermo)

      return
      end
!==============================================================  

      subroutine PARTICLEDATAINPUT(Ypart,Tpart,Ppart)

      use filesmod
      use MONITORPOSMOD
      use TMOD
      use PRESSmod
      use dimensmod
      use phimod
      use PHIMINMAXMOD
      use MOLWTMOD
      use DENMOD
      use logicmod
      
      real :: Ypart(nqsp1:nqlastm),Tpart,Ppart                            
       
      do nq = nqsp1,nqlastm
        phi(nq) = Ypart(nq)
        phihold(nq) = phi(nq)
        if(nq<=nqlast) phimagnitude(nq) = max(phimagnitude(nq),phi(nq))           
      enddo 
      Tcell = Tpart
      Pressure = Ppart
      
!      ! READ DATA
!      open(unit = filedataset,   file = 'soot/dataset.in',   status='old')
!      read(filedataset,*) Tcell     ! T(K)
!      read(filedataset,*) Pressure      ! P(kPa)
!      read(filedataset,*) dt 
!      do nq = nqsp1,nqlastm
!        read(filedataset,*) phi(nq)
!      enddo 
!      close(filedataset)
      

      call SpeciesSumCheck

      cellMW = WTMOLFN(phi(nqsp1:nqlast))
      dencell =  denfn(Tcell,cellMW)
      
      return
      end
!==============================================================  

      subroutine PARTICLEDATAOUTPUT(Ypart)

      use filesmod
      use MONITORPOSMOD
      use TMOD
      use PRESSmod
      use dimensmod
      use phimod
      use PHIMINMAXMOD
      use MOLWTMOD
      use DENMOD
      use logicmod
      
      real :: Ypart(nqsp1:nqlastm)                          
       
      do nq = nqsp1,nqlastm
        Ypart(nq) = phi(nq) 
      enddo 
      
      return
      end
!==============================================================  

      subroutine EACHPARTICLE(dt,Ypart,Tpart,Ppart)
      
!==============================================================  
      use bin1mod
      use bin2mod
      use dimensmod
      use elementsmod
      use filesmod
      use phimod
      use speciesmod
      use phitimemod
      use nitermod
      use isendmod
      use denmod
      use Tmod
      use pressmod
      use molwtmod
      use iterperstepmod                                      
      use logicmod
      
      logical :: ireducestep
      real :: Ypart(nqsp1:nqlastm),Tpart,Ppart                            


      deltimeuse = dt
DO      
      time = 0.0
      niter = 0
      ktime = iterperstep
      newstep = 0
      isend = .false.
     
      kount = kount + 1   ! counts particles           
      call TIMESTEPMANAGEMENT(dt,finaltime,ntimestep,maxit,iterend,deltimeuse,deltime,iterperstep)

     !----------------------------------------------     
     ! Time output 
     if(kount==1) then
     call SUM_FIELD(sum0_element)

      ! Allocate Time-step Output File 
      nt = ntimestep+2
      allocate (phitimeout(nt,nqlastm),fvtime(nt),sootNtime(nt), sootDtime(nt),sootHtoCtime(nt), stat=istat) 
      call STATCHECK(istat,'phitimeout')     
      phitimeout = 0.0

     if(itimeplot /= 0) call TIMEOUTPUT(newstep,time,Tcell,Pressure,dencell,cellMW)
     endif
     !----------------------------------------------     

     call ITERATIONS(kount,iterend,newstep, &
      time,deltime,finaltime,nzit,dencell,tcell,pressure,cellMW,ktime,ktimemax, &
      iterperstep,unsteadytol,resave,resavenorm,resmaxnorm,nqresmaxnorm,isend,ireducestep)
           
     if(.not.ireducestep) then
       if (deltimeuse < dt) then
         write(*,*) 'Reducing timestep has lead to successful convergence'
         write(file_err,*) 'Reducingtime step has lead to successful convergence'
       
       endif
       
       exit
       
     endif
     
     if(ireducestep) then
       if(deltimenom >= deltimeuse) then
         write(*,*) 'Reduce time step and still fails. CONVERGENCE FAILURE !!!!!!',dt,deltimenom
         write(file_err,*) 'Reduce time step and still fails. CONVERGENCE FAILURE !!!!!!',dt,deltimenom

         exit
       elseif(deltimenom > dt) then
         write(*,*) 'Cannot reduce timestep to get convergence because deltimenom >= dt. CONVERGENCE FAILURE !!!!!!',dt,deltimenom
         write(file_err,*) 'Cannot reduce timestep to get convergence because deltimenom >= dt. CONVERGENCE FAILURE !!!!!!',dt,deltimenom

         exit
         
       else
         write(*,*) 'Reducing timestep in attempt to get convergence'
         write(file_err,*) 'Reducing timestep in attempt to get convergence'
         
       endif
       
     endif
     
     deltimeuse = deltimenom
     
     call PARTICLEDATAINPUT(Ypart,Tpart,Ppart)
ENDDO 
      
     return
     end 

!==============================================================  
!===============================================================

      subroutine NEXTDATASET(finaltime,ntimestep,maxit,iterend,deltimenom,deltime,iterperstep)
      
      use filesmod
      use MONITORPOSMOD
      use TMOD
      use PRESSmod
      use dimensmod
      use phimod
      use PHIMINMAXMOD
      use MOLWTMOD
      use DENMOD
      use logicmod
            
! Read next data set     
      read(filedataset,*,end=999,err=999) iparticle
      read(filedataset,*) Tcell     ! T(K)
      read(filedataset,*) Pressure      ! P(kPa)
      read(filedataset,*,err=999) finaltime 

      do nq = nqsp1,nqlastm
      read(filedataset,*,err=999,end=999) phi(nq)
      if(nq<=nqlast) phimagnitude(nq) = max(phimagnitude(nq),phi(nq))           
      enddo 

      ! Time steps
      deltime = min(deltimenom,finaltime)
      ntimestep = nint( finaltime /deltime)
      maxit     = ntimestep * iterperstep    ! total iterations for whole time
      iterend = maxit

      call SpeciesSumCheck
      cellMW = WTMOLFN(phi(nqsp1:nqlast))
      dencell =  denfn(Tcell,cellMW)

!      write( fileRPT, '(/" Unsteady for this run." )')
!      write( fileRPT, '(" Initial time                    = ", 1pg12.3,1x,"s" )') time
!      write( fileRPT, '(" Final time                      = ", 1pg12.3,1x,"s" )') finaltime
!      write( fileRPT, '(" deltime                         = ", 1pg12.3,1x,"s" )') deltime
!      write( fileRPT, '(" No. of time steps               = ", i6 )')    ntimestep
!      write( fileRPT, '(" No. of iterations               = ", i6 )')    maxit

      return

!999   write(file_err,*)' Stop. Error reading species mass fractions. '
!      stop 

999   stop

      end
!=========================================================================================

      subroutine SHUTDOWN
      use filesmod
      use logicmod      

      close(filedataset)
      close(filephiout)
      close(fileRPT)
      close(filenumeric)
      close(filereact)
      close(file_err)
      if(itimeplot/=0) close(filephi)
      
      stop
      end
      
!=========================================================================================
      subroutine TIMESTEPMANAGEMENT(dt,finaltime,ntimestep,maxit,iterend,deltimenom,deltime,iterperstep)
      
      use filesmod
      use MONITORPOSMOD
      use TMOD
      use PRESSmod
      use dimensmod
      use phimod
      use PHIMINMAXMOD
      use MOLWTMOD
      use DENMOD
      use logicmod
            
      time = 0.0

      ! Time steps
      finaltime = dt
      deltime = min(deltimenom,finaltime)
      ntimestep = nint( finaltime /deltime)
      maxit     = ntimestep * iterperstep    ! total iterations for whole time
      iterend = maxit


!---------------------------------------------------
      return

      end
!=========================================================================================
                      
