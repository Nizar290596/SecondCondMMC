!==============================================================  

     subroutine TIMEOUTPUT(newstep,time,Tcell,Pressure,dencell,cellMW)

!==============================================================  
     use paramsmod
     use dimensmod
     use ELEMENTSmod
     use phitimemod
     use phimod
     
     call SUM_FIELD(sum1_element)

     istpout = newstep + 1
     phitimeout(istpout,1) = time     
     phitimeout(istpout,2) = Tcell
     phitimeout(istpout,3) = Pressure
     phitimeout(istpout,4) = dencell
     phitimeout(istpout,5) = cellMW

     if(sum0_element(1)>rtiny) phitimeout(istpout,6) = sum1_element(1)/sum0_element(1)  ! Cfinal/Cinitial
     if(sum0_element(2)>rtiny) phitimeout(istpout,7) = sum1_element(2)/sum0_element(2)  ! Hfinal/Hinitial
     phitimeout(istpout,nqsp1:nqLastm) = phi(nqsp1:nqlastm)
     Lastimestepout = istpout
     return
     end 

!==============================================================  

      subroutine PLOTIME

! for time-dependent Plot output
!=====================================================================      

      use paramsmod
      use dimensmod
      use logicmod
      use phitimemod
      use filesmod
      use speciesmod
      use SECTIONALMOD

      integer, parameter :: ngroup = 5 ! Grouped species sums No. 2 - ngroup

      character(Len=32) :: hedgroup(2:ngroup)
      character(Len=32) :: SectGphead(NsectoutGp), SectGpname 
      logical :: INSECTGROUP

      real :: scr(Lastimestepout,1,1), X(Lastimestepout), ydum(1)=0.0, zdum(1) = 0.0

!     phitimeout(newstep,1) = time     
!     phitimeout(newstep,2) = T
!     phitimeout(newstep,3) = Pressure
!     phitimeout(newstep,4) = den
!     phitimeout(newstep,5) = cellMW
!     phitimeout(newstep,6) = Cfinal/Cinitial
!     phitimeout(newstep,7) = Hfianl/Hinitial

      if(itimeplot == 0) return
      if(ntimestep == 0) return
           
      do n  = 2, Ngroup
      write(hedgroup(n),'("Group",i2.2)') n
      enddo

      do n  = 1, NsectoutGp
      write(SectGphead(n),'("SectGp",i2.2)') n
      enddo

      write(filephi,'(//20x," Species units are ",a10//)') SpUnits_out(iSpUnits_out)

      ifirst = 1
      X = phitimeout(1:Lastimestepout,1)/unit_tosec               

      Tmax = maxval( phitimeout(1:Lastimestepout,2) )          ! Temperature
      scr(1:Lastimestepout,1,1) = phitimeout(1:Lastimestepout,2) - (T0K*iCentigrade_out)
      call PLOTOUTPUT(filephi,ifirst,Lastimestepout,1,1,10,1,Tmax,X,ydum,zdum, &
      scr, 'Temperature',unitname_time,' ',CK_out)

      Pmax = maxval( phitimeout(1:Lastimestepout,3) )            ! Pressure
      scr(1:Lastimestepout,1,1) = phitimeout(1:Lastimestepout,3)
      call PLOTOUTPUT(filephi,ifirst,Lastimestepout,1,1,10,1,Pmax,X,ydum,zdum, &
      scr, 'Pressure',unitname_time,' ','kPa')

      phimax = maxval( phitimeout(1:Lastimestepout,6) )          ! Cfinal/Cinitial
      scr(1:Lastimestepout,1,1) = phitimeout(1:Lastimestepout,6)
      call PLOTOUTPUT(filephi,ifirst,Lastimestepout,1,1,10,1,phimax,X,ydum,zdum, &
      scr, 'C/Cinitial',unitname_time,' ',' ')

      phimax = maxval( phitimeout(1:Lastimestepout,7) )          ! Hfinal/Hinitial
      scr(1:Lastimestepout,1,1) = phitimeout(1:Lastimestepout,7)
      call PLOTOUTPUT(filephi,ifirst,Lastimestepout,1,1,10,1,phimax,X,ydum,zdum, &
      scr, 'H/Hinitial',unitname_time,' ',' ')

! Species: Single 
      do nq = nqsp1,nqlastm
      if(printphis(nq) == 1) then
       do i = 1,Lastimestepout
       scr(i,1,1) = Y_out(phitimeout(i,nq),iSpUnits_out,phitimeout(i,5),phitimeout(i,4),nq)
       enddo
      phimax = maxval( scr(1:Lastimestepout,1,1) )
      call PLOTOUTPUT(filephi,ifirst,Lastimestepout,1,1,10,1,phimax,X,ydum,zdum,scr, &
                      spname(nq),unitname_time,' ',SpUnits_out(iSpUnits_out))
      endif
      enddo


! Sum Gas-phase Species for Printphi = 2 - ngroup
      do igroup = 2,ngroup
      write(filephi,'(///" Plot SumSpecies ",i2," contains: ")') igroup              
      scr = 0.0
      isprintgroup= 0
            
       do nq = nqsp1,nqlastm
       if(abs(printphis(nq)) == igroup) then
       isprintgroup=1
       write(filephi,'(5x,a24)') spname(nq)
        do i = 1,Lastimestepout
        scr(i,1,1) = scr(i,1,1) + Y_out(phitimeout(i,nq),iSpUnits_out,phitimeout(i,5),phitimeout(i,4),nq)
        enddo
       endif
       enddo
       
      if(isprintgroup==1) then
      phimax = maxval( scr(1:Lastimestepout,1,1) )
      call PLOTOUTPUT(filephi,ifirst,Lastimestepout,1,1,10,1,phimax,X,ydum,zdum,scr, &
                      hedgroup(igroup),unitname_time,' ',SpUnits_out(iSpUnits_out))
      endif
      enddo


!  Sectional Groups
      write(filephi,'(//" Sectional Species Groups. ")')

      do igroup = 1,NsectoutGp
      scr = 0.0
      isprintgroup= 0
      write(filephi,'(///" SectGp ",i2.2," contains: ")') igroup
       do  nq = nqSectionmol1,nqSectionradn
       if(INSECTGROUP(nq,igroup)) then
       write(filephi,'(5x,a24)') spname(nq)
       isprintgroup=1
        do i = 1,Lastimestepout
        scr(i,1,1) = scr(i,1,1) + Y_out(phitimeout(i,nq),iSpUnits_out,phitimeout(i,5),phitimeout(i,4),nq)
        enddo
       endif
       enddo

      if(isprintgroup==1) then
      phimax = maxval( scr(1:Lastimestepout,1,1) )
      call PLOTOUTPUT(filephi,ifirst,Lastimestepout,1,1,10,1,phimax,X,ydum,zdum,scr, &
                      SectGphead(igroup),unitname_time,' ',SpUnits_out(iSpUnits_out))
      endif
      enddo

! converted to soot Fv, N, D   
      DO igroup = 1,NsectoutGp
      call SOOT_FV_D_N(fvtime,sootNtime,sootDtime,sootHtoCtime,phitimeout(1,2),phitimeout,phitimeout(1,4), &
              ntimestep,1,nqLastm,nqsp1,nqLastm,1,Lastimestepout,iany2,igroup)

      if(iany2==1) then
      write(SectGpname,'("SectGp",i2.2,"_FV")') igroup
      scr(:,1,1) = fvtime(:)
      phimax = maxval( scr(1:Lastimestepout,1,1) )
      call PLOTOUTPUT(filephi,ifirst,Lastimestepout,1,1,10,1,phimax,X,ydum,zdum,scr,SectGpName,unitname_time,' ','v/v')

      write(SectGpname,'("SectGp",i2.2"_N")')  igroup
      scr(:,1,1) = sootNtime(:)
      phimax = maxval( scr(1:Lastimestepout,1,1) )
      call PLOTOUTPUT(filephi,ifirst,Lastimestepout,1,1,10,1,phimax,X,ydum,zdum,scr,SectGpName,unitname_time,' ','/cm3')

      write(SectGpname,'("SectGp",i2.2"_D")')  igroup
      if(isD63) write(SectGpname,'("SectGp",i2.2"_D63")')  igroup
      scr(:,1,1) = sootDtime(:)
      phimax = maxval( scr(1:Lastimestepout,1,1) )
      call PLOTOUTPUT(filephi,ifirst,Lastimestepout,1,1,10,1,phimax,X,ydum,zdum,scr,SectGpName,unitname_time,' ','nm')

      write(SectGpname,'("SectGp",i2.2"_HoC")')  igroup
      scr(:,1,1) = sootHtoCtime(:)
      phimax = maxval( scr(1:Lastimestepout,1,1) )
      call PLOTOUTPUT(filephi,ifirst,Lastimestepout,1,1,10,1,phimax,X,ydum,zdum,scr,SectGpName,unitname_time,' ',' ')
      endif

      ENDDO

      return
      end
!=========================================================================           

      subroutine ITEROUT(niter,nzit,Tcell,resmaxnorm,nqresmaxnorm,phimax,time)

! ITERATION OUTPUT

!==============================================================
      use paramsmod
      use filesmod
      use dimensmod
      use speciesmod
      use logicmod 
      use monitormod 
      use numericmod 
      use PHIMINRESMOD
      use twodimmod
      use elementsmod
      use sectionalmod
                                     
      tempout = tcell - (T0K * iCentigrade_out)

!-----------------------------------------
! SCREEN each iteration
!@         write(*, '(i5,1pg9.1,0pf6.0, 3(1pg9.1),1x,a13)')  niter, time,tempout, &
!         resavenorm,resmaxnorm, Phimaxnq(nqresmaxnorm),spName(nqresmaxnorm)

        write(file_err, '(/"Warning: GENPOL unsteady solver has not converged")')
        write(file_err, '(/"See numerical.txt for details")')

        write(filenumeric, '(i5,1pg12.3,0pf6.0, 3(1pg10.2),1x,a14)')  niter,time,tempout, &
        resmaxnorm, phimax,spName(nqresmaxnorm)
        
        write(*, '(/"Warning: GENPOL unsteady solver has not converged")')
        write(*,'(/"  iter   time(s)     T    Rmax/max  Phimax ")')
        write(*, '(i5,1pg12.3,0pf6.0, 3(1pg10.2),1x,a14)') &
        niter,time,tempout,resmaxnorm, phimax,spName(nqresmaxnorm)


! SCREEN periodically

     IF(nzit == 0) then
     ! Cout/Cin   Hout/Hin    
     Coutoin = 0.0
     Houtoin = 0.0
     if(sum0_element(1) /= 0.0) Coutoin = sum1_element(1)/sum0_element(1)
     if(sum0_element(2) /= 0.0) Houtoin = sum1_element(2)/sum0_element(2)
!@     write(*,'(" Cin ",1pg12.2," Hin ",1pg12.2,5x," Cout/in ",0pf8.3,"  Hout/in ",0pf8.3)')Cin,Hin,Coutoin,Houtoin
     write(filenumeric,'(/" Cin ",1pg12.2," Hin ",1pg12.2,10x," Cout/in ",0pf8.3,"   Hout/in ",0pf8.3)')Cin,Hin,Coutoin,Houtoin
     ENDIF


! Screen headings
     if(nzit == 0) call ScreenHead1

      
      return
      end
!==============================================================

      subroutine SPECIESOUTPUTS(resave,resavenorm,resmaxnorm,nqresmaxnorm)

!======================================================================
      use dimensmod
      use filesmod
      use phimod
      use DENMOD
      use TMOD
      use MONITORMOD
      use logicmod
      use paramsmod
      
      real :: scr(0:nqLastm)

      ! Max of Mean Residual/phimax for screen output  
      scr(0:nqlast) = resmon(0:nqlast)
      call AVERESNORM(resavenorm,resave,Phimaxnq(1:nqlast),scr)
      call MAXRESNORM(resmaxnorm,nqresmaxnorm,Phimaxnq(1:nqlast),scr)

!---------------------      
      return
      end
!===================================================================                            


      function INSECTGROUP(nq,igroup)

      use LOGICMOD
      use sectionalmod
      use filesmod
            
      logical :: INSECTGROUP
      
      INSECTGROUP = .false.              
      if(nq < nqSectionmol1 .or. nq > nqSectionradn) return
       call SectionfromNQ(nq,iC,iH,iS,iT,m)
       if( m >= iSectOutGp(1,igroup) .and.  m <=iSectOutGp(2,igroup) .and. &    
          iC >= iSectOutGp(3,igroup) .and. iC <=iSectOutGp(4,igroup) .and. &    
          iH >= iSectOutGp(5,igroup) .and. iH <=iSectOutGp(6,igroup) .and. &    
          iS >= iSectOutGp(7,igroup) .and. iS <=iSectOutGp(8,igroup) .and. &    
          iT >= iSectOutGp(9,igroup) .and. iT <=iSectOutGp(10,igroup) ) INSECTGROUP = .true.

      return
      end

!=====================================================
      subroutine INITOUT (deltime,maxit)

!===========================================================
      use dimensmod
      use filesmod
      use speciesmod
      use TWODIMMOD
      use logicmod
      
!----------------------------------------------------------------


! Screen
!@      if(unsteady.and.maxit > 0) write(*,'("    time step      ",g10.3," sec ")') deltime
      if(unsteady.and.maxit > 0) write(filenumeric,'("    time step      ",g10.3," sec ")') deltime

       write(filenumeric,'(/"  iter   time(s)     T",a2,"    Rmax/max  Phimax ")') CK_out
!@       write(*          ,'(/"  iter  time(s)  T",a2,"    Rave/max  Rmax/max  Phimax     i  j ")') CK_out

      return
      end
!===========================================================

      subroutine ScreenHead1
      
!==============================================================

      use filesmod
      use logicmod

       write(filenumeric,'(/"  iter   time(s)     T",a2,"    Rave/max  Rmax/max  Phimax ")') CK_out
!@       write(*          ,'(/"  iter  time(s)  T",a2,"    Rave/max  Rmax/max  Phimax     i  j ")') CK_out

       return
       end
!==============================================================
      subroutine AVERESNORM(resnormave,resave,phimax,res)

!=====================================================================      

      use dimensmod
      use phiminresmod
      use filesmod

      real :: phimax(nqlast), res(0:nqlast)      
      
      resnormave = 0.0
      resave     = 0.0
      kount=0

      do nq =  nqsp1,nqLast         ! nqLast
      if(phimax(nq) > phiminres) then
      resnormave = resnormave + abs(res(nq))/phimax(nq)
      resave     = resave     + abs(res(nq))
      kount = kount + 1
      endif
      enddo

      resnormave = resnormave / max(float(kount), 1.0)
      resave     = resave     / max(float(kount), 1.0) 

      return
      end
!=====================================================================      

      subroutine MAXRESNORM(resnormax,nqresnormax,phimax,res)

!  max of res/phimax for nq range and nq for max
!=====================================================================      

      use dimensmod
      use phiminresmod
      use filesmod
       
      real :: phimax(nqlast), res(0:nqlast)      

      
      nqresnormax = nqsp1
      resnormax = 0.0

      do nq =  nqsp1,nqLast        
       if(phimax(nq) > phiminres) then
       resoverphi = res(nq)/phimax(nq)
        if(abs(resoverphi) > abs(resnormax)) then
        resnormax = resoverphi
        nqresnormax = nq
        endif
       endif
      enddo

      return
      end
!=====================================================================      
!=======================================================

      subroutine PLOTOUTPUT(ifile,ifirst,ni,nj,nk,ncol,irn,phimax,x,y,z,phi, &
                 title,Xaxisname,Raxisname,Unitsname)

!================================================================

      use dimensmod
      use logicmod
      use paramsmod
            
      real :: phi(ni,nj,nk)
      real :: x(ni), y(nj), z(nk) 
      character (LEN=*) :: title,Xaxisname,Raxisname,Unitsname
      character (len=24) :: form

             
! format
      form = "(1x,i4,1p50e12.3) " 
      
      if(ifirst == 1) then

      write(ifile, '(" Format  " )')
      write(ifile, '("          at  R1   R2    R3    R4 " )')
      write(ifile, '("  at  X1 " )')
      write(ifile, '("  at  X2 " )')
      write(ifile, '("  at  X3 " )')
      
      write(ifile, '(//"PLOTOUTPUT")')    
      write(ifile, '(i6,"     Max. No. of axials        " )') ni
      write(ifile, '(i6,"     Max. No. of radials       " )') nj
      write(ifile, '(i6,"     Data columns/line         " )') ncol
      write(ifile, '(i6,"     Region               " )') irn
      write(ifile, '((a),"    General species units" )') SpUnits_out(iSpUnits_out)
      write(ifile, '(" 2       format version" )')   
      write(ifile, '(/" X axis ", (a))') Xaxisname
      write(ifile, '(" Y axis ",(a))')  Raxisname
      ifirst = 0
      endif

      write(ifile, '(/"......................")')     
      write(ifile, '(a)') title
      write(ifile, '(a)') Unitsname
      write(ifile, '(1pg12.3,"     Max. value " )') phimax
      write(ifile, '(3i6,"     This variable has: total No. of cols(Rs), rows(Xs), nk " )') nj,ni,nk

      do k = 1,nk
      if(nk>1) write(ifile, '(/ 7x," PLANE k =  ",i3,"   Z = ",f12.5)') k, z(k)*rad2deg

      jsta =  -ncol + 1
      jend =  0

       do while(jend < nj)
       jsta = jsta + ncol
       jend = min(nj, (jsta + ncol-1))
       write(ifile, '(/1x,T18,100(i6,6x))')(j, j = jsta, jend)
       write(ifile, '("===============================================")')   
       write(ifile, '(1x,T18,1p50e12.3)') (y(j), j = jsta, jend)
       write(ifile, '("-----------------------------------------------")')   
        do   i = 1, ni
        write(ifile, form) i,x(i),(phi(i,j,k), j = jsta, jend)
        enddo
       enddo

      enddo

      return
      end

!================================================================


      subroutine DELAY(time)
                                            ! t0 = SECNDS(0.0)
      secs=0.0
      call system_clock(it1,icrate,icmax)

      do while (secs < time)
                                            !      secs = secnds(t0) 
      call system_clock(it2,icrate,icmax)
      if(it2 < it1) it2 = it2+icmax  
      secs =  real(it2-it1) / real(icrate)

      enddo

      return
      end
!=========================================================================           
