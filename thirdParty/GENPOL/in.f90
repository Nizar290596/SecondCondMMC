
!=============================================================

      subroutine INPUT(iterperstep,unsteadytol,deltimenom)

!       reads in values from fileIN and echo prints them to fileRPT

!====================================================================

      use paramsmod
      use dimensmod
      use filesmod
      use logicmod
      use fuelmod
      use TAGMOD
      use sootmod
      use speciesmod
      use Tresetmod
      use TWODIMMOD
      use phiminresmod
      use numericmod
      use sectionalmod
                                          
!------------------------------------------------------------------
      write( fileRPT, '(// 20x, " GENPOL "// )')

!  READ  INPUT FILE
      iterperstep  = 1

      unsteady = .true.        

! GENERAL SECTION
      call findheader(fileIN,isfound,'GENERAL_SECTION')
      if(isfound == 0) stop

      call findots(fileIN)
      read( fileIN,*,err=999) deltimenom
      
      ! unsteadytol is max phi fractional variation which terminates iterations per timestep, singl gridcell only.
      read(fileIN,*,err=999,end=999) iterperstep,unsteadytol

      write( fileRPT, '(/" Unsteady for all runs." )')
      write( fileRPT, '(" Iterations/time step              = ", i6 )')    iterperstep 
      write( fileRPT, '(" Skip tolerance (single grid cells only) = ", g12.3 )') unsteadytol
!------------------------------------------------------------

! PROGRAM METHODS

      EBUmixfrac = .false. 

      METHOD = 10 
      write(fileRPT, '(/1x,79(1h-)/)')
      write( fileRPT, '(" Laminar flame with full chemistry. ")')
      LamChemmethod = .true.
      UseTransportData = .false.
       
      isLamEnthEqn = .true.
      write(fileRPT, '(/1x,79(1h-)/)')

      solvenq = .false.
      istraj       = .false.
      isdispersion = .false.
      isradn       = .false.
      sootrad      = .false.

      solvenq(nqsp1:nqlast) = .true.
      isreacting   = .true.

     
      write( fileRPT, '( "    Program Control")')
      write(fileRPT,'(/"   Trajectories   Dispersion   Combustion      Radiation")')   
      write(fileRPT,'(6x,4(L1,14x))')istraj, isdispersion, isreacting, isradn


      write(fileRPT, '(/1x,79(1h-)/)')
!----------------------------------------------------------------
!--------------------------------------------------------------------
!  OTHER PARAMETERS

!-------------------------
! TEMPERATURE UNITS
!  0 = iCentigrade_in,  1 = Centigrade   input & output


      iCentigrade_in=0
      CK_in = ' K'
      write(fileRPT,*)' Input temperatures in DEGREES KELVIN     '


      iCentigrade_out = 0 
      CK_out = ' K'
      head_Tout =  'Temperature K   '
      write(fileRPT,*)' Output temperatures in DEGREES KELVIN     '

! ---------------------------------------

! OUTPUT SECTION
      call findheader(fileIN,isfound,'OUTPUT_SECTION')
      if(isfound == 0) stop
      call findots(fileIN)

!--------------------
! for PlotOut
      read(fileIN,*,err=999) timeunit
      call ILimitin(timeunit,1,5,"timeunit")
      call UNITofTIME

! Species Output Units 1 = Mole fraction; 0 = Mass fraction  -1 = Source term (kg/m^3 s) 

      read(fileIN,*,err=999)  iSpUnits_out
      call ILimitin(iSpUnits_out,0,3," OUTPUT UNITS needs to be 0 - 3 ")
      write(fileRPT,'(/"SPECIES OUTPUT UNITS are: ",a10)')SpUnits_out(iSpUnits_out) 

!-------------------------------------------
! Report normalised residuals 
      read(fileIN,*,err=999)  Phiminres

! Phimax threshold for Res/Phimax
      write(fileRPT,'(/" Phimax threshold for normalised residuals res/phimax is ",g12.3)') Phiminres

!-------------------------------------------

!  Particle sources printout (disabled)
      printspart(nqsp1:nqlast) = 0


!-------------------------------------------
! Contours for TAG Reaction Rates
      iTagcontour=0

!----------------------
! Special Outputs
! up to maxspecialout X's for special outputs
! up to maxspecialout R's for special outputs
      Xspecialout = 1.e30
      Rspecialout = 1.e30
      nXspecialout = 0
      nRspecialout = 0
!-------------------------------------------
! 0 or 1 Output for Chemkin input
   ichemkin=0  
!----------------------

! 63 = D63 for soot D; 0 = volume averaged soot 
      read(fileIN,*,err=999) iD63 
      isD63 = .true. 
      if(iD63 == 0) isD63 = .false.

! Time plot output
      read(fileIN,*,err=999) itimeplot 

!----------------------

! -------------------------------------------------------------

! Sectional Product Output summed Groups (up to maxsectoutGp) for graphics   (includes Fv D N)
!           Ai Aj     C1  C2   H1  H2    S1   S2     T1  T2                    
!   0 0 means full range
      call findots(fileIN)
      i = 0
      DO WHILE(.true.)
      call endlist(fileIN,iend)
      if(iend == 1) exit
      i = i+ 1
      read(fileIN,*,err=999) iSectOutGp(:,i)
      !  make 0  0  full range
       call ZEROisFULLRANGE(iSectOutGp(1,i),iSectOutGp(2,i),2)
       call ZEROisFULLRANGE(iSectOutGp(3,i),iSectOutGp(4,i),nCsections)
       call ZEROisFULLRANGE(iSectOutGp(5,i),iSectOutGp(6,i),nHsections)
       call ZEROisFULLRANGE(iSectOutGp(7,i),iSectOutGp(8,i),nSsections)
       call ZEROisFULLRANGE(iSectOutGp(9,i),iSectOutGp(10,i),nTsections)
      if(i == maxsectoutGp) exit 
      ENDDO
      NsectoutGp = i


! Reaction Rate summed for each Tag in list
! iTagsum(1,i)  = Tag No.

! iTagsum(2,i) is mode 
! 0 = rate 
! 1 = C-mole rate 
! 2 = H-mole rate 
! 3 = Sect net C-mole rate
! 4 = Sect net H-mole rate

! Tag   mode      R1 Sect Bin Range        R2 Sect Bin Range       Product section Bin Range               
!  No.   No.      C1-C2  H1-H2 T1 T2       C1-C2  H1-H2 T1 T2         C1-C2  H1-H2  T1 T2
!   0 0 means full range

      call findots(fileIN)

      i = 0
      DO WHILE(.true.)
      call endlist(fileIN,iend)
      if(iend == 1) exit
      i = i+ 1
      read(fileIN,*,err=999) iTagsum(:,i)
!      call ZEROisFULLRANGE(iTagsum(3,i),iTagsum(4,i),  nCsections)
!      call ZEROisFULLRANGE(iTagsum(5,i),iTagsum(6,i),  nHsections)
!      call ZEROisFULLRANGE(iTagsum(7,i),iTagsum(8,i),  nTsections)
!      call ZEROisFULLRANGE(iTagsum(9,i),iTagsum(10,i), nCsections)
!      call ZEROisFULLRANGE(iTagsum(11,i),iTagsum(12,i),nHsections)
!      call ZEROisFULLRANGE(iTagsum(13,i),iTagsum(14,i),nTsections)
!      call ZEROisFULLRANGE(iTagsum(15,i),iTagsum(16,i),nCsections)
!      call ZEROisFULLRANGE(iTagsum(17,i),iTagsum(18,i),nHsections)
!      call ZEROisFULLRANGE(iTagsum(19,i),iTagsum(20,i),nTsections)

       if(itagsum(1,i) < 1) then
       write(file_err,*)'Stop. Tag No. must be > 0 at line ',i
       stop
       endif
       if(itagsum(2,i)<0 .or. itagsum(2,i)>4) then
       write(file_err,*)'Stop. Tag mode must be 0 - 4, is ',itagsum(2,i)
       stop
       endif

      if(i == maxTagsum) exit 
      ENDDO
      NTaglist = i
      isTagsum = .false.
      if(NTaglist > 0) isTagsum = .true.

! -------------------------------------------------------------

      return
999   call errorline(filein,' 1    ')
      end

!=====================================================================

      subroutine ZEROisFULLRANGE(ibin1,ibin2,maxbin)
       
      if(ibin1 == 0 .or. ibin2 == 0) then
      ibin1 = 1
      ibin2 = maxbin
      endif
       
      return
      end
!=====================================================================

      subroutine DEFAULTBOUND(twall,uwallin)

!  Default wall temperature & default wall heat transfer coefficient

      use paramsmod
      use filesmod
      use dimensmod
      use logicmod
!-----------------------------------------------------------------
      call findots(fileIN)

      read(fileIN,*,err=999) Twall
      if(iCentigrade_in == 1)  twall = twall + T0K

      read(fileIN,*,err=999) Uwallin
      if(uwallin < 0.) uwallin = -1.e-20

! WALL TEMPERATURE FUNCTION
      read(fileIN,*,err=999) walltempfn

      if(walltempfn) then
      write( fileRPT, '( / " Wall conduction calculated. "/ )')
      else
      write( fileRPT, '( / " Wall temperature remains constant. "/ )')
      endif

!-------------------------

      return
999   call errorline(filein,' 2    ')
      end

!=====================================================================

      subroutine readline(nfile,str,iend,Linechange)

! if blank before ! then seeks next line

      use filesmod
      character(Len=*) :: str
      character*7 dash, dots
      data dash /'-------'/
      data dots  /'.......'/
      Linechange = 0

      str = ' '      
      DO while(.true.)
      read(unit=nfile, end=999, err=999, fmt = '((a))') str
      Linechange = Linechange + 1


! dots
      if(str(1:7) == dots(1:7)) then
       write(file_err,*)' Unexpected ............... at '
       backspace(nfile)
       backspace(nfile)
       Linechange = linechange -2
        do  n = 1,5
        read(nfile,'(a80)') str
        Linechange = linechange + 1
        write(file_err,'(1x,a80)') str
        enddo
      STOP
      endif

! end of list
      if(str(1:7) == dash(1:7)) then
      iend = 1
      RETURN
      endif

! blank line or ! following leading blanks
      do i =1, 80
      if(str(i:i) /= ' ' .and. str(i:i) /=  char(9)) exit
      enddo

      if(i < 80 .and. str(i:i) /= '!') exit  ! cycle for comment line or blank line for 80 cols
      ENDDO

! data
      iend = 0
      RETURN


999   iend = 1
      RETURN
      
     
      end
!=====================================================
 
       subroutine NOLEADINGBLANKS(Line)

! suppresses leading blanks in Line 
!========================================================

      use paramsmod
      character(Len=LineL) :: Line, temp
      
      L = Len_trim(Line)

      iblank = 0
      do i = 1,L
      if(Line(i:i) /= ' ' .and. Line(i:i) /= char(9)) exit
      iblank = i
      enddo 
      
      if(iblank == 0) return

      temp = Line
      Line = ' '
      Line(1:L-iblank) = temp(iblank+1:L) 

      return
      end

!========================================================
      subroutine findots(nfile)

! finds next row of dots
!===================================================================
      use filesmod
      character*7 chr,spec
      data spec /'.......'/


10    read(nfile,'(3x,a7)',end=20) chr
      if(chr /= spec) goto 10
      return

20    write(file_err,*)' String not found ',spec
      stop
      end

!=====================================================

      subroutine findheader(nfile,isfound,header)

! rewinds file and finds specified header
!=====================================================
      use filesmod
      use paramsmod
      character (Len = *) :: header
      character (Len = LineL) :: chr

      Lhead = LEN_TRIM(header)
      rewind(nfile)
      isfound = 1

10    read(nfile,'(a)',end=20) chr
      if(chr(1:Lhead) /= header(1:Lhead)) goto 10
      return

20    write(file_err,*)' Cannot find Section title: ',header
      write(file_err,*)' Title must start in column 1; is case sensitive '
      isfound = 0
      end

!=====================================================
      subroutine endlist(nfile,iend)

      use filesmod
      character*80 str
      character*7 dash, dots
      data dash /'-------'/
      data dots  /'.......'/

5     read(nfile,'(a80)') str

! comment line if first character is c or C or # or * or !
      if(str(1:1) == 'c'.or.str(1:1) == 'C' .or. str(1:1) == '#'.or.str(1:1) == '*'.or.str(1:1) == '!') then
      if(str(1:1) == '*') write(fileRPT,'(1x,a80)') str
      goto 5
      endif

! dots
      if(str(3:9) == dots(1:7)) then
      write(file_err,*)' End of list dashed line ----- missing before '
      backspace(nfile)
      backspace(nfile)
      read(nfile,'(a80)') str
      write(file_err,'(a80)') str
      stop
      endif

! end of list
      if(str(3:9) == dash(1:7)) then
      iend = 1
      return
      endif

! blank line
      do 10 i =1, 80
      if(str(i:i) /= ' ') goto 6
10    continue
      goto 5

! data
6     iend = 0
      backspace(nfile)
      return

      end
!=====================================================

      subroutine errorline(nfile,at)

      use filesmod
      character(Len=80):: str
      character(Len=*) :: at

      backspace(nfile)
      write(file_err,'(" ERROR in input file in first of following lines at ",(a))') at


      do n = 1,5
      read(nfile,'(a80)',end=20) str
      write(file_err,'(1x,a80)') str
      enddo

                    ! error trace
                    ! y=1.0; z=y/x; a=1/z

20    stop
      end

!=======================================================================

      subroutine ILimitin(idata,Limlow,Limhigh,name)
      
!==============================================================

! integer data must be equal to or within range of limits

      use filesmod
      character (LEN=*) ::  name

      if(idata < Limlow) then
      write(file_err,'(/1x,(a)," = ",i6," is lower than allowed limit ",i6  )') name,idata,Limlow
      stop
      elseif(idata > Limhigh) then
      write(file_err,'(/1x,(a)," = ",i6," is higher than allowed limit ",i6  )') name,idata,Limhigh
      stop
      endif

      return
      end
      
!==============================================================
      

      subroutine RLimitin(rdata,rLimlow,rLimhigh,name)
      
!==============================================================

! real data must be equal to or within range of limits

      use filesmod
      character (LEN=*) ::  name

      if(rdata < rLimlow) then
      write(file_err,'(/1x,(a)," = ",g12.5," is lower than allowed limit ",g12.5  )') name,rdata,rLimlow
      stop
      elseif(rdata > rLimhigh) then
      write(file_err,'(/1x,(a)," = ",g12.5," is higher than allowed limit ",g12.5 )') name,rdata,rLimhigh
      stop
      endif

      return
      end
      
!==============================================================
      

      subroutine Rorder(dirn,arr,n, reference)

!==============================================================
! check whether in order
! put dirn =1.0 for ascending order; -1.0 for descending order

      use filesmod
      real :: arr(n)
      character (Len=*) :: reference      
      do i = 2,n  
      if(dirn*arr(i) <= dirn*arr(i-1) ) then 
      write(file_err,*)' Input data is not in ascending/descending order for ',reference
      stop
      endif
      enddo
       
      return
      end
      
!==============================================================
      
     
      subroutine UNITofTIME

      use logicmod
      
       if(timeunit == 1) then
       unit_tosec = 1.0e-6
       unitname_time = 'usec'
       elseif(timeunit == 2) then
       unit_tosec = 1.0e-3
       unitname_time = 'msec'
       elseif(timeunit == 3) then
       unit_tosec = 1.0
       unitname_time = 's'
       elseif(timeunit == 4) then
       unit_tosec = 60.0
       unitname_time = 'min'
       elseif(timeunit == 5) then
       unit_tosec = 3600.0
       unitname_time = 'hr'
       endif

      return
      end

!==============================================================
      subroutine FINDINDEX(x,i,frac,a,n,iascend)

!=========================================

! a(1:n) is monotonically increasing or decreasing with i
! returns i, fr  so that   x = frac * a(i+1) + (1.-frac) * a(i)   

! if x < a(1) i = 1    frac = 0.0  
! if x > a(n) i = n-1  frac = 1.0  
!---------------------------
      use filesmod
      real :: a(1:n)
      
      if(n == 1) then
      i = 1
      frac = 0.0
      iascend = 0 
      return
      endif

      iascend = -1
      if(a(2)>a(1)) iascend = 1

      if( (x<=a(1) .and. iascend == 1) .or. (x>=a(1) .and. iascend == -1) ) then
      i = 1
      frac = 0.0
      return
      endif

      if( (x>=a(n) .and. iascend == 1) .or. (x<=a(n) .and. iascend == -1) ) then
      i = n-1
      frac=1.0
      return
      endif

      call LOCATE(a,1,n,x,i,1,n-1)

      frac = (x - a(i)) / (a(i+1) - a(i)) 
     

      return
      end

!======================================================
      subroutine LOCATE(a,n0,n,x,j,nmin,nmax)

!==================================================================
! returns index j where x lies between a(j) & a(j+1); 
! x >  a(j) for ascending 
! x <= a(j) for descending series 
! a must be monotonic increasing or decreasing      

! can specify nmin = n0, or n0-1 for minimum return index   
! can specify nmax = n-1, or n for maximum return index (make n-1 so that interpolation can be made)  

      real :: a(n0:n)

      jmin = max(n0-1,nmin); jmin = min(n0,jmin)
      jmax = min(n,nmax); jmax = max(n-1,jmax)        

      jL = jmin
      ju = jmax+1       
      
      do while(ju-jL > 1)
      jm = (ju+jL)/2
      
       if((a(n) > a(1)) .eqv. (x > a(jm))) then
       jL = jm
       else
       ju=jm
       endif

      enddo

      j=jL
      
      return
      end
      
!==================================================================

      subroutine LOCATE8(a,n0,n,x,j,nmin,nmax)

!==================================================================
! for real(8)

! returns index j where x lies between a(j) & a(j+1); 
! x >  a(j) for ascending 
! x <= a(j) for descending series 
! a must be monotonic increasing or decreasing      

! can specify nmin = n0, or n0-1 for minimum return index   
! can specify nmax = n-1, or n for maximum return index (make n-1 so that interpolation can be made)  

      real(8) :: a(n0:n), x

      jmin = max(n0-1,nmin); jmin = min(n0,jmin)
      jmax = min(n,nmax); jmax = max(n-1,jmax)        

      jL = jmin
      ju = jmax+1       
      
      do while(ju-jL > 1)
      jm = (ju+jL)/2
      
       if((a(n) > a(1)) .eqv. (x > a(jm))) then
       jL = jm
       else
       ju=jm
       endif

      enddo

      j=jL
      
      return
      end
      
!==================================================================
      subroutine AbsoluteFileName(directory,relfilename,absfilename,length)

      character (len=*) :: relfilename
      character (len=256) :: directory,absfilename

      absfilename = relfilename      
      Length =  len_trim(absfilename)     

      return
      end

!===================================================================

      function isTAGINRANGE(I,Imin,Imax)
      logical :: isTaginrange
      
      isTaginrange = .false. 

      if(Imax==0) then
      isTaginrange = .true.
      return
      
      elseif(I >= Imin .and. I <= Imax) then
      isTaginrange = .true.
      endif
       
      return
      end
!============================================                                                

      subroutine SPLITNAME(fullname,Lengthpath,Lengthfull)

! shows length of path of filename and length of full filename
!============================================                                                

! Lengthpath is character number of \ or :
 
      character(Len=256) fullname

      Lengthfull = LEN_TRIM(fullname)
             
      do n = Lengthfull,1,-1   
      if(fullname(n:n) == '\' .or. fullname(n:n) == ':') exit
      enddo

      Lengthpath = n 
      
      return
      end
!============================================                                                
               
      subroutine OPENFILES()
      
!======================================================                  

      use filesmod
      
      character (len=24) :: filename
      character (len=256) :: absfilename
      logical :: exists

! KEEP.oi absolute name
      keepfilename = 'keep.oi'
      cmcoifilename = 'cmc.oi' 
      flameletoifilename = 'flamelet.oi'

      open(unit = file_err, file = 'soot/error.txt', status='unknown')
      write(file_err,*,iostat=io) ' ' 
       if(io/=0) then
       write(*,'(///" Opening error: ERROR FILE IS OPEN  ")')
       call DELAY(3.0)
       stop 
       endif

       
      open(unit = fileRPT,        file = 'soot/report.txt',     status='unknown')
      open(unit = filereact,      file = 'soot/reaction.txt',   status='unknown')
      open(unit = filenumeric,    file = 'soot/numerical.txt',  status='unknown')
      !open(unit = filephiout,     file=  'soot/phiout.txt',     status='replace')
      !open(unit = filedataset,    file = 'soot/dataset.in',     status='old')



      call AbsoluteFileName(datadirectory,'soot/scheme.in',absfilename,length)
      inquire(file = absfilename, exist = exists)
       if(.not.exists) then
       write(file_err,*) 'Not available: ',absfilename
       stop
       endif
      open(unit=filescheme, file=absfilename, status='old')

      call AbsoluteFileName(datadirectory,'soot/thermo.in',absfilename,length)
      inquire(file = absfilename, exist = exists)
       if(.not.exists) then
       write(file_err,*) 'Not available: ',absfilename
       stop
       endif
      open(unit=filethermo, file=absfilename, status='old')


      !write(*, '( " GENPOL  -  J.H. KENT, Compuflow Solutions, Sydney, Australia."/)') 
      write(filenumeric,'( " GENPOL  -  J.H. KENT, Compuflow Solutions, Sydney, Australia."/)',err=999) 

!---------------------------------------------
      filename = 'sootSettings.in'
      absfilename = 'soot/sootSettings.in'

      !call AbsoluteFileName(datadirectory,filename,absfilename,L)
      open(unit=fileIN,err=10, file=absfilename, status='old')

!---------------------------------------------


      write( fileRPT, '(// "    GENPOL Version ",f6.1)',err=999) version
      write( fileRPT, '(   "    Copyright 2001 - 2013  J.H. KENT, Compuflow Solutions, Sydney, Australia.")')
      write( fileRPT, '(   "    All rights reserved. "//)') 


      write( fileRPT,      '(// "   INPUT FILE  ",256a//)') absfilename
      write( file_err,     '(// "   INPUT FILE  ",256a//)') absfilename
      write( filereact,    '(// "   INPUT FILE  ",256a//)',err=999) absfilename
      write( filenumeric,  '(// "   INPUT FILE  ",256a//)',err=999) absfilename

      return

999   write(file_err,*)' An output file is open. Close all files. '
      stop 
      
10    write(file_err,*) ' Input file soot/sootSettings.in is missing  '
      stop

      end
!======================================================  

      subroutine OPENFILES2()
      
!======================================================                  

      use filesmod
      use logicmod
      

! Phi file name for output
      if(iSpUnits_out==0) then
      phifilename = 'soot/phimass.txt'
      elseif(iSpUnits_out==1) then
      phifilename = 'soot/phimole.txt'
      else
      phifilename = 'soot/phi.txt'
      endif


      ! phimole, phimass
      if(itimeplot /= 0) open(unit = filephi,file = phifilename,status = 'unknown')
      
      return

999   write(file_err,*)' An output file is open. Close all files. '
      stop 

      end
!======================================================  


           
