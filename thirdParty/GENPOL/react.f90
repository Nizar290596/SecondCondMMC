
      subroutine ALLOCREACT
! Allocates with nreactarray instead of Maxreactlines
!====================================================================
      use dimensmod
      use filesmod
      use LOGICMOD
      use SPECIESMOD
      use MONITORMOD
      use REACTMOD
      use sizemod
      use THERMOMOD
      use sectionalmod
      use TAGmod

!--------------------------------        
      deallocate(nsignificantreact,nmaxreact,nminreact,nstorereact,storerate)
      allocate (nsignificantreact(nreactarray),nmaxreact(nreactarray),nminreact(nreactarray), &
                nstorereact(nreactarray),storerate(nreactarray),stat=istat) 
      call STATCHECK(istat,'nsignif   ')       
      mBmonitor = mBmonitor + nreactarray*7

      deallocate(iFORD)
      allocate (iFORD(nreactarray), stat=istat) 
      call STATCHECK(istat,'FORD      ')       
      iFord = 0
!---------------------------------------

      deallocate(chemrate, reactratesum, fwdratesum, phi_T) 
      allocate (chemrate(3,nreactarray), reactratesum(nreactarray), fwdratesum(nreactarray), &
                 phi_T(nreactarray),   stat=istat) 
      call STATCHECK(istat,'chemrate  ')        
      mbREACT = mbREACT + nreactarray*6 

      deallocate (reactmon, fwdreactmon) 
      allocate (reactmon(nreactarray), fwdreactmon(nreactarray), stat=istat) 
      call STATCHECK(istat,'reactmon  ')        
      mbMonitor = mbMonitor + 2*nreactarray  

      deallocate (Eqnimage, markreact) 
      allocate (Eqnimage(nreactarray), markreact(nreactarray), stat=istat) 
      call STATCHECK(istat,'eqnimage  ')        
      mbREACT = mbREACT + 2*nreactarray 

      deallocate (numbereact,Linereact, nforward, ireactype, nqorderSp, ordermols,reversible, iauxindex, reactantC) 
      allocate (numbereact(0:nreactarray+1),Linereact(nreactarray), nforward(nreactarray), ireactype(nreactarray),itag(nreactarray), &
              nqorderSp(maxtermsreaction,nreactarray), ordermols(maxtermsreaction,nreactarray),reversible(nreactarray), &
              iauxindex(nreactarray),reactantC(nreactarray), stat=istat) 
      call STATCHECK(istat,'numbereact')        
      mbREACT = mbREACT + nreactarray*7 + nreactarray*maxtermsreaction*2 

      numbereact = 0
      iTag = 0                
      markreact = '   '
      nqorderSp = 0
      ordermols = 0.0
      ireactype = 0    ! default for elementary reactions
      phi_T  = 0.0     ! default makes gamma coag = 1
             
      if(iallelementary /= 1)  then
      deallocate (expon)     
      allocate (expon(nqSp1:nqlastm,nreactarray),stat=istat)     
      call STATCHECK(istat,'expon     ')        
      mbREACT = mbREACT + (nqlastm-nqsp1+1)*nreactarray 
      else
      deallocate (expon)     
      allocate (expon(1,1),stat=istat)     
      call STATCHECK(istat,'expon     ')        
      endif
      
      deallocate (rLnKequil,reactmolnet) 
      allocate (rLnKequil(0:nreactarray),reactmolnet(nreactarray),stat=istat) 
      call STATCHECK(istat,'reactmol  ')        
      mbREACT = mbREACT + nreactarray*2 

      deallocate (specnetmass,nqspec, nq3rd, iterm3rd) 
      allocate (specnetmass(maxtermsreaction,nreactarray), &
                nqspec(maxtermsreaction,nreactarray), nq3rd(nreactarray), iterm3rd(2,nreactarray),stat=istat) 
      call STATCHECK(istat,'specmol   ')        
      mbREACT = mbREACT + nreactarray*2 + maxtermsreaction*nreactarray*2 


      reactratesum = 0.0
      fwdratesum   = 0.0
      reactmon     = 0.0
      fwdreactmon  = 0.0
       
!--------------------------------------

      return
      end

!================================================================================

      subroutine REACTIONCOUNT

! Initial count of number of reactions to get array sizes.
! Sizes are too large because omitted reactions are not yet counted.
!================================================================================
      use paramsmod
      use dimensmod 
      use filesmod 
      use logicmod
      use sectionalmod
      use sizemod
      use reactmod
                              
      character (LEN=LineL) ::  Line

!---------------------------------------------------------------
! Reactions
      call findheader(fileSCHEME,isfound,'REACTION_SECTION')
      if(isfound == 0) stop

      call findots(fileSCHEME)
      call findots(fileSCHEME)

      iallelementary = 1
      ireact = 0
      nGenericSectReact = 0 
      Linesreqd = 0
      maxnAiAj = 0
      naux   = 0
      maxtermsreaction = 0
      maxFORD = 0
      newiford = 1
      Linenumber = 0

      DO while (.true.)
      Line = ' '
      call READLINE(fileSCHEME,Line,iend,Linechange)
      Linenumber = Linenumber + Linechange
      if(iend == 1) exit
      call NOLEADINGBLANKS(Line)

      ! Count maximum No. of terms by "+" in product side
       kountnext = 1
       do i = 1,maxreactcols
       if(Line(i:i) == '=') kountnext = 1
       if(line(i:i) == '+') kountnext = kountnext + 1 
 
       call EXITREACTION(Line(i:i+1),i,iterminated)
       if(iterminated == 1) exit
       enddo

       maxtermsreaction = max(kountnext+3,maxtermsreaction)  

       ! Count non-consecutive FORD or RORD for special exponents        
       if(Line(1:4) == 'FORD' .or. Line(1:4) == 'RORD') then
       if(newiford==1) maxFORD = maxFORD + 1
       newiford = 0
       cycle
       else
       newiford = 1
       endif

       ! Auxiliary rate data: count by presence of LOW or HIGH
       if(Line(1:3) == 'LOW' .or. Line(1:4) == 'HIGH') then
       naux = naux + 1  
       cycle
       endif

      ! Count / for No. of 3rd body species in single line to be expanded as separate lines
      call RESERVEDWORDS(Line,isreserved)   !to avoid slashes in LOW etc.
      if(isreserved == 0) then
      nslash = 0
      do i = 1,LineL        ! search for /
      if(Line(i:i) == '!') exit            ! avoid any /'s after !
      if(Line(i:i) == '/') nslash = nslash + 1
      enddo
      N3rdbodyinline = (nslash+1)/2          ! allows for missing end /
      if(N3rdbodyinline > 0) ireact = ireact + N3rdbodyinline  - 1  ! -1 because line is counted anyway
      endif

       ! extra reactions if sectional method
       if(Sectionalmethod)  then           
       call EXTRAREACTIONS(Linenumber,ireact,maxreactcols,nAiAj,Linesreqd,Line)
       if(nAiAj > 0) nGenericSectReact = nGenericSectReact + 1 
       maxnAiAj = max(nAiAj,maxnAiAj)
       else
       ireact = ireact+1
       endif
      ENDDO

!      maxtermsnonsectreaction = maxtermsreaction

      if(Sectionalmethod .and. maxnAiAj == 0)  then  
      write(file_err,'( " There are no sectional reactions for Sectional Method ")')
      stop
      endif

       ! at least 6 for product matrix = 3 in REPEATCHECK
      maxtermsreaction = max(maxtermsreaction,6) 
      if(Sectionalmethod .and. nHsections == 1) maxtermsreaction = max(maxtermsreaction+1,6)  ! for single interpolation
      if(Sectionalmethod .and. nHsections > 1)  maxtermsreaction = max(maxtermsreaction+3,6)  ! for double interpolation
      if(Sectionalmethod .and. nSsections > 1)  maxtermsreaction = max(maxtermsreaction+4,6)  ! for triple interpolation

      maxReactLines = ireact 

!----------------------------------------

      ! Gamma data if present
      if(SectionalMethod) call READGAMMA
      
!--------------------------------------------
      return
999   call errorline(fileSCHEME,' 7    ')

      end
!=======================================================================

      subroutine REACTIN(icall,ierrsolve)

! 1st call just to determine array size nreactarray with omitted reactions to replace large maxreactlines
! 2nd call to read reactions
!=====================================================================


!          MOLE COEFFICIENTS   a.A + b.B + c.C => d.D + e.E + f.F

!Rate constant   k =  A  * T^n * exp (-E/RT)
!Molar reaction rate  r =  k  *  [A]^a * [B]^b  
!  if M present molar reaction rate  r =  k  *  [A]^a * [B]^b * [M]  

! Units to satisfy:
!          r   (mol/cm^3-s)
!         [A]  (mol/cm^3)


!Reaction type 0  Chem. rate only:  
!              1  EBU only (infinite Chem. rate)
!              2  min( EBU, Chem. rate)
                 

! Elementary reactions if molar coeffs. a=0 and b=0 
! Global reactions for non-zero a, b use Tabled a and b instead of molar coefficients (not for reversible). 

! Non-reversible       use =>  eg.   a.A + b.B  => c.C  + d.D
! Reversible reaction: use =   eg.   a.A + b.B  =  c.C  + d.D


! Equilibrium constant Table needed if Reversible.
!  T(K)  v's ln K data  


! For P-dependent reactions, brackets (M) are OK and are converted to blanks
! TROE needs 4 parameters; 4th is made a large +ve number when reading data if only 3 parameters supplied.

!==================================================================

      use paramsmod
      use dimensmod
      use filesmod
      use reactmod
      use propsmod
      use logicmod
      use speciesmod
      use fuelmod
      use sectionalmod
                  
      common / tomonout / nreaceq
      character (LEN=LineL) :: Line, Line1
      real    :: reactexp(3)
      real(8) :: A_3rd(max3rdbodyperline), A
      integer :: nq_3rd(max3rdbodyperline)
      integer :: iterm(9,3), nHSTlimited(3)
      logical :: is_SectionReaction              
      nreaceq = 0
!---------------------------------------------------------------------
! READ 
      ioutput = icall - 1  ! output only for 2nd call
      
      call findheader(fileSCHEME,isfound,'REACTION_SECTION')
      if(isfound == 0) stop

      call findots(fileSCHEME)
      !  Units for E: 1 = kJ/kmol;  2 = kJ/mol; 3 = kcal/mol; 4 = cal/mol      
      read(fileSCHEME,*,err=999)  iunitE
      call ILimitin(iunitE,1,4,"E units         ")

      read(fileSCHEME,*,err=999)  Tfreeze

      read(fileSCHEME,*,err=999) isReactionAnalysis   ! analysis output

      read(fileSCHEME,*,err=999) xsHcriterion  ! abs(xs H) criterion to omit reaction 
      xsHcriterion = max(xsHcriterion, 1.0e-6)  ! to avoid knocking out on rounding errors        

      read(fileSCHEME,*,err=999) omitcriterion   ! criterion omit reporting reaction in ranking if species production or consumption < criterion*|max rate|
      read(fileSCHEME,*,err=999) significantcrit !  significance criterion for reporting.  
      read(fileSCHEME,*,err=999) T_gamma_out     ! T(K) at which to report gamma in reaction.txt
      read(fileSCHEME,*,err=999) eltolerance     ! Element error tolerance for all reactions which will generate message but no Stop

!---------------------------------

      call findots(fileSCHEME)
      if(icall==2) then
      write(filereact,'(/20x," ECHO REACTION INPUT DATA  ")')
      write(filereact,'(/" Non-sectional echo is input image. Stored Data below has interpreter reactions. ")')
      write(filereact,'(/10x," INPUT UNITS for E are  ",a10)') E_units(iunitE)
      write(filereact,'(/10x," H mol error criterion to omit Sectional reaction is  ",1pg10.1)') xsHcriterion
      write(filereact,'( 10x," Excess H mols  xsH = Hproducts - Hreactants")')
      write(filereact,'(/10x," Element mol error tolerance all reactions except H in Sectional   ",1pg10.3)') eltolerance
      write(filereact,'( 10x," ### reactions exceed element mol error tolerance.   ")') 
      write(filereact,'(/"                 Rate mol/(cm^3 - s) = A * T^n * exp(-E/RT)")')   
      write(filereact,'(" For coagulation  Rate mol/(cm^3 - s) = gamma(T) * A * T^n * exp(-E/RT)")')

!@      write(*,*) 'Reading Reactions ' 

      endif

!-----------------------------------
      ierr_element = 0; ierrsolve = 0; interperr = 0 
      
      expon = 0.0;  nq3rd = 0;  iterm3rd = 0 
      Eqnimage   = ' '
      Mpresent   = 0

      nqspec = 0     ! nqspec(3) does not exist because 3rd body is put in nq3rd
      specnetmass = 0.0
                   
      iAuxIndex = 0;  iAuxType  = 0;   Auxrate = 0.0;  iaux = 0
      indxexpon = 0

      is_SectionReaction = .false.
      ifirstsectionline = 0
      iGenericSectn = 0 
      ierrAisearch = 0
      numbereact = 0
      linenumber = 0
      n = 0

      nHSTlimited = 0
      nlimitCexceeded = 0                                                  
      nlimitCunderflow = 0
      nlimitHexceeded = 0                                                  
      nsectreactaccepted = 0                                                  
      nelementerrexceeded = 0
!-------------------------------
! REACTION LOOP

      DO WHILE (.true.)

      Line = ' '
      call READLINE(fileSCHEME,line,iend,Linechange)   ! read reaction line
      Linenumber = Linenumber + Linechange
      if(iend == 1) exit
      call NOLEADINGBLANKS(Line)
      
      ! Special exponents - indxexpon is index of specexpon corresponding to Reaction Line n
      call READFORD(Line,n,indxexpon,iscycle)
      if(iscycle==1 .and. icall==2) write(filereact,'(4x,i8,5x,(a))') numbereact(n), Line         ! echo input
      if(iscycle==1) cycle                      ! no new line No. for FORD

      ! Auxiliary data  'LOW' or 'HIGH' or 'TROE' or 'SRI'
      call READAUXDATA(Line,n,iaux,iscycle)
      if(iscycle==1 .and. icall==2) write(filereact,'(4x,i8,5x,(a))') numbereact(n), Line         ! echo input

      !   insert /   / for Chemkin
      if(iChemkin==1 .and. iscycle==1 .and. icall==2) then
      Line1 = Line
      do i = 2,8
       if(Line1(i:i) == ' ') then
       Line1(i:i) = '/'
       exit 
       endif             
      enddo
      L = len_trim(Line1)
      Line1(L+1:L+2) = ' /'
      endif
         
      if(iscycle==1) cycle                      ! no new line No. for AUX


      !-------------------------------------------

      ! Look for 3rd body continuation line
      call READ3rdBODYLINE(N3bodyinline,A_3rd,nq_3rd,Line)

      if(N3bodyinline > 0) then 
       if(Mpresent == 0) then      ! 3rd body continuation must follow [M] reaction
       write(file_err,*)' 3rd body efficiency must follow [M] reaction ', numbereact(n)
       stop
       endif        
      
      do i = 1, N3bodyinline
      n = n + 1                       ! increment Reaction Line No. to get one 3rd body species per line       
       if(n > maxReactLines) then
       write(file_err,*)' ERROR: Need to Increase maxReactLines dimension ',maxReactLines
       stop
       endif
      numbereact(n) = numbereact(n-1)       ! 3rd body variants use same reaction number 
      ireactype(n)  = -1                    ! to mark 3rd body continuation 
      nq3rd(n)      = nq_3rd(i)                
      chemrate(1,n) = A_3rd(i)              ! contains multiplier for continuation 3rd body reactions
      reversible(n) = reversible(n-1)       ! no reversible for 3rd body reactions
      enddo
      
      if(icall==2) write(filereact,'(4x,i8,5x,(a))') numbereact(n), Line         ! echo input
      eqnimage(n) = Line
      cycle
      endif
      !-------------------------------------------

      ! n is reaction line 
      n = n + 1      ! increment Reaction Line No.       

       if(n > maxReactLines) then
       write(file_err,*)' Increase maxReactLines dimension '
       stop
       endif
      
      ! Check whether Sectional Reaction by presence of Ai or Aj
      if(Sectionalmethod) then 
      call Aisearch(Line,1,maxreactcols,nAi,iterm,ierrAisearch) 
       if(nAi > 0) then
       is_SectionReaction = .true.
       iGenericSectn = iGenericSectn + 1
       endif      
       
       if(is_SectionReaction .and. ifirstsectionline == 0) then
       LastgasphasereactNo = numbereact(n-1)  
       ifirstsectionline = n
        ! ALLOCATION       
        if(icall==1) then
        allocate (GenericSectionimage(nGenericSectReact),LineRangeSectn(nGenericSectReact,2), stat=istat) 
        call STATCHECK(istat,'Sectionim ')        
        endif
       GenericSectionimage = ' '
       endif
      
       ! To group all section reactions because reaction coeffs are different to elementary reactions.  
       if(.not. is_SectionReaction .and. ifirstsectionline > 0) then
       write(file_err,'(/" Stop: following non-sectional reaction must precede section reactions. ")')
       write(file_err,*) Line
       stop       
       endif

       ! cannot have Ai*n or Aj*n with n>nHsections 
       if(iterm(5,1) > nHsections .or. iterm(5,2) > nHsections) then
       write(file_err,'(/" Stop: following sectional refers to H section > No. of H sections = ",i3)') nHsections
       write(file_err,*) Line
       stop       
       endif
      endif

      ! Interpret Reaction Equation; does not see continuation lines
      call INTERPRETER(is_SectionReaction,Line,n,Reactmol(0:nqlastmp1),ordermols(:,n),nqorderSp(:,n), &
                        reversible(n),nq3rd(n), iterm3rd(1:2,n),Lstart,interperr)

      !UPDATE REACTION NUMBER
      numbereact(n) = numbereact(n-1) + 1     

      ! Associate Tag No. with Reaction Number and store
      if(icall==2) call TAGNUMBER(line,numbereact(n),iitag)
      
      if(is_SectionReaction) call SectionReactCheck(nq3rd(n),iterm3rd(1:2,n),reversible(n),Line)

      !------------------------------                 
      ! Type-dependant read rate data 
      if(interperr == 0) &
      call READRATECOEFFS(iallelementary,is_SectionReaction,Line,Lstart,A,Tn,AEin,reactexp,Cexp,Cfac,iCtype,S_X,ifragment)

      ! Check for -ve A in non-Section reaction
      if(.not. is_SectionReaction .and.  A <= 0.0) then
      write(file_err,*) ' Stop: negative rate constant in following non-section reaction. ', A
      write(file_err,*) Line
      stop
      endif 
      !------------------------------                 
      ! Save principal reaction equation [M] (sum of species) for any following specific 3rd body reactions. 
      Mpresent = 0
      if(nq3rd(n)==nqM) Mpresent  = 1
      
      !-----------------------------------                                    
      if(interperr == 1) cycle       ! error or species missing from specification list

      !-----------------------------------                                    
      ! Sectional Method Reaction Generator if Ai or Aj found anywhere in reaction 
      ! if Ai, Aj only in product they are converted to Ai_ Aj_ and there is only one reaction.   
      ! Section reactions use separate rate coeff. list to elementary reactions.
      
      if(is_SectionReaction) then
      LineRangeSectn(iGenericSectn,1) = n   ! first line No.

      call SECTIONREACTIONS(icall,ioutput,iitag,iGenericSectn,n,Line,reversible(n),A,Tn,AEin,Cexp,Cfac,iCtype,S_X,ifragment, &
                            nHSTlimited,nlimitCexceeded,nlimitCunderflow,nlimitHexceeded,nsectreactaccepted,iterm)
      LineRangeSectn(iGenericSectn,2) = n

      cycle    ! back to read next reaction
      endif
      
      !-----------------------------------                                    
      ! Continue for non-Sectional reactions

      IF(icall==2) THEN
      ! Elements balance Check 
      ierr_n = 0      
      call REACTIONELEMENTCHECK(n,ierr_n,xsC,xsH)
      ierr_element = max(ierr_n, ierr_element)  
      if(ierr_n > 0) nelementerrexceeded = nelementerrexceeded + 1
      call STOREREACTCOEFFS(n)

      ! echo input
      if(mod(n+19,20) == 0) write(filereact,'(/"Tag React.No.",40x,"A         n           E")')      

        if(ierr_n > 0) then      
        write(filereact,'("!###",i8,5x,(a))') numbereact(n), Line
        elseif(itag(numbereact(n)) > 0) then
        write(filereact,'(i3,1x,i8,5x,(a))') itag(numbereact(n)), numbereact(n), Line
        else
        write(filereact,'(4x,i8,5x,(a))') numbereact(n), Line
        endif
       ! OMIT ### reactions in file Chemkinout

      eqnimage(n) = Line     ! defines Eqnimage
     
      ! RATES storage: reactionline not reaction number to store multipliers for 3rd body
      chemrate(1,n) = A                    ! contains multiplier for continuation 3rd body reactions
      chemrate(2,n) = Tn
      chemrate(3,n) = Econversion(AEin,iunitE)    ! convert activation energy to kJ/kmol

      ! Species must be Solved if in a reaction.
      call isSOLVEDCHECK(n,ierrsolve)
      ENDIF

      ENDDO  ! END REACTION No. LOOP

!------------------------------------------------------------------------


      if(ierrAisearch /=0) then
      write(file_err,*)' Stop for Sectional Species specification error(s) in above list '
      stop
      endif

      if(icall==1) then
      nreactarray = n + 1
      else
      nreactionlines = n 
      endif
      
      if(icall==1) RETURN

      ! net moles for all species in each reaction  (for Kequil)     
      do n = 1, nreactionlines 
      reactmolnet(n) = sum( ordermols(1:maxtermsreaction,n) )   
      enddo
!= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
! OUTPUTS
      ! Sectional Ranges
      write(filereact, '(/1x,79(1h-))')
      write(filereact,'(/ " Reaction No. Range for Generic Sectional Reactions")') 
      do i = 1,nGenericSectReact
      Nr1 = numbereact(LineRangeSectn(i,1))
      Nr2 = numbereact(LineRangeSectn(i,2))
      write(filereact,'(5x,2i8,3x,i4,3x,(a))') Nr1,Nr2,i,GenericSectionimage(i)       
      enddo
      write(filereact, '(/1x,79(1h-))')
      
      call REACTION_DATA_OUTPUT(ifirstsectionline)   ! gas-phase

      ! Non sectional Repeated reactions check  ----------
      irepeaterr = 0
      nlast = nreactionlines
      if(ifirstsectionline > 0) nlast = ifirstsectionline - 1
      call REPEATCHECK(irepeaterr,nlast)
      
       if(irepeaterr /= 0) then      
       write(filereact,'(/" Non-Sectional Repeated reactions. ")')
       write(filereact,'(" First appearance: * ")')
       write(filereact,'(" = is same reaction ")')
       write(filereact,'(" + is same but different reaction rate or different 3rd body species ")')
       write(filereact,'(" S is Specific 3rd body, or M with all species"/)')
       do n = 1,nlast
       if(markreact(n) /= '   ') write(filereact,'(i6,2x,a3,5x,(a))') numbereact(n), markreact(n), eqnimage(n)
       enddo
       endif
                   
       if(irepeaterr==2) then
       write(file_err,*) ' STOP for Repeated Reaction error. See REACT.txt '      
       stop      
       endif
      ! ----------------------------------------------------
      
      write(filereact, '(/1x,79(1h-))')
      write(filereact,'(/ " =========== SUMMARY of REACTIONS ===========  ")')
      write(filereact,'(  " Separate Reaction No. for each gas-phase reaction, but not for additional 3rd body reactions. ")')
      write(filereact,'(  " Separate Reaction No. for each Sectional permutation. ")')
      if(nelementerrexceeded > 0) write(filereact,'(/" No. of reactions marked ### because element error > tolerance      ",i8)') nelementerrexceeded

      if(is_SectionReaction) then
      write(filereact,'(/ " Omitted Sectional Reactions ")')
      if(nlimitCexceeded > 0)    write(filereact,'(  " No. of reactions omitted because C > Last section C is             ",i8)') nlimitCexceeded
      if(nlimitCunderflow > 0)   write(filereact,'(  " No. of reactions omitted because C < First section C is            ",i8)') nlimitCunderflow
      if(nlimitHexceeded > 0)    write(filereact,'(  " No. of sectional reactions omitted because abs(xsH) > criterion is ",i8)') nlimitHexceeded
      if(nHSTlimited(1) > 0)     write(filereact,'(  " No. of reactions omitted because no H for section C is             ",i8)') nHSTlimited(1)
      if(nHSTlimited(2) > 0)     write(filereact,'(  " No. of reactions omitted because no S for section C is             ",i8)') nHSTlimited(2)
      if(nHSTlimited(3) > 0)     write(filereact,'(  " No. of reactions omitted because no T for section C is             ",i8)') nHSTlimited(3)
      nTot_omitted = sum(nHSTlimited) + nlimitCexceeded + nlimitCunderflow + nlimitHexceeded 
                              write(filereact,'(  " No. of sectional reactions omitted                                       ",i8)') ntot_omitted 
                              write(filereact,'(  " No. of sectional reactions accepted                                      ",i8)') nsectreactaccepted
      endif
       
      if(is_SectionReaction) write(filereact,'(/ " Last gas-phase reaction No.          ",i8)')  LastgasphasereactNo
      write(filereact,'(                         " Last Reaction No.                    ",i8)') numbereact(nreactionlines)

      write(filereact,'(/ " Total Number of Reactions is                              ",i8)') numbereact(nreactionlines)
      write(filereact,'(  " No. of all reactions lines accepted (array size)          ",i8)') nreactionlines
      write(filereact,'(  " Total Number of Reaction lines incl. 3rd body repeats is  ",i8)') nreactionlines
      write(filereact, '(/1x,79(1h-))')
        
      if(interperr /= 0) then
      write(file_err,*)' Stop for missing species or reaction error in above list '
      stop
      endif

       if(ierr_element /= 0) then
       write(file_err,*)' Stop for Element Balance error. See reactions marked ### '
       stop
       endif

!----------------------------------------

! Save Reaction Line No. from Reaction Number
      LineReact = 0
      do iline = 1, nreactionlines
      NR = numbereact(iline)   ! reaction No.
      Linereact(NR) = iline    ! Line No.   
      enddo

      call REACTIONELEMENTRATE   ! for C rate output
      
!----------------------------------------
!  save volatiles for VOLSPECIES
!        do i = nqsp1, nqLast
!        volsmol(i) = reactmol(i)      
!        enddo

!-----------------------------------------
! Check uni-molecular for elementary reactions
      ielementerr = 0
      if(iallelementary ==1) call ELEMENTARYMOLCHECK(ielementerr)
      if(ielementerr==1) stop      

!-----------------------------------------


      return
999   call errorline(fileSCHEME,' 6    ')

998   write(file_err,*) ' Error opening Reaction.in '
      stop

      end
!===================================================================      

      subroutine INTERPRETER(is_SectionReaction,Line,nreac,Reactmol,ordermol,nqorderSp,reversible, &
                              nq3rd,iterm3rd,Lstart,ierr)

!=========================================================================
      use paramsmod
      use dimensmod
      use filesmod
      use speciesmod

      character (LEN=LineL)  :: Line
      character (LEN=1)   :: chr1, chrp, chrsign
      character (LEN=80) ::  temp1
      character (LEN=24) ::  chrsymbol

      real(8) :: reactmol(0:nqlastmp1)           ! add 1 for M
      real(8) :: readmol(maxtermsreaction), ordermol(maxtermsreaction)
      integer :: nqSpecies(maxtermsreaction), nqorderSp(maxtermsreaction), iterm3rd(2)
      logical :: reversible
      logical :: is_SectionReaction              

!-----------------------------------
! Returns:
! Reactmol sums net mols for each species; -ve for reactants; +ve for products; 0 for 3rd body  eg. 0 0 0 0 -2 0 0 1 0 1
! Ordermol has mols for each of maxtermsreaction places: 3 LHS -ve or 0; 4 or more for RHS +ve or 0
! nqorderSp has nq's for each of maxtermsreaction places or 0
! reversible   true or false
! nq3rd is nq of 3rd body or 0
! iterm3rd(2) is term numbers of 3rd body
! Lstart is column to start reading subsequent reaction rate data
! ierr = 1 if species error


      reactmol   = 0.0
      readmol    = 0.0
      ordermol   = 0.0 
      nqSpecies    = 0       ! eg H + H = H2   20  20 35   0 0 0     no gaps            
      nqorderSp  = 0         ! eg H + H = H2   20  20  0   35 0 0    species for each place; 1-3 reactants; 4-7 products
      reversible = .true.

! Scan

      chr1 = ' '
      chrsign = '-'
      nn = 0 
      nt1 = 1
      temp1 = ' '
      nentry = 0      
      ns = 0
      itosymbol = 0
      inq = 0
      iorder = 0
      ichange=0
      nnchr = -1
                      ! number is No. of molecules of species
                      ! symbol is alpha and numerals in species formula
      inumber   = 0   ! =1 from start of number to finish of symbol
      isymbol   = 0   ! =1 from start of symbol to end of symbol
      kountnext = 1   ! count '+' to check for number of entries per side (species following + will be kountnext)
      iterminated = 0
      icolsymbol = 0

! brackets are made blanks
      do i = 1,maxreactcols       
      if(Line(i:i)  ==  '(' .or. Line(i:i)  ==  ')' ) Line(i:i)  =  ' ' 
      enddo 
               
      DO while (nn <= maxreactcols)
      nn = nn+1
      chrp = chr1
      chr1 = Line(nn:nn)
      if(Line(nn:nn)  /=  ' ' .and. Line(nn:nn)  /=  char(9) ) Lstart = nn+2   ! for reading subsequent rate data 

! exit 
      call EXITREACTION(Line(nn:nn+1),nn,iterminated)

      if(iterminated == 1)then  
       if(nn <= icolsymbol + 1) then ! to allow a complete cycle after last symbol is completed
       write(file_err,'(/" Reaction is truncated. Need extra space between reaction and rate data. ")')
       write(file_err,*) line
       stop
       endif
      exit
      endif

!--------------------------
! change sign after '='
      if(chr1 == '=') then
      chrsign = ' '
      kountnext = 1        ! reset counter
      ichange=1
      nnchr = nn
      elseif(chr1 == '+') then
      kountnext = kountnext + 1 
      endif
       
! Non-reversible if =>   otherwise > is part of sectional species name
      if(chr1 == '>') then
       if(nn == nnchr+1) then
       reversible = .false.
       endif
      endif


      istartsymbol = 0
      istartnumber = 0
! Number of moles    (? is for sectional product H or H2 for automatic H balance)
       if(isymbol == 0 .and. ((chr1 >= '0' .and. chr1 <= '9') .or. chr1 == '.' .or. chr1 == '?') ) then
       if(is_SectionReaction .and. ichange==1 .and. chr1 == '?') chr1 = '0'  ! change ? to 0 (for Sectional H or H2 product)  
       inumber = 1
       if((chrp < '0' .or. chrp > '9') .and. chrp /= '.') istartnumber = 1 

! end of species symbol
       elseif(isymbol == 1 .and. (chr1 == ' ' .or. chr1 == char(9) .or. chr1 == '+' .or. chr1 == '=')) then 
       isymbol = 0
       inumber = 0

        ! too many entries: 3 max for reactants 
        if(ichange==0 .and. kountnext > 3) then
        write(file_err,*)' Too many species listed in reactant side of reaction. '
        write(file_err,*) line
        ierr = 1 

        ! max products if irreversible
        elseif(ichange==1 .and. (.not. reversible) .and. kountnext > maxtermsreaction-3) then
        write(file_err,*)' Max. No. of species in product side of irreversible reaction is ',maxtermsreaction-3
        write(file_err,*) line
        ierr = 1 

        ! 3 max products if reversible
        elseif(ichange==1 .and. reversible .and. kountnext > 3) then
        write(file_err,*)' Max. 3 species in product side of reversible reaction. '
        write(file_err,*) line
        ierr = 1 
        endif
       
! begin species symbol A - Z or a - z
       elseif( (chr1 >= 'A' .and. chr1 <= 'Z') .or. (chr1 >= 'a' .and. chr1 <= 'z') ) then
       if(isymbol == 0) istartsymbol=1
       isymbol = 1
       endif


!--------------------------
! write mole numbers to temp1

!    one mole inserted
      if(istartsymbol == 1 .and. inumber == 0) then
      nt1 = nt1+3
      temp1(nt1:nt1) = chrsign
      nt1=nt1+1
      temp1(nt1:nt1) = '1'
      nentry = nentry +1

      elseif(inumber == 1.and.isymbol == 0) then
       if(istartnumber == 1) then
       nt1=nt1+3
       temp1(nt1:nt1) = chrsign
       nentry = nentry +1
       endif
      nt1=nt1+1
      temp1(nt1:nt1) = chr1
      endif

!--------------------------

! write to symbol

      if(isymbol == 1) then
       if(istartsymbol == 1) then
       ns = 0
       chrsymbol = ' '
       endif
      ns = ns+1
      chrsymbol(ns:ns) = chr1
      itosymbol = 1      
      icolsymbol = nn   ! to mark column of end of symbol


!--------------------------
! Find Species index No.

      elseif(itosymbol == 1) then
      itosymbol = 0      
      nq = IndexSpecies(chrsymbol)

! nqorderSp
      iorder = iorder + 1
      nqorderSp(iorder) = nq               ! output
              
       if(nq > 0) then
       inq=inq+1
       nqSpecies(inq) = nq
       else
       write(file_err,*) line(1:maxreactcols),'   has not-listed species  ',chrsymbol
       ierr = 1 
       endif 

      endif
      
      if(iorder < 3 .and. ichange==1) iorder = 3    ! after '=' or blank, start products
      ENDDO

!--------------------------
! End of reaction not found by maxreactcols
      if(iterminated == 0) then
      write(file_err,*)' End of reaction has not been found by column ', maxreactcols
      write(file_err,*) Line
      ierr = 1 
      return
!      stop
      endif
       
!--------------------------

      read(temp1,fmt=*,err=999) readmol(1:nentry)

! reactmol    output
      do i = 1,nentry
      n = nqSpecies(i) 
      reactmol(n) = reactmol(n) + readmol(i)  ! accumulate mole numbers for each species
      enddo

! ordermol output into maxtermsreaction number of places
      j = 0 
      do i = 1,nentry
       do while(j < maxtermsreaction)
       j = j + 1
        if(nqorderSp(j) > 0) then
        ordermol(j) = readmol(i)           ! output
        exit
        endif
       enddo
      enddo

! nq3rd output:   if any reactant species moles = same species moles in product
       nq3rd = 0
outer: do i = 1,3
       nqR = nqorderSp(i)            
       if(nqR == 0) cycle
        do j = 4,maxtermsreaction      ! 3 terms reserved for reactants side
        nqP = nqorderSp(j)           
        if(nqP == 0) cycle
         if(nqP == nqR .and. abs(ordermol(i)+ordermol(j)) < rsmall) then
         nq3rd = nqR
         iterm3rd(1) = i
         iterm3rd(2) = j
         exit outer
         endif
        enddo        
       enddo outer

! check if M is not on both reactants and products side
       Mreact= 0
       Mprod = 0
       do i = 1,3
       if(nqorderSp(i) == nqM) Mreact = 1            
       enddo  
       do i = 4,maxtermsreaction     
       if(nqorderSp(i) == nqM) Mprod = 1            
       enddo  
   
       if(Mreact /= Mprod) then
       write(file_err,*)' M does not appear on both sides of the reaction '
       write(file_err,*) Line
       ierr = 1 
       endif


      return

999   write(file_err,*)' Error reading reaction: Illegal character, or following data intrudes. '
      write(file_err,*) line
      stop

      end

!=================================================================

      subroutine READ3rdBODYLINE(N3bodyinline,A_3rd,nq_3rd,Line)

! For continuation line with specific 3rd bodies and rate multipliers
! eg.      H2/ 2.40/H2O/15.40/CH4/ 2.00/CO/ 1.75/CO2/ 3.60/C2H6/ 3.00/Ar/  .83/

!==================================================================================================

      use paramsmod
      use REACTMOD
      use filesmod
      use dimensmod

      real(8) :: A_3rd(max3rdbodyperline)
      integer :: nq_3rd(max3rdbodyperline)
      character (LEN=LineL) :: Line, nameline, numberline
      character (LEN=24)  :: speciesnames(max3rdbodyperline)

      N3bodyinline = 0 

      ! Count  / for No. of species in line
      nslash = 0
      do i = 1,LineL        ! search for /
      if(Line(i:i) == '!') exit            ! avoid any /'s after !
      if(Line(i:i) == '/') nslash = nslash + 1
      enddo
      N3bodyinline = (nslash+1)/2     ! allows for missing end /

        if(N3bodyinline > max3rdbodyperline) then
        write(file_err,*)' STOP: No. of 3rd body species ',N3bodyinline,' exceeds max3rdbodyperline ',max3rdbodyperline
        write(file_err,*) Line
        stop
        endif

      if(N3bodyinline == 0) return 
!-----------

      ! Read names and values
        ! make nameline only the species names
        ! make numberline only the values
      nameline   = Line
      numberline = Line 

      isname = 1
      do i = 1,LineL       
      if(Line(i:i) == '!') exit            ! avoid any /'s after !
      if(Line(i:i) == '/') isname = -1 * isname
      if(nameLine(i:i)   == '/')  nameLine(i:i)   = ' '       ! / becomes blank 
      if(numberline(i:i) == '/')  numberline(i:i) = ' '       ! / becomes blank 
      if(isname ==  1)  numberline(i:i) = ' '                 ! name becomes blank
      if(isname == -1)  nameLine(i:i) = ' '                   ! number becomes blank
      enddo


      call READMULTISPECIES(nameLine,Nspeciesread,speciesnames,N3bodyinline,24,lastchar)            
      read(numberline,fmt=*,err=999,end=999) A_3rd(1:N3bodyinline)

      do i = 1, N3bodyinline
      nq_3rd(i) = Indexspecies(speciesnames(i))
        if(nq_3rd(i) == 0) then
        write(file_err,*)' 3rd body species is not in list ',speciesnames(i)
        write(file_err,*) Line
        stop
        endif
      
      A_3rd(i) = A_3rd(i) - 1.0    ! -1 because have already included M = all species previously
      enddo
             
      return
999   call errorline(fileSCHEME,' 61  ')

      end
!==================================================================================================

      subroutine READRATECOEFFS(iallelementary,is_SectionReaction,Line,Lstart,A,Tn,Eact,reactexp,Cexp,Cfac,iCtype,S_X,ifragment)

! Read reaction rates
!==================================================================================================

      use paramsmod
      use filesmod
      use sectionalmod
      
      character (LEN=LineL) :: Line
      real(8) :: A
      real :: reactexp(3)
      logical :: is_SectionReaction              


      Cexp=0.0; Cfac=0.0; iCtype = 0; S_X=0 
      reactexp = 0.0

       if(is_SectionReaction) then       ! Sectionmethod reaction
       read(Line(Lstart:LineL),fmt=*,err=999,end=999) A, Tn, Eact, Cexp, Cfac, iCtype,S_X       

! Fragmentation if S_X < 0.0; then put S_X = 0
        iFragment = 0  
        if(S_X < 0.0) then
        iFragment= abs(nint(S_X))  ! =1 if S_X = -1.0,  =2 if S_X = -2.0
        S_X = 0.0 
        endif  
     
        if( (S_X<100.0 .and. S_X > real(nSsections+0.0001)) .or. (S_X>100.0 .and. S_X-100.0 > real(nSsections+0.0001)) .or. S_X < 0.0) then
        write(file_err,*) ' Need S_X or S_X-100 >= 0 and <= No. S sections but is ', S_X
        write(file_err,*) Line
        stop        
        endif

       elseif(iallelementary == 1) then     ! all elementary reactions
       read(Line(Lstart:LineL),fmt=*,err=999,end=999) A, Tn, Eact

       else                                 ! Global reactions: Exponents are read in order of reactants in equation. 
       read(Line(Lstart:LineL),fmt=*,err=999) A, Tn, Eact, reactexp(1:3)
       endif

      return
999   write(file_err,*) ' Wrong Reaction rate data format or need >= 5 spaces after reaction. '
      call errorline(fileSCHEME,' 6a   ')

      end
      
!==================================================================================================

      subroutine REACTIONELEMENTCHECK(n,ierr,xsC,xsH)

!==================================================================================================

      use dimensmod
      use paramsmod
      use reactmod
      use speciesmod
      use filesmod
      use sectionalmod
      use logicmod
      use elementsmod
            


      real(8) :: Rcarbon, Rhydrogen

                              
! Elements Mole Check 
       eltol2 = eltolerance / 1000.0
       xs_element = 0.0
       Rcarbon   = 0.0
       Rhydrogen = 0.0
       xsC = 0.0
       xsH = 0.0
       ierrC = 0
       ierrH = 0
       ierrX = 0

       if(n==1) write(file_err,'(/" Element mole error > ",1pg12.3," or element mol error/reactant mols > ",1pg12.3," in following"/)') eltolerance,eltol2 
       
       isSectional = 0               
       do i = 1,maxtermsreaction         
       if(nqorderSp(i,n) == 0) cycle
       if(nqorderSp(i,n) > nqLast) cycle     ! up to  & including Ai, Aj
       if(nqorderSp(i,n) >= nqSectionmol1  .and. nqorderSp(i,n) <= nqSectionN) isSectional = 1 
       
        ! net sum: +ve means Cprod > Creactants
        do ie = 1, maxelement
        xs_element(ie) = xs_element(ie) + ordermols(i,n) * spCHONSI(ie,nqorderSp(i,n))
        enddo

        if(i <= 3) then     ! Reactants sum
        Rcarbon   = Rcarbon   + ordermols(i,n) * spCHONSI(1,nqorderSp(i,n))
        Rhydrogen = Rhydrogen + ordermols(i,n) * spCHONSI(2,nqorderSp(i,n))
        endif
       enddo

       ! +ve xsC = Cproducts - Creactants   
       xsC = xs_element(1)
       xsH = xs_element(2)

       ! Carbon
       if( (abs(Rcarbon) < 1000.0 .and. abs(xsC) > eltolerance) .or. (abs(Rcarbon) >= 1000.0 .and. abs(xsC)/(abs(Rcarbon)+rtiny) > eltol2) ) then
       write(file_err,'(" React.",i5,"    C react. mols ",g12.3," Excess C mols ",g12.3)') numbereact(n),-Rcarbon,xsC
       ierrC = 1
       endif

       ! Hydrogen
       if( (abs(Rhydrogen) < 1000.0 .and. abs(xsH) > eltolerance) .or. (abs(Rhydrogen) >= 1000.0 .and. abs(xsH)/(abs(Rhydrogen)+rtiny) > eltol2) ) then
       ! for Sectional reactions xsH is shown in Reaction.txt output; no Stop-error for H        
        if(.not. (Sectionalmethod .and. isSectional == 1)) then
        write(file_err,'(" React.",i5,"    H react. mols ",g12.3," Excess H mols ",g12.3)')numbereact(n),-Rhydrogen,xsH
        ierrH = 1
        endif
       endif

      ! Other elements      
      do ie = 3, maxelement
       if(abs(xs_element(ie)) > eltolerance) then
       write(file_err,'(" Reaction ",i5, "    Excess ",a2," mols ",g12.3)') numbereact(n),ElementName(ie),xs_element(ie)
       ierrX = 1
       endif
      enddo 

      ierr = max(ierrC, ierrH, ierrX)      
            
      return
      end
!===============================================================                              


      subroutine REPEATCHECK(ierr,nlast)

      use filesmod
      use reactmod
      
      logical neqal8
      integer :: ir1in(3), ir2in(3)  
                        
! Check for repeated reaction 
      write(filereact,'(" Non-sectional Repeated Reactions ")')

      do n = 2,nlast
       if(numbereact(n) == numbereact(n-1)) cycle       ! auxiliary not considered
       do i = 1,n-1
       match=1
       
       ir1in = nqorderSp(1:3,n)
       ir2in = nqorderSp(1:3,i)
       call REPEATREACTION(match, ir1in, ir2in,iterm3rd(1,n),  nq3rd(n),iterm3rd(1,i),  nq3rd(i),isM)
       ir1in = nqorderSp(4:6,n)
       ir2in = nqorderSp(4:6,i)
       call REPEATREACTION(match, ir1in, ir2in,iterm3rd(2,n)-3,nq3rd(n),iterm3rd(2,i)-3,nq3rd(i),idum)

       ! check rates: if same reaction with different rates match = 2
       if(match==1) then
       if(neqal8(chemrate(1,n),chemrate(1,i)) ) match = 2
       if(neqal8(chemrate(2,n),chemrate(2,i)) ) match = 2
       if(neqal8(chemrate(3,n),chemrate(3,i)) ) match = 2
       if(nq3rd(n) /= nq3rd(i)) match = 2
       endif 

       if(match==1) then        ! error
       write(filereact,'(/" Reaction ",i6, " is repeat of ",i6)')numbereact(n),numbereact(i)
       write(filereact,'(" Reaction ",i6,5x,(a))')numbereact(n),Eqnimage(n)
       write(filereact,'(" Reaction ",i6,5x,(a))')numbereact(i),Eqnimage(i)
       ierr = 2
       markreact(i)(1:1) = '*'    ! same reaction: error
       markreact(n)(2:2) = '='    
       
       elseif(match==2) then    ! not error, but reported
       write(filereact,'(/" Reaction ",i6, " same reaction different rate to ",i6)')numbereact(n),numbereact(i)
       write(filereact,'(" Reaction ",i6,5x,(a))')numbereact(n),Eqnimage(n)
       ierr = max(ierr,1) 
       markreact(i)(1:1) = '*'    ! different rate
       markreact(n)(2:2) = '+'    
       endif
       if(match /= 0) then
       if(isM == 1) markreact(n)(3:3) = 'M'    
       if(isM == 1) markreact(i)(3:3) = 'S'    
       if(isM == 2) markreact(n)(3:3) = 'S'    
       if(isM == 2) markreact(i)(3:3) = 'M'    
       endif 
       enddo
      enddo

! If reverse reaction is repeated (for reversible reactions)
      do n = 2,nlast
      if(numbereact(n) == numbereact(n-1)) cycle               ! auxiliary not considered
       do i = 1,n-1
       if(.not.reversible(n) .and. .not.reversible(i)) cycle   ! allowed if both irreversible 
       match=1

       ir1in = nqorderSp(1:3,n)
       ir2in = nqorderSp(4:6,i)
       call REPEATREACTION(match, ir1in, ir2in, iterm3rd(1,n), nq3rd(n), iterm3rd(2,i)-3,nq3rd(i),isM)
       ir1in = nqorderSp(4:6,n)
       ir2in = nqorderSp(1:3,i)
       call REPEATREACTION(match, ir1in,ir2in, iterm3rd(2,n)-3,nq3rd(n), iterm3rd(1,i),  nq3rd(i),idum)

       ! check rates: same reaction with different rates allowed
       if(match==1) then
       if(neqal8(chemrate(1,n),chemrate(1,i)) ) match = 2
       if(neqal8(chemrate(2,n),chemrate(2,i)) ) match = 2
       if(neqal8(chemrate(3,n),chemrate(3,i)) ) match = 2
       endif 

       if(match==1) then
       write(filereact,'(" Reaction ",i6, " is reverse of ",i6)')numbereact(n),numbereact(i)
       write(filereact,'(" Reaction ",i6,5x,(a))')numbereact(n),Eqnimage(n)
       write(filereact,'(" Reaction ",i6,5x,(a))')numbereact(i),Eqnimage(i)
       ierr = 2
       markreact(i)(1:1) = '*'  ! same reaction: error
       markreact(n)(2:2) = '='    
       elseif(match==2) then    ! not error, but reported
       write(filereact,'(/" Reaction ",i6, " same reaction different rate to ",i6)')numbereact(n),numbereact(i)
       write(filereact,'(" Reaction ",i6,5x,(a))')numbereact(n),Eqnimage(n)
       ierr = max(ierr,1) 
       markreact(i)(1:1) = '*'  ! different rate
       markreact(n)(2:2) = '+'    
       endif
       if(match /= 0) then
       if(isM == 1) markreact(n)(3:3) = 'M'    
       if(isM == 1) markreact(i)(3:3) = 'S'    
       if(isM == 2) markreact(n)(3:3) = 'S'    
       if(isM == 2) markreact(i)(3:3) = 'M'    
       endif 
       enddo


      enddo


      return
      end
!====================================================================      

      subroutine REPEATREACTION(match,ir1in,ir2in,iterm3rd1,nq3rd1,iterm3rd2,nq3rd2,isM)

!==================================================================================================
! checks whether ir1(1:3) = ir2(1:3) in any order

      use dimensmod      
      integer :: ir1(3), ir2(3), ir1in(3), ir2in(3)  

      if(match==0) return
      ir1 = ir1in
      ir2 = ir2in
      
      ! Check for reaction with specific 3rd body against same reaction with M
      ! (put 3rd body = M to check if match is made) 
      isM = 0
      if(iterm3rd1 > 0 .and. iterm3rd2 > 0) then    
       if(nq3rd1 == nqM .and. nq3rd2 /= nqM) then
       ir2(iterm3rd2) = nqM
       isM = 1
       elseif(nq3rd1 /= nqM .and. nq3rd2 == nqM) then
       ir1(iterm3rd1) = nqM
       isM = 2
       endif
      endif
      
      do k = 1,3       
        if(ir1(k)==0) cycle
        if(ir1(k) /= ir2(1) .and. ir1(k) /= ir2(2) .and. ir1(k) /= ir2(3) ) then
        match = 0
        endif
      enddo

      do k = 1,3       
        if(ir2(k)==0) cycle
        if(ir2(k) /= ir1(1) .and. ir2(k) /= ir1(2) .and. ir2(k) /= ir1(3) ) then
        match = 0
        endif
      enddo

      return
      end
!==================================================================================================

      subroutine isSOLVEDCHECK(n,ierr)

!====================================================================      
      use dimensmod
      use reactmod
      use logicmod
      use filesmod
      use speciesmod
      
       do i = 1, maxtermsreaction
       nq = nqorderSp(i,n)
        if(nq>0 .and. nq<=nqlast) then
         if(.not.solvenq(nq)) then
         write(file_err,'(" ERROR:  Reaction ",i5," needs species to be solved: ",a24)')numbereact(n), spname(nq)
         ierr = 1
         endif
        endif
       enddo

      return
      end
!=====================================================================            

      subroutine ELEMENTARYMOLCHECK(ierr)

! Check for uni-molecular elementary reactions.
! Check on reactant side and if reversible reaction also on product side.
! Requirement is so that exponents are not needed for elementary reactions
!=====================================================================            
      use paramsmod
      use filesmod
      use reactmod
                  

      do n = 1,nreactionlines
       do i = 1,maxtermsreaction
       if(i > 3 .and. .not.reversible(n)) exit   ! skip products if not reversible  

        if((abs(ordermols(i,n)) -1.0) > rsmall) then    ! moles /= 1
        indx =  iFORD(n)                                ! if special exponent = 0 then OK      
         if(indx > 0) then                
          if(itermexpon(i,indx) == 1) then
          if(specexpon(i,indx)==0.0) cycle
          endif
         endif

        write(file_err,'(" Write 2CxHy as CxHy + CxHy for reaction ",i6,5x,(a))')numbereact(n),Eqnimage(n)
        ierr = 1
        endif

       enddo
      enddo

      return
      end
!====================================================================      

      subroutine READAUXDATA(Line,n,iaux,iscycle)

      ! Auxiliary data recognised by 'LOW' or 'HIGH' or 'TROE' or 'SRI'
      ! LOW is low P limit rate data with high P limit in usual data location
      ! TROE is Troe interpolation method

!============================================================      
      use REACTMOD
      use filesmod
      use paramsmod            

      character (LEN=LineL) :: Line

      iscycle = 0

      ! remove any slashes
      if(Line(1:3)=='LOW' .or. Line(1:4)=='HIGH' .or. Line(1:4)=='TROE' .or. Line(1:3)=='SRI') then
      do i = 1,LineL    
      if(Line(i:i) == '/' ) Line(i:i) = ' '
      enddo
      endif
      
                   
      if(Line(1:3) == 'LOW' .or. Line(1:4) == 'HIGH') then
       if(numbereact(n) == 0) then
       write(file_err,*)' LOW or HIGH or TROE data line before first reaction data '
       write(file_err,*) Line
       stop
       endif

             
      iaux = iaux + 1                      ! index counter for LOW or HIGH 
      iAuxIndex( numbereact(n) ) = iaux    ! numbereact points at index of Auxrate
      
       if(Line(1:3) == 'LOW')  then
       iauxtype(iaux) = -1                     ! 3 parameters
       read(Line(5:LineL),fmt=*,err=999,end=999) Auxrate(1:2,iaux), AEin    ! read E as real4

       elseif(Line(1:4) == 'HIGH') then
       iauxtype(iaux) = 1                      ! 3 parameters
       read(Line(6:LineL),fmt=*,err=999,end=999) Auxrate(1:2,iaux), AEin
       endif

       if(Auxrate(1,iaux) == 0.0 ) then
       write(file_err,*)' LOW or HIGH data line is zero '
       write(file_err,*) Line
       stop
       endif

      Auxrate(3,iaux) = Econversion(AEin,iunitE)    ! convert activation energy to kJ/kmol & store as real8

      iscycle = 1
      endif



!-----------------------------------
      if(Line(1:4) == 'TROE')  then          ! TROE
       if(abs(iauxtype(iaux)) /= 1 ) then
       write(file_err,*)' LOW or HIGH data line needed before TROE data '
       write(file_err,*) Line
       stop
       endif

      read(Line(6:LineL),fmt=*,iostat=ios) Auxrate(4:7,iaux)       ! try to read 4 parameters 
      if(ios/=0) then  
      read(Line(6:LineL),fmt=*,err=999,end=999) Auxrate(4:6,iaux)  ! else read 3 parameters
      Auxrate(7,iaux) = 1.0e12                                   ! makes exp(-T**/T) = 0                         
      endif

      iauxtype(iaux) = 2 * iauxtype(iaux)      ! -2 or +2

       if(sum(Auxrate(4:7,iaux))==0.0) then
       write(file_err,*)' TROE data line is zero '
       write(file_err,*) Line
       stop
       endif
            
      iscycle = 1
      
      elseif(Line(1:3) == 'SRI')  then       ! SRI
       if(abs(iauxtype(iaux)) /= 1 ) then
       write(file_err,*)' LOW or HIGH data line needed before SRI data '
       write(file_err,*) Line
       stop
       endif

      read(Line(5:LineL),fmt=*,err=999,end=999) Auxrate(4:8,iaux)
      iauxtype(iaux) = 3 * iauxtype(iaux)      ! -3 or +3

       if(sum(Auxrate(5:8,iaux))==0.0) then
       write(file_err,*)' SRI data line is zero '
       write(file_err,*) Line
       stop
       endif
            
      iscycle = 1
      endif

      return

999   call errorline(fileSCHEME,' aux  ')

      end
!===========================================================                    

      subroutine EXITREACTION(chr,icol,isend)

!=========================================================================

! reaction equation has ended:
! if one space or tab followed by . or 0-9  
! if no preceding + or =  and column icol > 1
! symbol A - Z cancels + =  status

      use dimensmod
      character (LEN=2)  :: chr

      integer :: iconnect = 0
      save :: iconnect
      
      isend = 0
      
! new line
      if(icol == 1) then
      iconnect = 1     ! 1 to avoid end if col 1 is blank
      
! connector
      elseif(chr(1:1) == '+' .or. chr(1:1) == '=') then
      iconnect = 1        

! symbol cancels connector
      elseif( (chr(1:1) >= 'A' .and. chr(1:1) <= 'Z') .or. (chr(1:1) >= 'a' .and. chr(1:1) <= 'z')) then
      iconnect = 0        

! end if iconnect = 0 and space followed by . or 0-9 or -
      elseif( (chr(1:1) == ' ' .or. chr(1:1) == char(9)) .and. (iconnect == 0) &
         .and.  (chr(2:2)=='.' .or. chr(2:2)=='-' .or. (chr(2:2) >= '0' .and. chr(2:2) <= '9'))) then 
      isend = 1

! end if iconnect = 0 and space and last col to be scanned
      elseif((chr(1:1) == ' ' .or. chr(1:1) == char(9)) .and. iconnect==0 .and. icol==maxreactcols) then
      isend = 1

      endif        

      return
      end
!====================================================================== 

      subroutine READFORD(Line,nreac,indxexpon,iscycle)    

! Read FORD or RORD data; more than one FORD line accepted
!==================================================================================================

      use paramsmod
      use REACTMOD
      use filesmod
            
      character (LEN=LineL) :: Line
      character (LEN=24)  ::  formula, formulae(2)

      iscycle = 0
      if(Line(1:4) /= 'FORD' .and. Line(1:4) /= 'RORD') return

       if(numbereact(nreac) == 0) then
       write(file_err,*)' FORD or RORD data line before first reaction data '
       write(file_err,*) Line
       stop
       endif

! eg.   FORD /H 0.0/

     ! Read Species; is 2nd formula because FORD is read as the first
      call READMULTISPECIES(Line,Nspeciesread,formulae,2,24,lastchar)            
      formula = formulae(2) 
      nq = Indexspecies(formula)

      ! read Exponent
      islash = 0
      do i = 5,60
      if(Line(i:i) == '/') islash = 1 
      if(islash == 1 .and. Line(i:i) == ' ') exit 
      enddo
      read(Line(i:60),fmt=*,err=999,end=999) expn


      ! Term number
      iterm = 0
      do i = 1,maxtermsreaction
       if(nq == nqorderSp(i,nreac)) then
       iterm = i
       endif 
      enddo

      if(iterm == 0) then
      write(file_err,*)' FORD or RORD species is not in reaction No. ',numbereact(nreac)
      write(file_err,*)' or FORD, RORD only immediately after principal reaction (not for special 3rd body species)'
      write(file_err,*) Line
      write(file_err,*) ' nq species ',nq, formula 
      stop
      endif

      ! Store
      ! Increment index if no FORD already for this reaction
      if(iFORD(nreac) == 0) indxexpon = indxexpon + 1
      iFORD(nreac) = indxexpon         

       if(itermexpon(iterm,indxexpon) /= 0) then
       write(file_err,*)' FORD or RORD already specified for term, reaction. No. ',iterm,numbereact(nreac)
       write(file_err,*) Line
       stop
       endif

      itermexpon(iterm,indxexpon) = 1        ! indicates which term is involved    
      specexpon(iterm,indxexpon) = expn
      write(filereact,'(" Special exponent = ",f7.3," for ",a24)') expn,formula 

      iscycle = 1
      return

999   call errorline(fileSCHEME,' 77   ')
      return
      end

!==================================================================================================

      subroutine STOREREACTCOEFFS(n)

!===============================================================      
      use dimensmod
      use reactmod
      use speciesmod      

! Reactants/products nq for species in reaction (not 3rd body) not in order
! 1 - 3 reactants  (3 not used)
! 4 - 7 products

      ir = 1;  ip = 4
      do nq = nqsp1,nqLast       

       if(reactmol(nq)  <  0.) then
       nqspec(ir,n) = nq
       specnetmass(ir,n) = reactmol(nq) * spMW(nq) 
       ir = ir+1

       elseif(reactmol(nq)  >  0.) then
       nqspec(ip,n) = nq
       specnetmass(ip,n) = reactmol(nq) * spMW(nq) 
       ip = ip+1
       endif

      enddo

      return
      end

!===============================================================      

      subroutine REACTION_DATA_OUTPUT(ifirstsectionline)
      
!===============================================================                              
      use filesmod
      use reactmod
      use speciesmod      
      use Limitsmod
      use Logicmod
      use paramsmod
            
      real(8) :: specnetmols(maxtermsreaction)
      character (LEN=8)  :: spM
      character (LEN=5)  :: reverse
      character (LEN=4)  :: chrtag
      character (LEN=20)  :: term(7)

      write(filereact,'(/" NON-Sectional REACTIONS Stored Data.")')
      write(filereact,'(/10x," INPUT UNITS for E are  ",a10)') E_units(iunitE)
      write(filereact,'("                  Rate mol/(cm^3 - s) = A * T^n * exp(-E/RT)   ")')
      write(filereact,'(" Reactants: -ve mols;     Products: +ve mols.  ")')
      write(filereact,'(/ "  Tag    Line No.   Reaction No.      A          n         E kJ/kmol",T80,"    Reaction ")')

      do n = 1, nreactionlines
      if(n >= ifirstsectionline) exit  ! only do gas-phase
      if(mod(n,30)==0) write(filereact,'("  Tag    Line No.   Reaction No.      A          n         E kJ/kmol",T80,"    Reaction ")')      

       reverse = ' =>  '
       if(reversible(n)) reverse = ' <=> '

       itagn = itag(numbereact(n))
       chrtag = ' ' 
       if(itagn > 0) write(chrtag,'(i4)') itagn

       if(ireactype(n) == -1) then       ! 3rd body continuation
       write(filereact,'(a4,1x,i8," 3rd = ",i8,4x,1pd12.3,3x,a24)') chrtag,n,numbereact(n),(chemrate(1,n)+1.0),Spname(nq3rd(n))

       else
       do i = 1,7
       term(i) = ' '
       if(ordermols(i,n) /= 0.0) write(term(i),'(0pf5.1,1x,a14)') ordermols(i,n), Spname(nqorderSp(i,n))   
       enddo
       write(filereact,'(a4,1x,i8,7x,i8,4x,1p3d12.3,T80,3(a20),a5,4(a20) )') &
       chrtag,n,numbereact(n),(chemrate(i,n),i=1,3),(term(i),i=1,3), reverse, (term(i),i=4,7)
       endif

       iaux = iAuxIndex(numbereact(n))     ! Auxiliary data
       if(iaux > 0) then            
       if(numbereact(n) > numbereact(n-1)) then   ! to write once only
        if(iauxtype(iAux) < 0) then   
        write(filereact,'(15x,i8," LOW P ",4x,1p3d12.3)') numbereact(n),(auxrate(i,iaux),i=1,3)
        elseif(iauxtype(iAux) > 0) then   
        write(filereact,'(15x,i8," HI P  ",4x,1p3d12.3)') numbereact(n),(auxrate(i,iaux),i=1,3)
        endif

        if(abs(iauxtype(iAux)) == 2) then  
        write(filereact,'(15x,i8," TROE  ",4x,1p4d12.3)') numbereact(n),(auxrate(i,iaux),i=4,7)
        elseif(abs(iauxtype(iAux)) == 3) then
        write(filereact,'(15x,i8," SRI   ",4x,1p5d12.3)') numbereact(n),(auxrate(i,iaux),i=4,8)
        endif
       endif
       endif
              
      enddo

!----------------------------------------------------
      write(filereact, '(/1x,79(1h-))')

!      write(filereact,'(/20x," REACTIONS ")')

      write(filereact,'(/"   MAX. LIMIT of Chemical rate k (mol,cm,s) ",g12.3)') chemratemax
      write(filereact,'(/"   LIMIT of LN(K equilibrium)               ",f12.1)') rLnmax
      write(filereact,'(/"   Freeze Temperature (K) for reactions without T dependence  ",f7.0)') Tfreeze

      write(filereact, '(/1x,79(1h-)/)')

!----------------------------------------------------
      ! show net moles in reaction (not 3rd body)
      if(netmoleslist == 1) then 
      write(filereact,'(//19x,"     NET REACTING MOLES")')
      write(filereact,'(5x," Type 0 = Chemistry only; -1 = 3rd body continuation;    1=EBU only; 2=min(EBU,Chemistry)")')

      do  n = 1, nreactionlines
      if(n >= ifirstsectionline) exit  ! only do gas-phase
      spM = '        '
      if(mod(n+19,20) == 0) write(filereact,'(/16x,"React. No.     type  net moles 3rd body     Reactant",14x,"Reactant",10(15x,"Product"))')           
      
      specnetmols = 0.0
      do i = 1, maxtermsreaction      
      nq = nqspec(i,n)
      if(nq > 0) specnetmols(i) = specnetmass(i,n)/spMW(nq)
      enddo
      
      if(nq3rd(n) > 0) spM = spName(nq3rd(n))
      if(ireactype(n) == -1) then
      write(filereact,'(6x,2i14,9x,6x,a8,20(f7.3,1x,a14))')numbereact(n),ireactype(n),spM
      else
      write(filereact,'(6x,2i14,f9.1,6x,a8,20(f7.3,1x,a14))')numbereact(n),ireactype(n), &
      reactmolnet(n),spM,(specnetmols(i), spName(nqspec(i,n)), i=1,2),(specnetmols(i), spName(nqspec(i,n)), i=4,maxtermsreaction)
      endif               
      enddo

      write(filereact, '(/1x,79(1h-)/)')
      endif
!----------------------------------------------------

      write(filereact,'(/20x," SPECIES LACKING SOURCES OR SINKS ")')
      nosource = 0

      do i = 1,4
      
       if(i == 1) then
       write(filereact,'(/8x," Species with Sinks Only                      ")')
       elseif(i == 2) then
       write(filereact,'(/8x," Species with Sources Only                    ")')
       elseif(i == 3) then
       write(filereact,'(/8x," Species with NO Sources or Sinks             ")')
       else
       write(filereact,'(/8x," Species occurs in ONE-only Reversible Reaction. ")')
       endif


       do nq = nqsp1, nqLast
       isink   = 0
       isource = 0
       kount   = 0 
        do  ir = 1, nreactionlines
         do ii = 1,maxtermsreaction
         if(nqspec(ii,ir) /= nq) cycle 
          if(specnetmass(ii,ir) < 0.0) then
          isink = isink+1
          if(reversible(ir)) isource = isource+1
          endif

          if(specnetmass(ii,ir) > 0.0) then
          isource = isource+1
          if(reversible(ir)) isink = isink+1
          endif
         if(specnetmass(ii,ir) /= 0.0 .and. reversible(ir)) kount = kount + 1
        
         enddo
        enddo

       if(i == 1 .and. isink > 0 .and. isource == 0) write(filereact,'(25x,a24)') spName(nq)
       if(i == 2 .and. isink == 0 .and. isource > 0) write(filereact,'(25x,a24)') spName(nq)
       if(i == 4 .and. kount == 1 .and. isource == 1 .and. isink == 1 ) write(filereact,'(25x,a24)') spName(nq)
       if(i == 3 .and. isink == 0 .and. isource == 0) then
       write(filereact,'(25x,a24)') spName(nq)
       nosource(nq) = 1
       endif 
       enddo

      enddo

      write(filereact, '(/1x,79(1h-)/)')
!----------------------
      write(filereact, '(/1x,79(1h-)/)')

      return
      end

!===============================================================      

       subroutine ALLREACTSOURCE(isJacob,mode_rateout,Tlet,philet,dencell,ssu,ssp,vol)

! One location, sources (kg/m^3-s):  all species, all reactions  (including Polymer)
! with Jacobian
! ssu ssp dfdQ are zeroed here
!=========================================
      use paramsmod
      use dimensmod
      use filesmod
      use REACTMOD
      use JACMOD
      use logicmod 
      use pressmod
      use PHIMINMAXMOD
      use speciesmod
      use TAGMOD
      use MONITORPOSMOD
                
      integer :: nqs(maxtermsreaction), mMr(2)
      real    :: conc7(0:maxtermsreaction)           
      real    :: philet(nqsp1:nqlastm)
      real(8) :: rateout, frwdchrate, bkwdchrate
      real(8) :: ssu(nqSp1:nqlast), ssp(nqSp1:nqlast) 

!-----------------------------------------
! Non Reaction-Dependent Section

      dfdQ = 0.0
      ssu  = 0.0
      ssp  = 0.0  

! Equilibrium constant for all reversible reactions at Tlet
      call LnEquilibriumConst(Tlet)

! Conversion to concentration
      phitoConc(nqSp1:nqLast) =  dencell / spMW(nqSp1:nqLast) * 1.0e-3     ! g/cm^3 * mol /g -> mol/cm^3 
      phitoConc(nqBalance)    =  dencell / spMW(nqBalance)    * 1.0e-3        
      phitoConc(nqM)          =  pressure/(Rgas * Tlet)   * 1.0e-3     ! M 3rd body sum of all species

! Concentration  mol/cm^3
      conc(nqSp1:nqLast) =  philet(nqSp1:nqLast) * phitoConc(nqSp1:nqLast)    ! g/cm^3 * mol /g -> mol/cm^3 
      conc(nqBalance)    =  philet(nqBalance)    * phitoConc(nqBalance)        
      conc(nqM)          =  phitoConc(nqM)           
      conc(0) = 1.0

!-------------
! REACTION LOOP     
      nreac =  1 

      DO WHILE(nreac <= nreactionlines)
      nqs(1:maxtermsreaction) = nqspec (1:maxtermsreaction,nreac)      ! by species not term order; no repeated species  

      conc7 = 1.0                      ! conc7 is ordered with repeats
      mMr = 0                          ! the two term numbers which are 3rd body species else 0
       
       ! each Term.       
       do m = 1,maxtermsreaction              
       nq = nqorderSp(m,nreac)                            
       conc7(m) = conc(nq)             ! Reactant concentrations in term order for kinetic rate, not nqs order

       if(maxFORD > 0) call FORD(nreac,m,conc7(m))   ! Special exponents (if FORD specification)
       enddo

      ! 3rd body 
      conc3rdwtd =  conc(nq3rd(nreac))      ! if no 3rd body conc(0) = 1.0           

      ! Sum of 3rd body weighted concs. for all continuation lines with this reaction number .
      Ncontinue = 0                                         ! number of continuation lines to increment nreac    
      if(numbereact(nreac+1) == numbereact(nreac))  &       ! continuation if same reaction number following                  
      call WEIGHTED3rdBODYCONC(nreac,conc3rdwtd,Ncontinue)  ! nreac is incremented below by Ncontinue

      mMr(:) = iterm3rd(:,nreac)                ! term No.'s of 3rd body (or 0)
      conc7(mMr(1:2)) = conc3rdwtd              ! 3rd body conc to designated terms; to conc7(0) if not used


      !  For LOW P-dependent reaction, [M] is contained in k and is not in high P reaction.
      indexaux = iAuxIndex(numbereact(nreac))
      if(indexaux /= 0) then
      if(iauxtype(indexaux) < 0) conc7(mMr(1:2)) = 1.0
      endif
      call EACHELEMENTARY &
      (numbereact(nreac),nreac,reversible(nreac),nqs,nqorderSp(1,nreac),conc7,conc3rdwtd,Tlet,frwdchrate,bkwdchrate,ssu,ssp)

      ! Jacobian  dfdQ   kg/(m^3-s)
      if(isJacob==1) call makeJACOBIAN(nreac,nqs,bkwdchrate)

      ! Output for each reaction after final iteration
      if(after_isend) call RATES_FOR_OUTPUT
      
      ! increment nreac
      nreac = nreac + 1 + Ncontinue
      ENDDO  ! end reaction loop 

!-------------------------------------------------
! Avoid all-zero Jacobian for non-solved species
      if(isJacob==1) then
      do nq = nqsp1,nqlast
      if(.not.solvenq(nq)) dfdq(nq,nq) =  - 1.0
      enddo
      endif

      return

! + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +  
CONTAINS
      subroutine RATES_FOR_OUTPUT

      real(8) :: Creact,Hreact,Creactsect,Hreactsect,Cprodsect,Hprodsect, Xmols
            
      rateout = frwdchrate - bkwdchrate                                      !  mol/cm^3-s
      
      if(mode_rateout == 0) then    ! Reaction Analysis of reaction rates   
      reactratesum(nreac) = reactratesum(nreac) + rateout    * (vol*1.0e3)   !  mol/s  for Analysis
      fwdratesum(nreac)   = fwdratesum(nreac)   + frwdchrate * (vol*1.0e3)   !  mol/s  for Analysis 
       if(at_monitor) then  
       reactmon(nreac)     = rateout                 
       fwdreactmon(nreac)  = frwdchrate 
       endif
      endif
      
      ! Tagged Reactions summed for PLOT and CONTOUR. 
      irnumber = numbereact(nreac)
      if(includeReactSumNo(irnumber)) then
      if(mode_rateout > 0) call CHmols(nreac,Creact,Hreact,Creactsect,Hreactsect,Cprodsect,Hprodsect)

        if(mode_rateout == 0) then       ! rate
        Xmols = 1.0
        elseif(mode_rateout == 1) then   ! C rate
        Xmols = Creact
        elseif(mode_rateout == 2) then   ! H rate
        Xmols = Hreact
        elseif(mode_rateout == 3) then   ! C sectional net rate
        Xmols = Cprodsect - Creactsect
        elseif(mode_rateout == 4) then   ! H sectional net rate
        Xmols = Hprodsect - Hreactsect
        endif
      endif

       
      return
      end subroutine RATES_FOR_OUTPUT
! + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +  


      end subroutine ALLREACTSOURCE
!=====================================================================

      subroutine EACHELEMENTARY(number,nreac,isreversible,nqs,nqorder,conc7,conc3rdwtd,T,frwdchrate,bkwdchrate,tsu,tsp)

! One reaction, all species sources (kg/m^3-s)
!=========================================================================================
      use paramsmod
      use dimensmod
      use reactmod
      use filesmod
      use speciesmod
      use logicmod
      use Limitsmod
      use PHIMINMAXMOD
            
      real(8) :: Bconst, equilK, chemk, Fconst, chrate,  frwdchrate, bkwdchrate, DrlnK
      real(8) :: chemkaux, chemkinter
      real(8) :: tsu(nqSp1:nqlast), tsp(nqSp1:nqlast), tsx
      real(8) :: chemrateN(3),  TroeSRIdata(5)
      integer :: nqs(maxtermsreaction), nqorder(maxtermsreaction)
      real :: conc7(0:maxtermsreaction)
      logical :: isreversible
!=========================================================================================
      chrate     = 0.0
      frwdchrate = 0.0
      bkwdchrate = 0.0
      
! --------------------------------

! CHEMICAL

! Chem. moles   a.A  + b.B  + c.C  => d.D  + e.E  + f.F
! or            a.A  + b.B  + m.M  => d.D  + e.E  + m.M

! Units used are:
! concentration:      [A] = Ya/Wa * den  * 1.e-3           ->   mol/cm^3
! Rate constant        k  =  A  * T^n * exp (-E/RT)
! Molar reaction rate  r  =  k  *  [A]^a * [B]^b           ->   mol/cm^3-s
!    source term rate     =  - a * Wa * r * vol * 1.e3     ->   kg/s

! Data needs to satisfy:
!         r   (mol/cm^3-s)
!        [A]  (mol/cm^3)
!-------------------------

! Chemical rate Constant

      chemrateN = chemrate(1:3,nreac)
      call RATECONSTANT(chemk,T,Tfreeze,chemrateN)
      
      ! Coagulation and kinetic rate controlled for sectionals 
      if(phi_T(nreac) > 0.0) then                
      call COAGRATE(phi_T(nreac),chemk,T,dum)        

      else                                     ! if coagulation no need to check aux. data
      ! Auxiliary rate data
      indexaux = iAuxIndex(number)
       if(indexaux > 0 .and. T > Tfreeze) then

       chemrateN = auxrate(1:3,indexaux)
       call RATECONSTANT(chemkaux,T,Tfreeze,chemrateN)

        TroeSRIdata = auxrate(4:8,indexaux)
        if(iauxtype(indexaux) < 0) then       ! LOW  
        call RATEINTERPOLATE(iauxtype(indexaux), chemkinter, chemk, chemkaux, conc3rdwtd,T,TroeSRIData)
        else                                  ! HIGH
        call RATEINTERPOLATE(iauxtype(indexaux), chemkinter, chemkaux, chemk, conc3rdwtd,T,TroeSRIData)
        endif
       
       chemk = chemkinter
       endif
      endif

      Fconst = min(chemk,chemratemax)  ! Chem. k limit

!--------------------------------------------

     !  Elementary - must be uni-molecular expression

      chrate             = Fconst * conc7(1) * conc7(2) * conc7(3)   ! mol/cm^3-s    
      ratenq(nqorder(1)) = Fconst * conc7(2) * conc7(3)    ! Omit each reactant in turn for jacobian  
      ratenq(nqorder(2)) = Fconst * conc7(1) * conc7(3)    ! For repeat reactant only one will be omitted; 
      ratenq(nqorder(3)) = Fconst * conc7(1) * conc7(2)    !     is fine for linearised jacobian. 

!--------------------------------------

! Forward
      ! Chem rate: reactants have sp term; products have su term

      do m = 1,3
      nqrm = nqs(m)
      if(nqrm > 0) then !  g/mol * mol/cm^3-s *1e3= g/cm^3-s *1e3 =  kg/m^3-s 
      tsx = -specnetmass(m,nreac) * ratenq(nqrm) * phitoConc(nqrm) *1.e3    ! tsx is +ve
      tsp(nqrm) = tsp(nqrm) - tsx
      if(tsx < 0.0) call elementaryerror(1,2)
      endif
      enddo     

      do m = 4,maxtermsreaction
      nqpm = nqs(m)
      if(nqpm > 0) then    !     g/mol * mol/cm^3-s *1e3= g/cm^3-s *1e3 =  kg/m^3-s 
      tsx =  specnetmass(m,nreac) * chrate *1.e3
      tsu(nqpm) = tsu(nqpm) + tsx
      if(tsx < 0.0) call elementaryerror(1,1)
      endif
      enddo   

      frwdchrate = chrate      ! mol/cm^3-s      
      bkwdchrate = 0.0

      if(.not. isreversible) return

! Reversible

       rlnK = rLnKequil(nreac)

       if(rlnK > rLnmax) then         
       Bconst = 0.0

       else
       rlnK = max(rlnK, -rLnmax)         

       DrlnK  =  rlnK 
       equilK = Dexp (DrlnK) 

       Bconst  = chemk / equilK
       Bconst = min(Bconst,chemratemax)
       endif

! only 3 products for backward reaction allowed  
! Reaction rate     [ ] mol/cm^3 = Yi/Wi * den * 1e-3
      bkwdchrate = Bconst * conc7(4) * conc7(5) * conc7(6)          ! +ve    Uni-molecular expression needed.  
      ratenq(nqorder(4)) = Bconst * conc7(5) * conc7(6)     
      ratenq(nqorder(5)) = Bconst * conc7(4) * conc7(6) 
      ratenq(nqorder(6)) = Bconst * conc7(4) * conc7(5) 


      do m = 1,3                                   
      nqrm = nqs(m)
      nqpm = nqs(m+3)                        
      if(nqrm > 0) then ! -ve        +ve 
      tsx = -specnetmass(m,nreac) * bkwdchrate *1.e3          ! tsx is +ve    kg/m^3-s
      tsu(nqrm) = tsu(nqrm) + tsx                             ! kg/m^3-s
      if(tsx < 0.0) call elementaryerror(2,1)
      endif
                                            
      if(nqpm > 0) then ! +ve          +ve
      tsx = specnetmass(m+3,nreac) * ratenq(nqpm) * phitoConc(nqpm) *1.e3 ! kg/m^3-s
      tsp(nqpm) = tsp(nqpm) - tsx ! kg/m^3-s
      if(tsx < 0.0) call elementaryerror(2,2)
      endif
      enddo
     

      return

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
CONTAINS
      subroutine ELEMENTARYERROR(idirn,is_susp)

      if(idirn==1) write(file_err,*)' Error in Forward reaction: number nreac m ',number,nreac,m 
      if(idirn==2) write(file_err,*)' Error in Backwardreaction: number nreac m ',number,nreac,m 
      write(file_err,*)'All quantities should be +ve.'
      if(is_susp==1) write(file_err,*)' TSU ',tsx
      if(is_susp==2) write(file_err,*)' TSP ',tsx
      if(idirn==1 .and.  is_susp==1) write(file_err,*)'specnetmass(m,nreac),chrate ',specnetmass(m,nreac),chrate 
      if(idirn==1 .and.  is_susp==2) write(file_err,*)'nqrm,-specnetmass(m,nreac),ratenq(nqrm),phitoConc(nqrm) ', &
                                     nqrm,-specnetmass(m,nreac),ratenq(nqrm),phitoConc(nqrm)

      if(idirn==2 .and.  is_susp==1) write(file_err,*)'-specnetmass(m,nreac),bkwdchrate ',-specnetmass(m,nreac),bkwdchrate
      if(idirn==2 .and.  is_susp==2) write(file_err,*)'nqpm,specnetmass(m+3,nreac),ratenq(nqpm),phitoConc(nqpm) ', &
                                     nqpm,specnetmass(m+3,nreac),ratenq(nqpm),phitoConc(nqpm) 

      write(file_err,*)' conc3rdwtd,T ', conc3rdwtd,T 
      write(file_err,*)' conc7(0:maxtermsreaction) ',conc7
      write(file_err,*)' Fconst, Bconst',Fconst, Bconst
      write(file_err,*)' '
      write(file_err,*)'STOP for this reaction error.'
      stop
!      return 
      end subroutine ELEMENTARYERROR

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      end

!===================================================================

      subroutine RATECONSTANT(chemk,T,Tfreeze,chemrate)

! output chemk =  A* T^n * exp(-E/RT)
!=========================================================================
      use paramsmod
      use limitsmod
      use filesmod
      use dimensmod

      real(8) :: chemk, act, T_n, chemrate(3), EoRT

!-------------------------------------------------------

      
      if(T <= Tfreeze) then
      chemk = 0.0
      return
      endif

!  Expresssion
      if(chemrate(2)  /=  0.0) then
      T_n =  real(T,8) ** chemrate(2)
      else
      T_n = 1.0
      endif

      if(chemrate(3)  /=  0.) then
      EoRT =  chemrate(3) / real((Rgas * T),8) 
      EoRT =  min(EoRT, DLnmax)
      act  =  Dexp(-EoRT)

      else
      act = 1.0
      endif

      chemk = chemrate(1) * T_n * act

!---------------------

      return
      end
!=========================================================================

      subroutine makeJACOBIAN(nreac,nqs,bkwdchrate)

!============================================================      
      use reactmod
      use Jacmod
      
      real(8) :: bkwdchrate
      integer :: nqs(maxtermsreaction)
            
       do m = 1,maxtermsreaction            ! m is species not terms
       nqsm = nqs(m)
       if(nqsm == 0) cycle
       coefnqsm = specnetmass(m,nreac)

        do n = 1,3           ! n species on reaction side
        nqrn = nqs(n)
        if(nqrn > 0) dfdQ(nqsm,nqrn) = dfdQ(nqsm,nqrn) + coefnqsm * ratenq(nqrn) * phitoConc(nqrn) *1.e3 

        if(bkwdchrate > 0.0) then     ! only 3 products for backward reaction allowed  
        nqpn = nqs(n+3)      ! species on products side 
        if(nqpn > 0) dfdQ(nqsm,nqpn) = dfdQ(nqsm,nqpn) - coefnqsm * ratenq(nqpn) * phitoConc(nqpn) *1.e3 
        endif

        enddo                                     !     g/mol * mol/cm^3-s *1e3 = kg/m^3-s 
       enddo

      return
      end

!============================================================      

      subroutine WEIGHTED3rdBODYCONC(nreac,conc3rdwtd,Ncontinue)

! Sum of [3rd body concentration X rate (multiplier-1) factor]
! multiplier - 1  because conc3rwtd already contains sum of all species      
! increments Ncontinue; does all continuations for this NUMBEREACT
!===========================================================                       


      use dimensmod
      use REACTMOD

      if(nq3rd(nreac) == nqM) then          ! do if parent reaction 3rd body is M (sum of all species)     
      nextreac = nreac + 1
       do while (nextreac <= nreactionlines)                 
       if(numbereact(nextreac) /= numbereact(nreac)) exit    ! continue for same reaction number                  
       conc3rdwtd = conc3rdwtd + conc(nq3rd(nextreac)) * chemrate(1,nextreac)   ! (multiplier - 1) stored in chemrate(1
       Ncontinue = Ncontinue + 1
       nextreac = nextreac + 1
       enddo

      conc3rdwtd = max(conc3rdwtd,0.0)   ! prevent rounding error making conc3rdwtd -ve as conc3rdwtd can -> 0 
      endif

      return
      end
      
!===========================================================                       

      subroutine RATEINTERPOLATE(iauxtype, k, kinf, ko, concM, T, TroeSRIdata)

!===========================================================                    
! Lindemann form:  Lindemann, F.: Trans. Faraday Soc., 17 598 (1922).

! See p29  
!UC-405, SAND96-8216
!Printed May 1996
!CHEMKIN-III: A FORTRAN CHEMICAL KINETICS PACKAGE FOR THE ANALYSIS OF GASPHASE
!CHEMICAL AND PLASMA KINETICS
!Robert J. Kee, Fran M. Rupley, and Ellen Meeks, James A. Miller
!Combustion Chemistry Department, Sandia National Laboratories, Livermore, CA 94551-0969      

!   k     [ cm^3/mole-s ]
!  [M]    [ mole/cm^3 ]

! for LOW auxiliary data  (iauxtype < 0). Cannot have [M] in reaction. 
!   k = kinf { Pr/(1+Pr) }  kinf is P-independent high-P rate on reaction line      [cm^3/mole-s]   
!   Pr = ko[M] / kinf       ko is low-P rate on auxiliary line  [cm^6 / mole^2-s]

! if Pr >> 1 then k -> kinf. Since reaction is not fn of M at high P, M is not wanted in reaction.
! if Pr << 1 then k -> kinf * Pr = ko[M]. Therefore M dependence is introduced in reaction.


! for HIGH auxiliary data (iauxtype > 0) need M in reaction.  
!   k = ko { 1/(1+Pr) }     ko is P-independent low-P rate on reaction line        [cm^3/mole-s]  
!   Pr = ko[M] / kinf       kinf is high-P rate on auxiliary line                   [ 1/s ]

! if Pr >> 1 then k -> ko/Pr = kinf/[M] 
!      Since reaction is not fn of M at high P, M is needed in reaction eqn to cancel term.
! if Pr << 1 then k -> ko. Therefore M dependence in reaction is used.



! TROE for iauxtype = -2 or +2 

!-----------------------------------------------------

      real(8) :: k, kinf, ko, Pr, M, F, TroeSRIdata(5)
      
      M = concM 
      Pr = ko * M / kinf           


      if(iauxtype < 0) then                ! LOW
      k  = kinf * ( Pr / (1.0 + Pr) )  

      else                                 ! HIGH
      k  = ko * ( 1.0 / (1.0 + Pr) ) 
      endif

      if(abs(iauxtype) == 2) then          ! TROE
      call TROEfn(F,Pr,T,TroeSRIdata(1:4))
      k = k * F

      elseif(abs(iauxtype) == 3) then      ! SRI
      call SRIfn(F,Pr,T,TroeSRIdata(1:5))
      k = k * F
      endif
      
      return
      end
!===========================================================                    

      subroutine TROEfn(F,Pr,T,Troe)

!===========================================================                    

! from Chemkin, but noting errors: +0.75 not -0.75 and exp(-T not +T for all terms      

! Troedata read in order: alpha, T***,  T*  T**

! Always have 4 paramters but can make T** a large number to eliminate 4th parameter influence.
 

      use filesmod
      
      real(8) :: Pr, F, Troe(4), small=1.0d-60
      real(8) :: c, n, Fcent, TT, logFcent, logF, logPr, Prfn
      real(8) :: p2, p3, p4, dexpbound, xlim = 50.0
      
      TT = T
      p2 = -TT/(Troe(2) + 1.e-30)
      p3 = -TT/(Troe(3) + 1.e-30)
      p4 = -Troe(4)/TT
      Fcent = (1.0 - Troe(1)) * dexpbound(p2,xlim) + Troe(1) * dexpbound(p3,xlim) + dexpbound(p4,xlim)


       if(Fcent <= 0.0) then
       write(file_err,*)' Error in Troe data gives -ve Fcent. Data: ',Troe
       write(file_err,*)' T =  ', TT
       stop
       endif  
            
      logFcent = log10(max(Fcent,small)) 

      c = -0.4 - 0.67 * logFcent
      n = 0.75 - 1.27 * logFcent
      logPr = log10(max(Pr,small)) 

      Prfn =  (logPr + c) / ( n - 0.14*(logPr + c) )
      logF =  logFcent / (1.0 + Prfn**2)
      
      F = 10.0 ** logF

      return
      end
      
!===========================================================                    

      subroutine SRIfn(F,Pr,T,SRI)

!===========================================================                    

! SRI function

! See p29  
!UC-405, SAND96-8216
!Printed May 1996
!CHEMKIN-III: A FORTRAN CHEMICAL KINETICS PACKAGE FOR THE ANALYSIS OF GASPHASE

!      X = 1 / ( 1 + log**2(Pr) )
!      F = d * [a * exp(-b/T) + exp(-T/c)]**X  * T**e 

      use filesmod
      
      real(8) :: X, Pr, F, TT, SRI(5), PrLog, small=1.0d-60
      real(8) :: p2, p3, dexpbound, xlim = 50.0
      

      TT = T  
      PrLog = Log10(max(Pr,small))
      X = 1.0 / ( 1.0 + PrLog**2 )

      p2 = -SRI(2)/TT
      p3 = -TT/(SRI(3) + 1.e-30)

      F = SRI(4) * ( (SRI(1)*dexpbound(p2,xlim) + dexpbound(p3,xlim) )**X )  * TT**SRI(5)


      return
      end
      
!===========================================================                    

      function dexpbound(x,xlim)

      use paramsmod
      real(8) :: dexpbound, x, xlim 

      if(x < -xlim)    then
      dexpbound = 0.0

      elseif(x > xlim) then
      dexpbound = dexp(xlim)

      else
      dexpbound = dexp(x)
      endif             

      return
      end

!==============================================

      subroutine FORD(nreac,m,conc7m)

! exponent for species in term m if required
      use filesmod      
      use reactmod

       indx =  iFORD(nreac) 
        if(indx > 0) then                
         if(itermexpon(m,indx) == 1) then
         expn = specexpon(m,indx)
          if(expn == 0.0) then
          conc7m = 1.0
          else
          conc7m = conc7m ** expn              ! ignore Jacobian fix for this
          endif                     
         endif        
        endif        


      return
      end
!=================================================================            

       subroutine RESERVEDWORDS(Line,isreserved)

! test if reserved word is present in Line 
!========================================================

      use elementsmod

      character(Len=*) :: Line

           
      isreserved = 0       
      if(Line(1:4) == 'FORD') isreserved = 1
      if(Line(1:4) == 'RORD') isreserved = 1 
      if(Line(1:3) == 'LOW')  isreserved = 1
      if(Line(1:4) == 'HIGH') isreserved = 1 
      if(Line(1:3) == 'SRI')  isreserved = 1 
      if(Line(1:4) == 'TROE') isreserved = 1 

      return
      end

!========================================================

      subroutine REACTIONELEMENTRATE

! for output of rate of C consumed per reaction 
!==================================================================================================

      use reactmod
      use speciesmod
      use filesmod
            
      real(8) :: Rcarbon
      
      irprev = 0
      do n = 1,nreactionlines                              
      irnumber = numbereact(n)
      if(irnumber == irprev) cycle
      irprev = irnumber
            
      Rcarbon   = 0.0
      do i = 1,3
      nq = nqorderSp(i,n)
      if(nq==0 .or. nq == nq3rd(n)) cycle    ! 3rd body not included
      Rcarbon = Rcarbon + ordermols(i,n) * spCHONSI(1,nqorderSp(i,n))
      enddo

      reactantC(irnumber) = -Rcarbon    ! to make +ve

! write(file_err,*) n,irnumber,reactantC(irnumber)
      enddo
                        
      return
      end
!===============================================================                              

      subroutine TAGNUMBER(line,ireact,iitag)
      
      use paramsmod
      use reactmod
      use filesmod

      character (LEN=LineL) :: Line
      character (LEN=3) :: chr3

      iitag = 0
      do i = 1,125
      chr3 = '   '
       if(Line(i:i) == '#') then
       chr3 = Line(i+1:i+3)                        
       exit
       endif
      enddo

      if(chr3 /= '  ') then
      read(chr3,fmt='(i3)',err=999) iitag
      itag(ireact) = iitag
      endif

      return

999   write(file_err,'(" Skip bad # tag in following line at Reaction No. ",i6)')ireact
      write(file_err,*) line 
      return
      end
!===============================================================                              
      