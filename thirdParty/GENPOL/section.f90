
!========================================================

      subroutine MAKESECTIONRADICALS

!========================================================

      use dimensmod 
      use Logicmod
      use sectionalmod
      use speciesmod
      use numericmod
      use fuelmod
      use filesmod
      use sootdenmod
                                      
      real(8) :: Cn, Hn
      
      if(Sectionalmethod) then 
      nqSectionmoln = nqSectionmol1 + NsectionsAi -1
      nqSectionN = nqSectionmoln  

       if(isDoubleAi==2) then 
       nqSectionrad1 = nqSectionmoln + 1
       nqSectionradn = nqSectionrad1 + NsectionsAi -1
       nqSectionN    = nqSectionradn

       do n = 1,NsectionsAi
       nqmol = n + nqSectionmol1 - 1        
       nqrad = n + nqSectionrad1 - 1        

       Cn = spCHONSI(1,nqmol)
       Hn = spCHONSI(2,nqmol) - dHradical          

       spCHONSI(:,nqrad) = spCHONSI(:,nqmol)
       spCHONSI(2,nqrad) = Hn

       call SECTIONFORMULA(Cn,Hn,spFormula(nqrad),length,30) ! returns formula
       spName(nqrad)(3:) = spName(nqmol)(3:) 
       spName(nqrad)(1:2) = 'Aj'

       urfs(nqrad) = URFspeciesdflt

       spRelDiff(nqrad) = spRelDiff(nqmol) 
       printphis(nqrad) = printphis(nqmol)
       pltphis(nqrad)   = pltphis(nqmol)    

       if(solvenq(nqmol)) then
       solvenq(nqrad) = .true.
       else
       solvenq(nqrad) = .false.
       endif
       enddo             
      
       endif
      endif

      return
      end

!=========================================================================

      subroutine ExtraReactions(Linenumber,ireact,maxreactcols,nAiAj,Linesreqd,Line)

!=========================================================================
      use paramsmod
      use filesmod
      use sectionalmod
      
      character (LEN=LineL) ::  Line
      integer :: iterm(9,3)

      call Aisearch(Line,1,maxreactcols,nAiAj,iterm,ierrAisearch) 
       if(ierrAisearch /=0) then
       write(file_err,*)' Stop for Sectional Species specification error(s) in above list '
       stop
       endif
       
      rn = 1.0
      n = 1
            
      do i = 1,2 
      if(iterm(1,i) > 0) then    ! is Ai or Aj
      rn = rn *itermNsections(iterm(2,i),iterm(3,i),nCHSTsections(1)) !2 operator for C;  3 C section No.   
      rnH = itermNsections(iterm(4,i),iterm(5,i),nCHSTsections(2)) !4 operator for H;  5 H section No.
      rnH = min(rnH,HaverageN)
      rn = rn * rnH
      rnS = itermNsections(iterm(6,i),iterm(7,i),nCHSTsections(3)) !6 operator for S;  7 S section No.
      rnS = min(rnS,SaverageN)
      rn = rn * rnS
      rnT = itermNsections(iterm(8,i),iterm(9,i),nCHSTsections(4)) !8 operator for T;  9 T section No.
      rnT = min(rnT,TaverageN)
      rn = rn * rnT
      n = nint(rn*1.01 +1.0)  ! for rounding safety 
      endif
      if(i==1) n1 = n
      enddo

      isame = 1
      do i = 1,9
      if(iterm(i,1) /= iterm(i,2)) isame = 0
      enddo 
            
      
      if(isame == 1 .and. iterm(1,1) > 0) then
      n = (n+1)/2 + n1 ! half matrix only needed for Ai + Ai or Aj + Aj
      endif

      ireact = ireact + n
      Linesreqd = Linesreqd + n

      return
      end
!========================================================
      function itermNsections(iterm1,iterm2,NsectionsAi)

      use filesmod

!  Full specification operators: 
!  _ = only
! \ = except
! > = greater than   ( >00 means full range)
! < = less than  

! iterm(1    Ai=1  Aj=2
!       2    operators for C  0 = blank, 1= _  2= \  3= >   4= <
!       3    C section No. 
!       4    operators for H
!       5    H section No.
!       6    operator for S
!       7    S section No.
!       8    operator for T
!       9    T section No.    
!  eg.  if Ait2  then iterm =   1  3  0  3  0  3  0  1  2
!                               Ai                   _  T2

      
      if(iterm1==1) then           ! _
      itermNsections = 1  
      elseif(iterm1==2) then       ! \
      itermNsections = NsectionsAi - 1   
      elseif(iterm1==3) then       ! >
      itermNsections = NsectionsAi - iterm2  
      elseif(iterm1==4) then       ! <
      itermNsections = iterm2 - 1  
!       else
!      write(file_err,'(/" Error in ITERMNSECTIONS. iterm1 =  ",i5)') iterm1 
!      stop
      endif
        
      return
      end


!========================================================

      subroutine READGAMMA

!=================================================================            
      use Logicmod
      use filesmod
      use sectionalmod
            

! Hamaker constant for coagulation
      call findots(fileSCHEME) 

      read(fileSCHEME,*,err=999) Hamaker   
      call rLimitin(Hamaker,1.0e-25,1.0e-15," Hamaker constant ")
      read(fileSCHEME,*,err=999,end=999) Hamaker_jHfac(1:nHsections)  ! Hamaker(jH) / Hamaker constant
      read(fileSCHEME,*,err=999,end=999) Hamaker_iSfac(1:nSsections)  ! Hamaker(iS) / Hamaker constant
      read(fileSCHEME,*,err=999,end=999) Hamaker_iTfac(1:nTsections)  ! Hamaker(iS) / Hamaker constant
      read(fileSCHEME,*,err=999) HamakerPAH   
      read(fileSCHEME,*,err=999) nGammaTables

      nGammaTables = max( nGammaTables, 1)       
       if(nGammaTables > maxGammaTables) then  
       write(file_err,*)' No. of Gamma Tables is limited to ',maxGammaTables
       stop
       endif 

       write(filereact,'(/" No. of Gamma Tables is ",i4)') nGammaTables      


! Gamma 2-D Table data if present
! 1st column is C number to check correspondence with sectionals

      do iGammaTable = 1, nGammaTables

      call findots(fileSCHEME) 
      ii = 0 
      do i = 1, nCsections
      read(fileSCHEME,*,err=9,end=10) Cgammadata(i), Gamma2data(1:nCsections,i,iGammaTable)
      ii = ii + 1 
      enddo
      write(filereact,'(" Coagulation efficiency Gamma read from data Table. No. ",i4)') iGammaTable
      iGammadimension = 2
      
       do i = 1, nCsections
       do j = 1, nCsections
        if(Gamma2data(i,j,iGammaTable) <= 0.0 ) then     ! allow > 1
        write(file_err,*)' Error in GammaTable ',iGammaTable
        write(file_err,*)' Error at i,j, in gamma 2D data: ',i,j,Gamma2data(i,j,iGammaTable)
        stop
        endif 
       enddo
       enddo
      
      enddo

      return


10    if(ii == 0 .and. iGammaTable == 1) then      ! no data; Gamma = 1 (put data = 0)
      iGammadimension = 0
      Gamma2data = 0.0
      write(filereact,'(/" End of file at first line error. No data for Coagulation efficiency. Gamma = 1.0 ")') 

      else                      ! Data Error
      write(file_err,'(/" Error in reading Coagulation efficiency data. Stop.")') 
      write(file_err,'(" Searching for ",i6, " rows and columns ")') nCsections 
      write(file_err,'(" Error in row ",i6)') i
      write(file_err,'(" Gamma Table No.  ",i6)') iGammaTable
      stop
      endif
      return

9     write(file_err,'(/" Error in reading Coagulation efficiency data. Stop.")') 
      write(file_err,'(" Searching for ",i6, " rows and columns ")') nCsections 
      write(file_err,'(" Error in row ",i6)') i
      stop
!      return

999   call errorline(fileSCHEME,'gamma ')

      end
!=================================================================            

      subroutine GammaDataCheck

!=================================================================            
      use Logicmod
      use filesmod
      use sectionalmod
      use speciesmod
      integer(8) :: iC1, iC2       
! C number check correspondence with sectionals
      do i = 1, nCsections
       iC1 = sectionC(i)   + 0.5
       iC2 = Cgammadata(i) + 0.5
       if(iC1 /= iC2) then  
       write(file_err,*)' Gamma data C number does not correspond with sectional C number '
       write(file_err,*)' Section, Cdata, Csection ',i,Cgammadata(i),sectionC(i)
       stop
       endif 
      enddo

      return

      end
!=================================================================            

      subroutine SETCOAGFN(iabsA,nq1,nq2,Cn1,Cn2,iC1,iC2,chemrate,phi_T,Hamean) 

!=================================================================            
      use sectionalmod      
      use speciesmod
      use filesmod
      use paramsmod
      
      real(8) :: chemrate(3)

      if(nq1==0 .or. nq2==0) then
      write(file_err,*) ' Sectional: need nq1 & nq2 > 0 for SETCOAGFN ',nq1,nq2
      stop
      endif

       ! Betacoag - calculates min(free molecule, Stokes-Cunningham)
       call BETACOAGMIN(Cn1,Cn2,iC1,iC2,Acoag,Tncoag,Eactcoag)

       chemrate(1) = Acoag        
       chemrate(2) = Tncoag
       chemrate(3) = Eactcoag

! iC = 1 default for Gamma Table
! jH = 0 non-sectional: use HamakerPAH
       iC1 = 1;  jH1 = 0        
       iC2 = 1;  jH2 = 0        
        
       if(nq1 >= nqSectionmol1 .and. nq1 <= nqSectionradn) call SectionfromNQ(nq1,iC1,jH1,iS1,iT1,m1)
       if(nq2 >= nqSectionmol1 .and. nq2 <= nqSectionradn) call SectionfromNQ(nq2,iC2,jH2,iS2,iT2,m2)

! put phi_T(C1,C2,H1,H2) = phi_TonA/A * A = phi_T.  Later T is inserted at each location. 

! Hamaker interpolation across jH: cross-weighted by MW of reactants to emphasize section with low MW 
       if(jH1 == 0) then
       H1 = HamakerPAH * spMW(nq2)
       else
       H1 = Hamaker * Hamaker_jHfac(jH1) * Hamaker_iSfac(iS1) * Hamaker_iTfac(iT1) * spMW(nq2)
       endif
         
       if(jH2 == 0) then
       H2 = HamakerPAH * spMW(nq1)
       else
       H2 = Hamaker * Hamaker_jHfac(jH2) * Hamaker_iSfac(iS2) * Hamaker_iTfac(iT2) * spMW(nq1)
       endif
       
       Hamean =  (H1 + H2)/(spMW(nq1)+spMW(nq2))

!  gamma Table in terms of C's
       phi_T    = gamma2data(iC1,iC2,iAbsA) * Hamean     ! put into reaction-line array phi_T(n)

      return
      end
!=================================================================            

      subroutine SET_KINETIC_ARHENIUS_FN(Cn1,Cn2,iC1,iC2,AEin,chemrate,phi_T) 

!    kinetic rate * exp(-E/RT)
!=================================================================            
      use sectionalmod      
      use filesmod
      use paramsmod
      
      real(8) :: chemrate(3)
      
       ! Betacoag - calculates min(free molecule, Stokes-Cunningham)
       call BETACOAGMIN(Cn1,Cn2,iC1,iC2,Acoag,Tncoag,Eactcoag)

       chemrate(1)   = Acoag        
       chemrate(2)   = Tncoag
       chemrate(3) = Econversion(AEin,iunitE)      ! convert activation energy to kJ/kmol
       phi_T      = rhuge                  ! ensures that gamma = 1 in COAGRATE        

      return
      end
!=================================================================            

      subroutine COAGRATE(phi_T,chemk,T,gamma)        

! calculate gamma at T and multiply chemk by gamma
! output chemk

!=========================================================================
      use filesmod
      real(8) :: chemk

! Coagulation efficiency gamma = 1 - (1 + phi) * exp( -phi)
! where phi = phi0/A * A/T  where A is Hamaker const.
! phi0/A given by D'Anna spreadsheet (phi_TonA).

      
      gamma = 1.0
      phi = phi_T/T        
      if(phi < 20.0) gamma = 1.0 - (1.0 + phi) * exp(-phi)     ! phi = 20, gamma = 0.999999957
      chemk = chemk * real(gamma,8)
       
      return
      end
!=========================================================================

      subroutine MAKESECTIONNAME(nq,iT,iC,iH,iS,outname)

! makes names for section molecule
! radical name is copied from molecule name      
!==========================================================================      

      use logicmod      
      use sectionalmod      
      use filesmod
      
      character (LEN=24) :: outname
            
      outname = ' '

       if(iC > 99 .or. iH > 99 .or. iS > 99 .or. iT > 9) then
       write(file_err,*)' Stop: cannot have more than 99C, 99H, 99S, 9T Sections because of species names.',nq,iH,iS,iT
       stop
       endif

      write(outname(1:13),'("Ai_",i2.2,"_",i2.2,"_",i2.2,"_",i1)') iC,iH,iS,iT


      return
      end
!====================================================================== 

      function Cnumberfunction(iCtype,Cn1,Cn2,Cmin)

       if(iCtype ==1) then       ! mean
       Cnumberfunction = 0.5 * (Cn1 + Cn2)  
       Cnumberfunction = max(Cnumberfunction,Cmin)  
       elseif(iCtype ==2) then   ! max
       Cnumberfunction = max(Cn1,Cn2,Cmin)  
       elseif(iCtype ==3) then   ! sum
       Cnumberfunction = max( (Cn1 + Cn2), Cmin)  
       elseif(iCtype ==4) then   ! minimum
       Cnumberfunction = min(Cn1,Cn2)  
       Cnumberfunction = max(Cnumberfunction,Cmin)  
       else
       Cnumberfunction = 1
       endif

      return
      end
!===========================================================                       
      
      subroutine CsectionInterpolate(n,C8,frC,i)

! for C number returns Section number i below C and fraction frC between i & section i+1
! no tolerance error message if frc > 1 because reactions C8 > C(nCsections) should have been omitted

      use filesmod
      use Sectionalmod
      use reactmod
      
      real(8) :: C8
                                      
      if(C8 <= sectionC(1)) then
      frC =  C8/sectionC(1)
      i = 0
            
      else 
      call LOCATE8(sectionC, 1, nCsections, C8, i, 1, nCsections-1)
      frC =  (C8 - sectionC(i) ) / (sectionC(i+1) - sectionC(i))
      endif
      
        if(frC < 0.0 .or. frC > 1.0) then
        write(file_err,*)' Reaction ',numbereact(n),' Bad Sectional C number interpolation; frC = ', frC
        write(file_err,*)' nCsections = ',nCsections
        write(file_err,*)' Interpolating C = ',C8
        write(file_err,*)' Between i,sectionC(i) and i+1, sectionC(i+1) ',i,sectionC(i), i+1, sectionC(i+1)
        stop
        endif
      
      return
      end
!==========================================================            
      
      subroutine HsectionInterpolate(iC,frC,isAiorAj,particleH8,frH,j)

      use filesmod
      use Sectionalmod
      use reactmod
      
      real(8) :: particleH8, Hweighted(nHsections)

      if(nHsections == 1) then
      j = 1
      frh = 0.0
      return
      endif

! Use Hweighted interpolation at each H-section between C-sections 
      dH = dHradical*(isAiorAj - 1)   ! for Aj

      do j = 1, nHsections
       if(iC == 0) then

       Hweighted(j) = frC * (sectionH(iC+1,j) - dH)
       else
       Hweighted(j) = sectionH(iC,j) - dH + frC * (sectionH(iC+1,j) - sectionH(iC,j))    
       endif
      enddo

! j and frH for Hweighted to match particleH8
      call LOCATE8(Hweighted, 1, nHsections, particleH8, j, 1, nHsections-1)
      frH =  (particleH8 - Hweighted(j) ) / (Hweighted(j+1) - Hweighted(j))

! frH limited to 0 - 1      
      frH = min(frH,1.0)  ! LOCATE limits j to nHsections-1 and here limit frh to 1
      frH = max(frH,0.0)
      return
      end
!==========================================================            

      subroutine S_SECTIONINTERPOLATE(n,Line,iterm,S_X,iS1,iS2,frS1,frS2)

!  S is the specified sectional product species structure No. Products are distributed as S1, S2 etc. 
! If S = 0 then Sproduct = MW-weighted average of S of sectional reactants (non-sectional are ignored).  
! If S >= 1-99 (can be non-integer) then Sproducts = S
! If S >= 101 then assign to non-sectional reactant Sreactant = S - 100. Then same as S = 0 above. 
! If both reactants are non=sect then must specify S (>= 1) 
!=======================================================
      use paramsmod
      use dimensmod
      use sectionalmod
      use reactmod
      use speciesmod
      use filesmod
                   
      character (LEN=LineL) :: Line
      integer :: iterm(9,3)

     ! nSsections = 1 
     if(nSsections == 1) then
     S = 1.0
     iS1 = 1;    iS2 = 1
     frs1 = 1.0; frs2 = 0.0
     return
     endif

     ! never have 0 < S < 1   
     if(S_X > 0.0 .and. S_X < 1.0) then
     write(file_err,*) ' Error. Cannot have 0 < S < 1 '
     write(file_err,*) Line
     stop
     endif

     ! never have 100 < S < 101   
     if(S_X >= 100.0 .and. S_X < 101.0) then
     write(file_err,*) ' Error. Cannot have 100 < S < 101 '
     write(file_err,*) Line
     stop
     endif
                
     ! Both reactants non-sect need S >= 1 specified.
     if(iterm(1,1) == 0 .and. iterm(1,2) == 0 .and. S_X == 0.0) then
     write(file_err,*) ' Error. Need to specify S >= 1 for sectional product when both reactants are non-sectional. '
     write(file_err,*) Line
     stop
     endif

     IF(S_X > 0.0 .and. S_X < 100.0) THEN  ! Specified S = 1 - 99
     S = S_X 
     iS1 = int(S)
     iS2 = iS1 +1 

     ELSE            ! S = 0 or S >= 101 Interpolate. Must have at least one Sect if S<100.    
     nq1 = nqorderSp(1,n)
     nq2 = nqorderSp(2,n)
     S1 = 0.0;  S2 = 0.0
     iS1 = 1;   iS2 = 1
     spMW1 = 0.0; spMW2 = 0.0 
          
     if(iterm(1,1) > 0) then     ! is Ai or Aj with its S
     call SectionfromNQ(nq1,iC1,jH1,iS1,iT1,m1)
     S1 = iS1         
     spMW1 = spMW(nq1)  
     elseif(S_X > 100.0) then  
     S1 = S_X - 100.0         
     spMW1 = spMW(nq1)  
     endif

     if(iterm(1,2) > 0) then     ! is Ai or Aj   with its S
     call SectionfromNQ(nq2,iC2,jH2,iS2,iT2,m2)
     S2 = iS2
     spMW2 = spMW(nq2)  
     elseif(S_X > 100.0) then  
     S2 = S_X - 100.0         
     spMW2 = spMW(nq2)  
     endif
                 
     S = (S1*spMW1 + S2*spMW2)/(spMW1 + spMW2)

     if(spMW1 == 0.0) then 
     iS1 = iS2
     elseif(spMW2 == 0.0) then 
     iS2 = iS1
     else
     iS1 = S
     iS2 = iS1 + 1
     endif

     ENDIF
      
     frS2 = S - real(iS1)
     if(iS1 == iS2) frS2 = 0.0     ! make exact to reduce interpolations
     frS1 = 1.0 - frS2

      return
      end
!==============================================

      subroutine SECTIONREACTIONS(icall,ioutput,iitag,iGenericSectn,n,Line,isreversible,A,Tn,AEin,Cexp,Cfac,iCtype,S_X,ifragment, &
                 nHSTlimited,nlimitCexceeded,nlimitCunderflow,nlimitHexceeded,nsectreactaccepted,iterm)

! Sectional Method Reaction Generator 
!====================================================================      

! set 2 loops
! determine product
! update n
! check n dimension
! size dependent rate 

! Use Interpreter nq settings for Ai Aj to start sectional series.

! Use nqorderSP & Ordermol for operations (reactmol is 0 for repeated Ai or Aj)  

! iterm(1    Ai=1  Aj=2
!       2    operators for C  0 = blank, 1= _  2= \  3= >   4= <
!       3    C section No. 
!       4    operators for H
!       5    H section No.
!       6    operator for S
!       7    S section No.
!       8    operator for T
!       9    T section No.    
!  eg.  if Ait2  then iterm =   1  3  0  3  0  3  0  1  2
!                               Ai                   _  T2

!====================================================================      

      use paramsmod
      use dimensmod
      use speciesmod
      use reactmod
      use sectionalmod
      use filesmod
      use logicmod
            
      character (LEN=LineL) :: Line
      
      integer :: iterm(9,3), iCHST(4), nHSTlimited(3)
      real(8) :: ordermolkeep(maxtermsreaction)
      real(8) :: A
      real(8) :: particleC8, particleH8  
      integer :: nqorderSpkeep(maxtermsreaction)  

      logical :: isreversible
      
!      save :: Acoag,Tncoag,Eactcoag,Cexpcoag,Cfaccoag,iCtypecoag 
!------------------------------------------------------

      ierrmoles = 0 
      k1max = 1;  k2max = 1

      ! Output Generic Reaction Line
      if(ioutput==1) then
       GenericSectionimage(iGenericSectn) = Line   ! generic reaction image for output
      if(iGenericSectn == 1) write(filereact,'(//"SECTIONAL REACTIONS")')
      write(filereact,'(/"Tag   Section Generic No.",i8,9x,(a))')iGenericSectn,Line      
      if(ifragment>0)write(filereact,'(" Fragmentation ")')
      endif
      
! REACTANTS
      ! determine if Reactant 1 is Ai or Aj on nq number and set first nq
      if(iterm(1,1) > 0) k1max = NsectionsAi     ! term 1 is Ai or Aj
      nq01 = nqorderSp(1,n) - 1             ! sets start & whether molecule or radical


      ! determine if Reactant 2 is Ai or Aj on nq number and set first nq
      if(iterm(1,2) > 0) k2max = NsectionsAi    ! term 2 is Ai or Aj
      nq02 = nqorderSp(2,n) - 1            ! sets start & whether molecule or radical

! PRODUCTS
     ! determine which Product number is Ai or Aj and which number is first vacant.
     i_AiAj = 0; i_vacant = maxtermsreaction+1
     do i = 4, maxtermsreaction                          ! i = 4 + ...

      if(nqorderSp(i,n) /= 0 .and. i_vacant<maxtermsreaction) then ! terms must be vacant after first vacant term        
      write(file_err,*) ' Sectional: Vacant Product terms must be consecutive '
      write(file_err,*) Line
      stop
      endif
       
      if(nqorderSp(i,n) == nqSectionmol1) then       ! initially nq of Ai = nqSectionmol1    
      i_AiAj = i                                     ! the term No. which is Ai
      isAiorAj = 1
      elseif(nqorderSp(i,n) == nqSectionrad1) then   ! initially nq of Aj = nqSectionrad1
      i_AiAj = i                                     ! the term No. which is Aj
      isAiorAj = 2          
      elseif(nqorderSp(i,n) == 0 .and. i_vacant==maxtermsreaction+1) then         
      i_vacant = i                                   ! first term No. which is vacant
      endif
     enddo 

      ! Checks
      if(i_AiAj == 0) then
      write(file_err,*) ' Sectional: No Ai or Aj found in products '
      write(file_err,*) Line
      stop
      endif

      call CHECK_VACANT_TERMS(Line,n,i_vacant)

      call MOVE_TERM(n,i_AiAj,i_vacant)

!  assume uni-mole Ai Aj reactions
!  2 reactants only - no 3rd body

! preserve record of all species from Interpreter in nqorderSpkeep ordermolkeep except Ai Aj 
! Reacting moles in Ordermol does not change on reactant side, does change on product side

      nqorderSpkeep(:) = nqorderSp(:,n)     
      ordermolkeep = ordermols(:,n)      ! mols of all reactants and of non- Ai Aj products 

      do i = 1,maxtermsreaction
      if(nqorderSpkeep(i) ==  nqSectionmol1 .or. nqorderSpkeep(i) ==  nqSectionrad1) then
      nqorderSpkeep(i) = 0                 ! nq = 0 for all Ai, Aj reactants & products
      if(i > 3) ordermolkeep(i) = 0.0      ! product mols Ai or Aj  = 0
      endif
      enddo
            
!-------------------------------------            
! Reaction starting number (to reverse above)
      n = n - 1                          
        
      iout = 0            

! SECTIONAL LOOP
LP1:  DO k1= 1,k1max 
      nq1 = nq01 + k1 
      iC1 = 0 
       if(k1max > 1) then  ! means it is Ai or Aj
       call SectionfromNQ(nq1,iCHST(1),iCHST(2),iCHST(3),iCHST(4),m)
       call SELECTREACTION(1,iskip)
       if(iskip==1) cycle LP1
       iC1 = iCHST(1)
       endif
       
      k2min = 1
      if(nq01 == nq02) k2min = k1  ! both are Ai or both are Aj;  avoid repeating reactions with same species      

LP2:   DO k2 = k2min,k2max 
       nq2 = nq02 + k2 
       iC2 = 0
       if(k2max > 1) then  ! means it is Ai or Aj
       call SectionfromNQ(nq2,iCHST(1),iCHST(2),iCHST(3),iCHST(4),m)
       call SELECTREACTION(2,iskip)
       if(iskip==1) cycle LP2
       iC2 = iCHST(1)
       endif

       call PARTICLE_CHMOLES(nq1,nq2,ordermolkeep,nqorderSpkeep,particleC8,particleH8)

       if(ifragment > 0) then   ! Fragmentation
       particleC8 = 0.5 * particleC8 
       particleH8 = 0.5 * particleH8    ! for molecule (not radical) 
       call CHECK_C_LOWERLIMIT(icycle)
       if(icycle==1) cycle
       else 
       call CHECK_C_UPPERLIMIT(icycle)
       if(icycle==1) cycle
       endif 

       ! New reaction line number
       n = n + 1
        if(n > maxReactLines) then          ! check
        write(file_err,*)' Sectional reaction No. > maxReactLines dimension ', n,maxReactLines
        stop
        endif


       ! update for new n    
       numbereact(n) = numbereact(n-1) + 1   
       ireactype(n)  = 0       ! default for elementary
       reversible(n) = isreversible
       nq3rd(n) = 0
       if(icall==2) itag(numbereact(n)) = iitag    ! tag reaction
       
       ! Default reactants & products which are not Ai Aj in new nqorderSp, ordermol which are not Ai Aj
       nqorderSp(:,n) = nqorderSpkeep(:)
       ordermols(:,n) = ordermolkeep           ! all reactants preserved
       nqorderSp(1,n) = nq1 
       nqorderSp(2,n) = nq2 

      ! Product moles of Ai or Aj sections by interpolation for particleC, particleH 
      call SECTIONPRODUCTS(ioutput,Line,n,nq1,nq2,iterm,S_X,i_vacant,particleC8, particleH8,isAiorAj,ifragment,iomit)

       if(iomit > 0) then
       nHSTlimited(iomit) = nHSTlimited(iomit) + 1
       n = n-1
       cycle
       endif

       nsectreactaccepted = nsectreactaccepted + 1 

      !   Balance excess H
      call SectionHbalance(n)

       ! New Reactmol
       Reactmol(:)  = 0.0
        do i = 1,maxtermsreaction              
        nq = nqOrderSp(i,n) 
        Reactmol(nq) = Reactmol(nq) + ordermols(i,n)
        enddo               

       call STOREREACTCOEFFS(n)
       call REACTIONELEMENTCHECK(n,ierrmoles,xsC,xsH)
       if(ierrMoles /= 0) then     
       write(file_err,*)' Stop for C Moles error in sectional reaction '
       write(file_err,*) Line
       stop
       endif

       ! Output No. of reactions with xsH /= 0
       if(abs(xsH) > xsHcriterion) then
        nlimitHexceeded = nlimitHexceeded+1  ! multiple series reaction loses last reaction: is normal,no Warning
        if(ioutput==1) write(filereact,'(" Omitted reaction ",2a15," abs(xsH) > criterion ",1pg12.3)')spName(nq1),spName(nq2),xsH
        n = n - 1   ! backtrack
        nsectreactaccepted = nsectreactaccepted - 1 
        cycle
        endif


       ! RATE: Chemical or Coagulation Rate 
       call SECTIONRATE

       ! Sectional Output
       if(ioutput==1) call SECTIONOUTPUT(n,iout,xsH,Hamean)


       ENDDO LP2
      ENDDO LP1
!-------------------------------------            

      return

!# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
CONTAINS

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine SELECTREACTION(it,iskip)

!  Selection by Operator  
!  _  include single section
!    \ exclude single section
!  > include sections greater than  
!  < include sections less than

       iskip = 0      
       do i = 2,8,2 
       if(    iterm(i,it)==1 .and. iCHST(i/2)/=iterm(i+1,it)) then
       iskip = 1   ! single include
       elseif(iterm(i,it)==2 .and. iCHST(i/2)==iterm(i+1,it)) then
        iskip = 1   ! single exclude
       elseif(iterm(i,it)==3 .and. iCHST(i/2)<=iterm(i+1,it)) then
        iskip = 1   ! >
       elseif(iterm(i,it)==4 .and. iCHST(i/2)>=iterm(i+1,it)) then
        iskip = 1   ! <
       endif
       enddo

      return
      end subroutine SELECTREACTION
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine MOVE_TERM(n,i_AiAj,i_vacant)
      
!  Make Ai,Aj the last product so that vacant terms are contiguous
!     iterm              4      5     6     7       8
!     before             X      Ai    X    vacant  vacant  
!     after              X      X     Ai   vacant  vacant
      if(i_AiAj < i_vacant-1) then
      nqorderSp(i_vacant,n)   =  nqorderSp(i_AiAj,n)   ! AiAj to vacant 
      ordermols(i_vacant,n)   =  ordermols(i_AiAj,n)
      nqorderSp(i_AiAj,n)     = nqorderSp(i_vacant-1,n)  ! last term to previous AiAj
      ordermols(i_AiAj,n)     = ordermols(i_vacant-1,n) 
      nqorderSp(i_vacant-1,n) = nqorderSp(i_vacant,n)  ! Ai (in ivacant) to ivacant-1
      ordermols(i_vacant-1,n) = ordermols(i_vacant,n)
      nqorderSp(i_vacant,n) = 0
      ordermols(i_vacant,n) = 0
      i_AiAj = i_vacant-1
      endif
      i_vacant = i_vacant-1   ! 1st product term writes over generic Ai or Aj

      return
      end SUBROUTINE MOVE_TERM
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine CHECK_VACANT_TERMS(Line,n,i_vacant)

      use paramsmod
      character (LEN=LineL) ::  Line
      
      ! Need 1 vacant term on product side to accept 2nd Ai or Aj for C interpolation
       ineed = 1 
       if((maxtermsreaction - i_vacant+1) < ineed) then
       write(file_err,*) ' Sectional: No vacant space for C interpolation; max terms = ',maxtermsreaction
       write(file_err,'(i8,5x,(a))') numbereact(n), Line
       stop 
       endif

      ! Need 3 vacant terms on product side to accept 4th Ai or Aj  for H interpolation
       if(nHsections > 1) ineed = ineed + 2
       if(nHsections > 1 .and. (maxtermsreaction-i_vacant+1) < ineed) then
       write(file_err,*) ' Sectional: No vacant space H interpolation; max terms = ',maxtermsreaction
       write(file_err,'(i8,5x,(a))') numbereact(n), Line
       stop 
       endif

      ! Need 7 vacant terms on product side to accept 4th Ai or Aj  for H interpolation
       if(nSsections > 1) ineed = ineed + 4
       if(nSsections > 1 .and. (maxtermsreaction-i_vacant+1) < ineed) then
       write(file_err,*) ' Sectional: No vacant space for Structure interpolations; max terms = ',maxtermsreaction
       write(file_err,'(i8,5x,(a))') numbereact(n), Line
       stop 
       endif

      return
      end subroutine CHECK_VACANT_TERMS
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

       subroutine CHECK_C_UPPERLIMIT(icycle)

       ! Upper C limit check 
       icycle = 0       
       if(particleC8 > SpCHONSI(1,nqSectionmoln)) then
        if(k1max==1 .and. k2max==1) then
        ! single reaction exceeds C limit: Error
        if(ioutput==1) write(file_err,*)' Stop for following Section SINGLE reaction product C > Last Section C; See REACTION.txt '
        if(ioutput==1) write(filereact,*)' ERROR: SINGLE reaction omitted; Product C > last section C',particleC8,SpCHONSI(1,nqSectionmoln)
        if(ioutput==1) write(filereact,*) Line
        stop
        else
        nlimitCexceeded = nlimitCexceeded+1  ! multiple series reaction loses last reaction: is normal,no Warning
        sumnC = SpCHONSI(1,nq1) + SpCHONSI(1,nq2)
        if(ioutput==1) write(filereact,'(" Omitted reaction ",2a15," C ",1pg12.3," exceeds max.")') spName(nq1),spName(nq2),sumnC
        icycle = 1
        endif
       endif

       return
       end subroutine CHECK_C_UPPERLIMIT
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       subroutine CHECK_C_LOWERLIMIT(icycle)

       ! Upper C limit check 
       icycle = 0       
       if(particleC8 < SpCHONSI(1,nqSectionmol1)) then
        if(k1max==1 .and. k2max==1) then
        ! single reaction under Cmin: Error
        if(ioutput==1) write(file_err,*)' Stop for following Section SINGLE reaction product C < First Section C; See REACTION.txt '
        if(ioutput==1) write(filereact,*)' ERROR: SINGLE reaction omitted; Product C < First section C',particleC8,SpCHONSI(1,nqSectionmol1)
        if(ioutput==1) write(filereact,*) Line
        stop
        else
        nlimitCunderflow = nlimitCunderflow+1  ! multiple series fragmentation loses reaction: is normal,no Warning
        if(ioutput==1) write(filereact,'(" Omitted reaction ",2a15," Frag. C ",1pg12.3," below minimum.")')spName(nq1),spName(nq2),SNGL(particleC8)
        icycle = 1
        endif
       endif

       return
       end subroutine CHECK_C_LOWERLIMIT
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

       subroutine SECTIONRATE
       use paramsmod
       ! RATE: Chemical or Coagulation Rate 

       ! C number as function of both species
       Cn1 = spCHONSI(1,nq1)
       Cn2 = 0.0
       if(nq2 > 0) Cn2 = spCHONSI(1,nq2)     ! if there is a 2nd term

       ! Coagulation or Kinetic rate * Arhenius for A < 0 
       Hamean = 0.0 
       if(A < 0.0) then                 !  k1, k2 = 1 for non-sectional species in sectional equation
       iabsA = abs(nint(A))

        if(iabsA > nGammaTables) then
        write(file_err,*) ' Sectional: Requested Gamma Table No. ',iabsA,' exceeds actual ',nGammaTables
        write(file_err,'(a)') Line
        stop
        endif

        if(iabsA == 100) then   ! A = -100:   Use kinetic collision rate * exp(-E/RT) 
        call SET_KINETIC_ARHENIUS_FN(Cn1,Cn2,iC1,iC2,AEin,chemrate(1,n),phi_T(n)) 

        else       ! A = -1 -> -5: Use Gamma Tables 1-5.  k1 or k2 = 1 for gas-phase species in sectional equation
        call SETCOAGFN(iabsA,nq1,nq2,Cn1,Cn2,iC1,iC2,chemrate(1,n),phi_T(n),Hamean) 
        endif
       
       
       ! Sectional rate expression from Data 
       else
       Cfn = Cnumberfunction(iCtype,Cn1,Cn2,0.0)
       chemrate(1,n) = A * real((Cfac * Cfn),8) **Cexp         ! A * (m*C)^k      
       chemrate(2,n) = Tn
       chemrate(3,n) = Econversion(AEin,iunitE)                       ! convert activation energy to kJ/kmol
       endif

       ! Check
       if(chemrate(1,n) <= 0.0) then
       write(file_err,*) ' Sectional Reaction rate constant <= 0, reaction ',numbereact(n)
       write(file_err,*) Line
       write(file_err,*) ' A = ',chemrate(1,n)
       stop
       endif         

       return
       end subroutine SECTIONRATE

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

      end

!============================================================      

      subroutine Aisearch(Line,i1,i2,nAiAj,iterm,ier) 

!============================================================      
! searches for number of occurrences of string Ai optionally excluding Ai_ and Aj_  
! nAiAj is No. of Ai or Aj  LHS + RHS

! Section Terminology
! Reaction side
!  Ai  Aj      All  loops  summing T1 and T2   
!   Ait1    Ait2   Ajt1   Ajt2    All loops for types T1 or T2   

!  Full specification operators: 
!  _ = only
! \ = except
! > = greater than   ( >00 means full range)
! < = less than  
! TERM(9,3)
! (:,1) refers to 1st Ai or Aj term
! (:,2) refers to 2nd Ai or Aj term
! (:,3) refers to product Ai or Aj term
 
          
! iterm(1    Ai=1  Aj=2
!       2    operators for C  0 = blank, 1= _  2= \  3= >   4= <
!       3    C section No. 
!       4    operators for H
!       5    H section No.
!       6    operator for S
!       7    S section No.
!       8    operator for T
!       9    T section No.    

!  eg.  if Ait2  then iterm =   1  3  0  3  0  3  0  1  2
!                               Ai >  0  >  0  >  0  _  T2


!  eg. Full format is    Aj>15>00\04_2   means all loops for C>15, all H, exclude S=4, only for T=2 


! Product side:   
!    Ai  Aj    not permitted because   T1 or T2 not specified. OK if there is only T1
!    Ait1    Ait2   Ajt1   Ajt2   specify T1 or T2. C, H, S follow interpolation rules.

! eg.      Ait2   + Aj>15>00\04_2   =>     Ait2
!--------------------------------------------------------------------
      use filesmod
      use sectionalmod
            
      character (LEN=*) ::  Line
      character (LEN=2) :: c2
      character (LEN=1) :: c1
      character (LEN=11) :: c11

      integer :: iterm(9,3)

       iterm = 0
       ier  = 0
       ierR = 0
       ierP = 0


       nAiAj = 0       ! total No. of occurrences of Ai and Aj   
       LHS = 1          ! LHS of reaction
       numberterm = 1
       kountpluses = 0   ! to check that there are only 2 reactants for Sectional reaction
       nAiAjprod = 0
       
       do i = i1,i2                
        if(Line(i:i) == '=') then
        LHS = 0             ! is product side   
        elseif(Line(i:i) == '+') then
        if(LHS==1) kountpluses = kountpluses + 1
        numberterm = 2      ! is second term
        endif
       c2  = Line(i:i+1)
       c1  = Line(i+2:i+2)

       isAiAj = 0
       if(c2 == 'Ai') isAiAj = 1  !   found Ai  
       if(c2 == 'Aj') isAiAj = 2  !   found Aj

       if(isAiAj > 0) nAiAj = nAiAj + 1  ! No. of occurences of Ai, Aj in LHS + RHS              
       if(LHS == 0 .and. isAiAj > 0) nAiAjprod = nAiAjprod + 1  ! No. of occurences of Ai, Aj  RHS              
       
       ! Product side check if Ait1 Ajt1 etc instead of Ai Aj. Record Ai or Aj and T section No.
       if(LHS==0)  then            
        if(isAiAj > 0) then

         if(c1 /= 't' .and. c1 /= 'T') then
          if(nCHSTsections(4) > 1) then
          ierP = 1
          else
          it = 1   ! if Ai, Aj and nTsections = 1 then take as AiT1 AjT1
          endif 
         else
         read(Line(i+3:i+3),*,iostat=ios) iT 
         if(ios/=0) ierP = 1
         endif        
        
         if(ierP/=0) then
         write(file_err,*) ' Sectional product must be of form Ait1 or Ajt1 or Ajt2 etc.' !,linenumber
         write(file_err,*) line
         endif

         if(it > nCHSTsections(4)) then
         write(file_err,*) ' Sectional product Type No > No. T sections; is ',iT
         write(file_err,*) line
         ierP = 1
         endif

        iterm(1,3) = isAiAj  ! product Ai or Aj
        iterm(9,3) = it      ! product T section No.
        endif      

       cycle                      ! cycle if Product side
       endif     



! LHS only
       if(isAiAj == 0) cycle    
       
       iterm(1,numberterm) = isAiAj    ! is Ai or Aj
       if(c1 == ' ' .or. c1 == char(9) .or. c1 == '+') then
       iterm(2,numberterm) = 3   ! > 0
       iterm(4,numberterm) = 3
       iterm(6,numberterm) = 3
       iterm(8,numberterm) = 3
       cycle  
       endif
       
       if(c1 == 't' .or. c1 == 'T') then         ! Ait1 Ait2 Ajt1 Ajt2
       iterm(2,numberterm) = 3   ! > 0
       iterm(4,numberterm) = 3
       iterm(6,numberterm) = 3
       iterm(8,numberterm) = 1       ! means _
       read(Line(i+3:i+3),*,iostat=ios) iterm(9,numberterm)
       cycle
       endif

       c11  = Line(i+2:i+12)   !  Aj>15>00\04_2
                      ! c11         123456789 11
                      !                      10 12

      iterm(2,numberterm) = ioperator(C11(1:1))         ! C operator
      read(c11(2:3),*,iostat=ios) iterm(3,numberterm)   ! iC 
      iterm(4,numberterm) = ioperator(C11(4:4))         ! H operator 
      read(c11(5:6),*,iostat=ios) iterm(5,numberterm)   ! iH
      iterm(6,numberterm) = ioperator(C11(7:7))         ! S operator
      read(c11(8:9),*,iostat=ios) iterm(7,numberterm)   ! iS
      iterm(8,numberterm) = ioperator(C11(10:10))       ! T operator
      read(c11(11:11),*,iostat=ios) iterm(9,numberterm) ! iT

      ! error in reading section numbers or operator error
      if(ios/=0 .or. &
      iterm(2,numberterm)==0 .or. iterm(4,numberterm)==0 .or. iterm(6,numberterm)==0 .or. iterm(8,numberterm)==0) then
      write(file_err,*) ' Error in folowing sectional species specification for term ',numberterm
      write(file_err,*) line
      ierR = 1
      endif
 
      ! Bounds check
      if(ierR == 0) then
       do j = 2,8,2
       if(iterm(j,numberterm) == 1 .and. (iterm(j+1,numberterm) < 1 .or. iterm(j+1,numberterm) > nCHSTsections(j/2))) ierr = 1
       if(iterm(j,numberterm) == 2 .and. (iterm(j+1,numberterm) < 1 .or. iterm(j+1,numberterm) > nCHSTsections(j/2))) ierr = 1
       if(iterm(j,numberterm) == 3 .and. (iterm(j+1,numberterm) < 0 .or. iterm(j+1,numberterm) > nCHSTsections(j/2) -1)) ierr = 1
       if(iterm(j,numberterm) == 4 .and. (iterm(j+1,numberterm) < 2 .or. iterm(j+1,numberterm) > nCHSTsections(j/2) +1)) ierr = 1
       enddo
      
       if(ierR > 0) then
       write(file_err,*) ' Error in folowing sectional species specification for term No. ',numberterm
       write(file_err,*) ' Section number(s) used in species exceeds section range. '
       write(file_err,*) line
       endif 
      endif
      
      enddo

      ! Check only 2 reactants; no 3rd body
      if(kountpluses > 1 .and. nAiAj>0) then
      write(file_err,*) ' Sectional reaction can only have 2 reactants ' !,linenumber
      write(file_err,*) line
      stop
      endif 

      ! Check only 1 Ai or Aij product
      if(nAiAjprod > 1) then
      write(file_err,*) ' Sectional reaction can only have 1 Ai or Aj product ' !,linenumber
      write(file_err,*) line
      stop
      endif 

      if(ierR /= 0 .or. ierP /= 0) ier = 1
      return
      end
!============================================================      

      function ioperator(c1)
      
      character (LEN=1) :: c1

! operators  1= _  2= \  3= >   4= <

      if(c1 == '_') then
      ioperator = 1
      elseif(c1 == '\') then
      ioperator = 2
      elseif(c1 == '>') then
      ioperator = 3
      elseif(c1 == '<') then
      ioperator = 4
      else      
      ioperator = 0
      endif 

      return
      end 

!============================================================      

      subroutine SECTIONFORMULA(Cn,Hn,formula,Length,maxLength)
      
!=====================================================================
! read one species formatted as:  @  xxxx  yyyyy (real xxxx, yyyy)  to represent CxxxxHyyyy   
! returns formula CxxxxHyyyy  (integer(8) xxxx yyyyy), spCHONSI, length of formula 

      use paramsmod
      use filesmod
      use dimensmod
      
      character (LEN=LineL) ::  Linetemp !,Line
      character (LEN=maxLength)  ::  formula
      integer(8) :: iCn, iHn      
      real(8) :: Cn, Hn

      iCn = Cn + 0.5
      iHn = Hn + 0.5   
      
      write(Linetemp,*)'C',iCn,'H',iHn     ! write as C xxx H yyyy in Linetemp

      formula = ' '
      j = 0                               ! copy to formula without spaces
      do i = 1,maxLength
       do while(j < LineL) 
       j = j+1
       if(Linetemp(j:j) /= ' ') exit
       enddo 
      formula(i:i) = Linetemp(j:j)
      ii = i
      if(j == LineL) exit 
      enddo 

      if(ii == maxLength) then
      write(file_err,*)' Need longer formula format to fit following species: '      
      write(file_err,*) Linetemp
      stop
      endif

      Length = ii - 1  ! since last ii is a blank
            
      return
      end
      
!====================================================================



      subroutine SectionReactCheck(nq3rd,iterm3rd,reversible,Line)
      
      use paramsmod
      use filesmod
      character (LEN=LineL) :: Line
      logical :: reversible
      integer :: iterm3rd(2)
      
       nq3rd = 0   ! Ai on both sides is mistaken for 3rd body; assume no 3rd body reactions for Sectional
       iterm3rd(:) = 0
       
       if(reversible) then   
       write(file_err,*) ' No reversible reactions allowed for Ai Aj; not unimolecular and no real thermo properties '
       write(file_err,*) Line
       stop
       endif       

      return
      end
!=========================================================================


!====================================================================== 
! NOT USED
      subroutine XYINTERP(xy,xin,yout,ierr)

!   X,Y Table interpolation
!   ierr = -1 if X underflow   yout = y(1)
!   ierr =  1 if X overflow    yout = y(n)
!   ierr =  0 if in range
!   ierr =  2 if x values equal or out of sequence   
!   No. of entries are in xy(1,0)
!   searches from previous access ii
!===============================================================
      use dimensmod
      data ii / 1 /

dimension xy(2,0:2)  ! dummy

      n = nint(xy(1,0))
      ierr = 0


!---------------------------
      if(xin < xy(1,1) ) then
      yout = xy(2,1)
      ierr = -1
      return
      endif

      if(xin > xy(1,n) ) then
      yout = xy(2,n)
      ierr = 1
      return
      endif

      if(xy(1,ii) < xin) then

       do i = ii+1, n
       if(xy(1,i-1) < Xin .and.  xy(1,i) >= Xin) exit
       enddo

      else

       do i = ii, 2, -1
       if(xy(1,i-1) < Xin .and.  xy(1,i) >= Xin) exit
       enddo

      endif

      ii = i
      frac = (xin - xy(1,i-1)) / (xy(1,i) - xy(1,i-1))
      yout = xy(2,i-1) + frac * (xy(2,i) - xy(2,i-1))
       

      return
      end
!=========================================================================

      subroutine SECTIONOUTPUT(n,iout,xsH,Hamean)

!=========================================================================

      use filesmod
      use speciesmod
      use dimensmod
      use reactmod
      use paramsmod
      use logicmod
                                 
      real(8) :: rate_gamma, rate_nogamma
      real :: gamma_at_T
      character(LEN=16)  :: Spnameout(0:maxtermsreaction)
      character(LEN=4)   :: arrow
      character(LEN=LineL) :: Line1
      character(LEN=20) :: chr20
      character(LEN=24) :: chr24
      character(LEN=32) :: chr32
      character(LEN=32) :: chr4='    '
      
      iout = iout + 1

       Spnameout = ' '
       maxterm = 0   ! No. of terms in equation
        do i = 1,maxtermsreaction
        Spnameout(i) = spName(nqOrderSp(i,n))
        if(ordermols(i,n) > 0.0) maxterm = i
        enddo 

       arrow = ' => '
       if(abs(xsH) < 1.0e-6) xsH = 0.0    ! xsH = 0 

! T_gamma_out
        ! heading
       iT_gamma_out = T_gamma_out
       write(chr4,'(i4)') iT_gamma_out

       chr24 = '  '
       chr32 =                    'k = A * T^n * exp(-E/RT)        ' 
       if(phi_T(n) > 0.0) then
       chr24 = '   Gamma@    K   Hamaker'
       chr32 = 'k = gamma*A * T^n * exp(-E/RT)  '
       endif

       ! E in input units 
       Ekjperkmol = chemrate(3,n)
       Einputunits = Ebackconversion(Ekjperkmol,iunitE)        
         
       if(iout==1) then      
       if(phi_T(n) > 0.0) write(chr24(10:13),'(i4)') iT_gamma_out
       write(filereact,'("      React. No.  Generic data:  ^n = ",1pg11.3,"   E(kJ/kmol) = ",1pg11.3," = ",1pg11.3,1x,a10,a32,  &
           & T290," A        k",a4,a24,"     xsH")') chemrate(2:3,n),Einputunits,E_units(iunitE),chr32,chr4,chr24
       
       elseif(mod(iout,20)==0) then 
       if(phi_T(n) > 0.0) write(chr24(10:13),'(i4)') iT_gamma_out
       write(filereact,'(T290," A        k",a4,a24,"     xsH")') chr4,chr24
       endif


       ! Gamma
       call RATECONSTANT(rate_nogamma,T_gamma_out,0.0,chemrate(1,n))    ! rate without coag
       rate_gamma = rate_nogamma
       gamma_at_T = 1.0
       if(phi_T(n) > 0.0) call COAGRATE(phi_T(n), rate_gamma, T_gamma_out, gamma_at_T)  ! rate with gamma

       
       Line1 = ' '
       if(itag(numbereact(n)) > 0) then
       write(Line1,'(i3,1x,i8,3x,a13," + ",a13,a4, 20(3x,1pg9.2,1x,a13) )') itag(numbereact(n)), & 
       numbereact(n), (spNameout(i)(1:13),i=1,2),arrow, (ordermols(i,n),spNameout(i)(1:13),i=4,maxterm)  
       else
       write(Line1,'(4x,i8,3x,a13," + ",a13,a4, 20(3x,1pg9.2,1x,a13) )') & 
       numbereact(n), (spNameout(i)(1:13),i=1,2),arrow, (ordermols(i,n),spNameout(i)(1:13),i=4,maxterm)  
       endif
       Length = Len_trim(Line1)

       chr20 = '((a),T283,1p5g12.3)' 
       if(phi_T(n) > 0.0 .and. gamma_at_T <= 1.0) then   ! <= 1.0 show gamma
       if(abs(xsH)> 0.0) write(filereact,chr20) Line1(1:Length), chemrate(1,n), rate_gamma, gamma_at_T, Hamean,xsH     
       if(abs(xsH)==0.0) write(filereact,chr20) Line1(1:Length), chemrate(1,n), rate_gamma, gamma_at_T, Hamean
       else
       if(abs(xsH)> 0.0) write(filereact,'((a),T283,1p2g12.3,24x,g12.3)') Line1(1:Length),chemrate(1,n),rate_gamma,xsH 
       if(abs(xsH)==0.0) write(filereact,'((a),T283,1p2g12.3)')           Line1(1:Length),chemrate(1,n),rate_gamma
       endif

! Eqnimage for Sectional reactions.  Omit reaction No. & make compact 
      j = 12           
      do i = 1, L_eqnimage
       do while(j < LineL-1) 
       j = j+1
       if(Line1(j:j+1) == '  ') j = j+1
       exit
       enddo 
      Eqnimage(n)(i:i) = Line1(j:j)
      enddo 

      return
      end
!=========================================================================
      
      subroutine SOOT_FV_D_N(fv,sootN,sootD,sootHtoC,Tcell,phi,dencell,idim1,idim2,idim3,nq1,nq2,ix1,ix2,ido,igroup)

! sums all species in Sectional Group
! R1Cmax is Range 1 maximum C number.
! if R1Cmax is -ve sum is for Cnumber <= R1Cmax
! if R1Cmax is +ve sum is for Cnumber >  R1Cmax
! if R1Cmax = 0  sum is for all Cnumber

! Treat particles as gas molecules to get number density
!  kmol/m^3 = P/(RT)
!  N/kmol = Avogadro = 6.022e26
!  N/m^3 =  X Avogadro P/(RT)  where X is mole fraction of species

!====================================================================
      use paramsmod
      use dimensmod
      use sootdenmod
      use speciesmod
      use logicmod
      use filesmod
      use PRESSMOD      
      use sectionalmod
            
      real :: fvmin = 1.0e-12    !  threshhold for showing H/C N D                          
      real :: fv(idim1), phi(idim2:idim3), SootN(idim1), SootD(idim1), sootHtoC(idim1) 
      real :: scr(nqsp1:nqLast) 
      logical :: INSECTGROUP
            
      APoR =  (Avog*1.0e3) * Pressure /Rgas

      ido = 0
      fv    = 0.0
      SootN = 0.0
      SootD = 0.0
      sootHtoC = 0.0  

      
      DO i = ix1, ix2
      if(Tcell == 0.0) cycle      ! not in flow
      scr = phi(nqsp1:nqLast)      
      wtmol = WTMOLFN(scr)
      sfv = 0.0
      sNi = 0.0
      sNiDi6 = 0.0

      sumCmole = 0.0
      sumHmole = 0.0
                         
       do nq = nq1, nq2
       if(.not. INSECTGROUP(nq,igroup)) cycle
       
       ido = 1
       iC = nqCHST(nq,1)
       fvi = Fv_(phi(nq),iC,dencell)
       Xi  = Xmolfrac_fn(wtmol, nq, phi(nq))
       rNi = Xi * APoR / Tcell 
       sfv = sfv+ fvi
       sNi = sNi + rNi

       rNiDi3 = fvi * 6.0/pi
       rNiDi6 = 0.0
       if(rNi > rtiny) rNiDi6 = rNi*(rNiDi3/rNi)**2
       sNiDi6 = sNiDi6 + rNiDi6         
       
       sumCmole = sumCmole + Xi * spCHONSI(1,nq)
       sumHmole = sumHmole + Xi * spCHONSI(2,nq)
       enddo 
      
      fv(i)    = sfv
      if(sfv > fvmin) SootN(i) = sNi * 1.e-6     ! N/cm^3    

      
       ! H/C    
      if(sfv > fvmin .and. sumCmole > 0.0) SootHtoC(i) = sumHmole/sumCmole       
      
      ! D63  = [sum(ND^6)/sum(ND^3)]^0.333         
      if(isD63) then
      sNiDi3 = sfv * 6.0/pi  
      D63 = 0.0
!      if(sNiDi3 > rtiny) D63 =   (sNiDi6 / sNiDi3) ** 0.333
      if(sfv > fvmin) D63 =   (sNiDi6 / sNiDi3) ** 0.333
      SootD(i) = D63 * 1.e9    ! nm

      ! Volume-averaged D from  fv = (pi/6)*ND^3 
      else
!      if(sNi > rtiny) SootD(i) = ( 6.0*sfv/(sNi*pi) )** 0.333  * 1.0e9   ! nm    
      if(sfv > fvmin) SootD(i) = ( 6.0*sfv/(sNi*pi) )** 0.333  * 1.0e9   ! nm    
      endif 

      ENDDO

      return
      end

!====================================================================
      
      subroutine SOOTMEAN_C(sootC,Tcell,phi,idim1,idim2,idim3,nq1,nq2,ix1,ix2,ido,igroup)

! sums all species in sectional group
! R1Cmax is Range 1 maximum C number.
! if R1Cmax is -ve sum is for Cnumber <= R1Cmax
! if R1Cmax is +ve sum is for Cnumber >  R1Cmax
! if R1Cmax = 0  sum is for all Cnumber

! Treat particles as gas molecules to get number density
!  kmol/m^3 = P/(RT)
!  N/kmol = Avogadro = 6.022e26
!  N/m^3 =  X Avo P/(RT)  where X is mole fraction of species

!====================================================================
      use paramsmod
      use dimensmod
      use sootdenmod
      use speciesmod
      use logicmod
                   
      real :: phi(idim2:idim3), SootC(idim1) 
      logical :: INSECTGROUP
      real :: Ysootmin = 1.0e-9   ! threshhold for display  = Fvmin = 1.0e-12            

      ido = 0
      SootC = 0.0
      
      DO i = ix1, ix2
      if(tcell == 0.0) cycle      ! not in flow

      Ysoot = 0.0
      YC = 0.0

       do nq = nq1, nq2
       if(.not. INSECTGROUP(nq,igroup)) cycle
       ido = 1
       Cnumber = spCHONSI(1,nq) 
       Ysoot = Ysoot + phi(nq)                        ! total soot mass 
       YC = YC + phi(nq) * Cnumber
       enddo 
      
      if(Ysoot > Ysootmin) sootC(i) = YC / Ysoot     
      ENDDO

      return
      end

!====================================================================

      function NQfromSection(iC,iH,iS,iT,m)

! gives the Species nq index for Section iC,iH,iS,iT,m
! m = 1 Section molecule;  m = 2 Section radical

      use sectionalmod
      use filesmod
             
      if(iC == 0 .or. iH==0) then
      NQfromSection = 0

      elseif(iC > nCsections .or. iH > nHsections .or. iS > nSsections .or. iT > nTsections ) then
      write(file_err,*)' Error in NQfromSection  iC,iH,iS,iT ',iC,iH,iS,iT
      write(file_err,*)'  nCsections,nHsections,nSsections,nTsections ', nCsections,nHsections,nSsections,nTsections
      stop

      else     

      NQfromSection = 0

      if(iH_table(iC,iH) == 0) return
      if(iS_table(iC,iS,iT) == 0) return
      if(iT_table(iC,iT) == 0) return

      if(m==1) then
      nq1 = nqSectionMol1
      nq2 = nqSectionRad1 - 1 
      elseif(m==2) then
      nq1 = nqSectionRad1
      nq2 = nqSectionRadn 
      else
      write(file_err,*)' m error in NQfromSECTION ',m 
      stop
      endif
              
      do nq = nq1,nq2
       if(nqCHST(nq,1)==iC .and. nqCHST(nq,2)==iH .and. nqCHST(nq,3)==iS .and. nqCHST(nq,4)==iT) then
       NQfromSection = nq
       exit
       endif
      enddo

      endif
      
      return
      end

!================================================================

      subroutine SectionfromNQ(nq,iC,iH,iS,iT,m)

! gives Section iC,iH,iS,iT,m  from nq
! m = 1 Section molecule;  m = 2 Section radical

      use sectionalmod
      use filesmod
             
      if(nq >= nqSectionrad1 .and. nq <= nqSectionradn) then
      m=2
      elseif(nq >= nqSectionmol1 .and. nq <= nqSectionmoln) then
      m=1
      else
      write(file_err,*)' Stop: nq out of section range in SectionfromNQ ',nq
      write(file_err,*)' nqsectionmol1,nqsectionmoln ', nqsectionmol1,nqsectionmoln 
      write(file_err,*)' nqsectionrad1,nqsectionradn ',  nqsectionrad1,nqsectionradn 
      stop
      endif

      iC = nqCHST(nq,1)
      iH = nqCHST(nq,2)
      iS = nqCHST(nq,3)
      iT = nqCHST(nq,4)
      return
      end

!================================================================

      subroutine PARTICLE_CHMOLES(nq1,nq2,ordermolkeep,nqorderSpkeep,particleC8,particleH8)

! particle C,H moles from C, H balance
 
      use speciesmod
      use dimensmod
            
      real(8) :: ordermolkeep(maxtermsreaction)
      real(8) :: particleC8, particleH8
      real(8) :: ReactantC, ProductC, ReactantH, ProductH
      integer :: nqorderSpkeep(maxtermsreaction)  



! Total C number in reactants 
       ReactantC =  -ordermolkeep(1) * SpCHONSI(1,nq1) 
       ReactantH =  -ordermolkeep(1) * SpCHONSI(2,nq1) 
       if(nq2 > 0) then
       ReactantC = ReactantC - ordermolkeep(2) * SpCHONSI(1,nq2)           
       ReactantH = ReactantH - ordermolkeep(2) * SpCHONSI(2,nq2)           
       endif
       
! Total C number in products without Ai, Aj
       ProductC  =  0.0
       ProductH  =  0.0
        do i = 4,maxtermsreaction
        if(nqorderSpkeep(i) > 0) then
        ProductC = ProductC + ordermolkeep(i) * SpCHONSI(1,nqorderSpkeep(i))
        ProductH = ProductH + ordermolkeep(i) * SpCHONSI(2,nqorderSpkeep(i))
        endif
        enddo     

       particleC8 = ReactantC - ProductC
       particleH8 = ReactantH - ProductH

      return
      end

!================================================================

      subroutine SectionHbalance(n)

! for Sectional reaction with ?H or ?H2 in products
! excess (+ve only)  H moles are assigned to moles H or moles H2 in products

      use dimensmod
      use reactmod
      use speciesmod
            
       xshydrogen = 0.0
       do i = 1,maxtermsreaction         
       if(nqorderSp(i,n) == 0) cycle
       xshydrogen = xshydrogen - ordermols(i,n) * spCHONSI(2,nqorderSp(i,n))
       enddo
       
       if(xshydrogen > 0.0) then
       do i = 4,maxtermsreaction    ! products
       nq = nqOrderSp(i,n) 
        if(nq == nqH2 .or. nq == nqH) then   ! search for H or H2 with 0 (?) moles 
         if(ordermols(i,n) /= 0.0) exit          
         if(nq==nqH)  ordermols(i,n) = xshydrogen         
         if(nq==nqH2) ordermols(i,n) = xshydrogen/2.0         
        exit
        endif 
       enddo
       endif

      return
      end
!================================================================

      subroutine SECTIONPRODUCTS(ioutput,Line,n,nq1,nq2,iterm,S_X,i_vacant,particleC8, particleH8,isAiorAj,ifragment,iomit)

!  interpolates particle C H S to get section Nos iC, jH and fractions frC, frH.
!================================================================
      use paramsmod
      use dimensmod
      use sectionalmod
      use reactmod
      use speciesmod
      use filesmod
      use logicmod
                          
      character (LEN=LineL) :: Line
      real(8) :: particleC8, particleH8
      integer :: iterm(9,3)

      fragfac = 1.0
      if(ifragment>0) fragfac = 2.0 

      ! Product Ai or Aj 
      ! Interpolate sections for C;  iC is section number below C 
      
      iT = iterm(9,3)  ! Product T section No. 

      iomit = 0

      call CsectionInterpolate(n,particleC8,frC,iC)      ! get iC, frC. iC = 0 if Cprod < C(section 1)          
      frCm = 1.0 - frC

       if(frCM > 0.0 .and. iT_table(iC,iT)   == 0) then
       if(ioutput==1)write(filereact,'(" Omitted reaction ",2a15," for T at product section C,T ",2i5)')spName(nq1),spName(nq2),iC,iT
       iomit = 3
       return
       endif
       if(frC  > 0.0 .and. iT_table(iC+1,iT) == 0) then
       if(ioutput==1)write(filereact,'(" Omitted reaction ",2a15," for T at product section C+1,T ",2i5)')spName(nq1),spName(nq2),iC+1,iT
       iomit = 3
       return
       endif

!     Structure S interpolation from specification S_X or reactant weighting 
      call S_SECTIONINTERPOLATE(n,Line,iterm,S_X,iS1,iS2,frS1,frS2)   ! S = iS1*frs1 + iS2*frs2;  frs2 = 1-frs1

      ii = 0
      DO k = 1,2   ! for sections S1, S2    
      
      if(k == 1) then
      frS = frS1
      iS  = iS1
      else
      frS = frS2
      iS  = iS2
      endif
      if(frS == 0.0) cycle

       if(iS_table(iC,iS,iT)   == 0) then
       if(ioutput==1)write(filereact,'(" Omitted reaction ",2a15," for S at product section C,S ",2i5)')spName(nq1),spName(nq2),iC,iS
       iomit = 2
       return
       endif
       if(iS_table(iC+1,iS,iT) == 0) then
       if(ioutput==1)write(filereact,'(" Omitted reaction ",2a15," for S at product section C+1,S ",2i5)')spName(nq1),spName(nq2),iC+1,iS
       iomit = 2
       return
       endif

      call HsectionInterpolate(iC,frC,isAiorAj,particleH8,frH,jH)   ! get jH, frH
      frHm = 1.0 - frH

       if(frHM > 0.0) then
       if(iH_table(iC,jH)     == 0) then
       if(ioutput==1)write(filereact,'(" Omitted reaction ",2a15," for H at product section C,H ",2i5)')spName(nq1),spName(nq2),iC,jH
       iomit = 1
       return
       endif 
       endif 

       if(frH  > 0.0) then
       if(iH_table(iC,jH+1)   == 0) then
       if(ioutput==1)write(filereact,'(" Omitted reaction ",2a15," for H at product section C,H+1 ",2i5)')spName(nq1),spName(nq2),iC,jH+1
       iomit = 1
       return
       endif 
       endif 

       if(frHM > 0.0) then
       if(iH_table(iC+1,jH)   == 0) then
       if(ioutput==1)write(filereact,'(" Omitted reaction ",2a15," for H at product section C+1,H ",2i5)')spName(nq1),spName(nq2),iC+1,jH
       iomit = 1
       return
       endif 
       endif 

       if(frH  > 0.0) then
       if(iH_table(iC+1,jH+1) == 0) then
       if(ioutput==1)write(filereact,'(" Omitted reaction ",2a15," for H at product section C+1,H+1 ",2i5)')spName(nq1),spName(nq2),iC+1,jH+1
       iomit = 1
       return
       endif 
       endif 
       
!     Product interpolation
      jp1 = min(jH+1,nHsections)
      nqprodij     = NQfromSection(iC,  jH,   iS,iT, isAiorAj)
      nqprodipj    = NQfromSection(iC+1,jH,   iS,iT, isAiorAj)
      nqprodijp    = NQfromSection(iC,  jp1,  iS,iT, isAiorAj)
      nqprodipjp   = NQfromSection(iC+1,jp1,  iS,iT, isAiorAj)

       if(frC>0.0 .and. frH>0.0) then                           
       nqorderSp(i_vacant+ii,n) = nqprodipjp                          ! iC+1 jH+1   
       ordermols(i_vacant+ii,n) = frC*frH*frS * fragfac
       ii = ii+1
       endif 

       if(frC>0.0 .and. frHm>0.0) then            
       nqorderSp(i_vacant+ii,n) = nqprodipj                        !  iC+1 jH 
       ordermols(i_vacant+ii,n) = frC*frHm*frS * fragfac
       ii = ii+1
       endif

       if(iC>0 .and. frCm>0.0 .and. frH>0.0) then      
       nqorderSp(i_vacant+ii,n) = nqprodijp                          !  iC jH+1  
       ordermols(i_vacant+ii,n) = frCm*frH*frS * fragfac
       ii = ii+1
       endif

       if(iC>0 .and. frHm>0.0 .and. frCm>0.0) then      
       nqorderSp(i_vacant+ii,n) = nqprodij                         !  iC jH 
       ordermols(i_vacant+ii,n) = frCm*frHm*frS * fragfac
       ii = ii+1
       endif
      ENDDO

      return
      end
!================================================================

      subroutine GETSECTIONAL_CH(Cn,HtoC,Line)

! reads sectional C H data line       
      use paramsmod
      use Logicmod
      use filesmod
      use speciesmod
      use sectionalmod
      
      real(8) :: Cn, Cnprev = 0.0
      real :: HtoC(nHsections) 
      character (LEN=LineL) ::  Line
      
       if(.not. Sectionalmethod) then
       write(file_err,*) 'Stop: Sectional species input when not sectional method '
       write(file_err,*) Line
       stop 
       endif      

      read(line(2:LineL),*,err=999,end=999) Cn, HtoC(1:nHsections)   ! real to allow e format
       if(Cn <= Cnprev) then
       write(file_err,*) 'Stop: Sectional species must be in ascending C order. Cprev, C =  ',Cnprev,Cn
       stop 
       endif
      Cnprev = Cn 
      
      ! Minimum Hmols
      if(nHsections > 1) then
      do i = 1,nHsections       
      Hmol = Cn * HtoC(i)  
       if(Hmol < sect_Hmol_min) then
       Hmol = sect_Hmol_min
       HtoC(i) = Hmol / Cn
       endif 
      enddo 
      endif 
      
       if(nHsections > 1) then
       do i = 2,nHsections
        if(HtoC(i) <= HtoC(i-1)) then
        write(file_err,*) 'Stop: Sectional species H/C must be in ascending order. ', HtoC(:)
        stop 
        endif
       enddo
       endif 


      return
999   call errorline(fileSCHEME,' 9    ')

      end

!========================================================

      subroutine SELECT_TAG_REACTIONS(ntag,NactualTagsum)

! selects reaction NUMBERS for iTagsum(1,ntag) to include in Reaction summation
      
       use dimensmod
       use paramsmod
       use filesmod      
       use logicmod 
       use reactmod
       use TAGMOD
       use sectionalmod
       use sizemod

       logical :: isTaginrange
       
       includeReactSumNo = .false.
       NactualTagsum = 0

       iitag = iTagsum(1,ntag)   ! Tag No.

       do ireactnumber = 1, numbereact(nreactionlines)
       if(iTag(ireactnumber) /= iitag) cycle

       ! Limit Reactant bin ranges if sectional species; no 3rd body reactions for Sectional
       ireactline = LineReact(ireactnumber)
       nq1 = nqorderSp(1,ireactline)
       nq2 = nqorderSp(2,ireactline)

       if(nq1 >= nqSectionmol1 .and. nq1 <= nqSectionradn) then
       call SectionfromNQ(nq1,iC1,iH1,iS1,iT1,m1)
       if(.not. isTaginrange(iC1,iTagsum(3 ,ntag),iTagsum(4 ,ntag)) .or. & 
          .not. isTaginrange(iH1,iTagsum(5 ,ntag),iTagsum(6 ,ntag)) .or. & 
          .not. isTaginrange(iT1,iTagsum(7 ,ntag),iTagsum(8 ,ntag)) ) cycle
       endif

       if(nq2 >= nqSectionmol1 .and. nq1 <= nqSectionradn) then
       call SectionfromNQ(nq2,iC2,iH2,iS2,iT2,m2)
       if(.not. isTaginrange(iC2,iTagsum(9,ntag),iTagsum(10,ntag)) .or. & 
          .not. isTaginrange(iH2,iTagsum(11,ntag),iTagsum(12,ntag)) .or. &
          .not. isTaginrange(iT2,iTagsum(13,ntag),iTagsum(14,ntag)) ) cycle
       endif                

       ! Product bin ranges - include if one product is in range or if all non-sect prods   
       ! check in terms 4 - maxtermsreaction for sectional species
       inrange = 0
       kountsect = 0
        do i = 4, maxtermsreaction    
        nqp = nqorderSp(i,ireactline)
        if(nqp < nqSectionmol1 .or. nqp > nqSectionradn) cycle
        kountsect = kountsect +1
        call SectionfromNQ(nqp,iCp,iHp,iSp,iTp,mp)
        if(isTaginrange(iCp,iTagsum(15 ,ntag),iTagsum(16 ,ntag)) .and. & 
           isTaginrange(iHp,iTagsum(17 ,ntag),iTagsum(18 ,ntag)) .and. &
           isTaginrange(iTp,iTagsum(19 ,ntag),iTagsum(20 ,ntag)) ) then
        inrange = nqp
        exit
        endif
        enddo
       if(kountsect > 0 .and. inrange==0) cycle  ! if there are Sectional prods and none are in range

       includeReactSumNo(ireactnumber) = .true.  ! included if tagged and satisfying above conditions
       NactualTagsum = NactualTagsum + 1   ! actual number of reactions included in sum
       enddo      

      return
      end

!================================================================

      subroutine CHmols(n,Creact,Hreact,Creactsect,Hreactsect,Cprodsect,Hprodsect)

! counts C H mols. Does not exclude reaction
!==================================================================================================

      use dimensmod 
      use reactmod
      use speciesmod
      use filesmod
      use sectionalmod
      use Tagmod
                        
      real(8) :: Creact, Hreact,Cprod, Hprod, Cmol, Hmol 
      real(8) :: Creactsect, Hreactsect,Cprodsect, Hprodsect
      logical :: isTaginrange
      Creact=0.0; Hreact=0.0;Cprod=0.0; Hprod=0.0 
      Creactsect=0.0; Hreactsect=0.0;Cprodsect=0.0; Hprodsect=0.0


      irnumber = numbereact(n)
      do i = 1,3
      nq = nqorderSp(i,n)
      if(nq==0 .or. nq == nq3rd(n)) cycle    ! 3rd body not included
      Cmol = - ordermols(i,n) * spCHONSI(1,nqorderSp(i,n)) 
      Hmol = - ordermols(i,n) * spCHONSI(2,nqorderSp(i,n)) 
      Creact = Creact + Cmol
      Hreact = Hreact + Hmol
      if(nq>=nqsectionmol1 .and. nq<=nqSectionN) then 
      Creactsect = Creactsect + Cmol
      Hreactsect = Hreactsect + Hmol
      endif
      enddo

      do i = 4, maxtermsreaction
      nq = nqorderSp(i,n)
      if(nq==0 .or. nq == nq3rd(n)) cycle    ! 3rd body not included
      Cmol =  ordermols(i,n) * spCHONSI(1,nqorderSp(i,n)) 
      Hmol =  ordermols(i,n) * spCHONSI(2,nqorderSp(i,n)) 
      Cprod = Cprod + Cmol    ! redundant: Cprod = Creact
      Hprod = Hprod + Hmol    ! not redundant because of H imbalances
      if(nq<nqSectionmol1 .or. nq>nqSectionN) cycle
      
      ! only include Cprodsect Hprodsect for products in C,H,T specified bin range
      call SectionfromNQ(nq,iCp,iHp,iSp,iTp,mp)
        if(isTaginrange(iCp,iTagsum(15 ,Ltag),iTagsum(16 ,Ltag)) .and. & 
           isTaginrange(iHp,iTagsum(17 ,Ltag),iTagsum(18 ,Ltag)) .and. &
           isTaginrange(iTp,iTagsum(19 ,Ltag),iTagsum(20 ,Ltag)) ) then
        Cprodsect = Cprodsect + Cmol
        Hprodsect = Hprodsect + Hmol
        ! write(file_err,*)'ir,coef,nC ',irnumber,ordermols(i,n),spCHONSI(1,nqorderSp(i,n))
        endif
      enddo
                        
      return
      end
!===============================================================                              
