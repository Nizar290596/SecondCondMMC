
      subroutine READSPECIESNUMERICAL(iterperstep)

! Read solver numerical parameters
!=======================================================================
      use paramsmod
      use dimensmod
      use filesmod
      use logicmod
      use PHIMINMAXMOD
      use numericmod
      use focusmod
      use sizemod
      use twodimmod
!---------------------------------------------------------------

      call findheader(fileIN,isfound,'SOLVE_SECTION')
      if(isfound == 0) stop

      call findots(fileIN)

!   Solver                     Needs                            Working for
! 1  = ADI solver          large souces                     Flameletlib, laminar, CMC  
! 2  = Newton (NR) per Xi  Jacobian                         Flameletlib, laminar, CMC        
! 3  = VODE                Jacobian                         Flameletlib, CMC 
! 4  = Newton all Xi       Jacobian                         Flameletlib, CMC
! 12 = ADI/NR              large source + Jacobian 
! 13 = ADI/VODE            large source + Jacobian  

! Automatic


      ispecsolver = 2
 

! Focus criterion (residual/average residual) and max focus iterations for N-R 
      maxfocusit = 0
            

! Tolerance       
!      read(fileIN,*,err=999)  tolerancevode
      
      !  Tolerance for sum of species: stop for  Ysum > 1 + tolerance
      read(fileIN,*,err=999)  Ysumtolerance        
      Ysumtolerance = max(Ysumtolerance, 1.0e-3)  ! for rounding     
       
!    Default phimag
!    Phimagnitudemin = lowest permitted value of phimagnitude set       
      read(fileIN,*,err=999) phimagdefault, phimagnitudemin
      phimagnitude(nqsp1:nqLast) = max(phimagdefault,phimagnitudemin)    ! default

      !  Update Phimag 
      ! if iupdatephimag > 0 then update every iupdatephimag iteration or newstep.
      read(fileIN,*,err=999) iupdatephimag

!   Max.Phi DownMovement/Phi; Max.Phi Up Movement/Phimax
      read(fileIN,*,err=999)  fracmovedown, fracmoveup
!@      if(ispecsolver == 3 .or. unsteady) then
!      fracmovedown = 0.0
!      fracmoveup = 0.0
!      endif
      
!      if(unsteady .or. ispecsolver==3) &
!      write(fileRPT,'(/" No Restrictions on phi movement for Unsteady or Flamelet Library or Vode solver.")') 
       
      write(fileRPT,'(/" Fraction of phi movement down, fraction of phimagnitude movement up ",2g12.3)') fracmovedown, fracmoveup
      if(fracmovedown <=0.0)  write(fileRPT,'(/" NO restriction on phi movement down.")') 
      if(fracmoveup <=0.0)    write(fileRPT,'(/" NO restriction on phi movement up.")') 

! URF for species (all equal at present) 
      read(fileIN,*,err=999)  URFspeciesdflt

      if((unsteady .and. iterperstep==1) .or. ispecsolver==3) then
      URFspeciesdflt = 1.0
      write(fileRPT,'(/" Species URF default = 1 for Unsteady with iter/step=1 or for Vode ")') 
      endif

      write(fileRPT,'(/" Species URF default ",f12.5)') URFspeciesdflt
      call RLimitin(URFspeciesdflt,rsmall,1.0,"URFspeciesdflt ")

      return
999   call errorline(filein,' 5a   ')

      end
      
!====================================================================

      subroutine ALLOCSPEC

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
      use elementsmod

      write(fileReact,'(/" Species & Reaction Allocations ")')
      write(fileReact,'(" nqlast       ",i6)') nqlast
      write(fileReact,'(" maxReactLines ",i8)') maxReactLines
      write(fileReact,'(" nAux         ",i6)') nAux
      write(fileReact,'(" maxtermsreaction ",i6)') maxtermsreaction
      write(fileReact,'(" max products ",i6)') maxtermsreaction-3
      write(fileReact,'(" maxFORD      ",i6)') maxFORD

      allocate (solvenq(nqlastm), printphis(nqlastm), printspart(nqsp1:nqlastm), &
      pltphis(nqlastm), iprnsoot(nqsp1:nqLastm), stat=istat) 
      call STATCHECK(istat,'solvenq   ')     
      mbSpecies = mbSpecies + nqlast*3 + (nqlast-nqsp1+1)*2
!--------------------------------        

      allocate (spName(0:nqlastmp1),spFormula(0:nqlastmp1), stat=istat)          
      call STATCHECK(istat,'          ')      

      allocate (sumsource(nqsp1:nqLastm), stat=istat) 
      call STATCHECK(istat,'          ')       

      allocate (spMW(nqsp1:nqlastm), spRelDiff(nqlastm), spEnth(0:2,nqsp1:nqlastm,2), stat=istat) 
      call STATCHECK(istat,'          ')       

      allocate (spCHONSI(maxelement,nqsp1:nqlastm),stat=istat) 
      call STATCHECK(istat,'          ')       
      mbSpecies = mbSpecies + (nqlastm-nqsp1+1)*maxelement 

      spName = ' '
      solvenq = .false.
      sumsource = 0.0
      spCHONSI  = 0.0 
      spFormula = '      '
!--------------------------------------------

      allocate (resmon(0:nqLastm), phimaxnq(0:nqLastm),stat=istat)
      call STATCHECK(istat,'phimon    ')       
      mBmonitor = mBmonitor + maxmon*nqlastm*4 + nqLastm*4 

      allocate (noverNq(nqsp1:nqLastm), nunderNq(nqsp1:nqLastm), nmoveNq(nqsp1:nqLastm),stat=istat) 
      call STATCHECK(istat,'novernq   ')       

      allocate (nsignificantreact(maxReactLines),nmaxreact(maxReactLines),nminreact(maxReactLines), &
                nstorereact(maxReactLines),storerate(maxReactLines),stat=istat) 
      call STATCHECK(istat,'nsignif   ')       

      allocate (star(nqsp1:nqLastm), temp(nqsp1:nqLastm), stat=istat) 
      call STATCHECK(istat,'star      ')       

      allocate (iFORD(maxReactLines), specexpon(maxtermsreaction,maxFORD),itermexpon(maxtermsreaction,maxFORD), stat=istat) 
      call STATCHECK(istat,'FORD      ')       

      Phimaxnq    = 0.0
      resmon      = 0.0
      nmovenq  = 0
      novernq  = 0
      nundernq = 0

      iFord = 0
      specexpon = 1.0                   
      itermexpon = 0
!---------------------------------------

      allocate (chemrate(3,maxReactLines), reactratesum(maxReactLines), fwdratesum(maxReactLines), &
                auxrate(maxauxparams,nAux), iauxtype(nAux), phi_T(maxReactLines),   stat=istat) 
      call STATCHECK(istat,'chemrate  ')        

      allocate (reactmon(maxReactLines), fwdreactmon(maxReactLines), stat=istat) 
      call STATCHECK(istat,'reactmon  ')        

      allocate (Eqnimage(maxReactLines), markreact(maxReactLines), stat=istat) 
      call STATCHECK(istat,'eqnimage  ')        

      allocate (numbereact(0:maxReactLines+1),Linereact(maxReactLines), nforward(maxReactLines), ireactype(maxReactLines),  &
              nqorderSp(maxtermsreaction,maxReactLines), ordermols(maxtermsreaction,maxReactLines),reversible(maxReactLines), &
              iauxindex(maxReactLines),reactantC(maxReactLines), stat=istat) 
      call STATCHECK(istat,'numbereact')        

      markreact = ' '
      nqorderSp = 0
      ordermols = 0.0
      ireactype = 0    ! default for elementary reactions
      phi_T  = 0.0  ! default makes gamma coag = 1
             
      if(iallelementary /= 1)  then
      allocate (expon(nqSp1:nqlastm,maxReactLines),stat=istat)     
      call STATCHECK(istat,'expon     ')        
      else
      allocate (expon(1,1),stat=istat)     
      call STATCHECK(istat,'expon     ')        
      endif
      
      allocate (reactmol(0:nqlastmp1),rLnKequil(0:maxReactLines),reactmolnet(maxReactLines),stat=istat) 
      call STATCHECK(istat,'reactmol  ')        

      allocate (specnetmass(maxtermsreaction,maxReactLines), &
                nqspec(maxtermsreaction,maxReactLines), nq3rd(maxReactLines), iterm3rd(2,maxReactLines),stat=istat) 
      call STATCHECK(istat,'specmol   ')        

      allocate (conc(0:nqlastmp1), PhitoConc(nqsp1:nqlastmp1), concout(nqsp1:nqlastmp1), Ratenq(0:nqlastmp1), stat=istat) 
      call STATCHECK(istat,'conc      ')        
      mbSpecies = mbSpecies + nqlast*3 

      allocate (nosource(nqsp1:nqLastm), stat=istat) 
      call STATCHECK(istat,'nosource')         
      mbSpecies = mbSpecies + (nqLastm-nqsp1+1) 

      concout = 0.0 
      reactratesum = 0.0
      fwdratesum   = 0.0
      reactmon     = 0.0
      fwdreactmon  = 0.0
       
! Thermodynamic data
      allocate (thermcoeff(7,2,nqsp1:nqlastm), T_thermorange(3,nqsp1:nqlastm), stat=istat) 
      call STATCHECK(istat,'thermcoeff')                
      mbSpecies = mbSpecies + (nqlastm-nqsp1+1) * 15

!--------------------------------------
! Data
      call ALLOCATABLEDATA

      return
      end

!================================================================================

      subroutine SPECIESCOUNT

! Count number of species
!================================================================================
      use paramsmod
      use dimensmod 
      use filesmod 
      use logicmod
      use sectionalmod
      use sizemod
      use elementsmod
      use SOOTDENMOD
                        
      character (LEN=LineL) ::  Line
      real(8) :: Cn = 0.0, Cnprev = 0.0 


!---------------------------------------------------------------

! Species
      call findheader(fileSCHEME,isfound,'SPECIES_SECTION')
      if(isfound == 0) stop

      call findots(fileSCHEME)
      call READ_ELEMENTS(Linechange)

      call findots(fileSCHEME)
      read(fileSCHEME,*,err=999) iSectionalmethod,isradical       ! 1 or 0  sectional method, radicals present 
      if(iSectionalmethod > 0) Sectionalmethod   = .true.  ! needed for allocation of sectionC
      isdoubleAi = 2 * isradical

      read(fileSCHEME,*,err=999) nHsections, nSsections, nTsections   ! No. of H, Structure, Type sections
      nHsections = max(1,nHsections)
      nSsections = max(1,nSsections)
      nTsections = max(1,nTsections)

      if(sectionalmethod) call MAKE_HST_RANGE()
      
      call findots(fileSCHEME)
      call findots(fileSCHEME)
      nq = 0
      nCsections = 0
      NsectionsAi  = 0
      
      DO WHILE (.true.)
      Line = ' '
      call READLINE(fileSCHEME,Line,iend,Linechange)
      if(iend == 1) exit
      call NOLEADINGBLANKS(Line)

      ! SECTION
      ! Each Cn must have at least one H & S & T section
      ! For each Cn No. of sections (species) is nH X nS X nT where nH is No. of H sections at this Cn etc
      IF(Sectionalmethod .and. Line(1:1) == '@') THEN
      read(line(2:LineL),*) Cn
       if(Cn <= Cnprev) then   ! C must be in ascending order (also test for Aj not listed any more)
       write(file_err,'(/" STOP: Section ",i4," C = ",f12.0," <= Cprevious = ",f12.0)')nCsections, Cn,Cnprev
       stop 
       endif
      Cnprev = Cn

      nCsections = nCsections + 1
      if(nCsections == 1) nqSectionmol1 = nq+nqsp1-1   ! Define nqSectionmol1 (nq+1-1  +1 for next nq  -1 because of balance species)

      ! Check how many T sections for this Csection
      nTatC = 0
       do iT = 1,nTsections
       if(nCsections >= iT_Ranges(1,iT) .and. nCsections <= iT_Ranges(2,iT) ) nTatC = nTatC + 1 
       enddo 
       if(nTatC == 0) then
       write(file_err,'(/" STOP: There is no T section range covering C Section ",i4)')nCsections
       stop 
       endif

      ! Check how many H sections for this Csection
      nHatC = 0
       do iH = 1,nHsections
       if(nCsections >= iH_Ranges(1,iH) .and. nCsections <= iH_Ranges(2,iH) ) nHatC = nHatC + 1 
       enddo 
       if(nHatC == 0) then
       write(file_err,'(/" STOP: There is no H section range covering C Section ",i4)')nCsections
       stop 
       endif

      ! Check how many S sections for this Csection
      nSatC = 0
       do iS = 1,nSsections
       if(nCsections >= iS_Ranges(1,iS) .and. nCsections <= iS_Ranges(2,iS) ) nSatC = nSatC + 1 
       enddo 
       if(nSatC == 0) then
       write(file_err,'(/" STOP: There is no S section range covering C Section ",i4)')nCsections
       stop 
       endif

      ! Count
      do iH = 1,nHsections
      if(nCsections < iH_Ranges(1,iH) .or. nCsections > iH_Ranges(2,iH) ) cycle
       do iS = 1,nSsections
       if(nCsections < iS_Ranges(1,iS) .or. nCsections > iS_Ranges(2,iS) ) cycle
        do iT = 1,nTsections 
        if(nCsections < iT_Ranges(1,iT) .or. nCsections > iT_Ranges(2,iT) ) cycle
        if(maxTforSrange < iT .and. iS>1) cycle
        NsectionsAi = NsectionsAi + 1
        nq = nq + 1
        enddo         
       enddo 
      enddo 
       
      ELSE 
      nq = nq + 1
      ENDIF
      ENDDO

      if(Sectionalmethod .and. nCsections < 2) then
      write(file_err,'( " Sectional method indicated by Cn > 0: need nCsections > 1 ")')
      stop
      endif

      if(Sectionalmethod) then
      if(isdoubleAi==2) nq = nq + NsectionsAi

      write(filereact,'(/ " No. of C sections is      ",i4)') nCsections
      write(filereact,'(  " Max. No. of H sections is ",i4)') nHsections
      write(filereact,'(  " Max. No. of S sections is ",i4)') nSsections
      write(filereact,'(  " Max. No. of T sections is ",i4)') nTsections
      write(filereact,'(  " Actual No. of Sections Ai, Aj each is ",i4)') NsectionsAi

      nCHSTsections(1) = nCsections 
      nCHSTsections(2) = nHsections 
      nCHSTsections(3) = nSsections 
      nCHSTsections(4) = nTsections 
      endif
         
      nq = nq - 1     ! because balance species is included in list
      
      maxspecies = nq                           ! max. No. of variables in list
      nqlast     = maxspecies + nqsp1 - 1       ! definition here
      nqlastm    = nqlast + mextra 
      nqlastmp1  = nqlast + mextrap1 
      nqBalance  = nqlastm
      nqM        = nqLastmp1

      if(Sectionalmethod) then
      call MAKE_HST_TABLES()
       HaverageN =  float(sum(iH_Table))/nCsections
       SaverageN =  float(sum(iS_Table))/nCsections
       TaverageN =  float(sum(iT_Table))/nCsections
      endif
      
!----------------------------------
! Sectional Method - need to allocate here to read Gamma data if it exists

      if(SectionalMethod) then 
      allocate (sectionC(nCsections), sectionH(nCsections,nHsections),Hamaker_jHfac(nHsections), &
                 Hamaker_iSfac(nSsections),Hamaker_iTfac(nTsections),Cgammadata(nCsections),sectiondensity(0:nCsections), stat=istat)
      call STATCHECK(istat,'sectionC  ')                
      mbSpecies = mbSpecies + nCsections*2 + nCsections*nHsections
      
      sectiondensity = Pdensitydefault 
      
      allocate (Gamma2data(nCsections,nCsections,maxGammaTables), stat=istat)  ! allocate even if no data
      call STATCHECK(istat,'Gamma2    ')                
      mbSpecies = mbSpecies + nCsections*nCsections*maxGammaTables      
      else      
      allocate (sectiondensity(0:1), stat=istat)
      call STATCHECK(istat,'sectionden')                
      endif


      return
999   call errorline(fileSCHEME,' 7    ')

! + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + 
CONTAINS

! + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + 
       subroutine MAKE_HST_RANGE()

      ! iH_Ranges  and similarly for S and T
      ! Cmin   Cmax        
      !   2      5      H1 range
      !   3      4      H2 range
      !   4      5      H3 range

      ALLOCATE (iT_Ranges(2,nTsections),iH_Ranges(2,nHsections),iS_Ranges(2,nSsections), stat=istat) 
      call STATCHECK(istat,'iT_ranges ')     
      read(fileSCHEME,*,err=999) iH_ranges    ! section No. ranges for each H
      read(fileSCHEME,*,err=999) iS_ranges, maxTforSrange    ! section No. ranges for each S; T range (0 = all) 
      read(fileSCHEME,*,err=999) iT_ranges    ! section No. ranges for each T

      if(maxTforSrange < 1) then
      write(file_err,'(/" Error STOP: Need T Section max for S rang >= 1 ")')
      stop
      endif

       
      iHSTerror = 0 
      do iH = 1,nHsections
      if(iH_ranges(1,iH) > iH_ranges(2,iH)) iHSTerror = 1
      if(iH_ranges(1,iH) < 1) iHSTerror = 1
      enddo

      iHSTerror = 0 
      do iS = 1,nSsections
      if(iS_ranges(1,iS) > iS_ranges(2,iS)) iHSTerror = 2
      if(iS_ranges(1,iS) < 1) iHSTerror = 2
      enddo

      iHSTerror = 0 
      do iT = 1,nTsections
      if(iT_ranges(1,iT) > iT_ranges(2,iT)) iHSTerror = 3
      if(iT_ranges(1,iT) < 1) iHSTerror = 3
      enddo

       if(iHSTerror>0) then   ! range must be in ascending order 
       write(file_err,'(/" Error STOP: Need H,S,T Section ranges min <= max and min > 0")')
       if(iHSTerror==1) write(file_err,*) iH_ranges
       if(iHSTerror==2) write(file_err,*) iS_ranges
       if(iHSTerror==3) write(file_err,*) iT_ranges
       stop 
       endif

       return      
999   call errorline(fileSCHEME,'T-rang')
       end subroutine MAKE_HST_RANGE
! + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + 

      subroutine MAKE_HST_TABLES()

! Make H_table  shows which Csections have H sections
!  C    H1  H2  H3
!  1     0   1   0
!  2     1   1   0
!  3     0   0   1
! similarly for S and T

      ALLOCATE (iH_table(0:nCsections,nHsections),iS_table(0:nCsections,nSsections,nTsections),iT_table(0:nCsections,nTsections), stat=istat) 
      call STATCHECK(istat,'iT_ranges ')     
      iH_table = 0; iS_table = 0; iT_table = 0

      do iH = 1,nHsections
      iH_ranges(2,iH) = min(iH_ranges(2,iH) , nCsections)
      enddo

      do iS = 1,nSsections
      iS_ranges(2,iS) = min(iS_ranges(2,iS) , nCsections)
      enddo

      do iT = 1,nTsections
      iT_ranges(2,iT) = min(iT_ranges(2,iT) , nCsections)
      enddo

      do iH = 1,nHsections
      iC1 = iH_Ranges(1,iH)   
      iC2 = iH_Ranges(2,iH)   
      iH_table(iC1:iC2,iH) = 1
      enddo
      iH_table(0,:) = iH_table(1,:)  ! iH_table(0,:) = iH_table(1,:) to allow iC=0 interpolations

      do iT = 1,nTsections
      do iS = 1,nSsections
      iC1 = iS_Ranges(1,iS)   
      iC2 = iS_Ranges(2,iS)   
      if(iS > 1 .and. iT > maxTforSrange) cycle
      iS_table(iC1:iC2,iS,iT) = 1    ! = 1 even if iT_table(iC1:iC2,iT)=0 
      enddo
      enddo
      iS_table(0,:,:) = iS_table(1,:,:)  ! iS_table(0,:) = iS_table(1,:) to allow iC=0 interpolations

      do iT = 1,nTsections
      iC1 = iT_Ranges(1,iT)   
      iC2 = iT_Ranges(2,iT)   
      iT_table(iC1:iC2,iT) = 1
      enddo
      iT_table(0,:) = iT_table(1,:)  ! iT_table(0,:) = iT_table(1,:) to allow iC=0 interpolations

!
! Make nqCHST Table Sections against nq
!  nq    C   H   S   T 
!  136   1   1   1   2
!  137   2   1   1   3
      ALLOCATE (nqCHST(nqLast,4), stat=istat) ! allocate all nq to get 0 0 0 0 for non-sectional nq 
      call STATCHECK(istat,'iT_ranges ')     
      nqCHST = 0 
      nq = nqSectionmol1 - 1 

! Array Sequence iC iH iS iT
      do idouble = 1,isdoubleAi
       do iC = 1,nCsections
        do iH = 1,nHsections
        if(iH_table(iC,iH) == 0) cycle 
         do iS = 1,nSsections
          do iT = 1,nTsections
          if(iT_table(iC,iT) == 0) cycle 
          if(iS_table(iC,iS,iT) == 0) cycle 
          nq = nq + 1
          nqCHST(nq,1) = iC         
          nqCHST(nq,2) = iH         
          nqCHST(nq,3) = iS         
          nqCHST(nq,4) = iT         
          enddo
         enddo
        enddo
       enddo
      enddo

      return
      end subroutine MAKE_HST_TABLES

      end
!=======================================================================

      subroutine SPECIESIN

!=======================================================================
      use paramsmod
      use dimensmod
      use filesmod
      use speciesmod
      use logicmod
      use propsmod
      use fuelmod
      use phiminresmod
      use PHIMINMAXMOD
      use numericmod
      use sectionalmod
      use elementsmod
      use SOOTDENMOD
                  
      real :: HtoC(nHsections) 
      character (LEN=LineL) ::  Line
      character (LEN=30)  ::  formulae(1), formulaParticleEnth

      real(8) :: Cn, Hn
      logical :: isSection
      
! SPECIES SECTION
      call findheader(fileSCHEME,isfound,'SPECIES_SECTION')
      if(isfound == 0) stop
!@      write(*,*) 'Reading Species ' 
       
      linenumber = 0
      nq   = nqsp1-1


! Info already in
      call findots(fileSCHEME)
      read(filescheme,*) Line
      read(filescheme,*) Line

      call findots(fileSCHEME)
      read(fileSCHEME,*,err=999)  idum   ! sectionalmethod
      read(fileSCHEME,*,err=999)  idum   ! nH nS nT sections 
      read(fileSCHEME,*,err=999)  idum   ! read H ranges
      read(fileSCHEME,*,err=999)  idum   ! read S ranges
      read(fileSCHEME,*,err=999)  idum   ! read T ranges

! Species from which to use enthalpy data for polymers  (is copied on mass basis)
      call READLINE(fileSCHEME,line,iend,Linechange)
      call NOLEADINGBLANKS(Line)
      call READMULTISPECIES(Line,kspecies,formulae,1,30,lastchar)            
      formulaParticleEnth = formulae(1) 

!  Minimum H mols in Sectional
      read(fileSCHEME,*,err=999)  sect_Hmol_min
      write(filereact,'( " HtoC minimum ",f7.3)') sect_Hmol_min

! Mean-gas Laminar Dynamic Viscosity (kg/(m-s) = a + bT^2 +cT^3 +dT^4 + e for T(K)
      call findots(fileSCHEME)
      read(fileSCHEME,*,err=999)  ViscCoefBase
      if( sum(ViscCoefBase) <= 0.0) then
      write(file_err,*) 'Stop: Viscosity coefficients <= 0 '
      stop
      endif

! Schmidt No. for N2 in Mean-gas; Prandtl No. for mean-gas
      read(fileSCHEME,*,err=999) Schmidtbase, Prandtlbase 
      if( Schmidtbase <= 0.0 .or. Prandtlbase <= 0.0) then
      write(file_err,*) 'Stop: Schmidt or Prandtl are zero ',Schmidtbase, Prandtlbase
      stop
      endif

      write(filereact,'(/" Mean-gas is a major species which represents average in field." )') 
      write(filereact,'(" Mean-gas laminar dynamic viscosity (kg/(m-s) = a + bT^2 +cT^3 +dT^4 + e for T(K) ")') 
      write(filereact,'("                           a            b             c           d            e")')
      write(filereact,'(20x,1p5g13.3)') ViscCoefBase

      write(filereact,'( " Prandtl No. for mean-gas       ",f7.3)') Prandtlbase
      write(filereact,'( " Schmidt No. for N2 in Mean-gas ",f7.3)') Schmidtbase
      write(filereact,'(/" For each species i, D/D-N2 = [Di in mean-gas] / [D-N2 in mean-gas]")') 

! Balance species; Is not solved (but may react). Obtained as 1 - sum of other gas mass fractions. 
      call READLINE(fileSCHEME,line,iend,Linechange)
      call NOLEADINGBLANKS(Line)
      call READMULTISPECIES(Line,kspecies,formulae,1,30,lastchar)            
      Spbalance = formulae(1) 

      write(filereact,'(/" Balance Species is ",a24)') Spbalance
      write(filereact,'(/" Default particle density for non-sectional species ",f8.1," kg/m^3")') Pdensitydefault


! READ SPECIES LIST  for all species and for C sections        

      call findots(fileSCHEME)

      irepeated = 0
      ispbalancefound = 0 

      isect = 0
      iC =0
      isectionalfound=0
             
      DO WHILE (.true.)

      Line = ' '
      call READLINE(fileSCHEME,Line,iend,Linechange)
      Linenumber = Linenumber + Linechange
      if(iend == 1) exit
      call NOLEADINGBLANKS(Line)

      if(ispbalancefound == 1) then   ! continue nq count after Balance species found
      ispbalancefound = 2
      nq = nqhold
      endif
       
      nq = nq+1

      if(nq > nqlastm) then
      write(file_err,*)' No. of Species exceeds ',nqlastm
      stop
      endif

      ! SECTIONS  - data listed for each C section; @ in 1st col. 
      IF(Line(1:1) == '@') then
      call GETSECTIONAL_CH(Cn,HtoC,Line)
      if(isectionalfound == 0) isectionalfound = nq   ! just to check that sectional species continue not gas phase
      isSection = .true.
      iC = iC + 1
      sectionC(iC) = Cn
      call COPYtoHST(nHsections,nSsections,nTsections)       ! Data for this nq and copied to H, S, T sections if Sectional

      ! GAS PHASE
      ELSE
       if(isectionalfound > 0) then    ! because section list terminated by no further species
       write(file_err,*)' Stop: Sectional species list must be last in species list. '
       stop
       endif

      call READSPECIESandFORMULA(Line,spName(nq),spFormula(nq),length)  ! reads name and formula (if different from name)          
      call READFORMULA(spFormula(nq),spCHONSI(1,nq),30)  ! fills spCHONSI

      call WRONGNAMECHECK(nq,irepeated)                       ! Check for repeated Species Name
      call LOCATEBALANCESPECIES(nq,nqhold,ispbalancefound)    ! Hold Balance species nq for end of list
      call ZEROELEMENTCHECK(nq)                               ! Zero elements check (for species)
      isSection = .false.
      call COPYtoHST(1,1,1)  ! Data for this nq and copied to H, S, T sections if Sectional
      ENDIF


      ENDDO
      solvenq(nqlastm) = .false.      ! Balance species but solvenq(nqlastm) not used anyway
      call WRONGNAMECHECK(nqlastm,irepeated)      ! Check for repeated Balance Species Name

! END of READ SPECIES LIST         
!++++++++++++++++++++++++++++++++++++++++++++

!  Make Sectional Radicals list
      call MAKESECTIONRADICALS

! Checks
      if(irepeated > 0) then
      write(file_err,*)' STOP: Repeated Species. See above. '
      stop 
      endif

      if(ispbalancefound == 0) then
      write(file_err,*)' Need Balance Species specified in list ',spbalance
      stop
      endif

      nqLastnonSect = nqLast     ! nq last non-sectional species solved
      if(SectionalMethod) nqLastnonSect = nqSectionmol1 -1
      

! Sectional enthalpy and gamma data check
     if(SectionalMethod) then

      ! Identify species nq for particles enthalpy
      nqParticleEnth = Indexspecies(formulaParticleEnth)
      if(nqParticleEnth  == 0) then
      write(file_err,*)' Species for Section enthalpy is not in list ',formulaParticleEnth
      stop
      endif

     if(igammadimension == 2) call GammaDataCheck
     endif

!--------------------------------------
! Element balance: Solve for Fuel species & O2 only
!      if(EBUmixfrac) then
!      kount=0
!      do nq = nqsp1, nqLast
!      if(solvenq(nq)) kount=kount+1
!      enddo
!       if(kount /= 2) then
!       write(file_err,*)' Element balance: solve only 2 species: Fuel and O2 '
!       write(*,*)     ' Element balance: solve only 2 species: Fuel and O2 '
!       stop
!       endif
!      endif


! Molecular weight of species all treated as gases  
      do nq = nqsp1, nqLastm
      spMW(nq) = SPECIESMW_fn(nq)
      enddo

! Special index numbers
      nqCO2 = IndexSpecies("CO2 ")
      nqCO  = IndexSpecies("CO  ")
      nqO2  = IndexSpecies("O2  ")
      nqH2O = IndexSpecies("H2O ")
      nqH2  = IndexSpecies("H2  ")
      nqN2  = IndexSpecies("N2  ")
      nqH   = IndexSpecies("H   ")

!----------------------------------
! OUTPUTS
      write(fileReact, '(/1x,79(1h-))')
      write(fileReact,'(/" No. of species solved (incl. Sections, not Balance) ",i6)')  maxspecies
      write(fileReact,'( " First species index nqSp1                           ",i6)')  nqSp1
      write(fileReact,'( " Last non-sectional species index (excl Balance)     ",i6)')  nqLastnonSect
      write(fileReact,'( " Last species index (excl Balance) nqLast            ",i6)')  nqLast

! Sectional Method  (if iSectionalmethod) 
      if(Sectionalmethod) then
      write(fileReact,'(//" Sectional Method ")')
      write(fileReact,'(/" No. of C sections (molecule)  ",i8)') nCsections
      write(fileReact,'( " Index of 1st molecule         ",i8)') nqSectionmol1
      write(fileReact,'( " Index of 1st radical          ",i8)') nqSectionrad1
      write(fileReact,'( " Index of enthalpy species     ",i8)') nqParticleEnth

      write(fileReact,'(//20x," T-structure Composition.")')
      do iT = 1, nTsections  
      write(fileReact,'(/" T-",i2)')iT
      write(fileReact,'("   C   H   S")')
       do nq = nqSectionmol1, nqSectionmoln
       if(nqCHST(nq,4) == iT) write(filereact,'(3i4)') nqCHST(nq,1),nqCHST(nq,2),nqCHST(nq,3)
       enddo       
      enddo

      else
      write(fileReact,'(//" NOT Sectional Method ")')
      endif

      write(fileReact, '(/1x,79(1h-))')
      write(fileReact,'(/"    Required Species & Index Nos. ")')
      write(fileReact,'(" CO2      ",i8)') nqCO2
      write(fileReact,'(" CO       ",i8)') nqCO
      write(fileReact,'(" O2       ",i8)') nqO2
      write(fileReact,'(" H2O      ",i8)') nqH2O
      write(fileReact,'(" H2       ",i8)') nqH2
      write(fileReact,'(" H        ",i8)') nqH
      write(fileReact,'(" N2       ",i8)') nqN2

      if(nqCO2 == 0 .or. nqCO == 0 .or. nqO2 == 0 .or. nqH2O == 0 .or. nqH2 == 0) then
      write(file_err,'(" Missing species from base set CO2,CO,O2,H2O,H2 ")') 
      stop
      endif

      write(fileReact, '(/1x,79(1h-))')
      write(fileReact,'(/"    SPECIES ELEMENTS & MW. ")')

      do nq = nqsp1, nqlastm
      if(mod(nq-nqsp1+20,20)==0 .or. (sectionalmethod .and. (nq==nqSectionmol1 .or. nq==nqSectionrad1))) then
       if(sectionalmethod .and. nq>=nqSectionmol1) then
       write(fileReact,'(/"    C  H  S  T",14x,"Formula",27x,"Index     MW",2x,2a12,8x,a12,20a6)') elementname(1:maxelement)
       else
       write(fileReact,'(/28x,"Formula",27x,"Index     MW",2x,2a12,8x,a12,20a6)') elementname(1:maxelement)
       endif
      endif
      write(fileReact,'(1x,a24,5x,a30,i6,1p3g14.6,0pf12.0,20f6.0)') spName(nq),spFormula(nq),nq,spMW(nq),spCHONSI(1:maxelement,nq)
      enddo
      write(fileReact,'(1x,(a),35x,i6)') spName(nqM),nqM 

     ! Check species with same No. of C H O N elements
     call SAME_ELEMENTCHECK

!      call CombinedAiAjPrint  ! 
      return

!----------------------------------------------
999   call errorline(fileSCHEME,'specie')

998   write(file_err,*) ' Error in Species specification at line ',linenumber
      write(file_err,*) ' If Magnitude factor = 0 then magnitude data list is needed.  '
      write(file_err,*) line
      stop

! + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + 
CONTAINS

! + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +

       subroutine COPYtoHST(iHmax,iSmax,iTmax)

! Data copied and to H, S T sections if Sections exist

       nq = nq - 1

       ! Array Sequence iC iT iH iS must be same as nqCHST(nq...
       do iH = 1,iHmax  

       if(isSection) then
       Hn = HtoC(iH) * Cn
       Hn = max(Hn,2.0d0)  ! H is for the molecule; 2 to allow 1 less for Aj
       sectionH(iC,iH) = Hn
       if(iH_table(iC,iH) == 0) cycle
       endif

        do iS = 1,iSmax  
         do iT = 1,iTmax  
         if(isSection) then
         if(iT_table(iC,iT) == 0) cycle
         if(iS_table(iC,iS,iT) == 0) cycle
         endif
         nq = nq + 1

         ! Section ------------------------
         if(isSection) then 

         ! check sequence
         ierr_sequence = 0
         if(nqCHST(nq,1) /= iC) ierr_sequence = 1        
         if(nqCHST(nq,2) /= iH) ierr_sequence = 2         
         if(nqCHST(nq,3) /= iS) ierr_sequence = 3         
         if(nqCHST(nq,4) /= iT) ierr_sequence = 4         
          if(ierr_sequence > 0) then
          write(file_err,*) ' STOP: nq sequence error ',ierr_sequence 
          write(file_err,*) ' nq,iC,iH,iS,iT ', nq,iC,iH,iS,iT 
          stop
          endif 

         isect = isect + 1
         spCHONSI(:,nq) = 0.0        
         spCHONSI(1,nq) = Cn
         spCHONSI(2,nq) = Hn          

         call SECTIONFORMULA(Cn,Hn,spFormula(nq),length,30) ! returns formula
         call MAKESECTIONNAME(nq,iT,iC,iH,iS,spName(nq))
         endif
         !---------------------------------
              
         ! D/DN2 = (diffusivity of species in mean-gas) / (diffusivity of N2 in mean-gas)
         ! Read D/DN2  Table, ContourPlot
         call SpeciesDataLine(Line,length,linenumber,nHsections,spRelDiff(nq),sectiondensity(iC),printphis(nq))    

         ! URF for species same for all regions
         ! urfs(nq,maxregns) is only really used for non-species (nq < nqsp1) in SCALARCOEFF
         urfs(nq) = URFspeciesdflt

         ! Pr-Sc number for polymers; use -1 for unknown diffusivities
         if(abs(spRelDiff(nq)) < rsmall) then
         write(file_err,*)' No zero diffusivity for species; use -1 for unknown polymer diffusivities ',spName(nq)
         stop
         endif

         if(isect>0 .and. sectiondensity(iC) < rsmall) then
         write(file_err,*)'STOP: Section particle density too small ',sectiondensity(iC),' for section',isect
         stop
         endif

         ! Print
         call ILimitin(Printphis(nq),-5,5,"Printphi Group  ")
       
         ! Plot
!         call ILimitin(Pltphis(nq),-5,5,"Plot phi Group  ")

         solvenq(nq) = .true.    ! all species in list are solved
         enddo
        enddo
       enddo

       return      
       end subroutine COPYtoHST

! + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + 

      end

!=====================================================================

      subroutine READFORMULA(formula,CHONSI,maxlength)

! returns CHONSI numbers and formula length
!===================================================

      use dimensmod
      use filesmod
      use speciesmod
      use ELEMENTSmod

      character (LEN=maxlength)  :: word, formula

      character (LEN=2)   :: chr2, chrE
      character (LEN=80) ::  temp

      integer :: itemp(64)
      real    :: rtemp(64)
      real(8) :: CHONSI(maxelement)


!-----------------------------------
      CHONSI = 0      

! Read Formula into element numbers

! Length to first blank
      word = '                        '

      do n = 1, maxlength
      if(formula(n:n)  ==  ' ') exit
      word(n:n) = formula(n:n)
      Length = n
      enddo

      if(formula(n:n)  /=  ' ') then
      write(file_err,*) ' Formula too long for word length ', formula 
      stop
      endif      

      formula = '                        '
      formula = word      

! Regular      
      nn = 1
      isnumber = 1
      ientry = 0      
      temp = ' '

      n = 1
      DO WHILE(n <= Length)
      isnumberprev = isnumber
      isnumber = 0

      nchr = 2;  if(n==Length) nchr=1    ! No. of characters to read
      chr2 = '  '
      chr2(1:nchr) = word(n:n+nchr-1)
      call IDENTIFYELEMENT(chr2,m_element,nchr)   ! m_element>0 is element No. ; = 0 is a number; -1 fail to recognise
      if(nchr==2) n = n + 1
       

      if(m_element > 0) then
      write(chrE,'(i2)') m_element    ! element No.
      ientry = ientry + 1

      elseif(m_element == 0) then     ! is a number
      isnumber = 1
      else
      write(file_err,*) ' Error in READFORMULA. Formula is ',formula
      endif

      if(isnumber == 0) then
       if(isnumberprev == 0) then    !insert 1 for default 1 before inserting next element No.
       temp(nn:nn) = '1'
       nn=nn+2 
       temp(nn:nn+1) = chrE
       nn=nn+3       
       else
       nn=nn+2                       ! insert element No.
       temp(nn:nn+1) = chrE
       nn=nn+3       
       endif
       if(n == Length) temp(nn:nn) = '1'   !add 1 for last symbols 

      elseif(isnumber == 1) then   ! No. of atoms to temp (one digit at a time)
      nn=nn+1
      temp(nn:nn) = chr2(1:1)
      endif      
      
      n = n + 1
      ENDDO

! convert to numerals      
      itemp = 0; rtemp = 0.0
      ! temp is a series of (element number followed by No. of atoms) eq C2H4 = 1 2  2 4
      read(temp,*,err=998,end=998) (itemp(n),rtemp(n),n=1,ientry) 


      do n = 1,ientry
      i = itemp(n)
      CHONSI(i) = CHONSI(i) + rtemp(n) 
      enddo

      return


998   write(file_err,*)' Bad formula read '
      write(file_err,*) word
      write(file_err,*) temp
      stop

      end

!================================================
            
      subroutine READMULTISPECIES(Line,Nspeciesread,speciesnames,maxspecies,L_Spname,lastchar)            

! read up to maxspecies species in 512 columns each with max length = L_Spname.
! returns formulae in speciesnames and returns No. of species read = Nspeciesread
!  '!' terminates search
! Lastchar = Column No. of last character read in line
            
!=========================================================================================

      use paramsmod
      character (LEN=LineL) ::  Line
      character (LEN=L_Spname)  :: speciesnames(maxspecies)
      character (LEN=L_Spname)  ::  formula

!--------------------------------------------------------------------------------------

      n = 0
      i = 0
      ion=0
      k = 0
      formula = ' '      
      speciesnames = ' '

      do while(n < LineL)
      n = n+1
      if(Line(n:n)== '!') exit   ! comment
      
      if(Line(n:n)==' ' .or. line(n:n)==char(9) .or. Line(n:n)=='/') then   ! blank or slash terminates formula
       if(ion==0) then    ! continuing blank    
       cycle
       else              ! blank after character
       k = k+1
       speciesnames(k)=formula
       if(k == maxspecies) exit
       ion=0
       i = 0
       formula = ' '
       endif

      else                  ! character     
      ion=1  
      i=i+1 
      if(i > L_Spname) cycle
      formula(i:i) = line(n:n)
      endif
      enddo

      Nspeciesread = k
      Lastchar = n - 1      
      return
      end

!=====================================================================
            
      subroutine READSPECIESandFORMULA(Line,name,formula,length)            

! reads species name and formula if there is an alias
! must have no spaces for formula eg. A1C2H=C8H6
!  '!' terminates search 
            
!=========================================================================================

      use paramsmod
      use filesmod
      use elementsmod
            
      character (LEN=LineL) ::  Line
      character (LEN=2)   ::  chr2
      character (LEN=24)  ::  Name
      character (LEN=30)  ::  Formula

!--------------------------------------------------------------------------------------

! Species name - expected to start in col 1
! no leading number in name
       if(Line(1:1)  == ' ' .or. Line(1:1) == char(9) .or. (Line(1:1) >= '0' .and. Line(1:1) <= '9')) then
       write(file_err,*) ' Error: First column of Species name; no blanks and no numbers allowed '
       write(file_err,*) Line
       stop
       endif   

!---------------------------
       do n =1,54                            ! 54  = 24 + 30 character field beginning in col 1     
       if(Line(n:n)== '!') exit   ! comment

       if(Line(n:n) == ' ' .or. Line(n:n) == char(9) .or. Line(n:n) == '=') exit
       enddo
       n1 = n-1
       name = ' '
       name(1:n1) = Line(1:n1) 
       length = n1  ! length to know where to read rest of data line

       formula = ' '
       
       ! If formula is supplied; name is alias
       np1 = n1+1
       if(Line(np1:np1) == '=') then
       do n = np1+1, 54
       if(Line(n:n) == ' ' .or. Line(n:n) == char(9)) exit
       enddo
       n2 = n-1
       formula(1:(n2-np1)) = Line(np1+1:n2) 
       Lenformula = n2-np1
       length = n2  ! length to know where to read rest of data line

       else   ! formula = name
       formula(1:24) = name(1:24)
       Lenformula = Len_trim(formula)
       endif

!--------------------------------
       ! Check if formula is valid
              
       ! Element symbol or number              
       i = 1
       do while(i <= Lenformula)
       nchr = 2;  if(i==lenformula) nchr=1
       chr2 = '  '
       chr2(1:nchr) = formula(i:i+nchr-1)
       call IDENTIFYELEMENT(chr2,m_element,nchr)
       
        if(m_element >= 0) then
        continue
        else                       
        write(file_err,*) ' Stop. Error, or non-listed element in species formula ', formula
        stop
        endif   

       i = i+nchr
       enddo


       ! No brackets in name because they are turned into spaces for 3rd body eg +(M) -> + M
       ! No + in name
       do i = 1,24
       if(name(i:i) == '(' .or. name(i:i) == ')' .or. name(i:i) == '+') then
       write(file_err,*) ' Error: Remove ( ) or + from species name ',name
       stop
       endif 
       enddo

      return
      end

!=====================================================================

      subroutine SpeciesDataLine(Line, Length,linenumber,nHsections,DiffrelN2,sectiondensity, iprintphi) 

!====================================================================      
      use paramsmod
      use filesmod
      character (LEN=*) ::  Line
      real :: HtoC(nHsections) 
      real(8) :: Cn
            
      if(Line(1:1) == '@') then   ! alternative format. Cn * Hn not used here; read from formula later
      read(line(2:LineL),*,iostat=ios) Cn, HtoC(1:nHsections), DiffrelN2, sectiondensity, iprintphi
      else  
      read(line(length+1:LineL),*,iostat=ios) DiffrelN2, iprintphi
      endif

! not used
! List is for Port Species    rprintphi = 1.07 means only iH = 7;  rprintphi = 1.0 means all iH
!      iprint  =  int(rprintphi) 
!      iHprint =  nint(100.0*(rprintphi - iprint))
!      if(iHprint == 0 .or. abs(iHprint) == iH) iprintphi = iprint 
!      iplt   =  int(rpltphi) 
!      iHplt  =  nint(100.0*(rpltphi - iplt))
!      if(iHplt == 0 .or. abs(iHplt) == iH) ipltphi = iplt 


      if(ios/=0) then
      write(file_err,*) ' Error in Species specification at line ',linenumber
      write(file_err,*) line
      stop
      endif
      
      return
      end
!============================================================      

      subroutine SAME_ELEMENTCHECK
      
!========================================================
     use filesmod
     use dimensmod
     use logicmod
     use sectionalmod
     use speciesmod
     use elementsmod
     
     write(fileReact, '(/1x,79(1h-))')
     write(fileReact,'(//"    Species which have the same No. of C H O N elements ")')
     
     nq2 = nqLast
     if(Sectionalmethod) nq2 = nqSectionmol1 - 1 
     nosource = 0                  ! use as scr file
     do nq = nqSp1, nq2-1
      nn = 1
      do n = nq+1,nq2
        isallequal = 1
        do i = 1,maxelement-1
        if(spCHONSI(i,nq) /= spCHONSI(i,n)) isallequal =0
        enddo
      if(isallequal==1) then        
      if(nosource(nq) == 0) nosource(nq) = nq
      if(nosource(n)  == 0)nosource(n)  = nq
      nn = nn + 1
      endif
      enddo 
     
     if(nn > 1) then
     write(filereact,'(10x,"----------------")')
     do n = nqsp1,nq2
     if(nosource(n) == nq)  write(filereact,'(10x,a24)') spName(n)     
     enddo     
     endif 
     enddo

      return
      end

!========================================================

      subroutine LOCATEBALANCESPECIES(nq,nqhold,ispbalancefound) 
      
! Balance species is made nqBalance = nqlastm  

      use SPECIESMOD
      use dimensmod
      

      if(ispbalancefound == 0 .and. spName(nq) == spBalance) then
       if(nq < nqlastm) then
       spName(nqlastm)    = spBalance
       spName(nq) = ' '                     ! to avoid repeated species
       spFormula(nqlastm) = spFormula(nq)
       spCHONSI(:,nqlastm) = spCHONSI(:,nq)       
       nqhold = nq - 1
       nq = nqlastm
       endif
      ispbalancefound = 1
      endif

      return
      end

!========================================================

      subroutine ZEROELEMENTCHECK(nq)
      
      ! Zero elements check (for species) 

      use dimensmod
      use SPECIESMOD
      use filesmod
      use elementsmod

       rmoles = sum(spCHONSI(1:maxelement,nq))
        if(rmoles == 0.0) then
        write(file_err,*)' STOP: No Elements for nq,species ',nq,spName(nq),spFormula(nq)
        stop
        endif

      return
      end

!========================================================

      subroutine WRONGNAMECHECK(nq,irepeated)

! Check for repeated Species Name

      use dimensmod
      use SPECIESMOD
      use filesmod

       ! Ai and Aj reserved names for sectional.
       if(spName(nq)(1:2) == 'Ai' .or. spName(nq)(1:2) == 'Aj') then
       write(file_err,*) 'Stop: Species names Ai and Aj are reserved for Sectional. '
       write(file_err,*) spName(nq)
       stop 
       endif

       ! Repeated Name             
       do i = nqsp1,nq-1
       if(spName(nq)  ==  spName(i) ) then
       write(file_err,'(" Repeated species name in list    ",(a24))') spName(i)
       irepeated = irepeated + 1 
       endif
       enddo

      return
      end

!========================================================

      subroutine IDENTIFYELEMENT(chr2,m_element,nchr)

! read 1 or 2 characters and determines if an element symbol or a number.
! nchr as input gives how many characters to check  (only 1 if last character in formula)
! nchr as output gives how many characters were used for the symbol and therefore how many to advance next
! m_element = -1 for unrecognised 
!           = 0 if character is a number
!           > 0 is element number in elementname(m_element)

      use dimensmod
      use ELEMENTSmod
      character (LEN=2)   :: chr2

      m_element = -1   ! failed recognition

      if(chr2(1:1) >= 'A' .and. chr2(1:1) <= 'Z') then

       ! check if symbol has 2 characters (2nd must be lower case)
       if(nchr>1 .and. chr2(2:2) >= 'a' .and. chr2(2:2) <= 'z') then
       nchr = 2
       else
       nchr = 1
       endif
              
      do m = 1,maxelement  ! search for element symbol     
       if(chr2(1:nchr) == elementname(m)(1:nchr) ) then
       m_element = m 
       exit
       endif
      enddo

      elseif(chr2(1:1) >= '0' .and. chr2(1:1) <= '9') then
      m_element = 0        ! is a number
      nchr = 1 
      endif


      return
      end

!========================================================

      subroutine READ_ELEMENTS(Linechange)

!================================================================================
      use paramsmod
      use dimensmod 
      use filesmod 
      use elementsmod
      use logicmod
                              
      character (LEN=LineL) ::  Line
      character (LEN=2)  ::  formulae(20) = ' ', chr2
      real :: emw(20) = 0.0 
!---------------------------------------------------------------

! Species
       ! READLINE no good when no data, seeks next line

      read(filescheme,'((a))') Line
      call READMULTISPECIES(Line,n_elements,formulae,20,2,lastchar)            
      read(filescheme,'((a))') Line
      read(Line,fmt=*,err=999,end=999) emw(1:n_elements)

      ierr = 0      
      do i = 1, n_elements
      L = len_trim(formulae(i))
      if(.not. (formulae(i)(1:1) >= 'A' .and. formulae(i)(1:1) <= 'Z') ) ierr = i
      if(L==2 .and. .not. (formulae(i)(2:2) >= 'a' .and. formulae(i)(2:2) <= 'z') ) ierr = i 
      enddo

      if(ierr > 0) then
      write(file_err,'(/" STOP: Element specification must be of form X or Xx; is ",a2)') formulae(ierr)
      stop
      endif

      maxelement = 4 + n_elements 

      allocate (ElementName(maxelement),ELMW(maxelement),xs_element(maxelement), &
                sum0_element(maxelement),sum1_element(maxelement), &
                sumin_element(maxelement),sumout_element(maxelement),stat=istat) 
      call STATCHECK(istat,'elements  ')     

      Elementname(1) = 'C '
      Elementname(2) = 'H '
      Elementname(3) = 'O '
      Elementname(4) = 'N '   

      ELMW(1) = 12.0
      ELMW(2) =  1.0
      ELMW(3) = 16.0
      ELMW(4) = 14.0
      
      
      do i = 5,maxelement
      Elementname(i) = formulae(i-4)   
      ELMW(i) = emw(i-4)
      enddo 

      ! check for repeats
      do i = 1,maxelement
       do ii = 1,maxelement
       if(ii == i) cycle 
       if(Elementname(i) == Elementname(ii)) then
       write(file_err,'(" STOP: Repeated element specification ",a2)') Elementname(i)
       stop
       endif
       enddo      
      enddo      

      write(fileReact,'(/"     Elements    MW")')
      do i = 1, maxelement
      write(fileReact,'(i4,3x,a2,f11.1)') i, Elementname(i), ELMW(i)
      enddo

          
      return
999   call errorline(fileSCHEME,'elemen')

      end
!================================================================================
      function IndexSpecies(species)

! nq from spName(nq)

      use dimensmod
      use speciesmod
      use logicmod
      use sectionalmod
      use filesmod
      
      character (LEN=*) species      

      IndexSpecies = 0
      L = len_trim(species) 
      if(L==0) return
      
      ! put all Ai generics into nqSectionmol1; Aj generics into nqSectionrad1 

      if(L>=13) then    
       if(Sectionalmethod .and. species(1:2) == 'Aj' & 
        .and. (species(3:3)/='_' .or. species(6:6)/='_' .or. species(9:9)/='_' .or. species(12:12)/='_') ) then   
       IndexSpecies = nqSectionrad1  ! nq for first radical
       return
       elseif(Sectionalmethod .and. species(1:2) == 'Ai' &
        .and. (species(3:3)/='_' .or. species(6:6)/='_' .or. species(9:9)/='_' .or. species(12:12)/='_') ) then   
       IndexSpecies = nqSectionmol1   ! nq for first molecule     
       return
       endif
      endif

      if(L==2) then    
       if(Sectionalmethod .and. species(1:2) == 'Aj') then     ! Aj
       IndexSpecies = nqSectionrad1  ! nq for first radical
       return
       elseif(Sectionalmethod .and. species(1:2) == 'Ai') then   !Ai
       IndexSpecies = nqSectionmol1   ! nq for first molecule     
       return
       endif
      endif
      
      if(L > 3) then    
       if(Sectionalmethod .and. (species(1:3) == 'Ajt' .or. species(1:3) == 'AjT') ) then     ! Aj
       IndexSpecies = nqSectionrad1  ! nq for first radical
       return
       elseif(Sectionalmethod .and. (species(1:3) == 'Ait' .or. species(1:3) == 'AiT') ) then   !Ai
       IndexSpecies = nqSectionmol1   ! nq for first molecule     
       return
       endif
      endif

        
       do nq = nqsp1,nqlastmp1                                   ! for all incl Ai_xx_xx_xx_x Aj_xx_xx_xx_x   
        if(trim(species) == trim(spName(nq)) ) then
        IndexSpecies = nq
        exit
        endif
       enddo


      return
      end

!==============================================
      subroutine SpeciesSumCheck

! check if sumspecies  > 1 + Ysumtolerance

      use dimensmod
      use phimod
      use numericmod
      use speciesmod                   
      use filesmod
      
      sumY = sum(Phi(nqsp1:nqlastm))
      sumerror = abs(1.0 - sumY)

      if(sumerror > Ysumtolerance) then      
      write(file_err,'(/"STOP: Species input sum ",1pg12.3," exceeds tolerance ",1pg12.3)')sumY,Ysumtolerance
      write(file_err,*) Phi(nqsp1:nqlastm)
      stop
      endif

      
      return
      end

!================================================================
