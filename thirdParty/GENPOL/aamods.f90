!===================================================================
      MODULE DIMENSMOD
 
!===================================================================

! SPECIES

! The phi index of the first scalar
      integer, parameter :: nqmin  = 5  

! The phi index of the first species and first variable in Jacobian 
      integer, parameter :: nqSp1  = 10  

! The maximum number of additional species which are not solved eg. Balance 
!   Must have MEXTRA > 0   
      integer, parameter :: mextra    = 1        ! one for Balance
      integer, parameter :: mextrap1  = 2        ! additional for [M]


!  Maximum No. of parameters for auxiliary rate data TROE, SRI
      integer, parameter :: maxauxparams  = 8    

!  Maximum No. of species in single line of 3rd body auxiliary data (can have multiple lines too)
      integer, parameter :: max3rdbodyperline  = 16    


! (Phi index)
!  1  2  3  4  5  6  7  8  9   10   11  12  13  14  15   ..............
!  P  u  v  w  te ed h  f  g  fuel  CO2 CO O2  H2O  H2    ............

! These names for convenience. Keep these index numbers fixed.
      integer, parameter :: nqte   = 5  
      integer, parameter :: nqed   = 6  
      integer, parameter :: nqEnth = 7  
      integer, parameter :: nqf    = 8  
      integer, parameter :: nqg    = 9  

!  Maximum No. of columns to scan for reaction
      integer, parameter :: maxreactcols = 300

!  Maximum No. of terms in a reaction including 3rd body M
      integer :: maxtermsreaction    

!  Maximum No. of non-consecutive FORD or RORD lines
      integer :: maxFORD    

      integer, parameter :: L_eqnimage  = 256  ! character length

      integer, parameter :: L_spName = 24      ! species name length

!--------------------------------- 
! Maximum number of fluid composition specifications for Ports.
        integer :: maxstreams=1 

!---------------------------------------------------
      integer :: maxcells 

!  maxlist  >= total no. of port cells including overlapping cells
      integer :: maxlist  

!-----------------------------------
! Number of elements
      integer :: maxelement  = 4 

! Species
      integer :: maxReactLines  ! No. of reaction lines including 3rd body continuations and sectional specific reactions
 
!      CH4 .................  A4        Sectionmol1....SectionRadN      Balance      M 
!      nqSp1             nqLastnonSect  nqSectionmol1      nqLast       nqLastm    nqLastmp1
!       |<            maxspecies                             >|

      integer :: maxspecies      ! Maximum number of species solved (not Balance )
      integer :: nqLast          ! Index of last mass fraction species (not Balance)    
      integer :: nqlastm         ! = nqlast + mextra 
      integer :: nqlastmp1       ! = nqlast + mextrap1 

      integer :: nqLastnonSect   ! = nqSectionmol1 -1 if sections, else = nqLast
      integer :: nqBalance       ! Balance species = nqlastm 
      integer :: nqM             ! Sum of all species for 3 body reactions.(nqM = nqLastmp1)
      integer :: nAux            ! number of reactions which have auxiliary rate data 

!  Phi index of base set.
      integer :: nqCO2, nqCO, nqO2, nqH2O, nqH2, nqN2, nqH  
      integer :: nqDroplet, nqHC_EBU

      end module
!--------------------------------------------------------------------

      MODULE PARAMSMOD

      real, parameter :: pi       =   3.1415927 
      real, parameter :: twopi    =   2. * pi   
      real, parameter :: deg2rad  =   pi / 180. 
      real, parameter :: rad2deg  =   180. / pi 


      real, parameter :: rtiny    =   1.0e-24   
      real, parameter :: rsmall   =   1.0e-10   
      real, parameter :: rbig     =   1.0e10    
      real, parameter :: rhuge    =   1.0e28   
      real, parameter :: r4min    =   2.0e-38   
      real, parameter :: r4max    =   2.0e38   

      integer, parameter :: footprnt =  -100
      integer, parameter :: overlap  =  -200

      integer, parameter :: LineL = 512

      real, parameter :: grav       = 9.81        
      real, parameter :: Rgas       = 8.3144      ! kJ/kmol-K
      real, parameter :: Avog       = 6.022e23    ! molecules per gmole   ( = 6.022e26 per kgmole) 
      real, parameter :: SBconst    = 5.67e-11    ! Stefan-Boltzmann  kW/(m^2 K^4)
      real, parameter :: Boltzconst = 1.38e-23    ! Boltzmann constant  J/(molecule - K)
                      !  k = Ru/Navog [J/kmolK / (molecules/kmol) ] = 8314 / 6.023e26 = 1.38e-23 (J/molecule.K)

      real, parameter :: Tstd = 298.15
      real, parameter :: T0K  = 273.15
      real, parameter :: Patm = 101.3  ! kPa  

      integer, parameter :: Keqperline =  10

      real    :: gravx, gravy, gravz
      integer :: iaxisgrav, IunitE, IunitEgasification

      character(Len=10), parameter :: E_units(4) = (/'kJ/kmol   ','kJ/mol    ','kcal/mol  ','cal/mol   '/)
      real, parameter :: CaltoJoule  = 4.186   ! calories * CaltoJoule -> Joules
      real, parameter :: Econversionfactor(4) = (/ 1.0, 1000.0, 4186.0, 4.186 /)
      
      end module
!----------------------------------------

      MODULE ELEMENTSmod

      character(Len=2),allocatable :: ElementName(:)

!     Element molecular weights for  C H O N +        
      real,allocatable :: ELMW(:)
      real,allocatable :: xs_element(:), sumin_element(:), sumout_element(:) 
      real,allocatable :: sum0_element(:), sum1_element(:)
      end module
!----------------------------------------

      MODULE Limitsmod

      real(8), parameter :: chemratemax = 1.0d100
      real(8), parameter :: DLnmax      =  230.0  !  1e-99 -> 1e99, (could have 460 -> 1e198)
      real,    parameter :: RLnmax      =  230.0

      real, parameter :: Punstable = 1000.0      ! absolute P as criterion for instability

      end module
!------------------------------------------

      MODULE TAGMOD

      integer, parameter :: maxTagsum  = 100        ! max No. of tags to select
      real(8),allocatable :: Tagsum(:)              ! Tagsum(ijkp2) sum of rates of selected-reactions
      logical :: isTagsum = .false.                 ! true if any specification
      integer :: iTagsum(20,maxTagsum)              ! specification of included reactions  1 is Tag No.
      integer :: NTaglist = 0                       ! No. of tags in list
      integer :: iTagContour = 0                    ! Contours for tagged reactions 
      integer :: Ltag = 0                           ! current Tag No. in list  
      logical, allocatable :: includeReactSumNo(:)  ! includeReactSumNo(numbereact(nreac)) = .true. if in list
      end module

!------------------------------------------


      MODULE LOGICMOD

      integer, parameter :: maxspecialout = 20   ! max specialout locations in each X and R directions 
      integer, parameter :: maxsectoutGp  = 100   ! max sectional output groups

      
      allocatable :: solvenq(:),  printphis(:), printspart(:), pltphis(:), iprnsoot(:)

      integer :: iSectOutGp(10,maxsectoutGp), NsectoutGp = 0, maxTforSrange = 1          
      logical :: after_isend = .false. 
                   
      logical :: solvenq, walltempfn, EBUmixfrac, isLamEnthEqn
      logical :: istraj= .false., isdispersion = .false., isreacting = .false., isradn = .false.
      logical :: isSOLID = .false. , isLIQUID = .true.
      logical :: isTransportData = .false., UseTransportData = .false. , isMixtureData = .false.
      integer ::  iSpUnits_out, nsootgp, & 
        printvis, printtemp, printden, printblock, printsup, printsvp, printswp, printsmap,  & 
        printshp, printsrh, printvo2, printcpm, printcpt,printchi,  & 
        printraj, printrajend, printUwall,printQwall, printphis, printspart, printFv,  printspeciesmon

      real :: Fv1Cmax, unit_tosec
      character(len=10) :: unitname_time
      integer :: pltphis, pltvisc, pltT, pltO2, pltlgPt , pltraj , pltprod, pltucoef, pltQwall
      integer :: pltchi, pltden

      integer :: iCentigrade_in, iCentigrade_out, timeunit     
      character (Len=2) :: CK_in, CK_out
      character (Len=16) :: head_Tout
      
      integer :: Largesource, isReactionAnalysis = 0
      integer :: itimeplot   ! time plot output      
      integer :: ispecsolver=1, nitercycle1,nitercycle2
      integer :: method
      logical :: LAMINARmethod     = .false.   ! solve 2-D laminar-flow flame field with chemistry
      logical :: FLAMELETLIBmethod = .false.   ! make flamelet library
      logical :: OPPOSEDIFFmethod  = .false.   ! opposed diffusive-reactive layer in physical space 
      logical :: EBUmethod         = .false.   ! 2-D turbulent flame field
      logical :: FLAMELETmethod    = .false.   ! 2-D turbulent flame field using flamelet library
      logical :: CMCmethod         = .false.   ! 2-D turbulent flame field using CMC
      logical :: Sectionalmethod   = .false.   ! Sectional reactions present
      logical :: LamChemmethod     = .false.   ! full chemistry is used for laminar flame (versus flamelet library)
      logical :: isTURBULENT       = .false.   ! flame is turbulent

      logical :: unsteady = .false.
      logical :: correctelement = .false.

      logical :: use_bicgstab, use_bicgstab_p, use_bicgstab_v    !mk

      real :: Xspecialout(maxspecialout), Rspecialout(maxspecialout)       
      integer :: nXspecialout=0, nRspecialout=0, iChemkin = 0
      logical :: isD63
      real ::  Xfilterwidth, Rfilterwidth

      character(Len=10) :: SpUnits_out(0:3) = (/'g/g       ','mol/mol   ','mol/cm^3  ','gm/cm^3   '/)
      end module
!----------------------------------------

      MODULE PHIMOD

      allocatable :: phi(:), phihold(:)

      end module
!----------------------------------------

      MODULE DENMOD
      real :: dencell, denocell
      end module
!----------------------------------------
      MODULE MOLWTMOD

      real :: cellMW
      end module
!----------------------------------------
      MODULE TMOD

      real :: Tcell
      end module
!----------------------------------------
      MODULE PHIOMOD

!  save phis in time stepping
      allocatable :: phio(:)

      end module
!----------------------------------------

      MODULE SIMPLERMOD

      real(8) :: su, ap
      
      end module
!----------------------------------------
      MODULE NUMERICMOD

      integer :: Isweepdirn, Jsweepdirn      
      real :: URFspeciesdflt, URFsourcedef, URFenthsource, delEnthmax
      allocatable :: resid(:), residmax(:), urfs(:), SorsURF(:)
    
      allocatable :: phiresmax(:)
      real :: Ysumtolerance = 0.0      

      end module
!----------------------------------------

      MODULE iterperstepmod
      integer :: iterperstep
      end module
!----------------------------------------

      MODULE TWODIMMOD

! twodim is for axial direction
! isaxi  is axisymmetric

      logical :: twodim=.false., isaxi=.false. , onedim=.false. , zerodim=.false.

      end module
!----------------------------------------
      MODULE SOOTMOD

!  data for the soot model

!  sootabs - constant such that abs = sootabs *fv*T
!  carbtosoot - fraction of carbon in fuel to form soot
!  sootrad - true for soot absorption coefficient in radiation
!  nqsoot - index of species representing soot

      logical :: sootrad= .false.
      real    :: sootabs=0.0, carbtosoot=0.0
      integer :: nqsoot 
      
      end module
!----------------------------------------
      MODULE SOOTDENMOD

!  Pdensitydefault - default density of particle or molecule

      real :: Pdensitydefault = 1000.0  ! kg/m^3
      real, allocatable :: sectiondensity(:)

      end module
!----------------------------------------

      module FOOTMOD

      integer :: footreg, isdone
      allocatable  footreg(:,:), isdone(:,:)

      end module footmod

!----------------------------------------

      MODULE PROPSMOD
      use dimensmod
       
! Properties
      real :: ViscCoefBase(0:4) = 0.0, Schmidtbase=0.0, Prandtlbase=0.0
      real :: airmols, stoichmols, fuelmols, mws, mwh2o, mwco2, mwso2, mwco, mwh2, & 
              hofco2, hofh2o, hofso2, hofco

      real :: Cp_particle, Hfgparticle


      data  mws,  mwh2o, mwco2, mwso2, mwco,  mwh2  & 
        /   32.0, 18.0,  44.0,  64.0,  28.0,  2.0 /

!     enthalpy of formation  (kJ/kmol)
      data hofco2 / -393522. /
      data hofh2o / -241827. /
      data hofco  / -110529. /

      end module
!----------------------------------------

      MODULE PRESSMOD      
! kPa
      real :: pressure ! = 0.0

      end module
      
!---------------------------------------------------      
      MODULE FUELMOD

!  mass fractions of elements in FUEL species
      real :: specenHC

!  mass fractions of elements in fuelstream mixture
      real :: Cinfuelstream, Hinfuelstream, Oinfuelstream
      real :: Cinairstream,  Hinairstream,  Oinairstream
                
      real :: PsatHC(2,50)
      integer :: maxpsat 

      real :: H_particleout = 0.0, dHf = 0.0
      integer :: ispf

      allocatable :: Y_stream(:,:), reservelist(:) 
      integer :: nstream

!  for FLT
! mass fractions of HC in fuel and O2 in Air
! Enthalpy and T of fuel and air 

      real    :: HfuelFLT=0.0, HairFLT=0.0, TfuelFLT, TairFLT    !=0 because Props before Portprops?
      integer :: ifuelstream=0, iairstream=0

      end module
!----------------------------------------

      MODULE STOICHMOD

      real :: fstoich, Ostoich, CO2stoich, H2Ostoich
      end module

!=======================================================================

      MODULE ENTHSOURCEMOD

      allocatable :: Enthsource(:)
      end module

!=======================================================================

      MODULE MONITORMOD
      use dimensmod
      
!      real(8), allocatable :: sumon(:,:), spmon(:,:)
      real(8), allocatable :: storerate(:) 

      integer :: nqMON(10) 

      real, allocatable :: resmon(:)
                           

      real, allocatable :: phimaxnq(:)

      integer, allocatable :: noverNq(:), nMoveNq(:), nunderNq(:)

      integer,allocatable :: nsignificantreact(:),nmaxreact(:),nminreact(:),nstorereact(:)


      character (LEN=3), allocatable :: star(:)
      real, allocatable :: temp(:)
      real ::  dQTijkmax=0.0
      integer :: iXidQTijkmax=0, icdQTmax=0, jcdQTmax=0,  kcdQTmax=0

      end module
!----------------------------------------
      module MONITORPOSMOD
      real    :: xmon=0.0,ymon=0.0,zmon=0.0
      integer :: iregnmon=1,imon=1,jmon=1,kmon=1,mx2mon
      integer :: iregnmon1,imon1,jmon1,kmon1,mx2mon1
      logical :: at_monitor = .false.
      end module
!----------------------------------------

      MODULE SPECIESMOD

! spCHONSI:  C H O N +   numbers in 1 mole of species   

      character(LEN=24), allocatable :: spName(:)
      character(LEN=24) :: Spbalance  ! name of species which is balance gas for diffusion
      character(LEN=30), allocatable :: spFormula(:)

      real(8), allocatable :: sumsource(:), spCHONSI(:,:) 
      real(8) :: Cin=0.0, Cout=0.0, Hin=0.0, Hout=0.0, Csource=0.0, Hsource=0.0, Osource=0.0

      real, allocatable :: spRelDiff(:),  spEnth(:,:,:), spMW(:)
      integer,allocatable :: nosource(:)
      integer :: nqParticleEnth                               ! nq from which to copy enthalpy for particles

      end module
!----------------------------------------


      MODULE REACTMOD
      use dimensmod

!  reactmolnet(maxReactLines)  = net sum of molar coeffs of all species per reaction (for Kequil)
!------------------------
!  nqorderSp(maxtermsreaction,maxReactLines)  = nq of species in term order of occurrence in each equation incl. 3rd body
!  ordermols(maxtermsreaction,maxReactLines)  = moles of species corresponding to nqorderSp; -ve for reactants; +ve for products 

!     terms               1   2  3  |  4   5   6  7     NOTE  3 terms reserved for reactants   
!     eg  nqorderSp      16  13  12  | 25  12  20  0   
!         ordermols      -1  -1  -1  |  1   1   1  0    

!------------------------

!  nqspec(maxtermsreaction,maxReactLines) = lists nq of species appearing (net moles /= 0) in order of nq; else 0 
!  species are not repeated in list; 3rd body does not appear (because net moles =0) 
!   no position order (in nq order); 0 fills out end of lists
!   1 - 2   = reactants (3 not used because is 3rd body)  have -ve net mols
!   4 - max = products                                    have +ve net mols

!  specnetmass(maxtermsreaction,maxReactLines) in nqspec order for each species appearing  
!              = net sum of mass coefficients (= netmols*MW)    -ve for reactants; +ve for products

!     nq order       1   2  3  |  4   5   6  7     NOTE: 3 terms reserved for reactants   

!     eg  nqspec    13  16  0  | 20  25   0  0   
!       net mols    -1  -1  0  |  1   1   0  0    

!------------------------
!  nq3rd(maxReactLines) = nq of 3rd body (0 = no 3rd body)   eg = 12 
!  iterm3rd(2,maxReactLines) = term No.'s of 3rd body species (0 = no 3rd body)    eg = 3 5


!  numbereact(0:maxReactLines) = reaction number (does not increase for 3rd body continuations) + sectional permutations
!  Linereact(numbereact)  = reaction line number of the reaction number
!  nreactionlines         = number of reactions incl 3rd body continuations + sectional permutations

!  chemrate(3,maxReactLines)   = A,n,E in  k =  A  * T^n * exp E/RT
!
!  ireactype    0 elementary reaction
!              -1 3rd body continuation
!   expon     = exponents for species conc. in global reactions.
!   auxrate   = auxiliary reaction rate parameters for P dependent reactions
!   iauxindex(numbereact(n)) =  numbereact points at index of Auxrate
!   iauxtype   = type of aux data 1 = LOW-P;  2 = HIGH-P
!   phi_T  contains particle size d1, d2 dependent coagulation data
!   iTAG(numbereact) = tag No. of reaction number for output
!   reactantC(numbereact(n)) is number of C atoms on reactant side not including 3rd body
 
      real(8), allocatable :: chemrate(:,:), reactratesum(:), fwdratesum(:), auxrate(:,:)
      real(8), allocatable :: reactmon(:), fwdreactmon(:)
      real(8), allocatable :: reactmol(:), ordermols(:,:), reactmolnet(:)
      real(8), allocatable :: specnetmass(:,:)
!      real(8), allocatable :: specmol(:,:)
         
      character(LEN = L_eqnimage), allocatable :: Eqnimage(:)
      character(LEN = 3),  allocatable :: markreact(:)

      integer,allocatable :: numbereact(:),Linereact(:), nqspec(:,:), nq3rd(:), iterm3rd(:,:), nforward(:), &
                             ireactype(:), nqorderSp(:,:), iauxindex(:), iauxtype(:), iTag(:)
      logical,allocatable :: reversible(:)

      real,allocatable ::  expon(:,:), rLnKequil(:), phi_T(:), reactantC(:)

      real, allocatable :: conc(:), PhitoConc(:), concout(:)
      real(8), allocatable :: Ratenq(:)
      
      real :: Tfreeze=0, xsHcriterion = 1000.0, omitcriterion = 1.0,  significantcrit=0.0, T_gamma_out=1500.0, eltolerance = 0.0
      integer :: nreactionlines=0, iallelementary, nreactarray=0, netmoleslist = 0


!  iFORD(ireactionline)       =  contains index of specexpon for ireactionline
!  specexpon(iterm,index)     = special exponents in elementary reactios; species in term order 
!  itermexpon(iterm,index)    = 1 for any terms for which specexpon(iterm,index) /= 1.0 else = 0 
 
      allocatable :: iFORD(:), specexpon(:,:), itermexpon(:,:)
      end module

!----------------------------------------

!----------------------------------------
      MODULE FILESMOD
      
      character (len=24) :: phifilename, specialname, filenamein
      character (len=256) :: datadirectory, keepfilename,cmcoifilename, keepRfilename, flameletoifilename

      logical :: fmtplo = .false. , fmtkeepin = .false. , fmtkeepout = .false.
      logical :: cmcoiexists = .false.
      
      real, parameter ::  version  =  6.3

        integer, parameter ::  fileIN        =  10 
        integer, parameter ::  filedataset   =  11 
        integer, parameter ::  filemontxt    =  12 
        integer, parameter ::  fileRPT       =  17 
        integer, parameter ::  filephi       =  18 
        integer, parameter ::  file_err      =  26 
        integer, parameter ::  filereact     =  35 
        integer, parameter ::  filescheme    =  37 
        integer, parameter ::  filethermo    =  38 
        integer, parameter ::  filenumeric   =  45 
        integer, parameter ::  filephiout    =  48 

        real :: versionkeep
        
      end module
!=======================================      
! memory requirements
      module SIZEMOD
     
      integer :: mbCFD=0, mbFLT=0, mbDFDQ=0, mbCMC=0, mbREACT=0, mbRAD=0, mbTRAJ=0, mbBdfdq=0
      integer :: mbUNSTEADY=0, mbSOURCE=0, mbVODE=0, mBmonitor=0, mbSpecies = 0      
      integer :: mbBICGSTAB = 0                 
      end module
!=======================================      
      module TRESETMOD
      
      real :: Tsetminimum, Treset, Tflameletmin, Tletmaxlimit, Taboveadiab
      integer :: ireset=0, iTlimit=0, iTlimitcmc=0
      integer :: iexceedTletmax=0            
      real ::  Tpartsetmin = 0.0,  Tpartsetmax = 1.0e6

      end module
!=======================================      
      module PHIMINRESMOD
      
      real :: Phiminres
            
      end module
!=======================================      
      module PHIMINMAXMOD
!  phimagnitude   = characteristic value for upward movement limitation
!  from default or phi maximums from previous run soln (is flamelet if soln started with flamelet library) 

!  phimagnitudemin = minimum value of phimagnitude       
!  phimagdefault   = default value initially except for species in fluid streams
!  iupdatephimag   = iteration interval for updating phi magnitude from previous iteration
!  fracmovedown     = max fraction of current value for downward movement
!  fracmoveup       = max fraction of phimax for upward movement
 
      real :: fracmovedown, fracmoveup, criterionmove = 1.0e-3
      real :: phimagnitudemin, phimagdefault    
      integer :: iupdatephimag
      allocatable :: phimagnitude(:)
      
      end module
!=======================================      
      module ISENDMOD
      
      logical :: isend    = .false.    ! last iteration for CFD
      logical :: isendflt = .false.    ! last iteration for Flamelet

      end module
!=======================================      
      module NEWLOOPMOD
      
      logical :: newijkloop = .false.

      end module

!=======================================================================

      MODULE FOCUSMOD

      real :: FocusCriterion, focusthreshold
      integer :: maxfocusit, ifocusolver=1

      end module
      
!=======================================================================
      MODULE THERMOMOD

      real, allocatable :: thermcoeff(:,:,:)       ! Chemkin coefficients for H, S
      real, allocatable :: T_thermorange(:,:)        ! Tmin, Tsplit, Tmax for coefficients
      real, parameter :: Tfit(5) = (/ 300.0, 600.0, 1000.0, 1700.0, 2400.0 /)
      real :: Tpolysplit = Tfit(3)
      end module
!----------------------------------------

      MODULE SECTIONALMOD

! Sectional method
! nCsections = No. of sections at given H with different C
! nHsections = No. of sections at given C with different H
! nSsections = No. of sections at given C,H with different Structure  
! nTsections = No. of sections at given C,H,S with different Type
! nCHSTsections(4) = above 4
! NsectionsAi = total No. of sections for each Ai and Aj 
! nqCHST(nq,4) = iC iH iS iT  section number for each nq
 
! sectionC(1:nCsections) contains C numbers of nqSectionmol1 -> nqSectionmoln
! sectionH(1:nCsections,1:nHsections) contains H numbers of nqSectionmol1 -> nqSectionmoln
! isdoubleAi = 1 for Ai only = 2 for Ai & Aj series
! Gamma2data = coagulation efficiency as function of section(i) & section(j)

! LineRangeSectn(:) are the line Nos. of the first & last specific reaction of each generic sectional reaction 
!                             index is the generic reaction No.
      use paramsmod
              
      real, parameter :: dHradical = 1.0  ! dHradical = Hsectionmolecule - Hsectionradical

      integer :: NsectionsAi, isdoubleAi, ifirstSectionLine, nGenericSectReact=0
      integer :: nCsections, nHsections,  nSsections, nTsections, nCHSTsections(4) = 0  
      integer :: nqSectionmol1=0,  nqSectionmoln=0   ! nq for first, last molecule
      integer :: nqSectionrad1= -1 , nqSectionradn= -1   ! nq for first, last radical if it exists
      integer :: nqSectionN = 0                ! nq for final section: may be nqSectionmoln or nqSectionradn 
      real(8), allocatable :: sectionC(:), sectionH(:,:), Gamma2data(:,:,:), Cgammadata(:)
      integer, allocatable ::  LineRangeSectn(:,:)
      integer, allocatable :: iT_Ranges(:,:), iT_table(:,:), nqCHST(:,:)
      integer, allocatable :: iH_Ranges(:,:), iH_table(:,:), iS_Ranges(:,:), iS_table(:,:,:)
 
      real :: HaverageN ,SaverageN ,TaverageN, sect_Hmol_min = 0.0 
      integer :: iGammadimension = 0  ! = 2 for 2-D data table  = 1 for function else = 0
      integer :: NgammaTables, maxGammaTables = 5 
      real :: Hamaker = 0.0, Hamakerexp = 1.0, HamakerPAH = 0.0
      real, allocatable :: Hamaker_jHfac(:),Hamaker_iSfac(:),Hamaker_iTfac(:)
      character(LEN = LineL), allocatable :: GenericSectionimage(:)    ! image of generic section reaction

      end module
!=======================================================================

      module NITERMOD
      integer :: niter=0
      end module

!=======================================================================

      module phitimemod
      allocatable :: phitimeout(:,:),fvtime(:),sootNtime(:),sootDtime(:), sootHtoCtime(:)
      integer :: Lastimestepout = 0, ntimestep = 0

      end module
!=======================================================================
      module bin1mod
      integer :: maxit=0, iterkeep, ktime
      real :: unsteadytol,deltime,deltimenom 
      end module
!=======================================================================
      module bin2mod
      integer :: kount=0,iterend,newstep=0,ktimemax=0,nqresmaxnorm=0,nzit
      real ::  time=0.0,finaltime,resave=0.0,resavenorm=0.0,resmaxnorm=0.0

      end module
!=======================================================================
      MODULE JACMOD

      real(8), allocatable :: dfdQ(:,:)

      end module

!=======================================================================

      MODULE FNQMOD

      real(8), allocatable :: fnq(:)

      end module

!=======================================================================
                                                          