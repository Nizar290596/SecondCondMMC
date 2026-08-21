
      function Fv_(Ysoot,iC,gasdensity)

! solid volume fraction (v/v) from mass fraction for section iC   
!==========================================================================      
      use sootdenmod
      use sectionalmod
      
      Fv_ = Ysoot * gasdensity / sectiondensity(iC)   ! sectiondensity(0) = Pdensitydefault

      return
      end
!=====================================================================      

      function Y_out(Ymass,iSpUnits_out,cellMW,dencell,nq)

!=====================================================

       if(iSpUnits_out == 0) then             ! mass fraction
       Y_out = Ymass
       elseif(iSpUnits_out == 1) then         ! mole fraction
       Y_out = Xmolfrac_fn(cellMW, nq, Ymass)
       elseif(iSpUnits_out == 2) then         ! mol/cm^3
       Y_out = Ymolespercm3_fn(dencell,nq, Ymass)
       elseif(iSpUnits_out == 3) then         ! gm/cm^3
       Y_out = Ygmspercm3_fn(dencell,Ymass)
       endif

      return
      end

!=====================================================

      function Xmolfrac_fn(wtmolmix, nq, Ymassfrac)

!====================================================================
! Converts mass fraction to mole fraction treating any species as a gas 
! uses mixture molecular weight input

      use SPECIESMOD

      Xmolfrac_fn = Ymassfrac * wtmolmix / spMW(nq)


      return
      end
!===================================================

      function Ymolespercm3_fn(dencell, nq, Ymassfrac)

!====================================================================
! converts mass fraction to  ! mol/cm^3 treating any species as gas

      use SPECIESMOD

      Ymolespercm3_fn = Ymassfrac * dencell / spMW(nq)  * 1.0e-3

      return
      end
!===================================================

      function Ygmspercm3_fn(dencell,Ymassfrac)

!====================================================================
! converts mass fraction to  ! gm/cm^3 

      use SPECIESMOD

      Ygmspercm3_fn = Ymassfrac * dencell * 1.0e-3

      return
      end

!====================================================================

      subroutine PROPERTIES(dencell,t,cellMW)

!====================================================================
      use paramsmod
      use dimensmod
      use filesmod
      use iterperstepmod
      use logicmod
      use fuelmod
      use speciesmod
      use phimod
      use Tresetmod
                     
      real :: Yspecies(nqsp1:nqlastm) 
!---------------------------------------------------------------------

      Yspecies(nqsp1:nqlast) = phi(nqsp1:nqlast)
      Yspecies(nqBalance)    = Ybalance(Yspecies) 
      phi(nqBalance)         = Yspecies(nqBalance) 

! Molecular weight.
      wtmol = WTMOLFN(Yspecies)
      cellMW = wtmol

! DENSITY
       dencell =  denfn(t,wtmol)


      return
      end
!===========================================================

      function VISCLAMFN(T)

!====================================================================
      use propsmod                                                    
      use filesmod
      
! Dynamic viscosity (kg/(m-s) for T(K) from Chemkin data for 300 - 3000K fitted range   

      T2 = T*T
      T3 = T2*T
      T4 = T3*T
      visclamfn = ViscCoefBase(0) + ViscCoefBase(1)*T + ViscCoefBase(2)*T2 + ViscCoefBase(3)*T3 + ViscCoefBase(4)*T4  

      if(visclamfn <= 0.0) then
      write(file_err,'(" STOP: Laminar viscosity <= 0.0 ",1pg12.3," at T(K) ",f7.1)') visclamfn, T
      stop
      endif
      
!      visclamfn = ( -3.15e-6 * T**2  +  3.42e-2 * T  +  9.46) * 1.0e-6    ! N2 
!      visclamfn = -5.81E-19 * T**4 + 5.08E-15 * T**3 - 1.78E-11 * T**2 + 5.01E-08 * T + 4.43E-06   ! N2
!      visclamfn = -9.38E-19 * T**4 + 8.19E-15 * T**3 - 2.84E-11 * T**2 + 7.50E-08 * T + 4.29E-06   ! Ar

      return
      end
!=====================================================================

      subroutine LnEquilibriumConst(T)

!==================================================================================================

! Returns Ln Equilibrium constant Kc = (C3*C4) / (C1*C2) where C's are mole/cm3

! Kp based on partial pressures from thermo: Ln(Kp)  = -d(G/RT) where dG = Gprod - Greact
! Kp = X3*X4 / (X1*X2) * (P/Po)^(n3+n4 -n1-n2)    X's are mole/mole
! X mol/mol = C mol/cm3 *  RT/P cm3/mol 

! Kp = C3*C4 / (C1*C2) * (P/Po)^(n3+n4 -n1-n2)  * {RT/P cm3/mol}^(n3+n4 -n1-n2)   

! Kc = Kp * {Po/(RT) [mole/cm3] } ^ (n3+n4 -n1-n2)   where n's are mole stoich coeffs
! Ln(Kc) = Ln(Kp) + (n3+n4 -n1-n2) * Ln{Po/(RT) [mole/cm3] }

! Ln(Kc) = -d(G/RT) + (n3+n4 -n1-n2) * Ln{Po/(RT) [mole/cm3] }
!  Po/RT = 101.3 /(8.314 T) *1.e-3  mole/cm3
!----------------------------------------------------------------------------------

      use dimensmod
      use filesmod
      use REACTMOD
      use speciesmod
      use paramsmod
      
      real :: GoRTnq(0:nqlastmp1)

      ! reduce calls to GoRT
      do nq = nqsp1, nqlastm
      call Gibbs(nq,T,GoRTnq(nq),H)
      enddo 
      GoRTnq(0) = 0.0                        ! for nq = 0
      GoRTnq(nqlastmp1) = 0.0          ! M
      
      rLnKequil = 0.0
      PoR = Patm/Rgas * 1.e-3        

       do  n = 1, nreactionlines             ! reaction sequence
       if(.not. reversible(n)) cycle           ! not reversible
       
       sGoRT = 0.0

        do i = 1,maxtermsreaction
        nq = nqorderSp(i,n) 
        sGoRT = sGoRT + ordermols(i,n) * GoRTnq(nq)                  
        enddo

       rLnKequil(n) = -sGorT

       if(abs(reactmolnet(n)) > 0.01) rLnKequil(n) = rLnKequil(n) + reactmolnet(n) * Log( PoR/T ) 
       enddo


      return
      end
!==================================================================================================

      function SOOTVOLFRAC(Yspecies,density)

! soot solid volume fraction  (v/v) from mass fraction of all sectionals
!==========================================================================      
      use paramsmod
      use dimensmod
      use logicmod
      use speciesmod
      use sootdenmod
      use sectionalmod
                  
      real ::  Yspecies(nqsp1:nqlastm)      

      SOOTVOLFRAC = 0.0
      
      ! Soot volume fraction: sum all sectionals

      Fvsoot = 0.0
      do  nq = nqSectionmol1,nqSectionradn
      iC = nqCHST(nq,1)
      Fvi = Fv_(Yspecies(nq),iC,density)
      Fvsoot = Fvsoot + Fvi
      enddo
      
      SOOTVOLFRAC = Fvsoot
     
      return
      end
      
      
!==========================================================================      

      subroutine BETACOAGMIN(Cn1,Cn2,iC1,iC2,Ab,Tncoag,Eact)

! minimum of free molecule regime and continuum
! returns  Ab = Avog*[betamin / (sqrtTnom)] 
! to get coagulation rate:  Acoag = Ab  * T^0.5 * gamma where gamma is sticking prob.

!=====================================================

      use filesmod
      use paramsmod
      use SOOTDENMOD

      real :: Boltzcgs = 1.381e-16    ! Boltzmann const (erg/molK)
      real :: Kn, Lambda
      
      Tncoag  = 0.5
      Eact    = 0.0

!--------------------------      
! Coagulation rate
! from S K Friedlander, Smoke, Dust and Haze Wiley, 1977 p179
! (same as Liang-Shih Fan and Chao Zhu, Principles of Gas-Solid Flows, Cambridge, 1998 

! cgs units
! kinetic rate [molecule collisions /(cm^3-s) ] = B * G * n1(molecules/cm^3) * n2(molecules/cm^3)
! where B beta[cm^3/s] is coagulation constant (obtained in cgs units erg,gm,cm,s)
! B = (8 pi kT)^0.5  (1/mi + 1/mj)^0.5  (ri + rj)^2  where k = Rgas/Avog = 1.38e-16 (erg/molecule.K); m[g]; r[cm]
! for monodispersion B = [192 k T r / rho ]^0.5

! G gamma is collision efficiency 


! Rij[mole coag-collisions/(s-cm^3]   =  G B[cm^3/s] Avog Xi(mol/cm3) Avog Xj(mol/cm3) / Avog
!                                     =  Avog G B[cm^3/s] Xi(mol/cm^3) Xj(mol/cm^3)
! where Avog = 6.023*10^23  molecules/ mol

! Particle mass (gm)       
      rM1 = 12.0 * Cn1 / Avog   
      rM2 = 12.0 * Cn2 / Avog

! Particle radius (cm) 
      Rp1 = 100.0 * 0.5 * Dparticle(Cn1,iC1)
      Rp2 = 100.0 * 0.5 * Dparticle(Cn2,iC2)

! beta = bfree * T^0.5
      Bfree   =  sqrt( 8.0*pi*Boltzcgs * (1.0/rM1 + 1.0/rM2) ) * (Rp1 + Rp2)**2


!-----------------------------------------
! Stokes continuum  p179 Friedlander for particles > 0.1 * gas mean free path  Kn > 10
! Cunningham correction factor for transition  0.2 < Kn < 10
!    see  Hidy and Brock, The Dynamics of Aerocolloidal Systems, Pergamon, 1970 pp305 - 311. 
! Use minimum of free molecule beta and Stokes-Cunningham beta to get basic matchup between regimes.

! beta = C * 2kT/3mu * (1/Ri + 1/Rj) * (Ri + Rj) 
! C is Cunningham correction factor = 1 + Kn(1.257 + 0.4 exp( -0.55/(Kn/2)))

! mu is gas dynamic viscosity; use nominal value at T = 1500K; mu = 5.4e-4 dyne-s/cm^2
      rmu = 5.4e-4   ! gas viscosity cgs
      Tnom = 1500.0; sqrtTnom = 38.7

! Mean free path (cm) 
      Lambda = 100.0 * rMeanfreepath(Tnom)

write(file_err,*)'Lambda ', Lambda 

! use mean molecule for Kn
      Kn = Lambda / (0.5*(Rp1+Rp2))
      C = Cunninghamfactor(Kn)
      
      betaSC = C * 2.0 * Boltzcgs * Tnom / (3.0 * rmu) * (1.0/Rp1 + 1.0/Rp2) * (Rp1 + Rp2) 
      Betafree = Bfree * sqrtTnom
      

! Take minimum beta and divide by sqrt(Tnom) in either case. beta/T^0.5 reasonably constant with varying T  
      Bmin = min(betaSC, Betafree) / sqrtTnom

      Ab   = Avog * Bmin
                                  write(file_err,*)' Cn1,Cn2,kn,C,Betafree,betaSC,Bmin,Ab,Rp1,Rp2,Lambda ', Cn1,Cn2,kn,C,Betafree,betaSC,Bmin,Ab,Rp1,Rp2,Lambda 

!  Acoag = Avog*beta * gamma = Ab  * T^0.5 * gamma
!  calculate gamma separately in order to limit to max = 1.
!  R (gm-moles/cm^3-s) = Acoag * X1(moles/cm^3) * X2(moles/cm^3)


        ! Check Output for T = 1500K 
              write(file_err,'(" C1 C2 Kn Bfree, Bsc Bmin",6g15.3)') Cn1,Cn2,Kn,betafree,betaSC,Bmin*sqrtTnom

            
      return
      end

!=====================================================

      function Cunninghamfactor(Kn)
      
! Cunningham correction factor
      real :: Kn, Knmin   ! Knudsen number 
      
      Knmin = max(Kn,0.02)        ! to avoid underflow in exp      
      expfac = exp(-1.1/Knmin)

      Cunninghamfactor = 1.0  +  (1.257 + 0.4 * expfac ) * Kn

      return
      end
!=====================================================

      function Gammacoagfn(Cn,c1,c2,c3,c4)
! not used
! coagulation collision efficiency
!=====================================================
! empirical result from D'Anna obtained from premixed flames
! Dp calculated using rho = 1.8g/cm3
! primary fit is gamma =  y = 4.707E-11x2 + 8.521E-06x  where x =  (D[nm])^6
! against given  gamma = 1.82e-4 *exp[1.32e7 * D(cm)]         

!      X = Cn**2 
!      Gammacoagfn = 1.11E-17 * X**2 + 4.157E-09 * X   

      if(Cn > 1.e9) then        ! to avoid very large Cn
      Gammacoagfn = 1.0  

      else
                   
      Gammacoagfn = c1* Cn + c2* Cn**2 + c3* Cn**3 + c4* Cn**4
      Gammacoagfn = min(Gammacoagfn, 1.0) 
      endif

                                
      return
      end

!=====================================================

      subroutine AllocatableData

! Fixed data for allocatable arrays      
!====================================================================== 

      use dimensmod
      use SPECIESMOD
      

!  Fixed Names
      spName(0)  =     '         ' 
      spName(1)  =     'P        ' 
      spName(2)  =     'U        ' 
      spName(3)  =     'V        ' 
      spName(4)  =     'W        ' 
      spName(nqte) =   'TE       ' 
      spName(nqed) =   'ED       ' 
      spName(nqEnth) = 'ENTH     ' 
      spName(nqf) =    'F        ' 
      spName(nqg) =    'G        ' 
      spName(nqM) =    'M        '     ! represents sum of all species    

      return
      end
!====================================================================== 

      function Reactionenthalpy(nq)

! std enthalpy of reaction kJ/kg (+ve) per species

      use PROPSMOD
      use SPECIESMOD

      reactionenthalpy = spenth(0,nq,1) - (spCHONSI(1,nq)*hofco2 + 0.5*spCHONSI(2,nq)*hofh2o )/spMW(nq)
      return
      end
!====================================================================== 

      function Ybalance(Y) 

      use dimensmod
                  
      real(8) :: sumY
      real :: Y(nqsp1:nqLast)
      
      sumY     = sum(Y(nqsp1:nqLast))
      Ybalance = max( (1.0 - sumY), 0.0) 

      return
      end
!================================================================

      function DENFN(t,wtmol)

!====================================================================

      use paramsmod
      use pressmod
!      use filesmod

      denfn =  pressure * wtmol/ ( Rgas * t )

      return
      end
!====================================================================

      function WTMOLFN(Y)

! MW of mixture; includes all particle species treated as gases
!====================================================================
      use DIMENSMOD
      use SPECIESMOD
      use logicmod
      use filesmod
      
      real Y(nqsp1:nqLast)

!-------------------------------------------------------------------

      wdenom = 0.

      do nq = nqsp1, nqLast  
      wdenom = wdenom + Y(nq)/spMW(nq)    ! kgi/kg /(kgi/kmoli) = kmoli/kg
      enddo


      YBal = Ybalance(Y)

      wdenom = wdenom + YBal/spMW(nqBalance)

      wtmolfn = 1.0 / wdenom

      return
      end
!=================================================================      

      function SPECIESMW_fn(nq)

! MW of species nq treated as gas from its elemental composition.
!====================================================================

      use DIMENSMOD
      use PROPSMOD
      use SPECIESMOD
      use elementsmod

      wmol = 0.0
       do n = 1, maxelement
       wmol = wmol + ELMW(n) * spCHONSI(n,nq)
       enddo

      SPECIESMW_fn = wmol       

      return
      end

!====================================================================

      subroutine YmasstoXmole(Y,X,wtmolmix)

!====================================================================
      use DIMENSMOD
      use SPECIESMOD
      use logicmod
      use filesmod
      
      real :: Y(nqsp1:nqLastm),X(nqsp1:nqLastm)

!-------------------------------------------------------------------

      wtmolmix = WTMOLFN(Y)

      do nq = nqsp1, nqLastm  
      X(nq) = Xmolfrac_fn(wtmolmix, nq, Y(nq))
      enddo

      return
      end
!=================================================================      


       subroutine ELEMENTSOURCEBAL

!===================================================================
!  Element sources due to chemical reaction. Should be zero.
       use filesmod
       use dimensmod
       use speciesmod 
       use propsmod


        Csource = 0.0
        do nq = nqsp1,nqLast
        Csource = Csource + sumsource(nq) * elfrac(1,nq)
        enddo

        Hsource = 0.0
        do nq = nqsp1,nqLast
        Hsource = Hsource + sumsource(nq) * elfrac(2,nq)
        enddo

        Osource = 0.0
        do nq = nqsp1,nqLast 
        Osource = Osource + sumsource(nq) * elfrac(3,nq)
        enddo

      return
      end

!=================================================                  

      function rMeanfreepath(T)

!====================================================================================                  
      use filesmod
      use paramsmod


!Mean free path of gas L = 1/[ sqrt(2) * pi * d^2 * n]
! from Bird Stewart Lightfoot p20    d=3angstroms for O2   (p31)
! n molecules/vol = Navog*P/RT 
! d = molecule diameter; 
! make d = 3.6e-10m to fit mean free path to Table2.1 Friedlander.

!   L = 1 / [sqrt(2) * pi * d^2 * Navog*P/RT] = T/([sqrt(2) * pi * d^2 * Navog*P/R]
!   L(m) = T/[ 1.4142 * pi * 3.6e-10 ^2 * 6.023e26 * 101/8.314 ]  = T / 4.21e9   = T * 2.37e-10

! Hinds,W.C. Aerosol Technolgy, Wiley, 1999  p21  L air = 0.066e-6 m  at 101kPa, 293K.  

      rMeanfreepath =  2.37e-10 * T 
write(file_err,*)' T,rMeanfreepath',T,rMeanfreepath 

      
      return
      end
      
!====================================================================================                  

      function elfrac(ielem,nq)
! mass fraction of element in species
! ielem: 1 = C  2 = H  3 = O etc.... N S I
!===================================================================

      use dimensmod
      use SPECIESMOD
      use PROPSMOD
      use filesmod
      use elementsmod

      elfrac = spCHONSI(ielem,nq) * ELMW(ielem) / spMW(nq)  

      return
      end 
!===================================================================

      real function PolymerDiffusivity(Cn,iC,T)

      use paramsmod
      use filesmod
      real :: Lambda, Kn

! Large polymer diffusivity based on particle using Stokes friction with Cunningham correction factor
! from S K Friedlander, Smoke, Dust and Haze Wiley, 1977 p32

! D = kT/f
! where Boltzmann k = Ru/Navog [J/kmolK / (molecules/kmol) ] = 8314 / 6.023e26 = 1.38e-23 (J/molecule.K)

!Corrected Stokes friction coeff  f = 3.pi.mu.d/C     
! therefore D = CkT/[3.pi.mu.d]

! mu is gas dynamic viscosity
! C is correction factor    
! C is Cunningham correction factor = 1 + Kn(1.257 + 0.4 exp( -0.55/(Kn/2)))
! Knudsen No.  Kn = L/rp
! L is mean free path of gas    


! meld diffusivity with that at C10 (otherwise there is a jump between this and Chemkin data)
      Cadd = 0.0       
      if(Cn < 200.0) Cadd = (200. - Cn)/(200. - 10.) * 90.   ! adds C90 at Cn=10 -> C0 at Cn = 200
      Cnmod = Cn + Cadd


! Mean free path (m) 
      Lambda = rMeanfreepath(T)

! Particle radius (m) 
      Rp = 0.5 * Dparticle(Cnmod,iC)

! Knudsen No. lambda/Rp       
      Kn = Lambda/Rp

! Cunningham correction factor
      Ccf = Cunninghamfactor(Kn)

      Stokes = Boltzconst * T / (6.0 * pi * visclamfn(T) * Rp)
      
      PolymerDiffusivity = Ccf * stokes        ! m^2/s
      return
      end
      
!===================================================

      function Dparticle(Cn,iC)

! particle diameter(m)
! from particle density(section) and Carbon number Cn per molecule -> MW(g/gmole or kg/kmol) 

!====================================================================================                  

      use filesmod
      use paramsmod
      use sootdenmod
      use sectionalmod

! Diameter                       ! units g, gmole, cm

      if(iC==0) then
      Pdensity = Pdensitydefault   ! default for molecules

      elseif(iC > 0 .and. iC <= nCsections) then
      Pdensity = sectiondensity(iC)

      else
      write(file_err,*)'STOP: error in Dparticle iC>nCsections ', iC,nCsections
      stop
      endif

      Dparticle = ( 12.0 * Cn / (Pdensity*1.e-3 * Avog * pi/6.0)  ) ** 0.333   ! (cm)
      Dparticle = Dparticle * 0.01                                             ! (m) 
write(file_err,*)'Cn,iC,Pdensity,Dparticle,Avog ', Cn,iC,Pdensity,Dparticle,Avog
      return
      end
      
!==============================================================

      function PsurfaceA(X,d)

!====================================================================================                  
! Particles surface area/volume gas   (m^2/m^3) = (m^-1)
! X particle concentration (mol/cm^3)
! d particle diameter      (m)

      use filesmod
      use paramsmod

                      ! units for calculation here are cm

               !molecules/mole*mole/cm^3     cm^2/molecule                
      PsurfaceA = Avog * X  * pi * (100.0 * d)**2            ! (cm^2/cm^3)
      PsurfaceA = PsurfaceA  * 100.0                         ! (m^2/m^3)
      
      return
      end
      
!====================================================================================                  

      function Umeanmolecule(T,Wmol)

!====================================================================================                  
! Molecule mean velocity for molecule with molecular weight Wmol

!Bird, Stewart and Lightfoot, Transport Phenomena p20

! u = sqrt[8KT/(pi*m)]  = sqrt[8RT/(pi*M)]  
!      R = 8314    J/kmol-K   M (kg/kmol)  -> u (m/s)
!      R = 8314e4 erg/mol-K   M (g/mol)    -> u (cm/s)

      use filesmod
      use paramsmod

      Umeanmolecule = sqrt( 8.0 * (Rgas * 1.0e3) * T / (pi * Wmol) )   ! (m/s)

!     Umeanmolecule = sqrt( 8.0 * (Rgas * 1.0e7) * T / (pi * Wmol) )   ! (cm/s)
      
      return
      end
      
!================================================

      subroutine Hquadratic

!==================================================================================================

! use Chemkin polynomials to obtain quadratic H = a + b(T-298) + c(T-298)^2
! fit at T = 300, 1000, 2000K or reduced range


      use paramsmod
      use dimensmod
      use filesmod
      use SPECIESMOD
      use logicmod
      use sectionalmod
      use thermomod
                  
      real :: a(0:2)


      write(fileReact, '(/1x,79(1h-))')
      write(fileReact,'(/"    SPECIES STANDARDISED ENTHALPY COEFFS. (kJ/kg) ")')
      write(fileReact,'(/"    H(kJ/kg)  = a0 + a1*(T-298) + a2*(T-298)^2  fitted for range Tmin - Tmax ")')

      write(fileReact,'(/27x,"                 LOW       RANGE                                HIGH       RANGE                     ")')
      write(fileReact,'( 27x,"Index     a0        a1           a2       Tmin  Tmax      a0         a1           a2       Tmin  Tmax")')


      DO mrange = 1,2
           
      do nq = nqsp1, nqLastm

      if(sectionalmethod .and. nq >= nqSectionmol1 .and. nq <= nqSectionN) then
      spEnth(0:2,nq,mrange) = spEnth(0:2,nqParticleEnth,mrange)      


      else
       if(mrange==1) then
       T1 = Tfit(1)
       T2 = Tfit(2)
       T3 = Tfit(3)
       else
       T1 = Tfit(3)
       T2 = Tfit(4)
       T3 = Tfit(5)
       endif
      
       if(T_thermorange(3,nq) < T3) then  ! for data with reduced T range
       T3 = T_thermorange(3,nq)
       T2 = 0.5*(T1 + T3)
       endif

      call Gibbs(nq, T1, GoRT, H1)   ! H(kJ/kmol)
      call Gibbs(nq, T2, GoRT, H2)
      call Gibbs(nq, T3, GoRT, H3)

      ! solve for a, b, c
      tt1 = T1 -Tstd
      tt2 = T2 -Tstd
      tt3 = T3 -Tstd

      d2  = (tt2  - tt1)
      d3  = (tt3  - tt1)
      d22 = (tt2**2 - tt1**2)
      d33 = (tt3**2 - tt1**2)
            
      denom = (d33*d2 - d22*d3) 
       if(denom == 0.0) then
       write(file_err,*)' error in Hquadratic nq ',nq
       stop
       endif

      a(2) = ( (H3 - H1) * d2 - (H2 - H1) * d3 ) / denom
      a(1) = (H2 - H1 - a(2) * d22) / d2
      a(0) = H1 - a(1) * tt1 - a(2) * tt1**2
      
      spEnth(0:2,nq,mrange) = a(0:2) / spMW(nq)     ! convert to kJ/kg and store
       if(mrange==2) then
       if(mod(nq,20)==0) write(fileReact,'(/27x,"Index     a0        a1           a2       Tmin  Tmax      a0         a1           a2       Tmin  Tmax")')
       write(fileReact,'(1x,a24,i6,1p3g12.3,1x,0p2f6.0,1p3g12.3,1x,0p2f6.0)') spName(nq),nq,spEnth(0:2,nq,1),Tfit(1),Tfit(3),spEnth(0:2,nq,2),T1,T3
       endif
      endif
      enddo 

      ENDDO
      
      write(fileReact, '(/1x,79(1h-))')
      
      return
      end
!=================================================                              

      subroutine Gibbs(nq,T,GoRT,H)

!==================================================================================================

! uses Chemkin polynomial data to return Gibbs energy for Kequil. and enthalpy (to set up quadratic H) 
! returns GoRT = H/RT - S/R  & H(kJ/kmol)

! G = H - TS
! lnK =  -dG/RT = -d(H/(RT) - S/R) = -d(GoRT)

      use paramsmod
      use dimensmod
      use filesmod
      use thermomod

      if(nq == 0 .or. nq == nqM) then
      GoRT = 0.0
      H = 0.0 
      return
      endif 

! Can use high T range coeffs for whole range
! cannot use low T range for high T's      

! Coefficients: 1 for T > Tcommon;  2 for T < Tcommon or if Tmax = Tcommon
      i = 1
      if(T <= T_thermorange(2,nq) .or. abs(T_thermorange(3,nq)-T_thermorange(2,nq)) < 1.0) i = 2 

! H/RT = a1 + a2/2 T + a3/3 T^2 + a4/4 T^3 + a5/5 T^4  + a6/T
      HoRT = thermcoeff(1,i,nq) + thermcoeff(6,i,nq)/T 
      do k = 2,5
      HoRT = HoRT + thermcoeff(k,i,nq)/k * T**(k-1)
      enddo

! S/R =  a1 lnT + a2T + a3/2 T^2  + a4/3 T^3 + a5/4 T^4 + a7
      SoR = thermcoeff(1,i,nq) * log(T) + thermcoeff(7,i,nq)
      do k = 2,5
      SoR = SoR + thermcoeff(k,i,nq)/(k-1) * T**(k-1)
      enddo

      GoRT = HoRT - SoR

      H = HoRT * Rgas * T    ! (kJ/kmol)  
      return
      end

!==================================================================================================

      subroutine THERMOIN

!=======================================================================
! Read Chemkin format thermodynamic data
! 7 coefficients for two temperature ranges

! Format is:

! Species                              Tlow    Thigh   Tcommon                  1   
! Coeffs 1 - 5  for  Tcommon -> Thigh                                           2 
! Coeffs 6 - 7  for  Tcommon -> Thigh     Coeffs 1 - 3  for  Tlow -> Tcommon    3 
! Coeffs 4 - 7  for  Tlow    -> Tcommon                                         4 
 

! Can use high T range coeffs for whole range
! cannot use low T range for high T's      

! H/RT = a1 + a2/2 *T + a3/3 *T^2 + a4/4 *T^3 + a5/5 *T^4  + a6/T
! S/R =  a1 *lnT + a2*T + a3/2 *T^2  + a4/3 *T^3 + a5/4 *T^4 + a7

      use paramsmod
      use filesmod
      use dimensmod
      use thermomod
      use speciesmod
      use logicmod
      use sectionalmod
                              
      character (LEN=LineL) ::  Line
      character (LEN=24)  ::  formula, alias

!-------------------------
      T_thermorange = 0.0
      Tdefault = 1000.0
      thermcoeff = rhuge
            
      DO WHILE (.true.)

      ! Read
      Line = ' '
      call READLINE(fileTHERMO,Line,iend,Linechange)
      if(iend == 1) exit

       ! Read default common T
       if(Line(1:6) == 'THERMO') then
       read(fileTHERMO,*,iostat=ios) T1dum,Tdefault,T2dum 
       if(ios/=0) write(file_err,'(/" Thermo data, no default Tcommon found ")') 
       cycle       
       endif

!-------------------------------
! Species name
       ! Word length for species formula
       do n =1,18                            ! 18 character field beginning in col 1     
       if(Line(n:n) == ' ' .or. Line(n:n) == char(9)  .or. Line(n:n) == '=') exit
       enddo
       n1 = n-1
       formula = ' '
       formula(1:n1) = Line(1:n1) 

       ! Search for any alias
       alias = ' '
       do n = n1+1, n1+18    
        if(Line(n:n) == '=') then      ! alias must start with '=' immediately preceding formula 
        n2 = n+1
        do nn =n2,n2+18                 ! search for blank 
        if(Line(nn:nn) == ' ' .or. Line(nn:nn) == char(9)) exit
        enddo
        alias(1:nn-n2+1) = Line(n2:nn) 
        exit        
        endif
       enddo
       
       !write(file_err,*)' formula,alias ',formula,alias

!-------------------------------
! Index if species is in scheme
      nq = IndexSpecies(formula)
      if(nq == 0) nq = IndexSpecies(alias)
      if(nq == 0) cycle

!-------------------------------
! Temperature range order  Tmin Tmax Tsplit
      read(Line(48:65),fmt=*,err=999,end=999) Tmin, Tmax
       ! check T range
       if(Tmin > 350.0 .or. Tmax < 3000.0) then
       write(file_err,'(" WARNING: Want Thermo data T range 350 - 3000, but is: ",2f7.0," for ",(a))') Tmin,Tmax,formula
       endif
       T_thermorange(1,nq) = Tmin
       T_thermorange(3,nq) = Tmax
       
      read(Line(66:76),fmt=*, iostat=io) T_thermorange(2,nq)   
       if(T_thermorange(2,nq) < Tmin) then
       T_thermorange(2,nq) = Tdefault    ! default in case no Tcommon specified
!       write(file_err,'(" Thermodynamic data common T set to Tdefault for ",(a))') formula
       endif


!-------------------------------
! Coefficients
      call READLINE(fileTHERMO,Line,iend,Linechange)
!      read(Line,fmt='(5e15.8)',err=999) thermcoeff(1:5,1,nq)
      read(Line,*,err=999) thermcoeff(1:5,1,nq)
      call READLINE(fileTHERMO,Line,iend,Linechange)
!      read(Line,fmt='(5e15.8)',err=999) thermcoeff(6:7,1,nq), thermcoeff(1:3,2,nq)
      read(Line,*,err=999) thermcoeff(6:7,1,nq), thermcoeff(1:3,2,nq)
      call READLINE(fileTHERMO,Line,iend,Linechange)
!      read(Line,fmt='(5e15.8)',err=999) thermcoeff(4:7,2,nq)
      read(Line,*,err=999) thermcoeff(4:7,2,nq)

!write(file_err,'(5e17.8)')thermcoeff(1:5,1,nq)
!write(file_err,'(5e17.8)')thermcoeff(6:7,1,nq), thermcoeff(1:3,2,nq)
!write(file_err,'(5e17.8)') thermcoeff(4:7,2,nq)
      ENDDO

!-------------------------------
! Check for Missing Species
      imiss = 0
      write(file_err,'(" ")')
      do nq = nqSp1, nqLastm
      if(sectionalmethod .and. nq >= nqSectionmol1 .and. nq <= nqSectionN) cycle     ! sectional polynomials to copy enthalpy 
       if(thermcoeff(1,1,nq)>=rhuge) then
       imiss = 1
       write(file_err,'(" Thermo data missing for species number, name ",i6,3x,(a))') nq, spname(nq)
       endif      
      enddo 

      if(imiss > 0)  stop
      
      write(fileRPT,'(/" Thermo data, default Tcommon is ",f7.0)') Tdefault
      return

999   write(file_err,'(/" Replace any tab in THERMO.IN with a single space. "/)')
      call errorline(fileTHERMO,'Thermoin')

      end
      
!==================================================================================================

      function TEMPFN_G(h, Yspecies,a0)

      use DIMENSMOD
      use SPECIESMOD
      use paramsmod
      use thermomod

      real :: Yspecies(nqsp1:nqlastm), a(0:2)

! Temperature of gas mixture from given H and Yspecies
!=====================================================================
! Enthalpy coefficients  H(kJ/kg) = a0 + a1(T-298) + a2(T-298)^2

! must try low range first to check if H < or > than Hsplit

      a = 0.0
       do nq = nqsp1, nqLastm  
       a(0) = a(0) + Yspecies(nq) * spenth(0,nq,1)
       a(1) = a(1) + Yspecies(nq) * spenth(1,nq,1)
       a(2) = a(2) + Yspecies(nq) * spenth(2,nq,1)
       enddo

      Hsplit = a(0) + a(1)*(Tpolysplit-Tstd) + a(2)*(Tpolysplit-Tstd)**2  

      if(H < Hsplit) then
      Tempfn_g = TEMPFN(h, a)

      else
      a = 0.0
       do nq = nqsp1, nqLastm  
       a(0) = a(0) + Yspecies(nq) * spenth(0,nq,2)
       a(1) = a(1) + Yspecies(nq) * spenth(1,nq,2)
       a(2) = a(2) + Yspecies(nq) * spenth(2,nq,2)
       enddo

      Tempfn_g = TEMPFN(h, a)
      endif
      
      a0 = a(0)      
       
      return
      end
!========================================================

      function TEMPFN(h, a)

      use paramsmod
      
      real :: a(0:2) 

! Temperature calculations     h = a0 + a1 T + a2 T^2
! NEEDS PRIOR  CALL  for a0 a1 a2
!=====================================================================
      
! Linear
      if(a(2) == 0.0) then

      Tempfn =  (h - a(0)) / a(1)  + Tstd

! Quadratic
      else
        del = a(1)**2 - 4.0 * a(2) * ( a(0) - h )

        if(del < 0.0) then
                               !    write(*,*) 'enthalpy quadratic wrong  del= ',del
                               !    write(*,*)' b1 b2 b3 h ',b1,b2,b3,h
        Tempfn = 0.0

        else
        Tempfn = ( -a(1) + sqrt(del) ) / (2.0 * a(2))  + Tstd
        endif
      endif

      return
      end
!========================================================

      function ENTHFN_G(T,Yspecies)

      use DIMENSMOD
      use SPECIESMOD
      use paramsmod
      use thermomod

      real :: Yspecies(nqsp1:nqlastm), a(0:2)

! Enthalpy for gas mixture kJ/kg
!-------------------------------------------------------------------
! Enthalpy coefficients  H(kJ/kg) = a0 + a1(T-298) + a2(T-298)^2

      a = 0.0

      m = 1
      if(T > Tpolysplit) m = 2

       do nq = nqsp1, nqLastm  
       a(0) = a(0) + Yspecies(nq) * spenth(0,nq,m)
       a(1) = a(1) + Yspecies(nq) * spenth(1,nq,m)
       a(2) = a(2) + Yspecies(nq) * spenth(2,nq,m)
       enddo
      
      enthfn_g =  a(0)  +  a(1) * (T - Tstd) + a(2) * (T - Tstd)**2 
      
      return
      end
!=================================================================      

      function ENTHFN_nq(T,nq)

      use DIMENSMOD
      use SPECIESMOD
      use paramsmod
      use thermomod

! H for single species kJ/kg
!-------------------------------------------------------------------
! Enthalpy coefficients  H(kJ/kg) = a0 + a1(T-298) + a2(T-298)^2



      m = 1
      if(T > Tpolysplit) m = 2

      enthfn_nq =  spenth(0,nq,m)  +  spenth(1,nq,m) * (T - Tstd) + spenth(2,nq,m) * (T - Tstd)**2 
      
      return
      end
!=================================================================      



      subroutine SUM_FIELD(sum_element)

!====================================================

      use dimensmod
      use phimod
             
      real :: sum_element(maxelement)

      sum_element = 0.0

      do nq = nqsp1,nqLast 
       do ie = 1,maxelement
       sum_element(ie)  = sum_element(ie)  + phi(nq) * elfrac(ie,nq)
       enddo
      enddo


      return
      end
!=====================================================

    
      function Econversion(Ein,iunit)

      use filesmod
      use paramsmod
      
! Activation Energy conversion to kJ/kmol
!  Units for E   1 = kJ/kmol;  2 = kJ/mol;   3 = kcal/mol  4 = cal/mol

        if(iunit >= 1 .and.  iunit <= 4) then
        Econversion = Ein * Econversionfactor(iunit)
        else
        write(file_err,*)' Activation Energy Units are wrong ', iunit
        stop
        endif 

      return
      end
!==============================================      
    
      function Ebackconversion(EkJperkmol,iunit)

      use filesmod
      use paramsmod
      
! Activation Energy conversion from kJ/kmol to iunit
!  Units for E   1 = kJ/kmol;  2 = kJ/mol;   3 = kcal/mol  4 = cal/mol

        if(iunit >= 1 .and.  iunit <= 4) then
        Ebackconversion = EkJperkmol / Econversionfactor(iunit)
        else
        write(file_err,*)' Activation Energy Units are wrong ', iunit
        stop
        endif 

      return
      end
!==============================================      

      function Hconversion(iunitH,Hin,spMWnq)

      use PARAMSMOD
      use filesmod

! Enthalpy conversion to kJ/kg
! Units for H input:  1 = kJ/kg;  2 = kJ/kmol;  3 = Kcal/mol; 4 = cal/mol !   

       if(iunith == 1) then
       Hconversion = Hin
       elseif(iunith == 2) then
       Hconversion = Hin / spMWnq
       elseif(iunith == 3) then
       Hconversion = (CaltoJoule * 1000.0) * Hin / spMWnq
       elseif(iunith == 4) then
       Hconversion = CaltoJoule  * Hin / spMWnq
       else
       write(file_err,*)' Enthalpy Units wrong ', iunitH
       stop
       endif


      return
      end
!==============================================      

      subroutine CycleTIMER(niter)

      use filesmod
      
      integer :: it1 = 0.0      

      if(niter == 2) call system_clock(it1,icrate,icmax)
      
      if(niter == 3) then
      call system_clock(it2,icrate,icmax)
      icmax=icmax
      
      secs = float(it2-it1)/ icrate
!      elapseh = secs * float(maxit) /3600.
      elapseh = secs/3600.
      ielapseh = elapseh
      elapsem = (elapseh - float(ielapseh)) * 60.

      ielapsem = elapsem
      elapses = (elapsem - float(ielapsem)) * 60.
      ielapses = elapses
      
!      write(*      ,'("  Time for",i7," iterations is",i4," hr  ",i4," min  ",i4," sec"/)') maxit,ielapseh,ielapsem,ielapses
!      write(fileRPT,'("  Time for",i7," iterations is",i4," hr  ",i4," min  ",i4," sec"/)') maxit,ielapseh,ielapsem,ielapses
!      write(*      ,'("  Time for 1 iteration is",i4," hr  ",i4," min  ",f7.2," sec"/)') ielapseh,ielapsem,elapses
!      write(fileRPT,'("  Time for 1 iteration is",i4," hr  ",i4," min  ",f7.2," sec"/)') ielapseh,ielapsem,elapses
      endif

      return
      end
!======================================================  
      subroutine elapsedTIME(iset,secs)

      integer :: it1 = 0.0      
       
      if(iset== 0) then
      secs = 0.0 
      call system_clock(it1,icrate,icmax)

      else      
      call system_clock(it2,icrate,icmax)
      if(it2 < it1) it1 = it1 - icmax  ! overflow
      secs = float(it2-it1)/ float(icrate)
      endif

      return
      end

!=========================================================================================
