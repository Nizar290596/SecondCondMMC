/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License

    This file is not part of OpenFOAM but an original routine developed 
    in the OpenFOAM techonology.

    You can redistribute it and/or modify it under the terms of the GNU General
    Public License as published by the Free Software Foundation; either version
    3 of the License, or (at your option) any later version. See GNU General
    Public License at <http://www.gnu.org/licenses/gpl.html>

    OpenFOAM is a trademark of OpenCFD.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.


Class
    Foam::syntheticTurbulence

Authors:
    A. Montorfano, F. Piscaglia
    Dip. di Energia, Politecnico di Milano

Contacts:
    andrea.montorfano@polimi.it , ph. +39 02 2399 3909
    federico.piscaglia@polimi.it, ph. +39 02 2399 8620    

\*---------------------------------------------------------------------------*/

#include "syntheticTurbulence.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "fvPatchFieldMapper.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
namespace Foam
{
    defineTypeNameAndDebug(syntheticTurbulence, 0);
}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


//- construct by components
Foam::syntheticTurbulence::syntheticTurbulence
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    coords_(mesh.C()),
    kName_(dict.lookup("kName")),
    omegaName_(dict.lookup("omegaName")),
    UrefName_(dict.lookup("UrefName")),
    k_(mesh.nCells(),0.0),
    omega_(mesh.nCells(),0.0),
    UMean_(mesh.nCells(),pTraits<vector>::zero),
    distType_(dict.lookupOrDefault<word>("distType","uniform")),
    mesh_(mesh),    
    runTime_(mesh.time()),
    betaStar_(dict.lookupOrDefault<scalar>("betaStar",0.09)),
    A_(dict.lookupOrDefault<scalar>("A",1.4526)),
    ft_(dict.lookupOrDefault<scalar>("ft",1.0)),
    fampl_(dict.lookupOrDefault<scalar>("fampl",1.0)),
    fk_(dict.lookupOrDefault<scalar>("fk",1.0)),
    fSpect_(dict.lookupOrDefault<scalar>("fSpect",1.0)),
    nModes_(readLabel(dict.lookup("nModes"))),
    minK_(readScalar(dict.lookup("minK"))),
    tolerance_(dict.lookupOrDefault<scalar>("tolerance",1e-6)),
    kwMin_(readScalar(dict.lookup("kwMin"))),
    kwMax_(0.0),
    UControl_(dict.lookupOrDefault<Switch>("UControl", false)),
    lamControl_(dict.lookupOrDefault<Switch>("lamControl", false)),
    Nx_(dict.lookupOrDefault<label>("Nx", 256)),
        tauAF_(dict.lookupOrDefault<label>("tauAF", 100)),
    sL_(dict.lookupOrDefault<scalar>("sL", 0.2)),
    Rt_(dict.lookupOrDefault<scalar>("Rt", 1.0)),
    Ttarg_(dict.lookupOrDefault<scalar>("Ttarg", 500.0)),
        Xtarg_(dict.lookupOrDefault<scalar>("Xtarg", 0.01)),
    dtRatio_(readLabel(dict.lookup("dtRatio"))),
    curTimeIndex_(-1),
    randNum_(label(0)),
    fluctuation_(mesh.nCells(),pTraits<vector>::zero),
    LPtr_(NULL),
    kw_(nModes_,0.0),
    kwZero_(coords_.size(),0.0),
    ke_(coords_.size(),0.0),
    kkol_(coords_.size(),0.0),
    y_(coords_.size(),0.0),
    taut_(coords_.size(),0.0),
    fbl_(coords_.size(),1.0),
        os(runTime_.path()/"UMeant.dat"),
    curTime_(0.0)
{
    initVolField(dict);
}


//- construct by components, plus fluctuations
Foam::syntheticTurbulence::syntheticTurbulence
(
    const dictionary& dict,
    const fvMesh& mesh,
    const vectorField& utf
)
:
    coords_(mesh.C()),
    kName_(dict.lookup("kName")),
    omegaName_(dict.lookup("omegaName")),
    UrefName_(dict.lookup("UrefName")),
    k_(mesh.nCells(),0.0),
    omega_(mesh.nCells(),0.0),
    UMean_(coords_.size(),pTraits<vector>::zero),        
    distType_(dict.lookupOrDefault<word>("distType","uniform")),
    mesh_(mesh),    
    runTime_(mesh.time()),
    betaStar_(dict.lookupOrDefault<scalar>("betaStar",0.09)),
    A_(dict.lookupOrDefault<scalar>("A",1.4526)),
    ft_(dict.lookupOrDefault<scalar>("ft",1.0)),
    fampl_(dict.lookupOrDefault<scalar>("fampl",1.0)),
    fk_(dict.lookupOrDefault<scalar>("fk",1.0)),
    fSpect_(dict.lookupOrDefault<scalar>("fSpect",1.0)),
    nModes_(readLabel(dict.lookup("nModes"))),
    minK_(readScalar(dict.lookup("minK"))),
    tolerance_(dict.lookupOrDefault<scalar>("tolerance",1e-6)),
    kwMin_(readScalar(dict.lookup("kwMin"))),
    kwMax_(0.0),
    UControl_(dict.lookupOrDefault<Switch>("UControl", false)),
    lamControl_(dict.lookupOrDefault<Switch>("lamControl", false)),
    Nx_(dict.lookupOrDefault<label>("Nx", 256)),
        tauAF_(dict.lookupOrDefault<label>("tauAF", 100)),
    sL_(dict.lookupOrDefault<scalar>("sL", 0.2)),
    Rt_(dict.lookupOrDefault<scalar>("Rt", 1.0)),
    Ttarg_(dict.lookupOrDefault<scalar>("Ttarg", 500.0)),
        Xtarg_(dict.lookupOrDefault<scalar>("Xtarg", 0.01)),
    dtRatio_(readLabel(dict.lookup("dtRatio"))),
    curTimeIndex_(-1),
    randNum_(label(0)),
    fluctuation_(utf),
    LPtr_(NULL),
    kw_(nModes_),
    kwZero_(coords_.size(),0.0),
    ke_(coords_.size(),0.0),
    kkol_(coords_.size(),0.0),
    y_(coords_.size(),0.0),
    taut_(coords_.size(),0.0),
    fbl_(coords_.size(),1.0),
        os(runTime_.path()/"UMeant.dat"),
        curTime_(0.0)
{
    initVolField(dict);
}


//- construct by components
Foam::syntheticTurbulence::syntheticTurbulence
(
    const dictionary& dict,
    const fvPatchVectorField& pvf,
    const scalarField& fbl
)
:
    coords_(pvf.patch().Cf()),
    kName_(dict.lookup("kName")),
    omegaName_(dict.lookup("omegaName")),
    UrefName_(dict.lookup("UrefName")),
    k_(coords_.size(),0.0),
    omega_(coords_.size(),0.0),
    UMean_(coords_.size(),pTraits<vector>::zero),
    distType_(dict.lookupOrDefault<word>("distType","uniform")),
    mesh_(pvf.patch().boundaryMesh().mesh()),    
    runTime_(pvf.db().time()),
    betaStar_(dict.lookupOrDefault<scalar>("betaStar",0.09)),
    A_(dict.lookupOrDefault<scalar>("A",1.4526)),
    ft_(dict.lookupOrDefault<scalar>("ft",1.0)),
    fampl_(dict.lookupOrDefault<scalar>("fampl",1.0)),
    fk_(dict.lookupOrDefault<scalar>("fk",1.0)),
    fSpect_(dict.lookupOrDefault<scalar>("fSpect",1.0)),
    nModes_(readLabel(dict.lookup("nModes"))),
    minK_(readScalar(dict.lookup("minK"))),
    tolerance_(dict.lookupOrDefault<scalar>("tolerance",1e-6)),
    kwMin_(readScalar(dict.lookup("kwMin"))),
    kwMax_(0.0),
    UControl_(dict.lookupOrDefault<Switch>("UControl", false)),
    lamControl_(dict.lookupOrDefault<Switch>("lamControl", false)),
    Nx_(dict.lookupOrDefault<label>("Nx", 256)),
        tauAF_(dict.lookupOrDefault<label>("tauAF", 100)),
    sL_(dict.lookupOrDefault<scalar>("sL", 0.2)),
    Rt_(dict.lookupOrDefault<scalar>("Rt", 1.0)),
    Ttarg_(dict.lookupOrDefault<scalar>("Ttarg", 500.0)),
        Xtarg_(dict.lookupOrDefault<scalar>("Xtarg", 0.01)),
    dtRatio_(readLabel(dict.lookup("dtRatio"))),
    curTimeIndex_(-1),
    randNum_(label(0)),
    fluctuation_(coords_.size(),pTraits<vector>::zero),
    LPtr_(NULL),
    kw_(nModes_),
    kwZero_(coords_.size(),0.0),
    ke_(coords_.size(),0.0),
    kkol_(coords_.size(),0.0),
    y_(coords_.size(),0.0),
    taut_(coords_.size(),0.0),
    fbl_(fbl),
        os(runTime_.path()/"UMeant.dat"),
        curTime_(0.0)    
{
    initFvPatchField(dict, pvf);    
}


//- construct by components
Foam::syntheticTurbulence::syntheticTurbulence
(
    const dictionary& dict,
    const fvPatchVectorField& pvf,
    const vectorField& utf,
    const scalarField& fbl
)
:
    coords_(pvf.patch().Cf()),
    kName_(dict.lookup("kName")),
    omegaName_(dict.lookup("omegaName")),
    UrefName_(dict.lookup("UrefName")),
    k_(coords_.size(),0.0),
    omega_(coords_.size(),0.0),
    UMean_(coords_.size(),pTraits<vector>::zero),
    distType_(dict.lookupOrDefault<word>("distType","uniform")),
    mesh_(pvf.patch().boundaryMesh().mesh()),    
    runTime_(pvf.db().time()),
    betaStar_(dict.lookupOrDefault<scalar>("betaStar",0.09)),
    A_(dict.lookupOrDefault<scalar>("A",1.4526)),
    ft_(dict.lookupOrDefault<scalar>("ft",1.0)),
    fampl_(dict.lookupOrDefault<scalar>("fampl",1.0)),
    fk_(dict.lookupOrDefault<scalar>("fk",1.0)),
    fSpect_(dict.lookupOrDefault<scalar>("fSpect",1.0)),
    nModes_(readLabel(dict.lookup("nModes"))),
    minK_(readScalar(dict.lookup("minK"))),
    tolerance_(dict.lookupOrDefault<scalar>("tolerance",1e-6)),
    kwMin_(readScalar(dict.lookup("kwMin"))),
    kwMax_(0.0),
    UControl_(dict.lookupOrDefault<Switch>("UControl", false)),
    lamControl_(dict.lookupOrDefault<Switch>("lamControl", false)),
    Nx_(dict.lookupOrDefault<label>("Nx", 256)),
        tauAF_(dict.lookupOrDefault<label>("tauAF", 100)),
    sL_(dict.lookupOrDefault<scalar>("sL", 0.2)),
    Rt_(dict.lookupOrDefault<scalar>("Rt", 1.0)),
    Ttarg_(dict.lookupOrDefault<scalar>("Ttarg", 500.0)),
        Xtarg_(dict.lookupOrDefault<scalar>("Xtarg", 0.01)),
    dtRatio_(readLabel(dict.lookup("dtRatio"))),
    curTimeIndex_(-1),
    randNum_(label(0)),
    fluctuation_(utf),
    LPtr_(NULL),
    kw_(nModes_),
    kwZero_(coords_.size()),
    ke_(coords_.size()),
    kkol_(coords_.size()),
    y_(coords_.size()),
    taut_(coords_.size()),
    fbl_(fbl),
        os(runTime_.path()/"UMeant.dat"),
        curTime_(0.0)
{
    initFvPatchField(dict, pvf);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


Foam::vectorField Foam::syntheticTurbulence::U()
{
    if (coords_.size()!=0)
    {
        if (UControl_ == true && runTime_.value() != curTime_)
        {
             const volScalarField& T = mesh_.thisDb().lookupObject<volScalarField>("T");
             //const volScalarField& dQ = mesh_.thisDb().lookupObject<volScalarField>("dQ");
            //const DimensionedField<scalar, volMesh>& dQ = mesh_.thisDb().lookupObject<DimensionedField<scalar, volMesh>>("chemistry::Sh");
            const volScalarField& dQ = mesh_.thisDb().lookupObject<volScalarField>("flameSens");
             const scalar deltT = runTime_.deltaTValue();
             scalar dtau = tauAF_*deltT;
             //Info<< "dtau = " << dtau << nl << endl;
             //scalar sT = sL_*Foam::pow(Rt_, 0.25);
             //Info<< "sT = " << sT << nl << endl;
             //Info<< "UMean_ = " << UMean_ << nl << endl;
            
            vectorField nInlet(UMean_/mag(UMean_));
             //Info<< "nInlet = " << nInlet << nl << endl;

            scalar dx = Foam::pow(mesh_.V()[0], 1.0/3.0);
             //Info<< "dx = " << dx << nl << endl;
             scalar Tx[Nx_];
             scalar dQx[Nx_];
            scalar Tox[Nx_];
            scalar count[Nx_];
            
             for (int j = 0; j < Nx_; j++)
             {
                Tx[j] = 0;
                dQx[j] = 0;
                Tox[j] = 0;
                count[j] = 0;
             }
             
             forAll(mesh_.cells(), cellI)
             {
                scalar X =  mesh_.cellCentres()[cellI].component(0);
                for (int i = 0; i < Nx_; i++)
                {
                    if (X < (i+1)*dx && X >= i*dx)
                    {
                        Tx[i] = Tx[i] + T[cellI];
                        dQx[i] = dQx[i] + dQ[cellI];
                        Tox[i] = Tox[i] + T.oldTime().oldTime()[cellI];
                        count[i] = count[i] + 1;
                    }
                }
             }

            for (int i = 0; i < Nx_; i++)
            {
                reduce(Tx[i], sumOp<scalar>());
                reduce(dQx[i], sumOp<scalar>());
                reduce(Tox[i], sumOp<scalar>());
                reduce(count[i], sumOp<scalar>());
            }
            
             scalar Tdev = GREAT;
             scalar dQmax = SMALL;
             scalar XdQmax = SMALL;
             scalar XTins = SMALL;
            scalar Tins = SMALL;
            scalar Tinsp1 = SMALL;
               scalar Tinsm1 = SMALL;

             scalar Todev = GREAT;
             scalar XToins = SMALL;
            scalar Toins = SMALL;
            scalar Toinsp1 = SMALL;
            scalar Toinsm1 = SMALL;

             for (int i = 0; i < Nx_; i++)
             {
                if(count[i] > 0)
                {
                       Tx[i] = Tx[i] / count[i];
                    dQx[i] = dQx[i] / count[i];
                       Tox[i] = Tox[i] / count[i];
                }
            }
                
            for (int j = 1; j < Nx_; j++)
            {
                if (mag(Tx[j] - Ttarg_) < Tdev)
                {
                    Tdev = mag(Tx[j] - Ttarg_);
                    XTins = (j+0.5)*dx;
                    Tins = Tx[j];
                    Tinsp1 = Tx[j+1];
                    Tinsm1 = Tx[j-1];
                }
                
                if (mag(Tox[j] - Ttarg_) < Todev)
                {
                    Todev = mag(Tox[j] - Ttarg_);
                    XToins = (j+0.5)*dx;
                    Toins = Tox[j];
                    Toinsp1 = Tox[j+1];
                    Toinsm1 = Tox[j-1];
                }
                
                if (dQx[j] > dQmax)
                {
                    dQmax = dQx[j];
                    XdQmax = (j+0.5)*dx;
                }
            }
            
            if (Tins > Ttarg_)
            {
                XTins = XTins - 2*dx*(Tins-Ttarg_)/mag(Tinsp1-Tinsm1); 
            }
            else
            {
                XTins = XTins + 2*dx*(Ttarg_-Tins)/mag(Tinsp1-Tinsm1);
            }
            
            if (Toins > Ttarg_)
            {
                XToins = XToins - 2*dx*(Toins-Ttarg_)/mag(Toinsp1-Toinsm1); 
            }
            else
            {
                XToins = XToins + 2*dx*(Ttarg_-Toins)/mag(Toinsp1-Toinsm1);
            }
            
              //Info<< "Xtarg = " << Xtarg_ << nl << endl;
            //Info<< "XTins = " << XTins << nl << endl;
            //Info<< "Xflame =" << XdQmax << nl << endl;
             //Info<< "XToins = " << XToins << nl << endl;
             Info << "Current flame position: " << XTins << "\t Target flame position: " << Xtarg_ << "\t Difference: " << (Xtarg_- XTins) << endl;
            
             scalar Uins = mag(UMean_[0]);
             //Info<< "Uins = " << Uins << nl << endl;

            IOdictionary transportDict
            (
                IOobject
                (
                    "transportProperties",
                    runTime_.constant(),
                    runTime_.time(),
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
            );

               //dimensionedScalar nu_(transportDict.lookup("nu"));
            dimensionedScalar nu_("nu", dimensionSet(0,2,-1,0,0,0,0), readScalar(transportDict.lookup("nu")));

            scalarField tkol = Foam::pow(nu_.value()/(betaStar_*k_*omega_), 0.5);
            //Info << "tkol[0] = " << tkol[0] << nl << endl;
            scalar deltX =  Xtarg_-XTins;
            if (lamControl_ == false)
            {
                if (mag(deltX) < tkol[0]*Uins || Uins <= sL_)
                {
                    deltX = 0.0;
                }
            }

            scalar Uvar = (XTins - XToins)/deltT;
            //scalar sTMean = max(sL_, Foam::pow(2.0/3.0*k_[0], 0.5));
            //Info<< "sTMean = " << sTMean << nl << endl;
            
             //scalar dUins = (deltX-dtau*(Uins-sTMean))/(0.5*dtau);
            scalar dUins = (deltX-dtau*Uvar)/(0.5*dtau);
             //Info<< "dUins = " << dUins << nl << endl;

            if(dUins < -0.1*max(Uins, sL_))
            {
                dUins = -0.1*max(Uins, sL_);
            }

            if(dUins > 0.1*max(Uins, sL_))
            {
                dUins = 0.1*max(Uins, sL_);        
            }

             //Info<< "dUins = " << dUins << nl << endl;

            vectorField Unxt = Uins*nInlet;  
            if (Uins >= sL_)
            {
                Unxt = (Uins + dUins*deltT/dtau)*nInlet;
            }
            if (Uins < sL_)
            {
                Unxt = (Uins + 0.1*sL_*deltT/dtau)*nInlet;
            }

             //Info << " Unxt = " << Unxt << nl << endl;

            /*scalarList magUMean;

            mesh_.thisDb().store
            (
                new IOList<scalar>
                (
                       IOobject
                    (
                        "magUMean",
                        mesh_,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    magUMean
                )
            );

            const IOList<scalar>& dbScalarList = mesh_.thisDb().lookupObject< IOList<scalar> >("magUMean");

            const_cast< IOList<scalar>& >( dbScalarList ) = mag(Unxt) ;*/

            taut_ = L()/mag(Unxt);

            //Info<< " taut = " << taut_ << nl << endl;    

            update();

            if (lamControl_ == false)
            {
                vectorField U = Unxt + fluctuation_;

                UMean_ = Unxt;
                
                Info << "Mean inlet velocity control" << endl;
                
                Info << "UMean = " << UMean_ << endl;
                
                Info << "Turbulent fluctuations are generated at the inlet" << endl;

                os << runTime_.value()  << token::TAB << Uins << token::TAB << XdQmax << token::TAB << XTins << endl;
            
                curTime_ = runTime_.value();

                return U;
            }
            else
            {
                vectorField U = Unxt;
                Info << "Laminar inlet velocity" << endl;

                UMean_ = Unxt;

                 Info << "UMean = " << UMean_ << endl;

                os << runTime_.value()  << token::TAB << Uins << token::TAB << XdQmax << token::TAB << XTins << endl;
            
                curTime_ = runTime_.value();

                return U;
            }
        }
        else
        {
            if (lamControl_ == false)
            {
                vectorField U = UMean_ + fluctuation_;

                Info << "Turbulent fluctuations are generated at the inlet" <<endl;

                //Info << "UMean = " << UMean_ << endl;

                return U;
            }
            else
            {
                vectorField U = UMean_;

                Info << "Laminar inlet velocity" <<endl;

                //Info << "UMean = " << UMean_ << endl;

                return U;                
            }
        }
    }
    else
    {
        return vectorField(0);
    }
}


void Foam::syntheticTurbulence::write(Ostream& os) const
{

    os.writeKeyword("kName")        << kName_ << token::END_STATEMENT << nl;
    os.writeKeyword("omegaName")  << omegaName_ << token::END_STATEMENT << nl;
    os.writeKeyword("UrefName")   << UrefName_ << token::END_STATEMENT << nl;
    
    os.writeKeyword("distType")   << distType_ << token::END_STATEMENT << nl;
    os.writeKeyword("betaStar")   << betaStar_ << token::END_STATEMENT << nl;
    os.writeKeyword("A")           << A_ << token::END_STATEMENT << nl;
    os.writeKeyword("ft")             << ft_ << token::END_STATEMENT << nl;
    os.writeKeyword("fampl")       << fampl_ << token::END_STATEMENT << nl;        
    os.writeKeyword("fk")           << fk_ << token::END_STATEMENT << nl;
    os.writeKeyword("fSpect")       << fSpect_ << token::END_STATEMENT << nl;    
    os.writeKeyword("nModes")       << nModes_ << token::END_STATEMENT << nl;
    os.writeKeyword("minK")       << minK_ << token::END_STATEMENT << nl;
    os.writeKeyword("tolerance")  << tolerance_ << token::END_STATEMENT << nl;
    os.writeKeyword("kwMin")       << kwMin_ << token::END_STATEMENT << nl;
    os.writeKeyword("kwMax")       << kwMax_ << token::END_STATEMENT << nl;
        os.writeKeyword("UControl")       << UControl_ << token::END_STATEMENT << nl;
        os.writeKeyword("lamControl")       << lamControl_ << token::END_STATEMENT << nl;
        os.writeKeyword("Nx")           << Nx_ << token::END_STATEMENT << nl;
        os.writeKeyword("tauAF")          << tauAF_ << token::END_STATEMENT << nl;
        os.writeKeyword("sL")           << sL_ << token::END_STATEMENT << nl;
        os.writeKeyword("Rt")           << Rt_ << token::END_STATEMENT << nl;
        os.writeKeyword("Ttarg")       << Ttarg_ << token::END_STATEMENT << nl;
        os.writeKeyword("Xtarg")       << Xtarg_ << token::END_STATEMENT << nl;
    os.writeKeyword("dtRatio")          << dtRatio_ << token::END_STATEMENT << nl;
    
    LPtr_->write(os);
    
}


Foam::syntheticTurbulence::~syntheticTurbulence()
{
    if(LPtr_)
    {
        delete LPtr_;
    }
}

// * * * * * * * * * * * * * * * Private Functions  * * * * * * * * * * * * * //

void Foam::syntheticTurbulence::initVolField(const dictionary& dict)
{

    volScalarField k
    (
        IOobject
        (
            kName_,
            runTime_.timeName(),
            runTime_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    );

    k_ = k.internalField();

    volScalarField omega
    (
        IOobject
        (
            omegaName_,
            runTime_.timeName(),
            runTime_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    );

    omega_  = omega.internalField();

    volVectorField URef
    (
        IOobject
        (
            UrefName_,
            runTime_.timeName(),
            runTime_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    );
    
    UMean_ = URef.internalField();

    if (dict.found("kwMax"))
    {
        kwMax_ = readScalar(dict.lookup("kwMax"));
    }
    else
    {
        kwMax_ = constant::mathematical::pi/
                (max(pow(mesh_.V(),1.0/3.0)).value());
    }

    const scalarField epsilon(betaStar_*k_*omega_);

    IOdictionary transportDict
    (
        IOobject
        (
            "transportProperties",
            runTime_.constant(),
            runTime_.time(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    //dimensionedScalar nu(transportDict.lookup("nu"));
    dimensionedScalar nu("nu", dimensionSet(0,2,-1,0,0,0,0), readScalar(transportDict.lookup("nu")));
    kkol_ = pow(epsilon,(1.0/4.0))*::pow(nu.value(),(-3.0/4.0));

    const dictionary& Ldict(dict.subDict("L"));

    LPtr_ = turbulentLengthScale::New(Ldict, runTime_, mesh_, Ldict.get<word>("type"));
    
    //ke_ = (9.0*constant::mathematical::pi/55.0)*(A_/LPtr_->L());
    ke_ = (9.0*constant::mathematical::pi/55.0)*(1.4526/LPtr_->L());
    
    kwZero_ = ke_/fk_;
    
    forAll (kwZero_, i)
    {
           if (kwZero_[i] < kwMin_)
        {
            kwZero_[i] = kwMin_;
        }


    }

    taut_ = ft_ * k_/epsilon + SMALL;

    update();
}


void Foam::syntheticTurbulence::initFvPatchField(const dictionary& dict, const fvPatchVectorField& pvf)
{
    k_ = readPatchField<scalar>(kName_, pvf);
    omega_ = readPatchField<scalar>(omegaName_, pvf);    
    UMean_ = readPatchField<vector>(UrefName_, pvf);

    if (dict.found("kwMax"))
    {
        kwMax_ = readScalar(dict.lookup("kwMax"));
    }
    else
    {
        kwMax_ = constant::mathematical::pi/
                (max(pow(pvf.patch().magSf(),0.5)));
    }

    const scalarField epsilon(betaStar_*k_*omega_);

    IOdictionary transportDict
    (
        IOobject
        (
            "transportProperties",
            runTime_.constant(),
            runTime_.time(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    //dimensionedScalar nu(transportDict.lookup("nu"));
    dimensionedScalar nu("nu", dimensionSet(0,2,-1,0,0,0,0), readScalar(transportDict.lookup("nu")));
    kkol_ = pow(epsilon,(1.0/4.0))*::pow(nu.value(),(-3.0/4.0));

    const dictionary& Ldict(dict.subDict("L"));

    LPtr_ = turbulentLengthScale::New(Ldict, runTime_, pvf, Ldict.get<word>("type"));
    
    //ke_ = (9.0*constant::mathematical::pi/55.0)*(A_/L());
    ke_ = (9.0*constant::mathematical::pi/55.0)*(1.4526/L());
    
    kwZero_ = ke_/fk_;
    
    forAll (kwZero_, i)
    {
           if (kwZero_[i] < kwMin_)
        {
            kwZero_[i] = kwMin_;
        }


    }

    taut_ = ft_ * k_/epsilon + SMALL;

}


const Foam::scalarField& Foam::syntheticTurbulence::L()
{
    return LPtr_->L();
}


void Foam::syntheticTurbulence::filter(vectorField& ut)
{

    dimensionedScalar dt = runTime_.deltaT();
    dimensionedScalar dtFilt = dt*dtRatio_;

    scalarField a = exp(-dtFilt.value() / taut_);

    scalarField b = fampl_ * Foam::sqrt((1.0-Foam::sqr(a)));

    if (curTimeIndex_>0)
    {
        fluctuation_ = a*fluctuation_ + b*blend(ut);
    }
    else // first timestep
    {
        fluctuation_ = blend(ut);
    }

}

//

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
