/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

\*---------------------------------------------------------------------------*/

// #include <cmath>
#include "TwoEquationModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::TwoEquationModel<CloudType>::TwoEquationModel
(
    const dictionary& dict,
    CloudType& owner
)
:
    SootModel<CloudType>(dict, owner, typeName),
    species_(this->coeffDict().template lookupOrDefault<wordList>("species", {"C2H2", "O2", "OH"})),
    dySootdt_(0.),
    dnSootdt_(0.),
    diaSoot_(0.),
    surfSoot_(0.),
    kineticsCoeffs_(this->coeffDict().template lookupOrDefault
            <Foam::List<Foam::Tuple2<string, Foam::scalarList> > > ("kinetics",
                {
                    {"inception_C2H2",  {1e+4, 0, 0,   21000, 0,   2}},
                    {"cogulation_Cs",   {9}},
                    {"surfGrowth_C2H2", {6000, 0, 0,   12100, 0.5, 2}},
                    {"oxidation_O2",    {1e+4, 0, 0.5, 19778, 1,   1}},
                    {"oxidation_OH",    {0.36, 0, 0.5, 0,     1,   1}}
                })),
    rhoSoot_(this->coeffDict().template lookupOrDefault<scalar>("rhoSoot", 1800)),
    nSootScaling_(this->coeffDict().template lookupOrDefault<scalar>("nSootScaling", 10e+15)),
    inceptionCarbon_(this->coeffDict().template lookupOrDefault<scalar>("inceptionCarbon", 60)),
    gasSource_(this->coeffDict().template lookupOrDefault<Switch>("gasSource", false)),
    gasSourceTable_(this->coeffDict().template lookupOrDefault
            <Foam::List<Foam::Tuple2<string, Foam::scalar> > > ("gasSourceTable",
                {
                    {"inception_C2H2",  -1.},
                    {"inception_H2",     1.},
                    {"surfGrowth_C2H2", -1.},
                    {"surfGrowth_H2",    1.},
                    {"oxidation1_O2",   -1.},
                    {"oxidation1_CO",    2.},
                    {"oxidation2_OH",   -1.},
                    {"oxidation2_CO",    1.},
                    {"oxidation2_H",     1.}
                })),
    YsourceTable_(gasSourceTable_),
    sootLimit_(this->coeffDict().template lookupOrDefault<scalar>("sootLimit", 0.03e-6)),
    turbMixing_(this->coeffDict().template lookupOrDefault<Switch>("turbMixing", false)),
    debug_(this->coeffDict().template lookupOrDefault<Switch>("debug", false)),
    sootContribs_(6, 0.)
{
    // string inceptionName = kineticsCoeffs_[0].first();
    // scalarList inceptionList = kineticsCoeffs_[0].second();
    // Info << "inceptionName is " << inceptionName << endl;
    // Info << "inceptionList is " << inceptionList << endl;

    Info << "-> transport equations for ySoot and nSoot are solved." << endl;

    Info << "-> turbulent mixing:            " << this->turbMixing_ <<endl;
    Info << "-> coupling back to gas-phase:  " << this->gasSource()<<endl;
    Info << "-> soot density:                " << this->rhoSoot() << " kg/m^3" << endl;
    Info << "-> carbon number for inception: " << this->inceptionCarbon() << endl;
    Info << "-> sootIm criterion:            " << sootLimit() * 1e6 << " ppm" << endl;


    // initial YsourceTable.second() (ie, a scalar variable) set the value to 0
    for (int i = 0; i < YsourceTable_.size(); i++)
    {
        YsourceTable_[i].second() = 0.;
    }

    this->speciesCheck();

    // print the species info and relevant rate constants
    this->speciesInfo();

}

template<class CloudType>
Foam::TwoEquationModel<CloudType>::TwoEquationModel
(
    const TwoEquationModel<CloudType>& cm
)
:
    SootModel<CloudType>(cm),
    species_(this->coeffDict().template lookupOrDefault<wordList>("species", {"C2H2", "O2", "OH"}))
{
    Info << "TwoEquationModel 2nd constructor is called here." << endl;
}


// * * * * * * * * * * * * * * *Member Functions * * * * * * * * * * * * * * //

template<class CloudType>
void Foam::TwoEquationModel<CloudType>::speciesCheck()
{
    // check the species in fuel mechanism
    for (int i=0; i < species().size(); i++)
    {
        if (this->owner().composition().slgThermo().carrier().contains(this->species()[i]))
        {
            continue;
        }
        else
        {
            FatalErrorInFunction
                << this->species()[i]
                << " is not included in current fuel mechanism."
                << abort(FatalError);
        }
    }

    // check the size of kineticsCoeffs, abort if it does not equal to 6
    for (int i = 0; i < kineticsCoeffs().size(); i++)
    {
        if (kineticsCoeffs()[i].first().find("cogulation") != string::npos)
        {
            continue;
        }
        else
        {
            if (kineticsCoeffs()[i].second().size() != 6)
            {
                FatalErrorInFunction
                    << "kineticsCoeffs list should contain following data: \n"
                    << "(0) A, the pre-factor,\n"
                    << "(1) pressureDependency,\n"
                    << "(2) TCoeff, the beta,\n"
                    << "(3) TAct, activation energy,\n"
                    << "(4) SootSurfaceCoeff,\n"
                    << "(5) Number of carbon in the chosen species."
                    << abort(FatalError);
            }
        }
    }

    // pressure dependence is not supported now
    for (int i = 0; i < kineticsCoeffs().size(); i++)
    {
        if (kineticsCoeffs()[i].first().find("cogulation") != string::npos)
        {
            continue;
        }
        else 
        {
            if (kineticsCoeffs()[i].second()[1] != 0)
            {
                FatalErrorInFunction
                    << "Note: the current version does not consider pressure impacts\n"
                    << "Therefore, please have a check of soot kinetics."
                    << abort(FatalError);
            }
            else
            {
                continue;
            }
        }
    }
}

template<class CloudType>
Foam::scalar Foam::TwoEquationModel<CloudType>::calRateCoeff
(
    const scalar& kineticsCoeff,
    const scalar& T,
    const scalar& ySoot,
    const scalar& nSoot,
    const scalar& rho
)
{
    Foam::scalar cgRateCoeff;
    const double PI = constant::mathematical::pi;
    const double Bm = 1.38e-23;

    if (ySoot <= ROOTVSMALL || nSoot <= ROOTVSMALL)
    {
        cgRateCoeff = 0;
    }
    else
    {
        cgRateCoeff = 2 * kineticsCoeff * std::sqrt(6.*Bm*T/rhoSoot()) *
            pow((6.*ySoot)/(PI*rhoSoot()*nSoot*nSootScaling()), 1./6.) *
                    pow(rho*nSoot*nSootScaling(), 2.);
    }

    cgRateCoeff = max(cgRateCoeff, 0);

    return cgRateCoeff;
}

template<class CloudType>
Foam::scalar Foam::TwoEquationModel<CloudType>::calRateCoeff
(
    const scalarList& kineticsCoeffs,
    const scalar& T
)
{

    Foam::scalar rateCoeff;

    if (kineticsCoeffs[1] == 0) // no pressure dependency
    {
        if (kineticsCoeffs[4] != 0) // soot surface dependency
        {
            rateCoeff = kineticsCoeffs[0] * pow(T, kineticsCoeffs[2]) *
                exp(-kineticsCoeffs[3]/T) * pow(this->surfSoot(), kineticsCoeffs[4]);
        }
        else // no soot surface dependency
        {
            rateCoeff = kineticsCoeffs[0] * pow(T, kineticsCoeffs[2]) * exp(-kineticsCoeffs[3]/T);
        }
    }
    else // considering the pressure dependency - not supported now (may work as a switch later)
    {
        rateCoeff = 0.;
    }

    return rateCoeff;
}

template<class CloudType>
void Foam::TwoEquationModel<CloudType>::fillInSource
(
    const string& reactspecies,
    const scalarField& Y,
    const scalar& T,
    const scalar& ySoot,
    const scalar& nSoot,
    const scalar& rho,
    scalar& source4ySoot,
    scalar& source4nSoot
)
{
    // Note to units:
    // Molar weight is given as kg/kmol
    // OpenFOAM provides the constant RR which is also in J/(kmol K) not
    // in J/(mol K), see also thermodynamicConstants.C


    // basic information [kg/kmol]
    const double MWCarbon = 12.01099968;
    // Multiply the with 1E+3 for unit consistency [1/kmol]
    const double NA = 1E+3*constant::physicoChemical::NA.value();
    const scalar& deltaT  = this->owner().time().deltaTValue();
    string gasSourceTag;

    // info from function argument
    const int delimiterLoc  = reactspecies.find("_");
    const string& reactName = reactspecies.substr(0, delimiterLoc);
    const string& species   = reactspecies.substr(delimiterLoc+1, reactspecies.size());
    const int index = indexInTable<scalarList>(reactspecies, this->kineticsCoeffs());

    if (index == -1)
    {
        FatalErrorInFunction
            << "The argument "
            << reactspecies
            << " is not in the list of kineticsCoeffs\n"
            << abort(FatalError);
    }
    else if (reactName == "cogulation") // cogulation
    {
        const scalar& kineticsCoeff = this->kineticsCoeffs()[index].second()[0];
        source4nSoot = calRateCoeff(kineticsCoeff, T, ySoot, nSoot, rho);
    }
    else
    {
        const scalarList& kineticsCoeffs = this->kineticsCoeffs()[index].second();
        const scalar& rateCoeff = this->calRateCoeff(kineticsCoeffs, T);

        if (species == "C2H2") // inception or surface growth -> C2H2 and H2
        {
            const double WC2H2 = 26.037;
            const double WH2   = 2.016;
            const label& IC2H2 = this->owner().composition().slgThermo().carrier().species()["C2H2"];
            const scalar YC2H2 = Y[IC2H2];
            const scalar cC2H2 = YC2H2*rho/WC2H2;

            const scalar source = rateCoeff * cC2H2;

            if (reactName == "inception")
            {
                source4ySoot = kineticsCoeffs[5] * source * MWCarbon;
                source4nSoot = kineticsCoeffs[5] * source * NA / inceptionCarbon();
            }
            else if (reactName == "surfGrowth")
            {
                source4ySoot = kineticsCoeffs[5] * source * MWCarbon;
            }

            if (gasSource() == true)
            {
                // here C2H2 and H2 are considered
                this->fillInGasSourceTable(reactName, "C2H2", WC2H2, source, deltaT, rho);
                this->fillInGasSourceTable(reactName, "H2", WH2, source, deltaT, rho);
            }
        }
        else if (reactName == "oxidation")
        {
            const double WCO = 28.010;

            if (species == "O2")
            {
                const double WO2 = 31.999;
                const scalar YO2 = Y[this->owner().composition().slgThermo().carrier().species()["O2"]];
                const scalar cO2 = YO2*rho/WO2;

                const scalar source = rateCoeff * cO2;

                source4ySoot = kineticsCoeffs[5] * source * MWCarbon;

                if (gasSource() == true)
                {
                    // here O2 and CO are considered
                    this->fillInGasSourceTable(reactName+"1", "O2", WO2, source, deltaT, rho);
                    this->fillInGasSourceTable(reactName+"1", "CO", WCO, source, deltaT, rho);
                }

            }
            else if (species == "OH")
            {
                const double WOH = 17.007;
                const double WH  =  1.008;
                const scalar YOH = Y[this->owner().composition().slgThermo().carrier().species()["OH"]];
                const scalar cOH = YOH*rho/WOH;
                const scalar source = rateCoeff * cOH;

                source4ySoot = kineticsCoeffs[5] * source * MWCarbon;

                if (gasSource() == true)
                {
                    // here OH, CO, and H are considered
                    this->fillInGasSourceTable(reactName+"2", "OH", WOH, source, deltaT, rho);
                    this->fillInGasSourceTable(reactName+"2", "CO", WCO, source, deltaT, rho);
                    this->fillInGasSourceTable(reactName+"2", "H",  WH,  source, deltaT, rho);
                }
            }
        }
    }
}

template<class CloudType>
void Foam::TwoEquationModel<CloudType>::updateSootVariables
(
    const scalarField& Y,
    const scalar& T,
    scalar& ySoot,
    scalar& nSoot,
    const scalar& rho
)
{
    const scalar& deltaT = this->owner().time().deltaTValue();
    const double PI = constant::mathematical::pi;
    scalar tempNoUse(0.);

    // update diaSoot and surfSoot
    if (ySoot <= ROOTVSMALL || nSoot <= SMALL)
    {
        diaSoot() = 0.;
    }
    else
    {
        diaSoot() = pow(6.* ySoot/PI/rhoSoot()/nSoot/nSootScaling(), 1./3.);
    }

    surfSoot() = PI * pow(diaSoot(), 2.) * rho * nSoot * nSootScaling();

    // inception
    scalar source4ySootInC2H2(0.);
    scalar source4nSootInC2H2(0.);
    this->fillInSource("inception_C2H2", Y, T, ySoot, nSoot, rho, source4ySootInC2H2, source4nSootInC2H2);

    // cogulation
    scalar source4nSootCg(0.);
    this->fillInSource("cogulation_Cs", Y, T, ySoot, nSoot, rho, tempNoUse, source4nSootCg);

    // surface growth
    scalar source4ySootSgC2H2(0.);
    this->fillInSource("surfGrowth_C2H2", Y, T, ySoot, nSoot, rho, source4ySootSgC2H2, tempNoUse);

    // oxidation
    scalar source4ySootOxO2(0.);
    scalar source4ySootOxOH(0.);
    this->fillInSource("oxidation_O2", Y, T, ySoot, nSoot, rho, source4ySootOxO2, tempNoUse);
    this->fillInSource("oxidation_OH", Y, T, ySoot, nSoot, rho, source4ySootOxOH, tempNoUse);

    // total source
    scalar source4ySoot(0.);
    scalar source4nSoot(0.);
    source4ySoot = source4ySootInC2H2 + source4ySootSgC2H2 - source4ySootOxO2 - source4ySootOxOH;
    source4nSoot = source4nSootInC2H2 - source4nSootCg;

    // write soot contributions (at current time step)
    sootContribs_[0] =  (source4ySootInC2H2);
    sootContribs_[1] =  (source4ySootSgC2H2);
    sootContribs_[2] = -(source4ySootOxO2);
    sootContribs_[3] = -(source4ySootOxOH);
    sootContribs_[4] =  (source4nSootInC2H2/nSootScaling());
    sootContribs_[5] = -(source4nSootCg/nSootScaling());

    // update ySoot and nSoot
    ySoot += (source4ySoot*deltaT/rho);
    nSoot += (source4nSoot*deltaT/rho/nSootScaling());

    // limiter for ySoot and nSoot
    ySoot = max(ySoot, 0);
    nSoot = max(nSoot, 0);

}

template<class CloudType>
template<class DataType>
int Foam::TwoEquationModel<CloudType>::indexInTable
(
    const string& reactspecies,
    const Foam::List<Tuple2<string, DataType> >& Table
)
{
    int index = -1;
    for (int i = 0; i < Table.size(); ++i)
    {
        if (reactspecies == Table[i].first())
        {
            index = i;
            break;
        }
        else
        {
            continue;
        }
    }

    return index;
}


template<class CloudType>
void Foam::TwoEquationModel<CloudType>::fillInGasSourceTable
(
    const string& reactName,
    const string& speciesTag,
    const scalar& W,
    const scalar& source,
    const scalar& deltaT,
    const scalar& rho
)
{

    const string tagName = reactName + "_" + speciesTag;

    const label indexTagName = indexInTable(tagName, gasSourceTable());

    if (indexTagName == -1)
    {
        FatalErrorInFunction
            << "The reaction name of "
            << tagName
            << " cannot be found in gasSourceTable()\n"
            << abort(FatalError);
    }

    // gasSourceTable()[indexTagName].second() is the
    // stoichiometric coefficients for each species in the reaction
    YsourceTable_[indexTagName].second() = (gasSourceTable()[indexTagName].second() * source * W) * deltaT / rho;
}


template<class CloudType>
void Foam::TwoEquationModel<CloudType>::updateYs(scalarField& Y)
{
    // get the species and corresponding index
    for(int i = 0; i < YsourceTable().size(); i++)
    {
        const int delimiterLoc   = YsourceTable()[i].first().find("_");
        const string speciesName = YsourceTable()[i].first().substr(delimiterLoc+1, YsourceTable()[i].first().size());
        const int indexInYField  = this->owner().composition().slgThermo().carrier().species()[speciesName];
        Y[indexInYField] += YsourceTable()[i].second();
    }
}

template<class CloudType>
void Foam::TwoEquationModel<CloudType>::speciesInfo()
{

    Info << "----------------------------------------------------------------" <<nl;
    Info << "step\t\t" << "Yi\t" << "A\t" << "n\t" << "TAc\t" << "surf\t" << "CsCoeff" << endl;
    Info << "----------------------------------------------------------------" <<nl;

    for (int i = 0; i < this->kineticsCoeffs().size(); i++)
    {
        const int delimiterLoci = kineticsCoeffs()[i].first().find("_");
        const string reactTagi  = kineticsCoeffs()[i].first().substr(0, delimiterLoci);
        const string speciesi   = kineticsCoeffs()[i].first().substr(delimiterLoci+1, kineticsCoeffs()[i].first().size());
        if (reactTagi == "cogulation")
        {
            Info << reactTagi << '\t' << speciesi << '\t' << kineticsCoeffs()[i].second()[0] << nl;
        }
        else
        {
            Info << reactTagi << '\t' << speciesi << '\t' << kineticsCoeffs()[i].second()[0]
                << '\t' << kineticsCoeffs()[i].second()[2]
                << '\t' << kineticsCoeffs()[i].second()[3]
                << '\t' << kineticsCoeffs()[i].second()[4]
                << '\t' << kineticsCoeffs()[i].second()[5] << endl;
        }
    }
    Info << "----------------------------------------------------------------" <<nl;

}

template<class CloudType>
void Foam::TwoEquationModel<CloudType>::calSoot
(
     scalarField& Y,
     const scalar& T,
     scalar& ySoot,
     scalar& nSoot,
     scalar& sootVf,
     scalar& sootIm,
     const scalar& rho,
     scalarList& sootContribs
)
{

    if(debug_)
    {
        if(sootContribs.size() != this->sootContribs().size())
        {
            FatalErrorInFunction
                << "size of sootContribs from ThermoPopeParticle"
                << "is different from sootContribs_ in TwoEquationModel"
                << abort(FatalError);
        }
    }

    //- calculate soot-related variables
    if (T > 400)
    {
        //- updates ySoot and nSoot
        this->updateSootVariables(Y, T, ySoot, nSoot, rho);

        //- update the Y fields, species feedback
        if (gasSource() == true)
        {
            this->updateYs(Y);
        }

        //- calculate soot volume fraction based on ySoot
        sootVf = ySoot * rho / this->rhoSoot();

        //- update sootIm
        if (sootVf > sootLimit())
        {
            sootIm = 0.0;
        }
        else
        {
            sootIm = 1.0;
        }

        //- update soot contributions
        for(int i = 0; i < this->sootContribs().size(); i++)
        {
            sootContribs[i] = this->sootContribs()[i];
        }

        // for debugging
        if (debug_)
        {
            if (gasSource() == false)
            {
                for (int i = 0; i < YsourceTable_.size(); i++)
                {
                    if (YsourceTable_[i].second() != 0)
                    {
                        Info << "Error here, since gasSource take effects.\n";
                    }
                }
            }
            if (T < 270 || T > 3000)
            {
                Pout << "Temperature here is weird: " << T << endl;
            }

            if (sootVf > 1e-5)
            {
                Pout << "soot volme fraction is larger than 10 ppm" << endl;
            }
        }
    }
    else
    {
        //- update soot contributions (set to 0 when T < 400)
        //  considering that soot contributions are not accumulated
        //  it suggests once an iteration is done, and before starting
        //  the next iteration, they all should be reset.
        for(int i = 0; i < this->sootContribs().size(); i++)
        {
            sootContribs[i] = 0;
        }
    }

}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::TwoEquationModel<CloudType>::~TwoEquationModel()
{}

// ************************************************************************* //
