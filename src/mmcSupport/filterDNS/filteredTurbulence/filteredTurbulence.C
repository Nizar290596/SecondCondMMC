/*---------------------------------------------------------------------------*\
                                       8888888888                              
                                       888                                     
                                       888                                     
  88888b.d88b.  88888b.d88b.   .d8888b 8888888  .d88b.   8888b.  88888b.d88b.  
  888 "888 "88b 888 "888 "88b d88P"    888     d88""88b     "88b 888 "888 "88b 
  888  888  888 888  888  888 888      888     888  888 .d888888 888  888  888 
  888  888  888 888  888  888 Y88b.    888     Y88..88P 888  888 888  888  888 
  888  888  888 888  888  888  "Y8888P 888      "Y88P"  "Y888888 888  888  888 
------------------------------------------------------------------------------- 

License
    This file is part of mmcFoam.

    mmcFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    mmcFoam is distributed in the hope that it will be useful, but 
    WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY 
    or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with mmcFoam. If not, see <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "filteredTurbulence.H"

Foam::filteredTurbulence::filteredTurbulence
(
    const fvMesh& mesh,
    const fvMesh& filterMesh,
    const volScalarField& rhoDNS,
    const volVectorField& UDNS,
    const bool filteredMMCCloudFields,
    const scalar Sc,
    const scalar Sct
)
:
    filter_(mesh,filterMesh),
    thermo_(psiReactionThermo::New(filterMesh)),
    Y_(thermo_.ref().composition().Y()),
    filteredMMCCloudFields_(filteredMMCCloudFields),
    Sc_(Sc),
    Sct_(Sct)
{
    // Note: The turbulence model is always generated from the filtered fields

    // Construct the turbulence model by first filtering the fields to the 
    // filter mesh
    UFilteredPtr_.reset
    (
        new volVectorField
        (
            IOobject
            (
                "UFiltered",
                filterMesh.time().timeName(),
                filterMesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            filter_.FavreAverage(UDNS,rhoDNS)
        )
    );

    rhoFilteredPtr_.reset
    (
        new volScalarField
        (
            IOobject
            (
                "rhoFiltered",
                filterMesh.time().timeName(),
                filterMesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            filter_.mean(rhoDNS)
        )
    );

    const volVectorField& UFiltered = UFilteredPtr_();

    phiFilteredPtr_.reset
    (
        new surfaceScalarField
        (
            IOobject
            (
                "phiFiltered",
                filterMesh.time().timeName(),
                filterMesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            fvc::flux(UFiltered)
        )
    );

    turbulence_.reset
    (
        compressible::turbulenceModel::New
        (
            rhoFilteredPtr_(),
            UFilteredPtr_(),
            phiFilteredPtr_(),
            thermo_()
        )
    );


    // The turbulent diffusivity fields must match the size of the field for 
    // the mmcCloud.
    // If filterDNSField is switched off, the D, Dt, and DEff field have the 
    // size of the mesh not filterMesh!
    if (filteredMMCCloudFields_)
    {
        D_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "D",
                    filterMesh.time().timeName(),
                    filterMesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                filterMesh,
                dimensionedScalar("0",dimViscosity,0)
            )
        );

        Dt_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "Dt",
                    filterMesh.time().timeName(),
                    filterMesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                filterMesh,
                dimensionedScalar("0",dimViscosity,0)
            )
        );

        DEff_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "DEff",
                    filterMesh.time().timeName(),
                    filterMesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                filterMesh,
                dimensionedScalar("0",dimViscosity,0)
            )
        );

        DeltaE_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "DeltaE",
                    filterMesh.time().timeName(),
                    filterMesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                filterMesh,
                dimensionedScalar("0",dimLength,0)
            )
        );
        if (filterMesh.objectRegistry::foundObject<volScalarField>("delta"))
        {
            const volScalarField& delta =
                filterMesh.objectRegistry::lookupObject<volScalarField>("delta");
            DeltaE_.ref() = delta;
        }
        else
        {
            Info << "The field 'delta' was not found in the object registry." << nl
                 << "Calculating DeltaE based on the cube root of the cell volume." 
                 << endl;
            volScalarField& DeltaE = DeltaE_.ref();
            forAll(DeltaE,i)
            {
                DeltaE[i] = Foam::cbrt(filterMesh.V()[i]);
            }
        }
    }
    else
    {
        D_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "D",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar("0",dimViscosity,0)
            )
        );

        Dt_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "Dt",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar("0",dimViscosity,0)
            )
        );

        DEff_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "DEff",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar("0",dimViscosity,0)
            )
        );

        DeltaE_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "DeltaE",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar("0",dimLength,SMALL)
            )
        );

        // We need the deltaE of the filtered mesh, thus construct that first
        // and than map it to the DNS mesh.
        volScalarField DeltaEFiltered
        (
            IOobject
            (
                "DeltaE",
                filterMesh.time().timeName(),
                filterMesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            filterMesh,
            dimensionedScalar("0",dimLength,SMALL)
        );
        forAll(DeltaEFiltered,i)
        {
            DeltaEFiltered[i] = Foam::cbrt(filterMesh.V()[i]);
        }

        DeltaE_.ref() = filter_.mapFilterFieldToDNS(DeltaEFiltered);
    }
}


void Foam::filteredTurbulence::correct
(
    const volVectorField& UDNS,
    const volScalarField& rhoDNS,
    const volScalarField& pDNS,
    const volScalarField& TDNS,
    const PtrList<volScalarField>& Y
)
{
    // Filter the UDNS, rhoDNS
    filter_.updateField(UDNS,rhoDNS,UFilteredPtr_.ref());
    filter_.updateField(rhoDNS,rhoFilteredPtr_.ref());
    const volVectorField& UFiltered = UFilteredPtr_();
    phiFilteredPtr_.ref() = fvc::flux(UFiltered);

    // First update the thermo model
    // =============================

    // Update all species mass fractions with the new filtered value
    forAll(Y,i)
    {
        const volScalarField& Yi = Y[i];
        volScalarField& YiFiltered = Y_[i];
        filter_.updateField(Yi,rhoDNS,YiFiltered);
    }

    auto TFiltered = filter_.FavreAverage(TDNS,rhoDNS);
    auto pFiltered = filter_.mean(pDNS);
    thermo_->he() = thermo_->he(pFiltered(),TFiltered());
    thermo_->correct();

    // Now update the turbulence model
    turbulence_->correct();

    // Calculate the new D, Dt, and DEff values
    if (filteredMMCCloudFields_)
    {
        D_.ref() = turbulence_->nu()/Sc_;
        Dt_.ref() = turbulence_->nut()/Sct_;
        DEff_.ref() = D_() + Dt_();
    }
    else
    {
        volScalarField temp1 = turbulence_->nu()/Sc_;
        D_.ref() = filter_.mapFilterFieldToDNS(temp1);

        temp1 = turbulence_->nut()/Sct_;
        Dt_.ref() = filter_.mapFilterFieldToDNS(temp1);

        DEff_.ref() = D_();
    }
}


// ************************************************************************* //
