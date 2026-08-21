#include "efficiencyFunction.H"

namespace Foam
{

    defineTypeNameAndDebug(efficiencyFunction, 0);
    defineRunTimeSelectionTable(efficiencyFunction, dictionary);

    autoPtr<efficiencyFunction> efficiencyFunction::New
    (
        const volVectorField& U,
        const volScalarField& TF,
        const volScalarField& eqRatio,
        volScalarField& gammafit,
        const volScalarField& uPrime,
        volScalarField& dynsl0,
        volScalarField& dyndeltal0
    )
    {
        word typeName;

        // Enclose the creation of the dictionary to ensure it is
        // deleted before the model is created otherwise the dictionary
        // is entered in the database twice
        {
            IOdictionary propertiesDict
                (
                    IOobject
                    (
                        "efficiencyFunctionProperties",
                        U.time().constant(),
                        U.db(),
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE
                    )
                );
            
            propertiesDict.lookup("efficiencyFunction") 
                >> typeName;
        }

        auto cstrIter =
            dictionaryConstructorTablePtr_->find(typeName); // Updated for the ofv2212
        
        if (cstrIter == dictionaryConstructorTablePtr_->end())
        {
            FatalErrorIn
                (
                    "efficiencyFunction::New()"
                )   << "Unknown efficiency function " << typeName
                    << endl << endl
                    << "Valid efficiency functions are :" << endl
                    << dictionaryConstructorTablePtr_->toc()
                    << exit(FatalError);
        }

        return autoPtr<efficiencyFunction>(cstrIter()(U, TF, eqRatio, gammafit, uPrime, dynsl0, dyndeltal0));
        
    }



    efficiencyFunction::efficiencyFunction
    (
        const volVectorField& U,
        const volScalarField& TF,
        const volScalarField& eqRatio,
        volScalarField& gammafit,
        const volScalarField& uPrime,
        volScalarField& dynsl0,
        volScalarField& dyndeltal0
    )
        : IOdictionary
    (
        IOobject
        (
            "efficiencyFunctionProperties",
            U.time().constant(),
            U.db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
          U_(U),
          TF_(TF),
          eqRatio_(eqRatio),
          gammafit_(gammafit),
          uPrime_(uPrime),
          dynsl0_(dynsl0),
          dyndeltal0_(dyndeltal0),
          efficiency_
        (
            IOobject
            (
                "efficiency",
                U.time().timeName(),
                U.mesh(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            U.mesh(),
            dimensionedScalar("", dimless, 1.0)
        )
    {
    }

    efficiencyFunction::~efficiencyFunction()
    {
    }


}
