#ifndef EFFICIENCYFUNCTION_H
#define EFFICIENCYFUNCTION_H

#include "fvCFD.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

class efficiencyFunction
: public IOdictionary
{

protected:
    const volVectorField& U_;
    const volScalarField& TF_;
    const volScalarField& eqRatio_;
    volScalarField& gammafit_;
    const volScalarField& uPrime_;
    volScalarField& dynsl0_;
    volScalarField& dyndeltal0_;

    volScalarField efficiency_;

public:
    TypeName("efficiencyFunction");

    declareRunTimeSelectionTable
    (
        autoPtr,
        efficiencyFunction,
        dictionary,
        (
            const volVectorField& U,
            const volScalarField& TF,
            const volScalarField& eqRatio,
            volScalarField& gammafit,
            const volScalarField& uPrime,
            volScalarField& dynsl0,
            volScalarField& dyndeltal0
        ),
        (U, TF, eqRatio, gammafit, uPrime, dynsl0, dyndeltal0)
    );

    //- Return a reference to the selected LES model
    static autoPtr<efficiencyFunction> New
    (
        const volVectorField& U,
        const volScalarField& TF,
        const volScalarField& eqRatio,
        volScalarField& gammafit,
        const volScalarField& uPrime,
        volScalarField& dynsl0,
        volScalarField& dyndeltal0
    );



    efficiencyFunction(const volVectorField& U,
                       const volScalarField& TF,
                       const volScalarField& eqRatio,
                       volScalarField& gammafit,
                       const volScalarField& uPrime,
                       volScalarField& dynsl0,
                       volScalarField& dyndeltal0);
    ~efficiencyFunction();

    inline const volScalarField& eff() const 
        { return efficiency_; }

    virtual void correct() =0;
};

}

#endif
