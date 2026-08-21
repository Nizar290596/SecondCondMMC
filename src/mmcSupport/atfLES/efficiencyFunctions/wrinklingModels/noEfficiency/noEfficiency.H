#ifndef NOEFFICIENCY_H
#define NOEFFICIENCY_H

#include "efficiencyFunction.H"

namespace Foam
{

    class noEfficiency
    : public efficiencyFunction
    {

    public:
        TypeName("none");

        noEfficiency(const volVectorField& U,
                     const volScalarField& TF,
                     const volScalarField& eqRatio,
                     volScalarField& gammafit,
                     const volScalarField& uPrime,
                     volScalarField& dynsl0,
                     volScalarField& dyndeltal0);
        ~noEfficiency();

        virtual void correct();
    };

}

#endif
