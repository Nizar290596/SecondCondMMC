#ifndef CHARLETTEWRINKLING_H
#define CHARLETTEWRINKLING_H

#define nRow_ 3 // Updated for flamelet
#define nCol_ 3

#include "efficiencyFunction.H"

namespace Foam
{

    class charletteWrinkling
    : public efficiencyFunction
    {

        scalar Ck_, b_, gamma_, uPrimeCoeff_;
        Switch constLamProperties_;
        double lamtable[nRow_][nCol_];

    public:
        TypeName("charletteWrinkling");

        charletteWrinkling(const volVectorField& U,
                           const volScalarField& TF,
                           const volScalarField& eqRatio,
                           volScalarField& gammafit,
                           const volScalarField& uPrime,
                           volScalarField& dynsl0,
                           volScalarField& dyndeltal0);
        ~charletteWrinkling();

        virtual void correct();
    };

}

#endif
