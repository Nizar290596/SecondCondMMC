#include "colinWrinkling.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

defineTypeNameAndDebug(colinWrinkling, 0);
addToRunTimeSelectionTable(efficiencyFunction, colinWrinkling, dictionary);


colinWrinkling::colinWrinkling
(
    const volVectorField& U,
    const volScalarField& TF,
    const volScalarField& eqRatio,
    volScalarField& gammafit,
    const volScalarField& uPrime
)
    : efficiencyFunction(U, TF, eqRatio, gammafit, uPrime),
      uPrime_
    (
        IOobject
        (
            "uPrime",
            U_.time().timeName(),
            U_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh(),
        dimensionedScalar("", dimLength/dimTime, 0.0)
    ),
      alpha_(readScalar(lookup("alpha"))),
      deltal0_(lookup("deltal0")),
      sl0_(lookup("sl0"))
{
}

colinWrinkling::~colinWrinkling()
{
}

void colinWrinkling::correct()
{
// check here if RL=TF_?
    volScalarField RL = TF_;
    volScalarField RV = uPrime_ / sl0_;

    efficiency_ = 
        (1.0 + alpha_ * 0.75*exp(-1.2/pow(RV, 0.3))*pow(RL, 2./3.) * RV)
        /
        (1.0 + alpha_ * 0.75*exp(-1.2/pow(RV, 0.3))*pow(RL/TF_, 2./3.) * RV);
}

}
