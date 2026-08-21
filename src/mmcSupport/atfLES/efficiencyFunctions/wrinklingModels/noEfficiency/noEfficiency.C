#include "noEfficiency.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

defineTypeNameAndDebug(noEfficiency, 0);
addToRunTimeSelectionTable(efficiencyFunction, noEfficiency, dictionary);


noEfficiency::noEfficiency
(
    const volVectorField& U,
    const volScalarField& TF,
    const volScalarField& eqRatio,
    volScalarField& gammafit,
    const volScalarField& uPrime,
    volScalarField& dynsl0,
    volScalarField& dyndeltal0
)
    : efficiencyFunction(U, TF, eqRatio, gammafit, uPrime, dynsl0, dyndeltal0)
{
}

noEfficiency::~noEfficiency()
{
}

void noEfficiency::correct()
{
    // leave default (1)
}

}
