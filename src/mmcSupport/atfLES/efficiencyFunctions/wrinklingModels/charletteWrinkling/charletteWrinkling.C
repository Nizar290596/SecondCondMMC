#include "charletteWrinkling.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

defineTypeNameAndDebug(charletteWrinkling, 0);
addToRunTimeSelectionTable(efficiencyFunction, charletteWrinkling, dictionary);


charletteWrinkling::charletteWrinkling
(
    const volVectorField& U,
    const volScalarField& TF,
    const volScalarField& eqRatio,
    volScalarField& gammafit,
    const volScalarField& uPrime,
    volScalarField& dynsl0,
    volScalarField& dyndeltal0
)
    : efficiencyFunction(U, TF, eqRatio, gammafit, uPrime, dynsl0, dyndeltal0),
// todo: make Ck and b not readable
      Ck_(readScalar(lookup("Ck"))),
      b_(readScalar(lookup("b"))),
      gamma_(readScalar(lookup("gamma"))),
      uPrimeCoeff_(readScalar(lookup("uPrimeCoeff"))),
      constLamProperties_(readBool(lookup("constLamProperties")))
{
      if(constLamProperties_) {
          scalar fixedsl0_ = readScalar(lookup("fixedsl0"));
          scalar fixeddeltal0_ = readScalar(lookup("fixeddeltal0"));
          forAll(dynsl0_,cellI)
          {
               dynsl0_[cellI] = fixedsl0_;
               dyndeltal0_[cellI] = fixeddeltal0_;
          }
      } else {
          fileName   Lamfile(fileName(lookup("LamFileName")).expand());
          #include "readLamTable.H"
      }
}

charletteWrinkling::~charletteWrinkling()
{
}

void charletteWrinkling::correct()
{
// nachdenken: put this into TF > 1
    if(!constLamProperties_) {
        int colEqRatio = 0;
        int coldeltal0 = 1;
        int colsl0 = 2;
        forAll(eqRatio_,cellI)
        {
        int rowI = 0;
//        if(eqRatio_[cellI] > lamtable[rowI][colEqRatio] && eqRatio_[cellI] < 1.0) {
        if(eqRatio_[cellI] > lamtable[0][colEqRatio]) {
            while(eqRatio_[cellI] > lamtable[rowI][colEqRatio]) {
                if(eqRatio_[cellI] <= lamtable[rowI + 1][colEqRatio]) {
                            //- interpolate delta_0 and s_0 from table
                dyndeltal0_[cellI] = (lamtable[rowI][coldeltal0] + (lamtable[rowI+1][coldeltal0]-lamtable[rowI][coldeltal0])*\
                                                     (eqRatio_[cellI]-lamtable[rowI][colEqRatio])/(lamtable[rowI+1][colEqRatio]-lamtable[rowI][colEqRatio]))/1000.0;
                dynsl0_[cellI] = (lamtable[rowI][colsl0] + (lamtable[rowI+1][colsl0]-lamtable[rowI][colsl0])*\
                                                 (eqRatio_[cellI]-lamtable[rowI][colEqRatio])/(lamtable[rowI+1][colEqRatio]-lamtable[rowI][colEqRatio]))/100.0;
                }
                rowI = rowI + 1;
            }
          }
        else {
                Info << "eqRatio is lower than the values presented in the table" << nl << endl;
             }
        }
    }
// todo: nullabfragen?
// boundary treatment is not correct? Call this function after U boundaries are updated
//    uPrime_ = 2.0 * mag(pow(delta(), 3) * fvc::laplacian(fvc::curl(U_)));
//    uPrime_ = uPrime2Mean(0);
//    uPrime_ = sqrt((2.0/3.0)*turbulence->k());

//    volScalarField RL = TF_ * deltal0_ / deltal0_;
    volScalarField RL = TF_;
    scalar eps = 1e-05;
    volScalarField RV = uPrimeCoeff_ * uPrime_ / dynsl0_;

    forAll(efficiency_,cellI) {
        if(TF_[cellI] < 1.001) {
            efficiency_[cellI] = 1.0;
        } else {
            scalar a = 0.6 + 0.2 * exp(-0.1*RV[cellI]) - 0.2 * exp(-0.01*RL[cellI]);
            scalar fu = max(eps,4.0*pow(27.0*Ck_/110.0,0.5) * (18.0*Ck_/55.0) * pow(RV[cellI],2.0));
// todo: access pi
            scalar fdelta = max(eps,pow(27.0*Ck_*pow(3.14,4.0/3.0)/110.0 * (pow(RL[cellI],4.0/3.0)-1.0),0.5));
            scalar Resgs = 4.0 * RL[cellI] * RV[cellI];
            scalar fre = max(eps,pow(9./55. * exp(-1.5*Ck_*pow(3.14,4./3.)/max(eps,Resgs)),0.5) * pow(Resgs,0.5));

            gammafit_[cellI] = pow(max(0.0,pow(max(0.0,pow(pow(fu,-a)+pow(fdelta,-a),-1/a)),-b_) + pow(fre,-b_)),-1/b_);
//            gammafit_[cellI] = pow(max(0.0,pow(max(0.0,pow(pow(fu,2)+pow(fdelta,2),2)),2) + pow(fre,2)),2);
            efficiency_[cellI] = pow(1.0 + min(RL[cellI],gammafit_[cellI]*RV[cellI]),gamma_);
        }   
    }   
    efficiency_.correctBoundaryConditions();
    gammafit_.correctBoundaryConditions();
}

}
